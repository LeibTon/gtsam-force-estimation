"""Joint calibration of Young's modulus and tendon-guide friction coefficient against experimental_data/calibration_data."""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import tendon_robot as tr
from camera_transform import apply_transform, camera_to_model_transform, load_arm_points_mm
from experiment_parameters import GRAVITY_DIRECTION, GRAVITY_M_S2, ROD_MASS_PER_LENGTH_KG_M, ROUTING_ANGLE_OFFSET_RAD

CALIBRATION_DATA_DIR = Path(__file__).resolve().parents[2] / "experimental_data" / "calibration_data"
INCLUDED_DRIVE_LENGTHS = {2.0, 2.5, 3.0, 3.5, 4.0}
INITIAL_GUESSES = [(50.0, 1.0), (150.0, 0.15), (20.0, 0.0)]  # (E_gpa, friction)


def load_drive_trial(drive_dir, R, t):
    records = sorted((drive_dir / "annotation_records").iterdir())

    tensions, shapes = [], []
    for record in records:
        force_path = record / "force_sensor.json"
        points_path = record / "arm_points.json"
        if not (force_path.exists() and points_path.exists()):
            continue

        force_data = json.loads(force_path.read_text())
        if not force_data.get("available", True):
            continue
        tensions.append(force_data["force_n"])
        shapes.append(apply_transform(load_arm_points_mm(points_path), R, t))

    tension = float(np.mean(tensions))
    shape = np.mean(np.stack(shapes, axis=0), axis=0)
    return tension, shape


def load_all_trials():
    R, t = camera_to_model_transform()

    drive_lengths, tensions, shapes = [], [], []
    for drive_dir in sorted(CALIBRATION_DATA_DIR.glob("lact_*")):
        drive_length = float(drive_dir.name.replace("lact_", ""))
        if drive_length not in INCLUDED_DRIVE_LENGTHS:
            continue

        tension, shape = load_drive_trial(drive_dir, R, t)
        drive_lengths.append(drive_length)
        tensions.append(tension)
        shapes.append(shape)

    return drive_lengths, tensions, shapes


def build_config(youngs_modulus_pa, friction_coefficient):
    config = tr.TendonRobotConfig()
    config.youngs_modulus = youngs_modulus_pa
    config.shear_modulus = youngs_modulus_pa / (2.0 * (1.0 + 0.30))
    config.tendon_friction_coefficient = friction_coefficient
    config.angle_params = [tr.RoutingFunctionParams(angle_offset=ROUTING_ANGLE_OFFSET_RAD, total_angle=0.0)]
    config.gravity_magnitude = GRAVITY_M_S2
    config.gravity_direction = GRAVITY_DIRECTION
    config.rod_mass_per_length = ROD_MASS_PER_LENGTH_KG_M
    return config


def predict_shape(config, tension):
    estimator = tr.GaussianBasisForceEstimator(config)
    gamma_zero = np.zeros(4 * config.num_basis_functions)
    solution = estimator.step_forward(np.array([tension]), gamma_zero)
    disc_pose_idx = estimator.tendon_disc_config.disc_pose_idx
    return np.array([solution.backbone_pose_mean[i][:3, 3] for i in disc_pose_idx])


def residuals(params, tensions, measured_shapes):
    config = build_config(params[0] * 1e9, params[1])
    errors = [(predict_shape(config, tension) - measured_shape).ravel()
             for tension, measured_shape in zip(tensions, measured_shapes)]
    return np.concatenate(errors)


def calibrate(tensions, measured_shapes):
    best = None
    for e_gpa_0, friction_0 in INITIAL_GUESSES:
        result = least_squares(
            residuals, x0=[e_gpa_0, friction_0],
            args=(tensions, measured_shapes),
            bounds=([1.0, 0.0], [500.0, 5.0]),
        )
        rmse_mm = float(np.sqrt(np.mean(residuals(result.x, tensions, measured_shapes) ** 2))) * 1000.0
        candidate = {
            "start_E_gpa": e_gpa_0,
            "start_friction": friction_0,
            "E_gpa": float(result.x[0]),
            "friction": float(result.x[1]),
            "shape_rmse_mm": rmse_mm,
            "cost": float(result.cost),
            "nfev": int(result.nfev),
            "status": int(result.status),
        }
        if best is None or candidate["cost"] < best["cost"]:
            best = candidate

    return best


def main():
    drive_lengths, tensions, measured_shapes = load_all_trials()
    print(f"Loaded {len(drive_lengths)} trials (drive lengths: {drive_lengths})")
    print(f"Tensions (N): {[round(t, 4) for t in tensions]}")

    best = calibrate(tensions, measured_shapes)

    print("\nBest joint fit:")
    print(f"  Young's modulus:      {best['E_gpa']:.6f} GPa")
    print(f"  Friction coefficient: {best['friction']:.6f}")
    print(f"  Pooled shape RMSE:    {best['shape_rmse_mm']:.4f} mm")
    print(f"  (started from E={best['start_E_gpa']} GPa, mu={best['start_friction']})")

    output_path = Path(__file__).parent / "calibration_experiment_result.json"
    output_path.write_text(json.dumps({
        "drive_lengths": drive_lengths,
        "tensions_n": tensions,
        "best_solution": best,
        "routing_angle_offset_rad": ROUTING_ANGLE_OFFSET_RAD,
    }, indent=2))
    print(f"\nSaved result to {output_path}")


if __name__ == "__main__":
    main()
