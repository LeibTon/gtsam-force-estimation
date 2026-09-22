"""Run estimate_experiment.py over every single-force case listed in BENCHMARK_EXPERIMENT_SCOPE.md."""
import csv
import json
import statistics
from pathlib import Path

import numpy as np

import tendon_robot as tr
from estimate_experiment import estimate_experiment
from experiment_parameters import (GRAVITY_DIRECTION, GRAVITY_M_S2, ROD_MASS_PER_LENGTH_KG_M,
                                   ROUTING_ANGLE_OFFSET_RAD, TENDON_FRICTION_COEFFICIENT,
                                   YOUNGS_MODULUS_PA)

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent


def build_config():
    config = tr.TendonRobotConfig()
    config.youngs_modulus = YOUNGS_MODULUS_PA
    config.shear_modulus = YOUNGS_MODULUS_PA / (2.0 * (1.0 + 0.30))
    config.tendon_friction_coefficient = TENDON_FRICTION_COEFFICIENT
    config.gravity_magnitude = GRAVITY_M_S2
    config.gravity_direction = GRAVITY_DIRECTION
    config.rod_mass_per_length = ROD_MASS_PER_LENGTH_KG_M
    config.angle_params = [tr.RoutingFunctionParams(angle_offset=ROUTING_ANGLE_OFFSET_RAD, total_angle=0.0)]
    config.num_basis_functions = 1
    config.basis_sigma = 0.012
    return config

DRIVE_LENGTHS = (3, 4)
REQUESTED_WEIGHTS = {
    6: (20, 30, 40, 50),
    8: (20, 30, 40, 50),
    10: (10, 20, 30),
}
DATASET_FOR_DISC = {
    6: "disc6",
    8: "disc8",
    10: "disc10",
}

CSV_COLUMNS = [
    "trial", "disc", "weight_g", "status",
    "tendon_tension_n", "actual_force_n", "actual_location_m",
    "estimated_force_n", "estimated_location_m",
    "estimated_force_var_n2", "estimated_location_var_m2", "shape_rmse_mm",
    "force_error_n", "force_percentage_error", "force_location_error_mm",
    "force_direction_error_deg",
]

SUMMARY_FIELDS = [
    ("force_error_n", "Force magnitude error", "N"),
    ("force_percentage_error", "Force % error", "%"),
    ("force_location_error_mm", "Location error", "mm"),
    ("force_direction_error_deg", "Direction error", "deg"),
    ("shape_rmse_mm", "Shape RMSE", "mm"),
]


def transverse_direction_error_deg(rotation, estimated_force, actual_force) -> float:
    estimated_xy = (rotation.T @ estimated_force)[:2]
    actual_xy = (rotation.T @ actual_force)[:2]
    if np.linalg.norm(estimated_xy) <= 1e-12 or np.linalg.norm(actual_xy) <= 1e-12:
        return float("nan")

    angle_difference = (np.arctan2(estimated_xy[1], estimated_xy[0])
                        - np.arctan2(actual_xy[1], actual_xy[0]))
    return abs(float(np.degrees(np.arctan2(np.sin(angle_difference), np.cos(angle_difference)))))


def shape_rmse_mm(estimated_shape, measured_shape) -> float:
    to_base_relative = lambda points: points - points[0]
    residual = to_base_relative(estimated_shape) - to_base_relative(measured_shape)
    return float(np.sqrt(np.mean(np.sum(residual ** 2, axis=1))) * 1e3)


def location_variance_m2(gamma_cov: np.ndarray, beta: float, rod_length: float) -> float:
    """Propagate beta's variance through location = (L/2)(1 + tanh(beta))."""
    jacobian = rod_length / 2.0 * (1.0 - np.tanh(beta) ** 2)
    return float(jacobian ** 2 * gamma_cov[3, 3])


def force_magnitude_variance(force_cov: np.ndarray, force: np.ndarray) -> float:
    """Propagate the 3x3 force covariance through magnitude = ||force||."""
    magnitude = np.linalg.norm(force)
    if magnitude <= 1e-12:
        return float("nan")
    jacobian = force / magnitude
    return float(jacobian @ force_cov @ jacobian)


def evaluate_single_force(config, folder: Path, disc: int, weight_g: float):
    result = estimate_experiment(config, folder, [(disc, weight_g)])
    estimated_force, estimated_location = result["bases"][0]
    actual_force, actual_location = result["actual_loads"][0]

    estimated_magnitude = float(np.linalg.norm(estimated_force))
    actual_magnitude = float(np.linalg.norm(actual_force))

    loaded_node_rotation = np.asarray(
        result["solution"].backbone_pose_mean[result["nodes"][disc - 1]])[:3, :3]
    direction_error_deg = transverse_direction_error_deg(
        loaded_node_rotation, estimated_force, actual_force)

    force_error_n = abs(estimated_magnitude - actual_magnitude)

    gamma_mean = np.asarray(result["solution"].gamma_mean)
    gamma_cov = np.asarray(result["solution"].gamma_cov)
    force_cov = gamma_cov[:3, :3]

    return {
        "tendon_tension_n": result["tension"],
        "actual_force_n": actual_magnitude,
        "actual_location_m": actual_location,
        "estimated_force_n": estimated_magnitude,
        "estimated_location_m": estimated_location,
        "estimated_force_var_n2": force_magnitude_variance(force_cov, estimated_force),
        "estimated_location_var_m2": location_variance_m2(gamma_cov, gamma_mean[3], config.rod_length),
        "estimated_force_cov": force_cov.tolist(),
        "shape_rmse_mm": shape_rmse_mm(result["estimated_shape"], result["measured_shape"]),
        "force_error_n": force_error_n,
        "force_percentage_error": 100.0 * force_error_n / actual_magnitude,
        "force_location_error_mm": abs(estimated_location - actual_location) * 1e3,
        "force_direction_error_deg": direction_error_deg,
    }


def run_trial(config, disc, drive, weight):
    name = f"lact{drive}_weight{weight}_disc{disc}"
    folder = ROOT / "experimental_data" / DATASET_FOR_DISC[disc] / name
    row = {"trial": name, "dataset": DATASET_FOR_DISC[disc],
           "drive_mm": drive, "disc": disc, "weight_g": weight, "status": "missing"}
    if not folder.is_dir():
        return row, None

    try:
        result = evaluate_single_force(config, folder, disc, weight)
    except ValueError:
        row["status"] = "incomplete"
        return row, None

    row["status"] = "ok"
    row.update({k: v for k, v in result.items() if k in CSV_COLUMNS})
    record = {**row, **result}
    return row, record


def print_summary(rows):
    ok_rows = [row for row in rows if row["status"] == "ok"]
    print(f"\nn = {len(ok_rows)} trials")
    for field, label, unit in SUMMARY_FIELDS:
        values = [row[field] for row in ok_rows if row.get(field) not in (None, "", "nan")]
        if not values:
            continue
        mean = statistics.mean(values)
        median = statistics.median(values)
        print(f"{label:24s} mean={mean:8.4f} {unit:<3s} median={median:8.4f} {unit}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    config = build_config()

    all_cases = [(disc, drive, weight)
                for disc, weights in REQUESTED_WEIGHTS.items()
                for drive in DRIVE_LENGTHS
                for weight in weights]
    total = len(all_cases)

    rows, records = [], []
    for i, (disc, drive, weight) in enumerate(all_cases, start=1):
        print(f"[{i}/{total}] lact{drive}_weight{weight}_disc{disc}")
        row, record = run_trial(config, disc, drive, weight)
        rows.append(row)
        if record is not None:
            records.append(record)

    json_output = OUT_DIR / "single_force_benchmark.json"
    json_output.write_text(json.dumps({
        "youngs_modulus_pa": YOUNGS_MODULUS_PA,
        "tendon_friction_coefficient": TENDON_FRICTION_COEFFICIENT,
        "records": records,
        "rows": rows,
    }, indent=2))

    csv_output = OUT_DIR / "single_force_benchmark.csv"
    with csv_output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(json_output)
    print(csv_output)
    completed = sum(row["status"] == "ok" for row in rows)
    print(f"completed={completed} missing_or_incomplete={len(rows) - completed}")

    print_summary(rows)


if __name__ == "__main__":
    main()
