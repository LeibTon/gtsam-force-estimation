"""Run estimate_experiment.py over every two-force case listed in BENCHMARK_EXPERIMENT_SCOPE.md."""
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
EXPERIMENTAL_DATA = ROOT / "experimental_data"
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
    config.num_basis_functions = 2
    config.basis_sigma = 0.006
    return config

# (drive, disc1, weight1, disc2, weight2)
LACT0_CASES = [
    (0, 6, 50, 8, 50), (0, 6, 20, 8, 50), (0, 6, 50, 8, 20),
    (0, 6, 30, 10, 30), (0, 6, 50, 10, 50), (0, 6, 10, 10, 20),
    (0, 6, 40, 10, 50), (0, 6, 50, 10, 40), (0, 6, 60, 10, 20),
]
LACT3_CASES = [
    (3, 6, 30, 8, 30), (3, 6, 30, 8, 10),
    (3, 6, 50, 10, 50), (3, 6, 40, 10, 50), (3, 6, 50, 10, 10),
]
ALL_CASES = LACT0_CASES + LACT3_CASES

CSV_COLUMNS = [
    "trial", "drive_mm", "shape_rmse_mm",
    "disc", "weight_g", "tendon_tension_n", "actual_force_n", "actual_location_m",
    "estimated_force_n", "estimated_location_m",
    "estimated_force_var_n2", "estimated_location_var_m2",
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


def trial_folder(drive: int, disc1: int, weight1: int, disc2: int, weight2: int) -> tuple[Path, str]:
    dataset = "disc6_disc8" if {disc1, disc2} == {6, 8} else "disc6_disc10"
    weight_disc6 = weight1 if disc1 == 6 else weight2
    other_disc = disc2 if disc1 == 6 else disc1
    weight_other = weight2 if disc1 == 6 else weight1
    name = f"lact{drive}_weight{weight_disc6}_disc6_weight{weight_other}_disc{other_disc}"
    return EXPERIMENTAL_DATA / dataset / name, name


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


def match_bases_to_loads(bases, actual_locations):
    remaining = list(range(len(bases)))
    matches = []
    for actual_location in actual_locations:
        best_index = min(remaining, key=lambda i: abs(bases[i][1] - actual_location))
        remaining.remove(best_index)
        matches.append(best_index)
    return matches


def location_variance_m2(gamma_cov: np.ndarray, beta: float, beta_idx: int, rod_length: float) -> float:
    """Propagate beta's variance through location = (L/2)(1 + tanh(beta))."""
    jacobian = rod_length / 2.0 * (1.0 - np.tanh(beta) ** 2)
    return float(jacobian ** 2 * gamma_cov[beta_idx, beta_idx])


def force_magnitude_variance(force_cov: np.ndarray, force: np.ndarray) -> float:
    """Propagate the 3x3 force covariance through magnitude = ||force||."""
    magnitude = np.linalg.norm(force)
    if magnitude <= 1e-12:
        return float("nan")
    jacobian = force / magnitude
    return float(jacobian @ force_cov @ jacobian)


def evaluate_two_force(config, folder: Path, disc1: int, weight1: int, disc2: int, weight2: int):
    loads = [(disc1, weight1), (disc2, weight2)]
    result = estimate_experiment(config, folder, loads)

    actual_locations = [location for _, location in result["actual_loads"]]
    matches = match_bases_to_loads(result["bases"], actual_locations)
    rmse = shape_rmse_mm(result["estimated_shape"], result["measured_shape"])

    gamma_mean = np.asarray(result["solution"].gamma_mean)
    gamma_cov = np.asarray(result["solution"].gamma_cov)

    per_force = []
    for (disc, weight_g), (actual_force, actual_location), basis_index in zip(
            loads, result["actual_loads"], matches):
        estimated_force, estimated_location = result["bases"][basis_index]
        estimated_magnitude = float(np.linalg.norm(estimated_force))
        actual_magnitude = float(np.linalg.norm(actual_force))

        loaded_node_rotation = np.asarray(
            result["solution"].backbone_pose_mean[result["nodes"][disc - 1]])[:3, :3]
        direction_error_deg = transverse_direction_error_deg(
            loaded_node_rotation, estimated_force, actual_force)

        force_error_n = abs(estimated_magnitude - actual_magnitude)

        alpha_idx = 4 * basis_index
        beta_idx = 4 * basis_index + 3
        force_cov = gamma_cov[alpha_idx:alpha_idx + 3, alpha_idx:alpha_idx + 3]

        per_force.append({
            "disc": disc,
            "weight_g": weight_g,
            "tendon_tension_n": result["tension"],
            "actual_force_n": actual_magnitude,
            "actual_location_m": actual_location,
            "estimated_force_n": estimated_magnitude,
            "estimated_location_m": estimated_location,
            "estimated_force_var_n2": force_magnitude_variance(force_cov, estimated_force),
            "estimated_location_var_m2": location_variance_m2(gamma_cov, gamma_mean[beta_idx], beta_idx, config.rod_length),
            "estimated_force_cov": force_cov.tolist(),
            "force_error_n": force_error_n,
            "force_percentage_error": 100.0 * force_error_n / actual_magnitude,
            "force_location_error_mm": abs(estimated_location - actual_location) * 1e3,
            "force_direction_error_deg": direction_error_deg,
        })

    return rmse, per_force


def run_trial(config, drive, disc1, weight1, disc2, weight2):
    folder, name = trial_folder(drive, disc1, weight1, disc2, weight2)
    base_row = {"trial": name, "drive_mm": drive}
    if not folder.is_dir():
        return [{**base_row, "status": "missing"}], []

    try:
        rmse, per_force = evaluate_two_force(config, folder, disc1, weight1, disc2, weight2)
    except ValueError:
        return [{**base_row, "status": "incomplete"}], []

    rows, records = [], []
    for force_result in per_force:
        row = {**base_row, "status": "ok", "shape_rmse_mm": rmse}
        row.update({k: v for k, v in force_result.items() if k in CSV_COLUMNS})
        rows.append(row)
        records.append({**row, **force_result})
    return rows, records


def print_summary(rows):
    ok_rows = [row for row in rows if row["status"] == "ok"]
    print(f"\nn = {len(ok_rows)} force estimates")
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

    total = len(ALL_CASES)
    rows, records = [], []
    for i, (drive, disc1, weight1, disc2, weight2) in enumerate(ALL_CASES, start=1):
        _, name = trial_folder(drive, disc1, weight1, disc2, weight2)
        print(f"[{i}/{total}] {name}")
        trial_rows, trial_records = run_trial(config, drive, disc1, weight1, disc2, weight2)
        rows.extend(trial_rows)
        records.extend(trial_records)

    json_output = OUT_DIR / "two_force_benchmark.json"
    json_output.write_text(json.dumps({
        "youngs_modulus_pa": YOUNGS_MODULUS_PA,
        "tendon_friction_coefficient": TENDON_FRICTION_COEFFICIENT,
        "records": records,
        "rows": rows,
    }, indent=2))

    csv_output = OUT_DIR / "two_force_benchmark.csv"
    with csv_output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(json_output)
    print(csv_output)
    trials = len(ALL_CASES)
    completed_trials = len({row["trial"] for row in rows if row["status"] == "ok"})
    print(f"completed={completed_trials} missing_or_incomplete={trials - completed_trials}")

    print_summary(rows)


if __name__ == "__main__":
    main()
