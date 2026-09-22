"""Run estimate_simulation.py over every sample in a SoRoSim single-force dataset,
writing an output .h5 with the estimated results alongside the ground truth."""
from pathlib import Path

import h5py
import numpy as np

import tendon_robot as tr
from estimate_simulation import estimate_simulation
from simulation_parameters import POSITION_NOISE_STD_M, GRAVITY_M_S2, POS_MES_STD, SHEAR_MODULUS_PA, TENDON_FRICTION_COEFFICIENT, TENSION_STD, YOUNGS_MODULUS_PA

ROOT = Path(__file__).resolve().parents[2]
IN_PATH = ROOT / "simulation_data" / "dataset_1forces.h5"

NOISY_RESULT = False

OUT_NAME = "single_force_batch_results_noise.h5" if NOISY_RESULT else "single_force_batch_results.h5"
OUT_PATH = Path(__file__).resolve().parent / OUT_NAME


def build_config():
    config = tr.TendonRobotConfig()
    config.youngs_modulus = YOUNGS_MODULUS_PA
    config.shear_modulus = SHEAR_MODULUS_PA
    config.gravity_magnitude = GRAVITY_M_S2
    config.tendon_friction_coefficient = TENDON_FRICTION_COEFFICIENT
    config.num_basis_functions = 1
    config.basis_sigma = 0.006
    config.tension_meas_std = TENSION_STD
    config.tip_position_meas_std = POSITION_NOISE_STD_M if NOISY_RESULT else POS_MES_STD
    return config


def transverse_direction_error_deg(rotation, estimated_force, actual_force) -> float:
    estimated_xy = (rotation.T @ estimated_force)[:2]
    actual_xy = (rotation.T @ actual_force)[:2]
    if np.linalg.norm(estimated_xy) <= 1e-12 or np.linalg.norm(actual_xy) <= 1e-12:
        return float("nan")

    angle_difference = (np.arctan2(estimated_xy[1], estimated_xy[0])
                        - np.arctan2(actual_xy[1], actual_xy[0]))
    return abs(float(np.degrees(np.arctan2(np.sin(angle_difference), np.cos(angle_difference)))))


def nearest_node_rotation(backbone_pose_mean, nodes, rod_length, location_m):
    disc_arc_lengths = np.linspace(0.0, rod_length, len(nodes))
    disc_index = int(np.argmin(np.abs(disc_arc_lengths - location_m)))
    return np.asarray(backbone_pose_mean[nodes[disc_index]])[:3, :3]


def location_variance_m2(gamma_cov: np.ndarray, beta: float, rod_length: float) -> float:
    """Propagate beta's variance through location = (L/2)(1 + tanh(beta))."""
    jacobian = rod_length / 2.0 * (1.0 - np.tanh(beta) ** 2)
    beta_var = gamma_cov[3, 3]
    return float(jacobian ** 2 * beta_var)


def full_backbone_ground_truth(dense_shape_transforms, num_nodes: int = 41):
    """Interpolate SoRoSim's dense 101-point shape (evenly spaced in s) onto
    GTSAM's 41 backbone node positions (k/40 fractions of the rod length)."""
    dense_positions = dense_shape_transforms.reshape(101, 4, 4)[:, 3, :3]  # translation row
    dense_s = np.linspace(0.0, 1.0, 101)
    node_s = np.linspace(0.0, 1.0, num_nodes)
    return np.column_stack([np.interp(node_s, dense_s, dense_positions[:, axis]) for axis in range(3)])


def run_sample(config, disc_positions):
    result = estimate_simulation(config, disc_positions)
    estimated_force, estimated_location = result["bases"][0]
    solution = result["solution"]
    nodes = result["nodes"]

    gamma_mean = np.asarray(solution.gamma_mean)
    gamma_cov = np.asarray(solution.gamma_cov)

    return {
        "estimated_force": estimated_force,
        "estimated_location": estimated_location,
        "estimated_location_var": location_variance_m2(gamma_cov, gamma_mean[3], config.rod_length),
        "estimated_shape": result["estimated_shape"],
        "full_backbone_shape": result["full_backbone_shape"],
        "gamma_mean": gamma_mean,
        "gamma_cov": gamma_cov,
        "shape_cov": np.asarray([solution.backbone_pose_cov[node] for node in nodes]),
        "backbone_pose_mean": solution.backbone_pose_mean,
        "nodes": nodes,
        "solve_time_ms": solution.solve_time_ms,
        "total_time_ms": solution.total_time_ms,
        "num_iterations": solution.num_iterations,
    }


def main():
    config = build_config()

    with h5py.File(IN_PATH, "r") as f:
        n_samples = int(f["metadata/n_samples"][0, 0])
        disc_positions_all = f["ground_truth/disc/positions"][:]
        force_vec_all = f["ground_truth/force/vec"][:, 0, :]
        force_position_all = f["ground_truth/force/position"][:, 0]
        dense_shape_all = f["ground_truth/shape/transforms"][:]

    n_nodes = disc_positions_all.shape[2]
    num_backbone_nodes = 41
    estimated_force = np.zeros((n_samples, 3))
    estimated_location = np.zeros(n_samples)
    estimated_location_var = np.zeros(n_samples)
    estimated_shape = np.zeros((n_samples, n_nodes, 3))
    full_estimated_shape = np.zeros((n_samples, num_backbone_nodes, 3))
    full_ground_truth_shape = np.zeros((n_samples, num_backbone_nodes, 3))
    shape_cov = np.zeros((n_samples, n_nodes, 6, 6))
    gamma_mean = np.zeros((n_samples, 4 * config.num_basis_functions))
    gamma_cov = np.zeros((n_samples, 4 * config.num_basis_functions, 4 * config.num_basis_functions))

    force_error_n = np.zeros(n_samples)
    force_location_error_mm = np.zeros(n_samples)
    force_direction_error_deg = np.zeros(n_samples)
    solve_time_ms = np.zeros(n_samples)
    total_time_ms = np.zeros(n_samples)
    num_iterations = np.zeros(n_samples, dtype=np.int64)

    rng = np.random.default_rng(0)
    for i in range(n_samples):
        print(f"[{i + 1}/{n_samples}]")
        disc_positions = disc_positions_all[i].T
        if NOISY_RESULT:
            disc_positions = disc_positions + rng.normal(0.0, POSITION_NOISE_STD_M, disc_positions.shape)
        actual_force = force_vec_all[i]
        actual_location = float(force_position_all[i]) * config.rod_length

        r = run_sample(config, disc_positions)
        estimated_force[i] = r["estimated_force"]
        estimated_location[i] = r["estimated_location"]
        estimated_location_var[i] = r["estimated_location_var"]
        estimated_shape[i] = r["estimated_shape"]
        full_estimated_shape[i] = r["full_backbone_shape"]
        full_ground_truth_shape[i] = full_backbone_ground_truth(dense_shape_all[i], num_backbone_nodes)
        shape_cov[i] = r["shape_cov"]
        gamma_mean[i] = r["gamma_mean"]
        gamma_cov[i] = r["gamma_cov"]
        solve_time_ms[i] = r["solve_time_ms"]
        total_time_ms[i] = r["total_time_ms"]
        num_iterations[i] = r["num_iterations"]

        rotation = nearest_node_rotation(r["backbone_pose_mean"], r["nodes"], config.rod_length, actual_location)
        force_error_n[i] = abs(np.linalg.norm(r["estimated_force"]) - np.linalg.norm(actual_force))
        force_location_error_mm[i] = abs(r["estimated_location"] - actual_location) * 1e3
        force_direction_error_deg[i] = transverse_direction_error_deg(rotation, r["estimated_force"], actual_force)

    with h5py.File(OUT_PATH, "w") as f:
        f.create_dataset("ground_truth/disc/positions", data=disc_positions_all)
        f.create_dataset("ground_truth/force/vec", data=force_vec_all)
        f.create_dataset("ground_truth/force/position", data=force_position_all)

        f.create_dataset("estimated/force/vec", data=estimated_force)
        f.create_dataset("estimated/force/location_m", data=estimated_location)
        f.create_dataset("estimated/force/location_var_m2", data=estimated_location_var)
        f.create_dataset("estimated/shape", data=estimated_shape)
        f.create_dataset("estimated/full_backbone_shape", data=full_estimated_shape)
        f.create_dataset("ground_truth/full_backbone_shape", data=full_ground_truth_shape)
        f.create_dataset("estimated/shape_cov", data=shape_cov)
        f.create_dataset("estimated/gamma_mean", data=gamma_mean)
        f.create_dataset("estimated/gamma_cov", data=gamma_cov)
        f.create_dataset("estimated/solve_time_ms", data=solve_time_ms)
        f.create_dataset("estimated/total_time_ms", data=total_time_ms)
        f.create_dataset("estimated/num_iterations", data=num_iterations)

        f.attrs["youngs_modulus_pa"] = config.youngs_modulus
        f.attrs["shear_modulus_pa"] = config.shear_modulus
        f.attrs["rod_length_m"] = config.rod_length
        f.attrs["basis_sigma_m"] = config.basis_sigma
        f.attrs["n_samples"] = n_samples

    print(OUT_PATH)
    print(f"force_error_n           mean={force_error_n.mean():.4f} N   median={np.median(force_error_n):.4f} N")
    print(f"force_location_error_mm mean={force_location_error_mm.mean():.4f} mm  median={np.median(force_location_error_mm):.4f} mm")
    print(f"force_direction_error_deg mean={np.nanmean(force_direction_error_deg):.4f} deg median={np.nanmedian(force_direction_error_deg):.4f} deg")
    print(f"solve_time_ms           mean={solve_time_ms.mean():.4f} ms   median={np.median(solve_time_ms):.4f} ms")
    print(f"num_iterations          mean={num_iterations.mean():.4f}     median={np.median(num_iterations):.4f}")


if __name__ == "__main__":
    main()
