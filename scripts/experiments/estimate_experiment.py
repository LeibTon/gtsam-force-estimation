"""Run the GTSAM Gaussian-basis force estimator on a single experimental trial."""
import json
from pathlib import Path

import numpy as np

from camera_transform import (apply_transform, apply_transform_direction,
                              camera_to_model_transform,
                              force_direction_camera_to_model_rotation,
                              load_arm_points_mm)
from experiment_parameters import EXTERNAL_FORCE_TRANSMISSION_SCALE, GRAVITY_M_S2
from tendon_robot import GaussianBasisForceEstimator

ROPE_KEYS_SINGLE = ("force_direction",)
ROPE_KEYS_MULTI = ("rope1", "rope2")


def read_trial(folder: Path, num_forces: int = 1):
    rope_keys = ROPE_KEYS_SINGLE if num_forces == 1 else ROPE_KEYS_MULTI[:num_forces]
    frames = sorted((folder / "annotation_records").glob("auto_*"))
    required_files = ("arm_points.json", "force_sensor.json", "rope_directions.json")

    shapes, tensions = [], []
    directions = [[] for _ in rope_keys]
    for frame in frames:
        if not all((frame / name).is_file() for name in required_files):
            continue

        force_record = json.loads((frame / "force_sensor.json").read_text())
        if not force_record.get("available", True):
                    continue
        tension = force_record.get("force_n")
        if tension is None:
            continue

        rope = json.loads((frame / "rope_directions.json").read_text())
        frame_directions = []
        for key in rope_keys:
            points = np.asarray(rope[key]["points_3d"], dtype=float)
            if points.shape != (2, 3):
                break
            frame_directions.append(points[1] - points[0])
        if len(frame_directions) != len(rope_keys):
            continue

        shapes.append(load_arm_points_mm(frame / "arm_points.json"))
        tensions.append(float(tension))
        for i, direction in enumerate(frame_directions):
            directions[i].append(direction)

    if not shapes:
        raise ValueError(f"Incomplete trial: {folder}")

    shape_camera_mm = np.mean(shapes, axis=0)
    tension = max(float(np.mean(tensions)), 0.0)
    directions_camera = []
    for rope_directions in directions:
        direction_camera = np.mean(rope_directions, axis=0)
        direction_camera /= np.linalg.norm(direction_camera)
        directions_camera.append(direction_camera)
    return shape_camera_mm, tension, directions_camera, len(shapes)


def known_applied_force(weight_g: float, direction_model: np.ndarray) -> np.ndarray:
    magnitude = EXTERNAL_FORCE_TRANSMISSION_SCALE * weight_g * GRAVITY_M_S2 / 1000.0
    return magnitude * direction_model


def gamma_bases(gamma: np.ndarray, rod_length: float):
    num_basis = gamma.shape[0] // 4
    bases = []
    for i in range(num_basis):
        alpha = gamma[4 * i: 4 * i + 3]
        beta = gamma[4 * i + 3]
        location = rod_length / 2.0 * (1.0 + np.tanh(beta))
        bases.append((alpha, float(location)))
    return bases


def estimate_experiment(config, folder: Path, loads: list[tuple[int, float]]):
    """loads: list of (disc, weight_g), one per known applied force.

    Returns solution, measured_shape (model frame), estimated_shape (model
    frame), bases (list of (force_vec, location_m) from gamma), and per-load
    (actual_force_vec, actual_location_m).
    """
    R, t = camera_to_model_transform()
    R_force = force_direction_camera_to_model_rotation()

    shape_camera_mm, tension, directions_camera, num_frames = read_trial(folder, len(loads))
    measured_shape = apply_transform(shape_camera_mm, R, t)

    actual_loads = []
    for (disc, weight_g), direction_camera in zip(loads, directions_camera):
        direction_model = apply_transform_direction(direction_camera, R_force)
        direction_model /= np.linalg.norm(direction_model)
        actual_force = known_applied_force(weight_g, direction_model)
        actual_location = (disc - 1) / 10.0 * config.rod_length
        actual_loads.append((actual_force, actual_location))

    estimator = GaussianBasisForceEstimator(config)
    nodes = list(estimator.tendon_disc_config.disc_pose_idx)
    position_meas = {nodes[index]: measured_shape[index] for index in range(11)}
    solution = estimator.step(np.asarray([tension]), position_meas, {}, 0)

    gamma = np.asarray(solution.gamma_mean)
    bases = gamma_bases(gamma, config.rod_length)
    estimated_shape = np.asarray([solution.backbone_pose_mean[node][:3, 3] for node in nodes])

    return {
        "solution": solution,
        "nodes": nodes,
        "tension": tension,
        "num_frames": num_frames,
        "measured_shape": measured_shape,
        "estimated_shape": estimated_shape,
        "bases": bases,
        "actual_loads": actual_loads,
    }


if __name__ == "__main__":
    import tendon_robot as tr
    from experiment_parameters import (GRAVITY_DIRECTION, ROD_MASS_PER_LENGTH_KG_M,
                                       ROUTING_ANGLE_OFFSET_RAD, TENDON_FRICTION_COEFFICIENT,
                                       YOUNGS_MODULUS_PA)

    config = tr.TendonRobotConfig()
    config.youngs_modulus = YOUNGS_MODULUS_PA
    config.shear_modulus = YOUNGS_MODULUS_PA / (2.0 * (1.0 + 0.30))
    config.tendon_friction_coefficient = TENDON_FRICTION_COEFFICIENT
    config.gravity_magnitude = GRAVITY_M_S2
    config.gravity_direction = GRAVITY_DIRECTION
    config.rod_mass_per_length = ROD_MASS_PER_LENGTH_KG_M
    config.angle_params = [tr.RoutingFunctionParams(angle_offset=ROUTING_ANGLE_OFFSET_RAD, total_angle=0.0)]
    config.num_basis_functions = 1
    config.basis_sigma = 0.006

    folder = Path(__file__).resolve().parents[2] / "experimental_data/disc6/lact3_weight20_disc6"
    result = estimate_experiment(config, folder, [(6, 20.0)])
    print("tension:", result["tension"])
    print("bases:", result["bases"])
    print("actual_loads:", result["actual_loads"])
