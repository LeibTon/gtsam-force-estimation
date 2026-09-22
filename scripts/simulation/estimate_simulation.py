"""Run the GTSAM Gaussian-basis force estimator on a single SoRoSim-generated sample."""
import numpy as np

from tendon_robot import GaussianBasisForceEstimator


def gamma_bases(gamma: np.ndarray, rod_length: float):
    num_basis = gamma.shape[0] // 4
    bases = []
    for i in range(num_basis):
        alpha = gamma[4 * i: 4 * i + 3]
        beta = gamma[4 * i + 3]
        location = rod_length / 2.0 * (1.0 + np.tanh(beta))
        bases.append((alpha, float(location)))
    return bases


def estimate_simulation(config, disc_positions: np.ndarray):
    """disc_positions: (11, 3) array, model/global frame, metres.

    Independent of how many external forces were actually applied in the
    sample; num_basis_functions on config controls how many the estimator
    fits. Returns solution, node indices, and the estimated shape/bases.
    """
    estimator = GaussianBasisForceEstimator(config)
    nodes = list(estimator.tendon_disc_config.disc_pose_idx)
    position_meas = {nodes[index]: disc_positions[index] for index in range(11)}
    tension = np.zeros(estimator.tendon_disc_config.num_tendons)
    solution = estimator.step(tension, position_meas, {}, 0)

    gamma = np.asarray(solution.gamma_mean)
    bases = gamma_bases(gamma, config.rod_length)
    estimated_shape = np.asarray([solution.backbone_pose_mean[node][:3, 3] for node in nodes])
    full_backbone_shape = np.asarray([pose[:3, 3] for pose in solution.backbone_pose_mean])

    return {
        "solution": solution,
        "nodes": nodes,
        "estimated_shape": estimated_shape,
        "full_backbone_shape": full_backbone_shape,
        "bases": bases,
    }


if __name__ == "__main__":
    import h5py
    import tendon_robot as tr
    from simulation_parameters import POS_MES_STD, TENSION_STD
    from scripts.simulation.simulation_parameters import SHEAR_MODULUS_PA, YOUNGS_MODULUS_PA


    config = tr.TendonRobotConfig()
    config.youngs_modulus = YOUNGS_MODULUS_PA
    config.shear_modulus = SHEAR_MODULUS_PA
    config.gravity_magnitude = 0.0
    config.tendon_friction_coefficient = 0.0
    config.num_basis_functions = 1
    config.basis_sigma = 0.006
    config.tension_meas_std = TENSION_STD
    config.tip_position_meas_std = POS_MES_STD

    dataset_path = "../../simulation_data/dataset_100samples_1forces.h5"
    with h5py.File(dataset_path, "r") as f:
        disc_positions = f["ground_truth/disc/positions"][0].T
        actual_force = f["ground_truth/force/vec"][0, 0]
        actual_location = float(f["ground_truth/force/position"][0, 0]) * config.rod_length

    result = estimate_simulation(config, disc_positions)
    print("bases:", result["bases"])
    print("actual_force:", actual_force, "actual_location:", actual_location)
