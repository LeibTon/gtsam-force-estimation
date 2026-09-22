namespace gtsam {

enum class RoutingAngleFunction {
    CONSTANT = 0,
    LINEAR = 1
};

struct RoutingFunctionParams {
    double angle_offset = 1.490626;  // Starting angle (radians), 85.40654 deg
    double total_angle = 0.0;   // For LINEAR: total angle change across the rod
};

struct TendonDiscConfig {
    int num_tendons;
    int num_discs;
    double routing_radius;
    std::vector<int> disc_pose_idx;
    std::vector<std::vector<Vector3>> local_holes;  // (disc, tendon)
};

struct TendonRobotSolution {
    std::vector<Matrix4> backbone_pose_mean;
    std::vector<Matrix6> backbone_pose_cov;
    std::vector<Matrix4> pose_samples;

    std::vector<Vector6> applied_wrench_mean;
    std::vector<Matrix6> applied_wrench_cov;

    Eigen::VectorXd tensions_mean;
    Eigen::MatrixXd tensions_cov;

    Eigen::MatrixXd J_pose_tensions;

    // Gaussian basis force parameterization (only populated by GaussianBasisForceEstimator)
    Eigen::VectorXd gamma_mean;
    Eigen::MatrixXd gamma_cov;

    double solve_time_ms = 0;
    double extract_time_ms = 0;
    double total_time_ms = 0;
    int num_iterations = 0;

    TendonDiscConfig tendon_disc_config;

    TendonRobotSolution() = default;

    TendonRobotSolution(size_t num_backbone_poses, size_t num_samples = 0) {
        backbone_pose_mean.resize(num_backbone_poses);
        backbone_pose_cov.resize(num_backbone_poses);
        applied_wrench_mean.resize(num_backbone_poses - 1);
        applied_wrench_cov.resize(num_backbone_poses - 1);

        pose_samples.resize(num_samples);
    }
};

struct TendonRobotConfig{
    // Defaults are this repository's experimental robot (see PARAMETERS.md).
    // Rod base is at the origin, extending along the world/local +z axis.

    // --- Rod geometry and discretisation ---
    int num_discs = 11;
    int poses_between_discs = 3;
    double rod_length = 0.270;      // m
    double rod_diameter = 0.001;    // m
    bool use_midpoint = true;

    // --- Material properties ---
    double youngs_modulus = 71.723167e9;  // Pa
    double shear_modulus = 27.585833e9;   // Pa

    // --- Tendon routing and friction ---
    double routing_radius = 0.008;  // m
    double tendon_friction_coefficient = 0.0;  // unitless
    double tendon_friction_std = 1e-5;  // N
    std::vector<RoutingAngleFunction> angle_functions = {RoutingAngleFunction::CONSTANT};
    std::vector<RoutingFunctionParams> angle_params = {RoutingFunctionParams{1.490626, 0.0}};
    std::vector<double> tendon_routing_radii;  // m, per tendon
    std::vector<double> tendon_reach_norm;     // normalised arc-length, per tendon

    // --- Gravity ---
    double gravity_magnitude = 0.0;      // m/s^2
    Vector3 gravity_direction = Vector3(0, 0, 1);  // unit vector, world frame
    double rod_mass_per_length = 0.002 / 0.027;    // kg/m
    double disc_mass = 0.0;              // kg

    // --- Mechanics-constraint noise (soft equality tolerances) ---
    double cosserat_twist_r_std = 1e-2;  // rad
    double small_force_std = 1e-4;       // N
    double small_moment_std = 1e-5;      // N.m
    double small_r_std = 1e-3;           // rad
    double small_p_std = 1e-5;           // m

    // --- Gaussian basis force parameterisation ---
    int num_basis_functions = 1;        // M
    double basis_sigma = 0.006;          // m
    double gamma_alpha_prior_std = 1e-1; // N
    double gamma_beta_prior_std  = 1; // unitless (tanh-mapped location)

    // --- Measurement noise ---
    double tension_meas_std = 0.010;         // N
    double tip_position_meas_std = 0.00272;    // m
    double angular_strain_meas_std = 1e-2;  // rad/m
};
}
