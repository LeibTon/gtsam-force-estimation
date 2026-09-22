#include <chrono>
#include <map>
#include <cmath>
#include <algorithm>

#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/LevenbergMarquardtParams.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/DoglegOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/linear/GaussianBayesNet.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/nonlinear/NonlinearEquality.h>

#include "factors.h"
#include "types.h"

namespace gtsam{

TendonDiscConfig generate_tendon_disc_config(
    int num_discs,
    int num_poses,
    double routing_radius,
    const std::vector<RoutingAngleFunction>& angle_functions,
    const std::vector<RoutingFunctionParams>& angle_params,
    const std::vector<double>& tendon_routing_radii = {},
    const std::vector<double>& tendon_reach_norm    = {})
{
    int num_tendons = static_cast<int>(angle_functions.size());

    // Per-tendon radii, falling back to global routing_radius if not provided.
    std::vector<double> radii(num_tendons, routing_radius);
    for (int i = 0; i < num_tendons; ++i)
        if (i < static_cast<int>(tendon_routing_radii.size()))
            radii[i] = tendon_routing_radii[i];

    // Per-tendon reach (normalised arc-length), default = full length.
    std::vector<double> reaches(num_tendons, 1.0);
    for (int i = 0; i < num_tendons; ++i)
        if (i < static_cast<int>(tendon_reach_norm.size()))
            reaches[i] = tendon_reach_norm[i];

    TendonDiscConfig config;
    config.num_tendons    = num_tendons;
    config.num_discs      = num_discs;
    config.routing_radius = routing_radius;
    config.disc_pose_idx.reserve(num_discs);
    config.local_holes.reserve(num_discs);

    std::vector<double> pose_s(num_poses);
    std::vector<double> disc_s(num_discs);

    for (int i = 0; i < num_poses; ++i)
        pose_s[i] = static_cast<double>(i) / (num_poses - 1);
    for (int i = 0; i < num_discs; ++i)
        disc_s[i] = static_cast<double>(i) / (num_discs - 1);

    for (int disc_idx = 0; disc_idx < num_discs; ++disc_idx) {
        double s = disc_s[disc_idx];

        int closest_pose_idx = 0;
        double min_dist = std::abs(s - pose_s[0]);
        for (int i = 1; i < num_poses; ++i) {
            double dist = std::abs(s - pose_s[i]);
            if (dist < min_dist) { min_dist = dist; closest_pose_idx = i; }
        }
        config.disc_pose_idx.push_back(closest_pose_idx);

        std::vector<Vector3> holes;
        holes.reserve(num_tendons);
        for (int tendon_idx = 0; tendon_idx < num_tendons; ++tendon_idx) {
            // Tendons that don't reach this disc get a zero placeholder.
            // get_single_tendon_wrench will zero the force via norm > 1e-3 guard.
            if (s > reaches[tendon_idx] + 1e-9) {
                holes.emplace_back(0.0, 0.0, 0.0);
                continue;
            }

            double theta;
            if (angle_functions[tendon_idx] == RoutingAngleFunction::CONSTANT) {
                theta = angle_params[tendon_idx].angle_offset;
            } else if (angle_functions[tendon_idx] == RoutingAngleFunction::LINEAR) {
                theta = angle_params[tendon_idx].angle_offset + s * angle_params[tendon_idx].total_angle;
            } else {
                theta = 0.0;
            }

            double r = radii[tendon_idx];
            holes.emplace_back(r * std::cos(theta), r * std::sin(theta), 0.0);
        }
        config.local_holes.push_back(holes);
    }

    return config;
}

using symbol_shorthand::T; // poses
using symbol_shorthand::D; // disc wrenches
using symbol_shorthand::S; // internal stresses
using symbol_shorthand::Q; // tendon tensions
using symbol_shorthand::G; // Gaussian basis parameters (gamma)

class TendonRobotGtsam {
public:
    TendonRobotGtsam(const TendonRobotConfig& config) {
        num_backbone_poses_ = config.num_discs + (config.num_discs - 1) * config.poses_between_discs;
        ds_ = config.rod_length / (num_backbone_poses_ - 1);
        rod_diameter_ = config.rod_diameter;
        use_midpoint_ = config.use_midpoint;

        // Gravity force per segment, constant in the world/spatial frame.
        gravity_force_per_segment_ = config.gravity_magnitude * config.rod_mass_per_length * ds_
            * config.gravity_direction.normalized();
        disc_gravity_force_ = config.gravity_magnitude * config.disc_mass
            * config.gravity_direction.normalized();

        // Build stiffness matrix
        double cross_section_area = M_PI * std::pow(config.rod_diameter, 2) / 4.0;
        double cross_section_moment = M_PI * std::pow(config.rod_diameter, 4) / 64.0;

        double k_bending_x = config.youngs_modulus * cross_section_moment;
        double k_bending_y = config.youngs_modulus * cross_section_moment;
        double k_torsion = 2.0 * config.shear_modulus * cross_section_moment;
        double k_shear = config.shear_modulus * cross_section_area;
        double k_extension = config.youngs_modulus * cross_section_area;

        // Rod along +z: ω_x/ω_y=bending, ω_z=torsion, v_x/v_y=shear, v_z=extension
        K_inv_ = Matrix6::Zero();
        K_inv_(0, 0) = 1 / k_bending_x; // du[0] = ω_x = bending
        K_inv_(1, 1) = 1 / k_bending_y; // du[1] = ω_y = bending
        K_inv_(2, 2) = 1 / k_torsion;   // du[2] = ω_z = torsion
        K_inv_(3, 3) = 1 / k_shear;     // dv[0] = v_x = shear
        K_inv_(4, 4) = 1 / k_shear;     // dv[1] = v_y = shear
        K_inv_(5, 5) = 1 / k_extension; // dv[2] = v_z = extension

        // Tendon/Disc config
        tendon_config_ = generate_tendon_disc_config(
            config.num_discs, num_backbone_poses_, config.routing_radius,
            config.angle_functions, config.angle_params,
            config.tendon_routing_radii, config.tendon_reach_norm);

        // Noise models
        const int n_t = static_cast<int>(config.angle_functions.size());
        tensions_cov_ = noiseModel::Isotropic::Sigma(n_t, config.tension_meas_std);
        tendon_friction_cov_ = noiseModel::Isotropic::Sigma(n_t, config.tendon_friction_std);
        tendon_friction_coefficient_ = config.tendon_friction_coefficient;

        small_wrench_cov_ = noiseModel::Diagonal::Sigmas((Vector(6) << 
            config.small_moment_std, config.small_moment_std, config.small_moment_std, 
            config.small_force_std, config.small_force_std, config.small_force_std).finished());
        
        base_frame_cov_ = noiseModel::Diagonal::Sigmas((Vector(6) << 
            config.small_r_std, config.small_r_std, config.small_r_std, 
            config.small_p_std, config.small_p_std, config.small_p_std).finished());
        
        cosserat_twist_cov_ = noiseModel::Diagonal::Sigmas((Vector(6) << 
            config.cosserat_twist_r_std, config.cosserat_twist_r_std, config.cosserat_twist_r_std, 
            config.small_p_std, config.small_p_std, config.small_p_std).finished());
        

        num_basis_functions_     = config.num_basis_functions;
        basis_sigma_             = config.basis_sigma;
        gamma_alpha_prior_std_   = config.gamma_alpha_prior_std;
        gamma_beta_prior_std_    = config.gamma_beta_prior_std;
    }

    int num_backbone_poses_;
    bool use_midpoint_;
    double ds_;
    double rod_diameter_;
    double tendon_friction_coefficient_;
    Matrix66 K_inv_;

    Vector3 gravity_force_per_segment_;
    Vector3 disc_gravity_force_;

    // Gaussian basis parameters (used by GaussianBasisForceEstimator)
    int num_basis_functions_;
    double basis_sigma_;
    double gamma_alpha_prior_std_;
    double gamma_beta_prior_std_;

    TendonDiscConfig tendon_config_;

    noiseModel::Diagonal::shared_ptr tensions_cov_;
    noiseModel::Diagonal::shared_ptr tendon_friction_cov_;
    noiseModel::Diagonal::shared_ptr small_wrench_cov_;
    noiseModel::Diagonal::shared_ptr base_frame_cov_;
    noiseModel::Diagonal::shared_ptr cosserat_twist_cov_;


    Ordering ordering_;
    bool is_first_solve_ = true;
    int last_num_iterations_ = 0;
    Values values_;
    NonlinearFactorGraph graph_;
    Marginals marginals_;

    void solve_graph()
    {
        // Reusing the variable ordering can save a few ms.
        if (is_first_solve_) {
            ordering_ = Ordering::Colamd(graph_);
            is_first_solve_ = false;
        }

        DoglegParams params;
        params.setVerbosity("SILENT");
        params.setOrdering(ordering_);
        params.setLinearSolverType("MULTIFRONTAL_QR");
        params.setMaxIterations(500);  // default is 100
        DoglegOptimizer optimizer(graph_, values_, params);

        values_ = optimizer.optimize();
        last_num_iterations_ = optimizer.iterations();
    }

    std::vector<Vector> sample_cov(const Matrix& cov, int num_samples) {
        const int dim = cov.rows();
        Eigen::LLT<Matrix> llt(cov);
        Matrix L = llt.matrixL();

        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::normal_distribution<> normal(0.0, 1.0);  // N(0,1)

        std::vector<Vector> samples;
        samples.reserve(num_samples);

        for (int n = 0; n < num_samples; ++n) {
            Vector z(dim);
            for (int i = 0; i < dim; ++i)
                z(i) = normal(gen);

            Vector delta = L * z;
            samples.push_back(delta);
        }

        return samples;
    }

    void sample_tip_pose(TendonRobotSolution& solution, int num_samples) {
        Pose3 tip_pose_mean = Pose3(solution.backbone_pose_mean.back());
        Matrix6 tip_pose_cov = solution.backbone_pose_cov.back();

        std::vector<Vector> d_tip_pose = sample_cov(tip_pose_cov, num_samples);
        d_tip_pose.reserve(num_samples);

        for (int i = 0; i < num_samples; i++) {
            solution.pose_samples[i] = tip_pose_mean.retract(d_tip_pose[i]).matrix();
        }
    }

    void sample_solution(TendonRobotSolution& solution, int num_samples)
    {
        sample_tip_pose(solution, num_samples);
    }
};


// GaussianBasisForceEstimator
class GaussianBasisForceEstimator : public TendonRobotGtsam {
public:
    GaussianBasisForceEstimator(const TendonRobotConfig& config)
        : TendonRobotGtsam(config)
    {
        angular_strain_meas_cov_ = noiseModel::Isotropic::Sigma(3, config.angular_strain_meas_std);
        tip_position_meas_cov_   = noiseModel::Isotropic::Sigma(3, config.tip_position_meas_std);

        initialize_values();
    }

private:
    noiseModel::Isotropic::shared_ptr angular_strain_meas_cov_;
    noiseModel::Isotropic::shared_ptr tip_position_meas_cov_;

    void initialize_values() {
        values_.clear();

        for (size_t tension_idx = 0; tension_idx < tendon_config_.disc_pose_idx.size(); ++tension_idx) {
            values_.insert(Q(tension_idx), Eigen::VectorXd(Eigen::VectorXd::Zero(tendon_config_.num_tendons)));
        }

        // D(i) for every disc node (needed by TendonDiscWrenchFactor)
        for (size_t disc_idx = 1; disc_idx < tendon_config_.disc_pose_idx.size(); ++disc_idx) {
            values_.insert(D(tendon_config_.disc_pose_idx[disc_idx]), Vector6(Vector6::Zero()));
        }

        for (int i = 0; i < num_backbone_poses_; ++i) {
            values_.insert(T(i), Pose3(Rot3::Identity(), Point3(0.0, 0.0, i * ds_)));
            values_.insert(S(i), Vector6(Vector6::Zero()));
        }

        // beta initialised so mu_i = (L/2)*(1+tanh(beta_i)) lands at the desired s fractions.
        // M=1: s=0.5 (midpoint); M=2: s=0.3,0.7; M=3: s=0.2,0.5,0.8; M=4: s=0.2,0.4,0.6,0.8
        static const std::vector<std::vector<double>> init_fracs = {
            {0.5},
            {0.3, 0.7},
            {0.2, 0.5, 0.8},
            {0.2, 0.4, 0.6, 0.8},
        };

        // Amplitudes: small non-zero to avoid flat Jacobian at initialisation.
        Eigen::VectorXd gamma_init = Eigen::VectorXd::Zero(4 * num_basis_functions_);
        for (int i = 0; i < num_basis_functions_; ++i) {
            gamma_init[4 * i + 0] = 0.001;  // alpha_x
            gamma_init[4 * i + 1] = 0.001;  // alpha_y
            gamma_init[4 * i + 2] = 0.001;  // alpha_z
            int M = num_basis_functions_;
            double frac = (M >= 1 && M <= 4)
                          ? init_fracs[M - 1][i]
                          : 0.05 + 0.9 * i / (M - 1);  // fallback: uniform with 5% inset
            // Invert tanh map: beta = atanh(2*frac - 1)
            gamma_init[4 * i + 3] = std::atanh(2.0 * frac - 1.0);  // beta
        }
        values_.insert(G(0), gamma_init);
    }

    void build_graph_basis(const Eigen::VectorXd& tensions) {
        graph_.resize(0);

        // --- Mechanics factors ---

        // Kinematics factor 
        for (int i = 0; i + 1 < num_backbone_poses_; ++i) {
            graph_.add(CosseratRodTwistFactor(
                T(i), T(i + 1), S(i), S(i + 1), ds_, K_inv_, use_midpoint_, cosserat_twist_cov_));
        }

        // Wrench balance factor :
        //   disc nodes  -> GaussianBasisStressWithDiscFactor with D(i+1) + G(0)
        //   non-disc    -> GaussianBasisStressFactor with G(0)
        const double rod_length = ds_ * (num_backbone_poses_ - 1);
        for (int i = 0; i + 1 < num_backbone_poses_; ++i) {
            bool is_disc_next = std::find(
                tendon_config_.disc_pose_idx.begin(),
                tendon_config_.disc_pose_idx.end(), i + 1)
                != tendon_config_.disc_pose_idx.end();

            if (is_disc_next) {
                double s_k   = i       * ds_;
                double s_kp1 = (i + 1) * ds_;
                graph_.add(GaussianBasisStressWithDiscFactor(
                    T(i), T(i + 1),
                    S(i), S(i + 1),
                    D(i + 1),
                    G(0),
                    s_k, s_kp1, basis_sigma_, rod_length, num_basis_functions_,
                    small_wrench_cov_,
                    gravity_force_per_segment_ + disc_gravity_force_));
            } else {
                double s_k   = i       * ds_;
                double s_kp1 = (i + 1) * ds_;
                graph_.add(GaussianBasisStressFactor(
                    T(i), T(i + 1),
                    S(i), S(i + 1),
                    G(0),
                    s_k, s_kp1, basis_sigma_, rod_length, num_basis_functions_,
                    small_wrench_cov_,
                    gravity_force_per_segment_));
            }
        }

        // Tendon actuation factor : disc wrench balance, plus
        // tendon-tension friction propagation between adjacent spans.
        for (size_t disc_idx = 1; disc_idx < tendon_config_.disc_pose_idx.size(); ++disc_idx) {
            int pose_idx      = tendon_config_.disc_pose_idx[disc_idx];
            int pose_idx_prev = tendon_config_.disc_pose_idx[disc_idx - 1];
            std::vector<Vector3> holes_prev = tendon_config_.local_holes[disc_idx - 1];
            std::vector<Vector3> holes      = tendon_config_.local_holes[disc_idx];

            bool is_tip;
            int pose_idx_next;
            std::vector<Vector3> holes_next;

            if (disc_idx == tendon_config_.disc_pose_idx.size() - 1) {
                is_tip = true;
                pose_idx_next = T(0);
                holes_next    = tendon_config_.local_holes[0];
            } else {
                is_tip = false;
                pose_idx_next = tendon_config_.disc_pose_idx[disc_idx + 1];
                holes_next    = tendon_config_.local_holes[disc_idx + 1];
            }

            if (is_tip) {
                graph_.add(PriorFactor<Eigen::VectorXd>(Q(disc_idx),
                    Eigen::VectorXd(Eigen::VectorXd::Zero(tendon_config_.num_tendons)), tendon_friction_cov_));
            } else {
                if (tendon_friction_coefficient_ > 0.0) {
                    graph_.add(CapstanTendonTensionFactor(
                        T(pose_idx_prev), T(pose_idx), T(pose_idx_next),
                        Q(disc_idx - 1), Q(disc_idx),
                        tendon_friction_coefficient_, holes_prev, holes, holes_next,
                        tendon_friction_cov_));
                } else {
                    graph_.add(TendonTensionEqualityFactor(
                        Q(disc_idx - 1), Q(disc_idx), tendon_friction_cov_));
                }
            }

            graph_.add(TendonDiscWrenchFactor(
                T(pose_idx_prev), T(pose_idx), T(pose_idx_next),
                D(pose_idx), Q(disc_idx - 1), Q(disc_idx),
                is_tip, holes_prev, holes, holes_next,
                small_wrench_cov_));
        }

        // --- Measurement factors ---

        // Tendon tension measurement
        graph_.add(PriorFactor<Eigen::VectorXd>(Q(0), tensions, tensions_cov_));

        // --- Prior and boundary condition factors ---

        // Base boundary condition (Eq. 35) — identity
        graph_.add(PriorFactor<Pose3>(T(0), Pose3(Rot3::Identity(), Point3()), base_frame_cov_));

        // Tip boundary condition (Eq. 37) — free-space, zero tip wrench
        graph_.add(PriorFactor<Vector6>(S(num_backbone_poses_ - 1), Vector6::Zero(), small_wrench_cov_));

        // Force parameter prior (Eq. 39): weak zero-mean prior on all gamma
        // amplitudes (alpha_x, alpha_y, alpha_z per basis). Penalises large
        // forces without constraining their direction.
        {
            Eigen::VectorXd gamma_prior = Eigen::VectorXd::Zero(4 * num_basis_functions_);
            Eigen::VectorXd gamma_sigmas = Eigen::VectorXd::Ones(4 * num_basis_functions_);
            for (int i = 0; i < num_basis_functions_; ++i) {
                gamma_sigmas[4 * i + 0] = gamma_alpha_prior_std_;
                gamma_sigmas[4 * i + 1] = gamma_alpha_prior_std_;
                gamma_sigmas[4 * i + 2] = gamma_alpha_prior_std_;
                gamma_sigmas[4 * i + 3] = gamma_beta_prior_std_;
            }
            graph_.add(PriorFactor<Eigen::VectorXd>(
                G(0), gamma_prior,
                noiseModel::Diagonal::Sigmas(gamma_sigmas)));
        }
    }

    void extract_solution_basis(TendonRobotSolution& solution) {
        marginals_ = Marginals(graph_, values_, Marginals::QR);

        for (int i = 0; i < num_backbone_poses_; ++i) {
            solution.backbone_pose_mean[i] = values_.at<Pose3>(T(i)).matrix();
            solution.backbone_pose_cov[i]  = marginals_.marginalCovariance(T(i));
            // applied_wrench fields are not meaningful in basis mode — zero them out
            if (i > 0) {
                solution.applied_wrench_mean[i - 1] = Vector6::Zero();
                solution.applied_wrench_cov[i - 1]  = Matrix6::Zero();
            }
        }

        solution.gamma_mean = values_.at<Eigen::VectorXd>(G(0));
        solution.gamma_cov  = marginals_.marginalCovariance(G(0));

        solution.tensions_mean = values_.at<Eigen::VectorXd>(Q(0));
        solution.tensions_cov  = marginals_.marginalCovariance(Q(0));

        solution.tendon_disc_config = TendonDiscConfig(tendon_config_);

        KeyVector keys;
        keys.push_back(Q(0));
        keys.push_back(T(num_backbone_poses_ - 1));
        JointMarginal tensions_pose_joint = marginals_.jointMarginalCovariance(keys);

        Eigen::MatrixXd sigma_tensions_tensions = tensions_pose_joint(Q(0), Q(0));
        Eigen::MatrixXd sigma_pose_tensions     = tensions_pose_joint(T(num_backbone_poses_ - 1), Q(0));

        const int n_t2 = sigma_tensions_tensions.rows();
        Eigen::LDLT<Eigen::MatrixXd> ldlt2(sigma_tensions_tensions);
        solution.J_pose_tensions = sigma_pose_tensions * ldlt2.solve(Eigen::MatrixXd::Identity(n_t2, n_t2));
    }

    TendonRobotSolution update_basis(int num_samples) {
        TendonRobotSolution solution = TendonRobotSolution(num_backbone_poses_, num_samples);

        auto t0 = std::chrono::high_resolution_clock::now();
        solve_graph();
        auto t1 = std::chrono::high_resolution_clock::now();
        extract_solution_basis(solution);
        sample_solution(solution, num_samples);
        auto t2 = std::chrono::high_resolution_clock::now();

        solution.solve_time_ms   = std::chrono::duration<double, std::milli>(t1 - t0).count();
        solution.extract_time_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
        solution.total_time_ms   = std::chrono::duration<double, std::milli>(t2 - t0).count();
        solution.num_iterations  = last_num_iterations_;

        return solution;
    }

public:
    TendonRobotSolution step(const Eigen::VectorXd& tensions_meas,
                             const std::map<int, Vector3>& position_meas,
                             const std::map<int, Vector3>& angular_strain_meas,
                             int num_samples)
    {
        build_graph_basis(tensions_meas);

        // Position measurements — one factor per measured pose
        for (const auto& [pose_idx, pos] : position_meas) {
            graph_.add(PositionMeasurementFactor(T(pose_idx), pos, tip_position_meas_cov_));
        }

        // Angular strain measurement factors — only at measured poses
        for (const auto& [pose_idx, du] : angular_strain_meas) {
            graph_.add(AngularStrainFactor(S(pose_idx), du, K_inv_, angular_strain_meas_cov_));
        }

        return update_basis(num_samples);
    }

    // Calibration: forward solve for a known load (gamma_true pinned by a
    // tight prior). Poses only, no Marginals.
    TendonRobotSolution step_forward(const Eigen::VectorXd& tensions_meas,
                                     const Eigen::VectorXd& gamma_true)
    {
        build_graph_basis(tensions_meas);
        graph_.add(PriorFactor<Eigen::VectorXd>(
            G(0), gamma_true,
            noiseModel::Isotropic::Sigma(gamma_true.size(), 1e-9)));
        values_.update(G(0), gamma_true);

        TendonRobotSolution solution = TendonRobotSolution(num_backbone_poses_);
        auto start = std::chrono::high_resolution_clock::now();
        solve_graph();
        auto end = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_backbone_poses_; ++i)
            solution.backbone_pose_mean[i] = values_.at<Pose3>(T(i)).matrix();
        solution.tensions_mean = values_.at<Eigen::VectorXd>(Q(0));
        solution.tendon_disc_config = TendonDiscConfig(tendon_config_);
        solution.solve_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        solution.total_time_ms = solution.solve_time_ms;
        solution.num_iterations = last_num_iterations_;
        return solution;
    }

};
}
