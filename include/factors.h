
namespace gtsam{

// ============================================================================
// Helper functions
// ============================================================================

// Computes Ad^T(T_k^{-1}T_{k+1})*Lambda_k, i.e. propagates a body-frame wrench from pose_k to pose_kp1.
Vector6 propagate_wrench_backward(
    const Pose3& pose,
    const Pose3& tip_pose,
    const Vector6& tip_wrench,
    OptionalJacobian<6, 6> H_pose = {},
    OptionalJacobian<6, 6> H_tip_pose = {},
    OptionalJacobian<6, 6> H_tip_wrench = {})
{
    Matrix66 d_tip_pose_inv_d_tip_pose;
    Matrix66 d_delta_d_tip_pose_inv, d_delta_d_pose;
    Matrix66 d_wrench_d_delta, d_wrench_d_tip_wrench_;

    Pose3 tip_pose_inv = tip_pose.inverse(
        (H_pose ? &d_tip_pose_inv_d_tip_pose : 0));

    Pose3 delta = tip_pose_inv.compose(pose,
        (H_tip_pose ? &d_delta_d_tip_pose_inv : 0),
        (H_pose ? &d_delta_d_pose : 0));

    Vector6 wrench = delta.AdjointTranspose(
        tip_wrench,
        (H_pose || H_tip_pose ? &d_wrench_d_delta : 0),
        (H_tip_wrench ? &d_wrench_d_tip_wrench_ : 0));

    // Assign Jacobians if needed
    if (H_pose) {
        *H_pose = d_wrench_d_delta * d_delta_d_pose;
    }
    if (H_tip_pose) {
        *H_tip_pose = d_wrench_d_delta * d_delta_d_tip_pose_inv * d_tip_pose_inv_d_tip_pose;
    }
    if (H_tip_wrench) {
        *H_tip_wrench = d_wrench_d_tip_wrench_;
    }

    return wrench;
}

// Computes Rot^T(T_{k+1})*F, i.e. rotates a spatial/world-frame wrench into pose's body frame.
Vector6 spatial_to_body_wrench(
    const Vector6& wrench_spatial, const Pose3& pose,
    OptionalJacobian<6, 6> H_wrench = {}, OptionalJacobian<6, 6> H_pose = {})
{
    Matrix3 d_moment_d_rotation, d_force_d_rotation, d_moment_d_moment, d_force_d_force;
    Matrix36 d_rotation_d_pose;

    Vector6 wrench_body;

    Rot3 rot = pose.rotation(d_rotation_d_pose);

    wrench_body.head<3>() = rot.unrotate(wrench_spatial.head<3>(),
        H_pose ? &d_moment_d_rotation : 0,
        H_wrench ? &d_moment_d_moment : 0);

    wrench_body.tail<3>() = rot.unrotate(wrench_spatial.tail<3>(),
        H_pose ? &d_force_d_rotation : 0,
        H_wrench ? &d_force_d_force : 0);

    if (H_pose) {
        H_pose->setZero();
        H_pose->block<3,3>(0,0) = d_moment_d_rotation;
        H_pose->block<3,3>(3,0) = d_force_d_rotation;
    }

    if (H_wrench) {
        H_wrench->setZero();
        H_wrench->block<3,3>(0,0) = d_moment_d_moment;
        H_wrench->block<3,3>(3,3) = d_force_d_force;
    }

    return wrench_body;
}

// One tendon segment's wrench on a disc, pulling toward the adjacent hole (paper's f_d,i^{-/+}).
Vector6 get_single_tendon_wrench(
    const double tension,
    const Pose3 pose,
    const Pose3 pose_other,
    const Point3 hole,
    const Point3 hole_other,
    OptionalJacobian<6, 1> H_tension = {},
    OptionalJacobian<6, 6> H_pose = {},
    OptionalJacobian<6, 6> H_pose_other = {})
{
    // Compute the position of the other hole in the world frame
    Matrix36 d_hole_other_world_d_pose_other;
    Point3 hole_other_world = pose_other.transformFrom(hole_other,
        H_pose_other ? &d_hole_other_world_d_pose_other : 0);

    Matrix36 d_hole_other_local_d_pose;
    Matrix3 d_hole_other_local_d_hole_other_world;
    Point3 hole_other_local = pose.transformTo(hole_other_world,
        H_pose? &d_hole_other_local_d_pose : 0,
        d_hole_other_local_d_hole_other_world);

    Vector3 hole_diff = hole_other_local - hole;
    double norm = hole_diff.norm();

    Vector3 force_dir;
    Matrix3 d_force_dir_d_hole_diff = Matrix3::Zero();

    bool valid = hole_diff.allFinite() && norm > 1e-3;

    if (valid) {
        force_dir = normalize(hole_diff, H_pose || H_pose_other ? &d_force_dir_d_hole_diff : 0);
    } else {
        force_dir = Vector3::Zero();
    }

    Vector3 force = tension * force_dir;
    Matrix31 d_force_d_tension = force_dir;
    Matrix33 d_force_d_force_dir = tension * Matrix3::Identity();

    Matrix33 d_moment_d_force;
    Vector3 moment = cross(hole, force, nullptr,
         H_tension || H_pose || H_pose_other ? &d_moment_d_force : 0);

    Vector6 wrench;
    wrench << moment, force;

    if (H_tension) {
        H_tension->head<3>() = d_moment_d_force * d_force_d_tension;
        H_tension->tail<3>() = d_force_d_tension;
    }

    if (H_pose) {
        Matrix36 d_force_dir_d_pose = d_force_dir_d_hole_diff * d_hole_other_local_d_pose;
        Matrix36 d_force_d_pose = d_force_d_force_dir * d_force_dir_d_pose;
        Matrix36 d_moment_d_pose = d_moment_d_force * d_force_d_pose;

        H_pose->block<3,6>(0,0) = d_moment_d_pose;
        H_pose->block<3,6>(3,0) = d_force_d_pose;
    }

    if (H_pose_other) {
        Matrix36 d_force_dir_d_pose_other =
            d_force_dir_d_hole_diff *
            d_hole_other_local_d_hole_other_world *
            d_hole_other_world_d_pose_other;

        Matrix36 d_force_d_pose_other = d_force_d_force_dir * d_force_dir_d_pose_other;
        Matrix36 d_moment_d_pose_other = d_moment_d_force * d_force_d_pose_other;

        H_pose_other->block<3,6>(0,0) = d_moment_d_pose_other;
        H_pose_other->block<3,6>(3,0) = d_force_d_pose_other;
    }

    return wrench;
}

// Integrated Gaussian-basis force over [s_k, s_kp1] (paper Eq. 14-15) and its
// Jacobian w.r.t. gamma (3 x 4M). gamma layout per basis i (block of 4):
// [ax, ay, az, beta], with mu_i = (L/2)*(1 + tanh(beta_i)) in (0, L).
// Shared by GaussianBasisStressFactor and GaussianBasisStressWithDiscFactor.
Vector3 compute_integrated_gaussian_basis_force(
    const Eigen::VectorXd& gamma,
    double s_k, double s_kp1, double sigma, double rod_length, int M,
    Eigen::Matrix<double, 3, Eigen::Dynamic>* J_gamma = nullptr)
{
    const double sqrt2   = std::sqrt(2.0);
    const double sqrt2pi = std::sqrt(2.0 * M_PI);

    Vector3 f_e = Vector3::Zero();

    if (J_gamma) {
        J_gamma->setZero(3, 4 * M);
    }

    for (int i = 0; i < M; ++i) {
        int base = 4 * i;
        double ax   = gamma[base + 0];
        double ay   = gamma[base + 1];
        double az   = gamma[base + 2];
        double beta = gamma[base + 3];

        // tanh reparameterisation: mu = (L/2) * (1 + tanh(beta))  =>  mu in (0, L)
        // At beta=0: mu = L/2 (midpoint). d_mu/d_beta = (L/2) * sech^2(beta)
        double tanh_beta   = std::tanh(beta);
        double mu          = (rod_length / 2.0) * (1.0 + tanh_beta);
        double d_mu_d_beta = (rod_length / 2.0) * (1.0 - tanh_beta * tanh_beta);

        double z_kp1    = (s_kp1 - mu) / (sqrt2 * sigma);
        double z_k      = (s_k   - mu) / (sqrt2 * sigma);
        double erf_diff = 0.5 * (std::erf(z_kp1) - std::erf(z_k));

        f_e[0] += ax * erf_diff;
        f_e[1] += ay * erf_diff;
        f_e[2] += az * erf_diff;

        if (J_gamma) {
            // df/d(alpha): erf_diff * I (linear in amplitudes)
            (*J_gamma)(0, base + 0) = erf_diff;
            (*J_gamma)(1, base + 1) = erf_diff;
            (*J_gamma)(2, base + 2) = erf_diff;

            // df/d(beta) = df/d(mu) * d_mu/d_beta
            // d/dmu [0.5*(erf(z_{k+1}) - erf(z_k))]
            //   = (1/(sqrt(2pi)*sigma)) * (exp(-z_k^2) - exp(-z_{k+1}^2))
            double d_erf_diff_d_mu   = (std::exp(-z_k * z_k) - std::exp(-z_kp1 * z_kp1))
                                       / (sqrt2pi * sigma);
            double d_erf_diff_d_beta = d_erf_diff_d_mu * d_mu_d_beta;
            (*J_gamma)(0, base + 3) = ax * d_erf_diff_d_beta;
            (*J_gamma)(1, base + 3) = ay * d_erf_diff_d_beta;
            (*J_gamma)(2, base + 3) = az * d_erf_diff_d_beta;
        }
    }
    return f_e;
}

// ============================================================================
// Mechanics factors
// ============================================================================

// Kinematics factor (paper Fig. 2: yellow, combined kinematics/constitutive
// factor): enforces strain-pose consistency between the discrete kinematic
// equation and the constitutive midpoint approximation.
class CosseratRodTwistFactor: public NoiseModelFactorN<Pose3, Pose3, Vector6, Vector6> {
    double ds_;  // segment length
    gtsam::Matrix66 K_inv_;  // Assuming constant stiffness inverse per factor
    bool use_midpoint_;

public:
    using Base = NoiseModelFactorN<Pose3, Pose3, Vector6, Vector6>;
    using Base::evaluateError;

    CosseratRodTwistFactor(Key pose_0_key,
                           Key pose_1_key,
                           Key stress_0_key,
                           Key stress_1_key,
                           double ds,
                           const Matrix66& K_inv,
                           bool use_midpoint,
                           const SharedNoiseModel& model):
        Base(model, pose_0_key, pose_1_key, stress_0_key, stress_1_key),
        ds_(ds), K_inv_(K_inv), use_midpoint_(use_midpoint) {}

    Vector evaluateError(
        const Pose3& pose_0,
        const Pose3& pose_1,
        const Vector6& stress_0,
        const Vector6& stress_1,
        OptionalMatrixType H1,
        OptionalMatrixType H2,
        OptionalMatrixType H3,
        OptionalMatrixType H4) const override {

        Matrix66 d_delta_d_pose_0, d_delta_d_pose_1;
        Pose3 delta = pose_0.between(pose_1,
            H1 ? &d_delta_d_pose_0 : 0,
            H2 ? &d_delta_d_pose_1 : 0);

        Matrix66 d_twist_d_delta;
        Vector6 twist = Pose3::Logmap(delta,
            H1 || H2 ? &d_twist_d_delta : 0);

        Vector6 stress_mid = 0.5 * (stress_0 + stress_1);
        Vector6 stress = use_midpoint_ ? stress_mid : stress_0;

        Vector6 nominal_strain = Vector6::Zero();
        nominal_strain[5] = 1.0;  // Straight rod: linear velocity in z direction only
        Vector6 twist_p = ds_ * (K_inv_ * stress + nominal_strain);

        Vector6 twist_error = twist_p - twist;

        if (H1) {
            *H1 = -d_twist_d_delta * d_delta_d_pose_0;
        }

        if (H2) {
            *H2 = -d_twist_d_delta * d_delta_d_pose_1;
        }

        if (H3) {
            *H3 = ds_ * K_inv_;
            if (use_midpoint_) {
                *H3 *= 0.5;
            }
        }

        if (H4) {
            if (use_midpoint_) {
                *H4 = 0.5 * ds_ * K_inv_;
            } else {
                *H4 = Matrix6::Zero();
            }
        }

        return twist_error;
    }
};

// Wrench balance factor at non-disc nodes (paper Eq. 25; Fig. 2: blue),
// driven by a shared Gaussian basis gamma in R^{4M} (amplitudes alpha_i,
// tanh-bounded centers mu_i via beta_i, erf-integrated over [s_k, s_{k+1}]).
class GaussianBasisStressFactor: public NoiseModelFactorN<Pose3, Pose3, Vector6, Vector6, Eigen::VectorXd> {
    double s_k_;         // arc-length at node k (metres)
    double s_kp1_;       // arc-length at node k+1
    double sigma_;       // fixed Gaussian width (metres)
    double rod_length_;  // L: total rod length, used for sigmoid reparameterisation
    int M_;              // number of basis functions
    Vector3 gravity_force_spatial_;  // constant per-segment gravity force, world frame (zero = disabled)

public:
    using Base = NoiseModelFactorN<Pose3, Pose3, Vector6, Vector6, Eigen::VectorXd>;
    using Base::evaluateError;

    GaussianBasisStressFactor(Key pose_k_key,
                              Key pose_kp1_key,
                              Key stress_k_key,
                              Key stress_kp1_key,
                              Key gamma_key,
                              double s_k,
                              double s_kp1,
                              double sigma,
                              double rod_length,
                              int M,
                              const SharedNoiseModel& model,
                              const Vector3& gravity_force_spatial = Vector3::Zero())
        : Base(model, pose_k_key, pose_kp1_key, stress_k_key, stress_kp1_key, gamma_key),
          s_k_(s_k), s_kp1_(s_kp1), sigma_(sigma), rod_length_(rod_length), M_(M),
          gravity_force_spatial_(gravity_force_spatial) {}

    Vector evaluateError(
        const Pose3& pose_k,
        const Pose3& pose_kp1,
        const Vector6& stress_k,
        const Vector6& stress_kp1,
        const Eigen::VectorXd& gamma,
        OptionalMatrixType H1,
        OptionalMatrixType H2,
        OptionalMatrixType H3,
        OptionalMatrixType H4,
        OptionalMatrixType H5) const override
    {
        // 1. Compute integrated external force and its Jacobian w.r.t. gamma
        Eigen::Matrix<double, 3, Eigen::Dynamic> d_fe_d_gamma;
        Vector3 f_e = compute_integrated_gaussian_basis_force(
            gamma, s_k_, s_kp1_, sigma_, rod_length_, M_, H5 ? &d_fe_d_gamma : nullptr);

        // 2. Build the 6D spatial wrench (moment = 0, force = f_e + gravity)
        Vector6 wrench_spatial = Vector6::Zero();
        wrench_spatial.tail<3>() = f_e + gravity_force_spatial_;

        // Paper: Lambda_{k+1} = Ad^T(T_k^{-1}T_{k+1})*Lambda_k - Rot^T(T_{k+1})*[0;v_e]
        // Residual: stress_kp1 - Ad^T(T_k^{-1}T_{k+1})*stress_k + Rot^T(T_{k+1})*[0;v_e] = 0
        // propagate_wrench_backward(pose_kp1, pose_k, stress_k) = Ad^T(T_k^{-1}T_{k+1})*stress_k
        Matrix6 d_wrench_body_d_pose_kp1, d_wrench_body_d_wrench_spatial;
        Vector6 wrench_body = spatial_to_body_wrench(
            wrench_spatial, pose_kp1,
            H2 || H5 ? &d_wrench_body_d_wrench_spatial : nullptr,
            H2       ? &d_wrench_body_d_pose_kp1       : nullptr);

        Matrix6 d_stress_p_d_pose_kp1, d_stress_p_d_pose_k, d_stress_p_d_stress_k;
        Vector6 stress_propagated = propagate_wrench_backward(
            pose_kp1, pose_k, stress_k,
            H2 ? &d_stress_p_d_pose_kp1 : nullptr,
            H1 ? &d_stress_p_d_pose_k   : nullptr,
            H3 ? &d_stress_p_d_stress_k : nullptr);

        Vector6 error = stress_kp1 - stress_propagated + wrench_body;

        if (H1) {
            *H1 = -d_stress_p_d_pose_k;
        }

        if (H2) {
            *H2 = d_wrench_body_d_pose_kp1 - d_stress_p_d_pose_kp1;
        }

        if (H3) {
            *H3 = -d_stress_p_d_stress_k;
        }

        if (H4) {
            *H4 = Matrix6::Identity();
        }

        if (H5) {
            Matrix63 d_wrench_spatial_d_fe;
            d_wrench_spatial_d_fe.topRows<3>().setZero();
            d_wrench_spatial_d_fe.bottomRows<3>().setIdentity();

            Eigen::MatrixXd d_wrench_body_d_fe =
                d_wrench_body_d_wrench_spatial * d_wrench_spatial_d_fe;  // 6 x 3

            *H5 = d_wrench_body_d_fe * d_fe_d_gamma;  // 6 x 4M
        }

        return error;
    }
};


// Like GaussianBasisStressFactor, but also adds a disc wrench D for segments
// where tendon actuation acts alongside the distributed gamma load (paper
// Fig. 2: cyan, wrench balance at disc/routing nodes).
class GaussianBasisStressWithDiscFactor
    : public NoiseModelFactorN<Pose3, Pose3, Vector6, Vector6, Vector6, Eigen::VectorXd> {
    double s_k_;         // arc-length at node k (metres)
    double s_kp1_;       // arc-length at node k+1
    double sigma_;       // fixed Gaussian width (metres)
    double rod_length_;  // L: total rod length, used for sigmoid reparameterisation
    int M_;              // number of basis functions
    Vector3 gravity_force_spatial_;  // constant per-segment gravity force, world frame (zero = disabled)

public:
    using Base = NoiseModelFactorN<Pose3, Pose3, Vector6, Vector6, Vector6, Eigen::VectorXd>;
    using Base::evaluateError;

    GaussianBasisStressWithDiscFactor(Key pose_k_key,
                                      Key pose_kp1_key,
                                      Key stress_k_key,
                                      Key stress_kp1_key,
                                      Key disc_wrench_key,
                                      Key gamma_key,
                                      double s_k,
                                      double s_kp1,
                                      double sigma,
                                      double rod_length,
                                      int M,
                                      const SharedNoiseModel& model,
                                      const Vector3& gravity_force_spatial = Vector3::Zero())
        : Base(model, pose_k_key, pose_kp1_key, stress_k_key, stress_kp1_key, disc_wrench_key, gamma_key),
          s_k_(s_k), s_kp1_(s_kp1), sigma_(sigma), rod_length_(rod_length), M_(M),
          gravity_force_spatial_(gravity_force_spatial) {}

    Vector evaluateError(
        const Pose3& pose_k,
        const Pose3& pose_kp1,
        const Vector6& stress_k,
        const Vector6& stress_kp1,
        const Vector6& disc_wrench,
        const Eigen::VectorXd& gamma,
        OptionalMatrixType H1,
        OptionalMatrixType H2,
        OptionalMatrixType H3,
        OptionalMatrixType H4,
        OptionalMatrixType H5,
        OptionalMatrixType H6) const override
    {
        Eigen::Matrix<double, 3, Eigen::Dynamic> d_fe_d_gamma;
        Vector3 f_e = compute_integrated_gaussian_basis_force(
            gamma, s_k_, s_kp1_, sigma_, rod_length_, M_, H6 ? &d_fe_d_gamma : nullptr);

        Vector6 wrench_spatial = Vector6::Zero();
        wrench_spatial.tail<3>() = f_e + gravity_force_spatial_;

        // Paper: Lambda_{k+1} = Ad^T(T_k^{-1}T_{k+1})*Lambda_k - Rot^T(T_{k+1})*[0;v_e] - F_ten_{k+1}
        // Residual: stress_kp1 - Ad^T(T_k^{-1}T_{k+1})*stress_k + Rot^T(T_{k+1})*[0;v_e] + disc_wrench = 0
        // propagate_wrench_backward(pose_kp1, pose_k, stress_k) = Ad^T(T_k^{-1}T_{k+1})*stress_k
        Matrix6 d_wrench_body_d_pose_kp1, d_wrench_body_d_wrench_spatial;
        Vector6 wrench_body = spatial_to_body_wrench(
            wrench_spatial, pose_kp1,
            H2 || H6 ? &d_wrench_body_d_wrench_spatial : nullptr,
            H2       ? &d_wrench_body_d_pose_kp1       : nullptr);

        Matrix6 d_stress_p_d_pose_kp1, d_stress_p_d_pose_k, d_stress_p_d_stress_k;
        Vector6 stress_propagated = propagate_wrench_backward(
            pose_kp1, pose_k, stress_k,
            H2 ? &d_stress_p_d_pose_kp1 : nullptr,
            H1 ? &d_stress_p_d_pose_k   : nullptr,
            H3 ? &d_stress_p_d_stress_k : nullptr);

        Vector6 error = stress_kp1 - stress_propagated + wrench_body + disc_wrench;

        if (H1) {
            *H1 = -d_stress_p_d_pose_k;
        }

        if (H2) {
            *H2 = d_wrench_body_d_pose_kp1 - d_stress_p_d_pose_kp1;
        }

        if (H3) {
            *H3 = -d_stress_p_d_stress_k;
        }

        if (H4) {
            *H4 = Matrix6::Identity();
        }

        if (H5) {
            *H5 = Matrix6::Identity();
        }

        if (H6) {
            Matrix63 d_wrench_spatial_d_fe;
            d_wrench_spatial_d_fe.topRows<3>().setZero();
            d_wrench_spatial_d_fe.bottomRows<3>().setIdentity();

            Eigen::MatrixXd d_wrench_body_d_fe =
                d_wrench_body_d_wrench_spatial * d_wrench_spatial_d_fe;

            *H6 = d_wrench_body_d_fe * d_fe_d_gamma;
        }

        return error;
    }
};

// Tendon actuation factor (paper Eq. 27; paper Fig. 2: violet): disc wrench
// balance, D equal to the sum of tendon-guide forces from the adjacent
// spans.
class TendonDiscWrenchFactor: public NoiseModelFactorN<Pose3, Pose3, Pose3, Vector6, Eigen::VectorXd, Eigen::VectorXd> {
    bool is_tip_;
    std::vector<Point3> holes_prev_;
    std::vector<Point3> holes_;
    std::vector<Point3> holes_next_;
public:
    using Base = NoiseModelFactorN<Pose3, Pose3, Pose3, Vector6, Eigen::VectorXd, Eigen::VectorXd>;
    using Base::evaluateError;

    TendonDiscWrenchFactor(Key pose_prev_key,
                           Key pose_key,
                           Key pose_next_key,
                           Key wrench_key,
                           Key tension_in_key,
                           Key tension_out_key,
                           const bool is_tip,
                           const std::vector<Point3>& holes_prev,
                           const std::vector<Point3>& holes,
                           const std::vector<Point3>& holes_next,
                           const SharedNoiseModel& model):
        Base(model, pose_prev_key, pose_key, pose_next_key, wrench_key, tension_in_key, tension_out_key),
        is_tip_(is_tip), holes_prev_(holes_prev), holes_(holes), holes_next_(holes_next) {}

    Vector evaluateError(
        const Pose3& pose_prev,
        const Pose3& pose,
        const Pose3& pose_next,
        const Vector6& wrench,
        const Eigen::VectorXd& tension_in,
        const Eigen::VectorXd& tension_out,
        OptionalMatrixType H1,
        OptionalMatrixType H2,
        OptionalMatrixType H3,
        OptionalMatrixType H4,
        OptionalMatrixType H5,
        OptionalMatrixType H6) const override
    {
        const int n_tendons = static_cast<int>(tension_in.size());
        Vector6 wrench_tendons = Vector6::Zero();
        Eigen::Matrix<double, 6, Eigen::Dynamic> d_wrench_d_tension_in =
            Eigen::Matrix<double, 6, Eigen::Dynamic>::Zero(6, n_tendons);
        Eigen::Matrix<double, 6, Eigen::Dynamic> d_wrench_d_tension_out =
            Eigen::Matrix<double, 6, Eigen::Dynamic>::Zero(6, n_tendons);
        Matrix66 d_wrench_d_pose = Matrix66::Zero();
        Matrix66 d_wrench_d_pose_prev = Matrix66::Zero();
        Matrix66 d_wrench_d_pose_next = Matrix66::Zero();

        for (int tendon_idx = 0; tendon_idx < n_tendons; ++tendon_idx) {
            Vector6 d_wrench_prev_d_tension = Vector6::Zero();
            Matrix6 d_wrench_prev_d_pose = Matrix6::Zero(), d_wrench_prev_d_pose_prev = Matrix6::Zero();

            Vector6 wrench_prev = get_single_tendon_wrench(
                tension_in[tendon_idx], pose, pose_prev,
                holes_[tendon_idx], holes_prev_[tendon_idx],
                H5 ? &d_wrench_prev_d_tension : 0,
                H2 ? &d_wrench_prev_d_pose    : 0,
                H1 ? &d_wrench_prev_d_pose_prev : 0);

            wrench_tendons += wrench_prev;
            d_wrench_d_pose      += d_wrench_prev_d_pose;
            d_wrench_d_pose_prev += d_wrench_prev_d_pose_prev;

            if (!is_tip_) {
                Vector6 d_wrench_next_d_tension = Vector6::Zero();
                Matrix6 d_wrench_next_d_pose = Matrix6::Zero(), d_wrench_next_d_pose_next = Matrix6::Zero();

                Vector6 wrench_next = get_single_tendon_wrench(
                    tension_out[tendon_idx], pose, pose_next,
                    holes_[tendon_idx], holes_next_[tendon_idx],
                    H6 ? &d_wrench_next_d_tension    : 0,
                    H2 ? &d_wrench_next_d_pose        : 0,
                    H3 ? &d_wrench_next_d_pose_next   : 0);

                wrench_tendons += wrench_next;
                d_wrench_d_tension_out.col(tendon_idx) = d_wrench_next_d_tension;
                d_wrench_d_pose     += d_wrench_next_d_pose;
                d_wrench_d_pose_next += d_wrench_next_d_pose_next;
            }

            d_wrench_d_tension_in.col(tendon_idx) = d_wrench_prev_d_tension;
        }

        if (H1) *H1 = -d_wrench_d_pose_prev;
        if (H2) *H2 = -d_wrench_d_pose;
        if (H3) *H3 = -d_wrench_d_pose_next;
        if (H4) *H4 = Matrix6::Identity();
        if (H5) *H5 = -d_wrench_d_tension_in;
        if (H6) *H6 = -d_wrench_d_tension_out;

        return wrench - wrench_tendons;
    }
};

// ============================================================================
// Tendon friction factors (extension)  these implement disc-to-disc Coulomb friction from a
// separate reference.
// ============================================================================

// Friction-free tendon: tension is identical in adjacent routing spans.
class TendonTensionEqualityFactor
    : public NoiseModelFactorN<Eigen::VectorXd, Eigen::VectorXd> {
public:
    using Base = NoiseModelFactorN<Eigen::VectorXd, Eigen::VectorXd>;
    using Base::evaluateError;

    TendonTensionEqualityFactor(Key tension_in_key, Key tension_out_key,
                                const SharedNoiseModel& model)
        : Base(model, tension_in_key, tension_out_key) {}

    Vector evaluateError(const Eigen::VectorXd& tension_in,
                         const Eigen::VectorXd& tension_out,
                         OptionalMatrixType H1,
                         OptionalMatrixType H2) const override {
        const int n = static_cast<int>(tension_in.size());
        if (H1) *H1 = -Matrix::Identity(n, n);
        if (H2) *H2 = Matrix::Identity(n, n);
        return tension_out - tension_in;
    }
};

// Disc-to-disc Coulomb friction using the capstan relation from
// Feliu-Talegon et al., IEEE/ASME T-Mech. 2025, Eq. (10):
// T_{i+1}=exp(-mu*phi_i)T_i, phi_i=acos(-v_L^T v_R).
class CapstanTendonTensionFactor
    : public NoiseModelFactorN<Pose3, Pose3, Pose3, Eigen::VectorXd, Eigen::VectorXd> {
    double mu_;
    std::vector<Point3> holes_prev_, holes_, holes_next_;
public:
    using Base = NoiseModelFactorN<Pose3, Pose3, Pose3, Eigen::VectorXd, Eigen::VectorXd>;
    using Base::evaluateError;

    CapstanTendonTensionFactor(Key pose_prev_key, Key pose_key, Key pose_next_key,
                               Key tension_in_key, Key tension_out_key,
                               double mu, const std::vector<Point3>& holes_prev,
                               const std::vector<Point3>& holes,
                               const std::vector<Point3>& holes_next,
                               const SharedNoiseModel& model)
        : Base(model, pose_prev_key, pose_key, pose_next_key, tension_in_key, tension_out_key),
          mu_(mu), holes_prev_(holes_prev), holes_(holes), holes_next_(holes_next) {}

    Vector evaluateError(const Pose3& pose_prev, const Pose3& pose, const Pose3& pose_next,
                         const Eigen::VectorXd& tension_in, const Eigen::VectorXd& tension_out,
                         OptionalMatrixType H1, OptionalMatrixType H2, OptionalMatrixType H3,
                         OptionalMatrixType H4, OptionalMatrixType H5) const override {
        const int n = static_cast<int>(tension_in.size());
        Vector error(n);
        Matrix dphi_prev = Matrix::Zero(n, 6);
        Matrix dphi_pose = Matrix::Zero(n, 6);
        Matrix dphi_next = Matrix::Zero(n, 6);
        Vector attenuation(n);

        for (int j = 0; j < n; ++j) {
            Matrix36 d_prev_world_d_prev, d_prev_local_d_pose;
            Matrix3 d_prev_local_d_prev_world;
            const Point3 prev_world = pose_prev.transformFrom(holes_prev_[j],
                H1 ? &d_prev_world_d_prev : 0);
            const Point3 prev_local = pose.transformTo(prev_world,
                H2 ? &d_prev_local_d_pose : 0, d_prev_local_d_prev_world);

            Matrix36 d_next_world_d_next, d_next_local_d_pose;
            Matrix3 d_next_local_d_next_world;
            const Point3 next_world = pose_next.transformFrom(holes_next_[j],
                H3 ? &d_next_world_d_next : 0);
            const Point3 next_local = pose.transformTo(next_world,
                H2 ? &d_next_local_d_pose : 0, d_next_local_d_next_world);

            const Vector3 h_left = prev_local - holes_[j];
            const Vector3 h_right = next_local - holes_[j];
            const double left_norm = h_left.norm();
            const double right_norm = h_right.norm();
            double phi = 0.0;
            if (left_norm > 1e-9 && right_norm > 1e-9) {
                Matrix3 dvl_dhl, dvr_dhr;
                const Vector3 v_left = normalize(h_left, (H1 || H2) ? &dvl_dhl : 0);
                const Vector3 v_right = normalize(h_right, (H2 || H3) ? &dvr_dhr : 0);
                const double c = std::clamp(-v_left.dot(v_right), -1.0, 1.0);
                phi = std::acos(c);
                const double sin_phi_sq = 1.0 - c * c;
                // The derivative of acos is singular in the exactly straight state.
                // The capstan loss itself is zero there, so use the limiting zero
                // linearisation and let the next nonlinear iteration update it.
                if (sin_phi_sq > 1e-10) {
                    const double inv_sin = 1.0 / std::sqrt(sin_phi_sq);
                    const Eigen::Matrix<double, 1, 3> dphi_dvl = v_right.transpose() * inv_sin;
                    const Eigen::Matrix<double, 1, 3> dphi_dvr = v_left.transpose() * inv_sin;
                    const Matrix36 dvl_prev = dvl_dhl * d_prev_local_d_prev_world * d_prev_world_d_prev;
                    const Matrix36 dvl_pose = dvl_dhl * d_prev_local_d_pose;
                    const Matrix36 dvr_pose = dvr_dhr * d_next_local_d_pose;
                    const Matrix36 dvr_next = dvr_dhr * d_next_local_d_next_world * d_next_world_d_next;
                    dphi_prev.row(j) = dphi_dvl * dvl_prev;
                    dphi_pose.row(j) = dphi_dvl * dvl_pose + dphi_dvr * dvr_pose;
                    dphi_next.row(j) = dphi_dvr * dvr_next;
                }
            }
            attenuation[j] = std::exp(-mu_ * phi);
            error[j] = tension_out[j] - attenuation[j] * tension_in[j];
        }

        if (H1) {
            *H1 = Matrix::Zero(n, 6);
            for (int j = 0; j < n; ++j) H1->row(j) = mu_ * attenuation[j] * tension_in[j] * dphi_prev.row(j);
        }
        if (H2) {
            *H2 = Matrix::Zero(n, 6);
            for (int j = 0; j < n; ++j) H2->row(j) = mu_ * attenuation[j] * tension_in[j] * dphi_pose.row(j);
        }
        if (H3) {
            *H3 = Matrix::Zero(n, 6);
            for (int j = 0; j < n; ++j) H3->row(j) = mu_ * attenuation[j] * tension_in[j] * dphi_next.row(j);
        }
        if (H4) *H4 = -attenuation.asDiagonal().toDenseMatrix();
        if (H5) *H5 = Matrix::Identity(n, n);
        return error;
    }
};

// ============================================================================
// Measurement factors (paper Section III-A2)
// ============================================================================

// Strain measurement factor (paper Eq. 33; paper Fig. 2: green, measurement
// factor).
class AngularStrainFactor: public NoiseModelFactorN<Vector6> {
    Vector3 du_meas_;
    Matrix6 K_inv_;
public:
    using Base = NoiseModelFactorN<Vector6>;
    using Base::evaluateError;

    AngularStrainFactor(Key stress_key,
                        const Vector3& du_meas,
                        const Matrix6& K_inv,
                        const SharedNoiseModel& model):
        Base(model, stress_key), du_meas_(du_meas), K_inv_(K_inv) {}

    Vector evaluateError(const Vector6& stress, OptionalMatrixType H1) const override {
        // strain = K_inv * stress; angular part = strain.head<3>()
        Vector3 du_pred = (K_inv_ * stress).head<3>();
        Vector3 error = du_pred - du_meas_;

        if (H1) {
            *H1 = K_inv_.topRows<3>();
        }

        return error;
    }
};

// Pose measurement factor (paper Eq. 31), position-only (paper Fig. 2:
// green, measurement factor).
class PositionMeasurementFactor: public NoiseModelFactorN<Pose3> {
    Vector3 position_meas_;

public:
    using Base = NoiseModelFactorN<Pose3>;
    using Base::evaluateError;

    PositionMeasurementFactor(Key pose_key,
                              Vector3 position_meas,
                              const SharedNoiseModel& model):
        Base(model, pose_key), position_meas_(position_meas) {}

    Vector evaluateError(
        const Pose3& pose,
        OptionalMatrixType H1) const override
    {
        Matrix36 d_position_d_pose;
        Vector3 error = pose.translation(d_position_d_pose) - position_meas_;

        if (H1) {
            *H1 = d_position_d_pose;
        }

        return error;
    }
};


}
