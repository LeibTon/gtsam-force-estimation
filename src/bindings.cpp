#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "tendon_robot.h"

namespace py = pybind11;
using namespace gtsam;


PYBIND11_MODULE(tendon_robot, m) {
    py::class_<TendonDiscConfig>(m, "TendonDiscConfig")
        .def(py::init<>())
        .def_readwrite("num_tendons", &TendonDiscConfig::num_tendons)
        .def_readwrite("num_discs", &TendonDiscConfig::num_discs)
        .def_readwrite("routing_radius", &TendonDiscConfig::routing_radius)
        .def_readwrite("disc_pose_idx", &TendonDiscConfig::disc_pose_idx)
        .def_readwrite("local_holes", &TendonDiscConfig::local_holes);

    py::class_<TendonRobotSolution>(m, "TendonRobotSolution")
        .def(py::init<>())
        .def_readwrite("backbone_pose_mean", &TendonRobotSolution::backbone_pose_mean)
        .def_readwrite("backbone_pose_cov", &TendonRobotSolution::backbone_pose_cov)
        .def_readwrite("pose_samples", &TendonRobotSolution::pose_samples)
        .def_readwrite("applied_wrench_mean", &TendonRobotSolution::applied_wrench_mean)
        .def_readwrite("applied_wrench_cov", &TendonRobotSolution::applied_wrench_cov)
        .def_readwrite("tensions_mean", &TendonRobotSolution::tensions_mean)
        .def_readwrite("tensions_cov", &TendonRobotSolution::tensions_cov)
        .def_readwrite("J_pose_tensions", &TendonRobotSolution::J_pose_tensions)
        .def_readwrite("gamma_mean", &TendonRobotSolution::gamma_mean)
        .def_readwrite("gamma_cov", &TendonRobotSolution::gamma_cov)
        .def_readwrite("solve_time_ms", &TendonRobotSolution::solve_time_ms)
        .def_readwrite("extract_time_ms", &TendonRobotSolution::extract_time_ms)
        .def_readwrite("total_time_ms", &TendonRobotSolution::total_time_ms)
        .def_readwrite("num_iterations", &TendonRobotSolution::num_iterations)
        .def_readwrite("tendon_disc_config", &TendonRobotSolution::tendon_disc_config);

    py::enum_<RoutingAngleFunction>(m, "RoutingAngleFunction")
        .value("CONSTANT", RoutingAngleFunction::CONSTANT)
        .value("LINEAR", RoutingAngleFunction::LINEAR)
        .export_values();

    py::class_<RoutingFunctionParams>(m, "RoutingFunctionParams")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("angle_offset"), py::arg("total_angle"))
        .def_readwrite("angle_offset", &RoutingFunctionParams::angle_offset)
        .def_readwrite("total_angle", &RoutingFunctionParams::total_angle);

    py::class_<TendonRobotConfig>(m, "TendonRobotConfig")
        .def(py::init<>())
        .def_readwrite("num_discs", &TendonRobotConfig::num_discs)
        .def_readwrite("poses_between_discs", &TendonRobotConfig::poses_between_discs)
        .def_readwrite("rod_length", &TendonRobotConfig::rod_length)
        .def_readwrite("rod_diameter", &TendonRobotConfig::rod_diameter)
        .def_readwrite("youngs_modulus", &TendonRobotConfig::youngs_modulus)
        .def_readwrite("shear_modulus", &TendonRobotConfig::shear_modulus)
        .def_readwrite("routing_radius", &TendonRobotConfig::routing_radius)
        .def_readwrite("tendon_friction_coefficient", &TendonRobotConfig::tendon_friction_coefficient)
        .def_readwrite("tendon_friction_std", &TendonRobotConfig::tendon_friction_std)
        .def_readwrite("use_midpoint", &TendonRobotConfig::use_midpoint)

        .def_readwrite("gravity_magnitude", &TendonRobotConfig::gravity_magnitude)
        .def_readwrite("gravity_direction", &TendonRobotConfig::gravity_direction)
        .def_readwrite("rod_mass_per_length", &TendonRobotConfig::rod_mass_per_length)
        .def_readwrite("disc_mass", &TendonRobotConfig::disc_mass)

        .def_readwrite("cosserat_twist_r_std", &TendonRobotConfig::cosserat_twist_r_std)
        .def_readwrite("small_force_std", &TendonRobotConfig::small_force_std)
        .def_readwrite("small_moment_std", &TendonRobotConfig::small_moment_std)
        .def_readwrite("small_r_std", &TendonRobotConfig::small_r_std)
        .def_readwrite("small_p_std", &TendonRobotConfig::small_p_std)

        .def_readwrite("tension_meas_std", &TendonRobotConfig::tension_meas_std)
        .def_readwrite("tip_position_meas_std", &TendonRobotConfig::tip_position_meas_std)
        .def_readwrite("angular_strain_meas_std", &TendonRobotConfig::angular_strain_meas_std)

        .def_readwrite("angle_functions", &TendonRobotConfig::angle_functions)
        .def_readwrite("angle_params", &TendonRobotConfig::angle_params)
        .def_readwrite("tendon_routing_radii", &TendonRobotConfig::tendon_routing_radii)
        .def_readwrite("tendon_reach_norm", &TendonRobotConfig::tendon_reach_norm)

        .def_readwrite("num_basis_functions", &TendonRobotConfig::num_basis_functions)
        .def_readwrite("basis_sigma", &TendonRobotConfig::basis_sigma)
        .def_readwrite("gamma_alpha_prior_std", &TendonRobotConfig::gamma_alpha_prior_std)
        .def_readwrite("gamma_beta_prior_std", &TendonRobotConfig::gamma_beta_prior_std);

    py::class_<GaussianBasisForceEstimator>(m, "GaussianBasisForceEstimator")
    .def(py::init<const TendonRobotConfig&>())
    .def_readonly("tendon_disc_config", &GaussianBasisForceEstimator::tendon_config_)
    .def_readonly("num_backbone_poses", &GaussianBasisForceEstimator::num_backbone_poses_)
    .def("step", &GaussianBasisForceEstimator::step,
         py::arg("tensions_meas"),
         py::arg("position_meas"),
         py::arg("angular_strain_meas"),
         py::arg("num_samples"),
         py::call_guard<py::gil_scoped_release>())
    .def("step_forward", &GaussianBasisForceEstimator::step_forward,
         py::arg("tensions_meas"),
         py::arg("gamma_true"),
         py::call_guard<py::gil_scoped_release>());
}
