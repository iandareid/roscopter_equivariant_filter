#ifndef ROSCOPTER_EQF_HPP
#define ROSCOPTER_EQF_HPP

#include <Eigen/Dense>
#include <geometry_msgs/msg/vector3_stamped.hpp>
#include <rclcpp/rclcpp.hpp>

#include "eqf/error_state.hpp"
#include "eqf/equivariant_state.hpp"
#include "eqf/state.hpp"
#include "ekf/estimator_ros.hpp"

class EstimatorEQF : public roscopter::EstimatorROS
{
public:
  EstimatorEQF();
  using EqfCovariance = Eigen::Matrix<double, 21, 21>;

private:
  void estimate(const Input & input, Output & output) override;
  void calculate_nees(const roscopter_msgs::msg::State & truth, Output & output) override;
  void declare_parameters();
  void initialize_state();
  void initialize_covariance();
  void publish_accel_bias();
  bool calc_mag_field_properties(const Input & input);
  Eigen::Vector3d calculate_magnetic_reference(
    double declination_rad,
    double inclination_rad) const;
  double declination_rad() const;
  double inclination_rad() const;

  eqf::State state_;
  EqfCovariance P_;
  eqf::State reference_state_;
  eqf::ErrorState error_state_;
  eqf::EquivariantState observer_state_;
  Eigen::Vector3d magnetic_reference_G_ = Eigen::Vector3d::UnitX();
  double declination_ = 0.0;
  double inclination_ = NOT_IN_USE;
  bool declination_fallback_warned_ = false;
  rclcpp::Publisher<geometry_msgs::msg::Vector3Stamped>::SharedPtr accel_bias_pub_;
};

#endif  // ROSCOPTER_EQF_HPP
