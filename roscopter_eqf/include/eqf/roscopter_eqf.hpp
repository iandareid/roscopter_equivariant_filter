#ifndef ROSCOPTER_EQF_HPP
#define ROSCOPTER_EQF_HPP

#include <Eigen/Dense>
#include <geometry_msgs/msg/vector3_stamped.hpp>
#include <rclcpp/rclcpp.hpp>

#include <cstdint>

#include "eqf/covariance.hpp"
#include "eqf/error_state.hpp"
#include "eqf/error_dynamics.hpp"
#include "eqf/equivariant_state.hpp"
#include "eqf/initialization.hpp"
#include "eqf/process_noise.hpp"
#include "eqf/robust_measurement.hpp"
#include "eqf/state.hpp"
#include "eqf/utils.hpp"
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
  bool initialize_heading(const Input & input);
  void propagation_step(const Input & input);
  void measurement_update(const Input & input);
  void gnss_position_measurement_update_step(const Input & input);
  void gnss_velocity_measurement_update_step(const Input & input);
  void magnetometer_measurement_update_step(const Input & input);
  eqf::RobustUpdateConfig robust_update_config(double alpha);
  void commit_measurement_update(
    const eqf::RobustUpdateResult & result,
    const char * sensor_name);
  void publish_accel_bias();
  bool calc_mag_field_properties(const Input & input);
  Eigen::Vector3d calculate_magnetic_reference(
    double declination_rad,
    double inclination_rad) const;
  double declination_rad() const;
  double inclination_rad() const;

  eqf::State state_;
  EqfCovariance P_;
  EqfCovariance A_begin_ = EqfCovariance::Zero();
  eqf::Matrix21x18d L_ = eqf::Matrix21x18d::Zero();
  eqf::Matrix18d Qc_ = eqf::Matrix18d::Zero();
  eqf::State reference_state_;
  eqf::ErrorState error_state_;
  eqf::EquivariantState observer_state_;
  Eigen::Vector3d magnetic_reference_G_ = Eigen::Vector3d::UnitX();
  double declination_ = 0.0;
  double inclination_ = NOT_IN_USE;
  bool declination_fallback_warned_ = false;
  bool state_initialized_ = false;
  bool has_last_time_ = false;
  int64_t last_imu_stamp_nanoseconds_ = 0;
  rclcpp::Publisher<geometry_msgs::msg::Vector3Stamped>::SharedPtr accel_bias_pub_;
};

#endif  // ROSCOPTER_EQF_HPP
