#include "eqf/roscopter_eqf.hpp"

#include "ekf/geomag.h"

#include <cmath>
#include <limits>
#include <memory>

EstimatorEQF::EstimatorEQF()
: roscopter::EstimatorROS()
{
  declare_parameters();
  initialize_state();
  initialize_covariance();
  if (load_magnetic_model() != 0) {
    RCLCPP_ERROR(this->get_logger(), "Failed to load magnetic model.");
  }

  accel_bias_pub_ =
    this->create_publisher<geometry_msgs::msg::Vector3Stamped>("/estimated_state/accel_bias", 10);
}

void EstimatorEQF::declare_parameters()
{
  params_.declare_double("gravity", 9.80665);
  params_.declare_double("declination", 10.6);
  params_.declare_double("inclination", NOT_IN_USE);
  params_.declare_bool("force_declination_param", false);
  params_.declare_bool("require_gnss_for_heading_init", false);
  params_.declare_double("gps_lever_arm_north", 0.0);
  params_.declare_double("gps_lever_arm_east", 0.0);
  params_.declare_double("gps_lever_arm_down", 0.0);
  params_.declare_double("init_covariance_att", std::pow(0.05, 2));
  params_.declare_double("init_covariance_pos", std::pow(0.5, 2));
  params_.declare_double("init_covariance_vel", std::pow(0.05, 2));
  params_.declare_double("init_covariance_gyro_bias", std::pow(0.001, 2));
  params_.declare_double("init_covariance_accel_bias", std::pow(0.001, 2));
  params_.declare_double("init_covariance_gnss_lever_arm", std::pow(0.001, 2));
  params_.declare_double("init_covariance_mag_rotation", std::pow(0.001, 2));
}

void EstimatorEQF::initialize_state()
{
  reference_state_ = eqf::EqfState::Identity();
  reference_state_.setGnssLeverArmImu(
    Eigen::Vector3d(
      this->get_parameter("gps_lever_arm_north").as_double(),
      this->get_parameter("gps_lever_arm_east").as_double(),
      this->get_parameter("gps_lever_arm_down").as_double()));
  observer_state_.reset();
  state_ = eqf::phi(observer_state_, reference_state_);
  has_last_time_ = false;
}

void EstimatorEQF::initialize_covariance()
{
  P_.setZero();
  P_.block<3, 3>(eqf::ErrorState::kAttitudeOffset, eqf::ErrorState::kAttitudeOffset) =
    this->get_parameter("init_covariance_att").as_double() * Eigen::Matrix3d::Identity();
  P_.block<3, 3>(eqf::ErrorState::kVelocityOffset, eqf::ErrorState::kVelocityOffset) =
    this->get_parameter("init_covariance_vel").as_double() * Eigen::Matrix3d::Identity();
  P_.block<3, 3>(eqf::ErrorState::kPositionOffset, eqf::ErrorState::kPositionOffset) =
    this->get_parameter("init_covariance_pos").as_double() * Eigen::Matrix3d::Identity();
  P_.block<3, 3>(eqf::ErrorState::kGammaGyroOffset, eqf::ErrorState::kGammaGyroOffset) =
    this->get_parameter("init_covariance_gyro_bias").as_double() * Eigen::Matrix3d::Identity();
  P_.block<3, 3>(eqf::ErrorState::kGammaAccelOffset, eqf::ErrorState::kGammaAccelOffset) =
    this->get_parameter("init_covariance_accel_bias").as_double() * Eigen::Matrix3d::Identity();
  P_.block<3, 3>(eqf::ErrorState::kLeverArmOffset, eqf::ErrorState::kLeverArmOffset) =
    this->get_parameter("init_covariance_gnss_lever_arm").as_double() *
    Eigen::Matrix3d::Identity();
  P_.block<3, 3>(eqf::ErrorState::kMagRotationOffset, eqf::ErrorState::kMagRotationOffset) =
    this->get_parameter("init_covariance_mag_rotation").as_double() *
    Eigen::Matrix3d::Identity();
}

void EstimatorEQF::estimate(const Input & input, Output & output)
{
  if (!mag_init_) {
    calc_mag_field_properties(input);
  }

  propagation_step(input);
  measurement_update(input);

  const Eigen::Vector3d euler = state_.eulerAngles();

  output.pn = state_.pn();
  output.pe = state_.pe();
  output.pd = state_.pd();

  output.vx = state_.v_x();
  output.vy = state_.v_y();
  output.vz = state_.v_z();

  output.phi = euler.x();
  output.theta = euler.y();
  output.psi = euler.z();
  output.quat = state_.quaternionGlobalFromImu().cast<float>();

  output.p = input.gyro_x - state_.b_gx();
  output.q = input.gyro_y - state_.b_gy();
  output.r = input.gyro_z - state_.b_gz();

  output.bx = state_.b_gx();
  output.by = state_.b_gy();
  output.bz = state_.b_gz();

  output.nees_valid = false;

  publish_accel_bias();
}

void EstimatorEQF::propagation_step(const Input & input)
{
  const eqf::ImuInput imu_input{
    Eigen::Vector3d(input.gyro_x, input.gyro_y, input.gyro_z),
    Eigen::Vector3d(input.accel_x, input.accel_y, input.accel_z)};

  if (!imu_input.angular_velocity.allFinite() || !imu_input.linear_acceleration.allFinite()) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(),
      *this->get_clock(),
      1000,
      "Skipping EqF propagation because the IMU input contains NaN or infinity.");
    return;
  }

  const rclcpp::Time current_time = this->get_clock()->now();
  if (!has_last_time_) {
    last_time_ = current_time;
    has_last_time_ = true;
    state_ = eqf::phi(observer_state_, reference_state_);
    return;
  }

  const rclcpp::Duration dt_duration = current_time - last_time_;
  last_time_ = current_time;
  const double dt = dt_duration.nanoseconds() / 1e9;

  if (!std::isfinite(dt) || dt <= 0.0) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(),
      *this->get_clock(),
      1000,
      "Skipping EqF propagation because the computed time step is not positive and finite.");
    state_ = eqf::phi(observer_state_, reference_state_);
    return;
  }

  const Eigen::Vector3d gravity_in_inertial(
    0.0,
    0.0,
    this->get_parameter("gravity").as_double());
  const eqf::ImuPropagationResult propagation = eqf::propagateImuBeforeA(
    observer_state_,
    reference_state_,
    imu_input,
    gravity_in_inertial,
    dt);

  if (!propagation.propagated) {
    state_ = eqf::phi(observer_state_, reference_state_);
    return;
  }

  observer_state_ = propagation.observer_after;
  state_ = propagation.state_after;
}

void EstimatorEQF::measurement_update(const Input &)
{
}

void EstimatorEQF::calculate_nees(const roscopter_msgs::msg::State &, Output & output)
{
  output.nees_valid = false;
}

void EstimatorEQF::publish_accel_bias()
{
  geometry_msgs::msg::Vector3Stamped accel_bias_msg;
  accel_bias_msg.header.stamp = this->get_clock()->now();
  accel_bias_msg.vector.x = state_.b_ax();
  accel_bias_msg.vector.y = state_.b_ay();
  accel_bias_msg.vector.z = state_.b_az();
  accel_bias_pub_->publish(accel_bias_msg);
}

double EstimatorEQF::declination_rad() const
{
  return declination_ * M_PI / 180.0;
}

double EstimatorEQF::inclination_rad() const
{
  if (std::abs(inclination_ - NOT_IN_USE) <= std::numeric_limits<double>::epsilon()) {
    return 0.0;
  }

  return inclination_ * M_PI / 180.0;
}

Eigen::Vector3d EstimatorEQF::calculate_magnetic_reference(
  const double declination_rad,
  const double inclination_rad) const
{
  const Eigen::Matrix3d inclination_rotation =
    Eigen::AngleAxisd(-inclination_rad, Eigen::Vector3d::UnitY()).toRotationMatrix();
  const Eigen::Matrix3d declination_rotation =
    Eigen::AngleAxisd(declination_rad, Eigen::Vector3d::UnitZ()).toRotationMatrix();

  const Eigen::Vector3d reference =
    declination_rotation * inclination_rotation * Eigen::Vector3d::UnitX();
  const double reference_norm = reference.norm();
  if (!std::isfinite(reference_norm) || reference_norm < 1e-12) {
    return Eigen::Vector3d::UnitX();
  }

  return reference / reference_norm;
}

bool EstimatorEQF::calc_mag_field_properties(const Input & input)
{
  const bool require_gnss = this->get_parameter("require_gnss_for_heading_init").as_bool();
  const bool force_param = this->get_parameter("force_declination_param").as_bool();
  const double configured_declination = this->get_parameter("declination").as_double();
  const double configured_inclination = this->get_parameter("inclination").as_double();

  auto use_configured_mag = [&]() {
      declination_ = configured_declination;
      inclination_ = configured_inclination;
      magnetic_reference_G_ = calculate_magnetic_reference(declination_rad(), inclination_rad());
      mag_init_ = true;
      return true;
    };

  if (require_gnss && !gps_init_) {
    return false;
  }

  if (force_param) {
    return use_configured_mag();
  }

  if (has_fix_) {
    if (!input.gps_new) {
      return false;
    }

    double total_intensity;
    double grid_variation;
    const double decimal_year =
      static_cast<double>(input.gps_year) + static_cast<double>(input.gps_yday) / 365.0;

    const int mag_success = geomag_calc(
      input.gps_alt / 1000.0,
      input.gps_lat,
      input.gps_lon,
      decimal_year,
      &declination_,
      &inclination_,
      &total_intensity,
      &grid_variation);

    if (mag_success == 0) {
      magnetic_reference_G_ = calculate_magnetic_reference(declination_rad(), inclination_rad());
      mag_init_ = true;
      return true;
    }

    RCLCPP_ERROR(
      this->get_logger(),
      "Something went wrong while calculating inclination and declination.");
    RCLCPP_ERROR(
      this->get_logger(),
      "Inclination and declination not set, estimation will likely be poor.");
    return false;
  }

  if (!require_gnss) {
    if (!declination_fallback_warned_) {
      RCLCPP_WARN(
        this->get_logger(),
        "GNSS unavailable for magnetic model lookup; using configured magnetic field parameters.");
      declination_fallback_warned_ = true;
    }
    return use_configured_mag();
  }

  return false;
}

int main(int argc, char ** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<EstimatorEQF>());
  rclcpp::shutdown();
  return 0;
}
