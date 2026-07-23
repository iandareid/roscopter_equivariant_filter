#include "eqf/roscopter_eqf.hpp"

#include "ekf/geomag.h"

#include <cmath>
#include <limits>
#include <memory>

EstimatorEQF::EstimatorEQF()
: roscopter::EstimatorROS()
{
  enable_event_driven_updates();
  declare_parameters();
  params_.set_parameters();
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
  params_.declare_double("declination", 10.6);
  params_.declare_double("inclination", NOT_IN_USE);
  params_.declare_bool("force_declination_param", false);
  params_.declare_bool("require_gnss_for_heading_init", true);
  params_.declare_double("gps_lever_arm_north", 0.0);
  params_.declare_double("gps_lever_arm_east", 0.0);
  params_.declare_double("gps_lever_arm_down", 0.0);
  params_.declare_double("init_covariance_att", std::pow(5.0 * 0.017, 2));
  params_.declare_double("init_covariance_pos", 0.0001);
  params_.declare_double("init_covariance_vel", 0.0001);
  params_.declare_double("init_covariance_gyro_bias", 0.0001);
  params_.declare_double("init_covariance_accel_bias", std::pow(0.001, 2));
  params_.declare_double("init_covariance_gnss_lever_arm", std::pow(0.001, 2));
  params_.declare_double("init_covariance_mag_rotation", std::pow(0.001, 2));
  params_.declare_bool("eqf_measurement_update_enabled", true);
  params_.declare_bool("eqf_robust_inflation_enabled", true);
  params_.declare_double("eqf_robust_alpha_magnetometer", 0.5);
  params_.declare_double("eqf_robust_alpha_gnss_position", 0.5);
  params_.declare_double("eqf_robust_alpha_gnss_velocity", 0.5);
  params_.declare_double("eqf_innovation_covariance_jitter", 1.0e-9);
  params_.declare_double("eqf_covariance_jitter", 1.0e-12);
  params_.declare_double("eqf_minimum_direction_norm", 1.0e-9);
  params_.declare_double("eqf_maximum_condition_number", 1.0e12);
  params_.declare_double("sigma_n_gps", 7.0);
  params_.declare_double("sigma_e_gps", 7.0);
  params_.declare_double("sigma_d_gps", 7.0);
  params_.declare_double("sigma_vn_gps", 0.08);
  params_.declare_double("sigma_ve_gps", 0.08);
  params_.declare_double("sigma_vd_gps", 0.08);
  params_.declare_double("sigma_mag", 0.03);
  params_.declare_double("eqf_gyro_noise_density", 0.4 * M_PI / 180.0);
  params_.declare_double("eqf_accel_noise_density", 0.6);
  params_.declare_double("eqf_gyro_bias_random_walk", 1.0e-7);
  params_.declare_double("eqf_accel_bias_random_walk", 0.001);
  params_.declare_double("eqf_gnss_lever_arm_random_walk", 0.0);
  params_.declare_double("eqf_mag_calibration_random_walk", 0.0);
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
  state_initialized_ = false;
  has_last_time_ = false;
  last_imu_stamp_nanoseconds_ = 0;
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

  eqf::ProcessNoise noise;
  noise.gyro_noise_density.setConstant(params_.get_double("eqf_gyro_noise_density"));
  noise.accel_noise_density.setConstant(params_.get_double("eqf_accel_noise_density"));
  noise.gyro_bias_random_walk.setConstant(params_.get_double("eqf_gyro_bias_random_walk"));
  noise.accel_bias_random_walk.setConstant(params_.get_double("eqf_accel_bias_random_walk"));
  noise.gnss_lever_arm_random_walk.setConstant(
    params_.get_double("eqf_gnss_lever_arm_random_walk"));
  noise.mag_calibration_random_walk.setConstant(
    params_.get_double("eqf_mag_calibration_random_walk"));
  L_ = eqf::makeProcessNoiseMapL();
  Qc_ = eqf::makeContinuousProcessNoiseCovariance(noise);
}

void EstimatorEQF::estimate(const Input & input, Output & output)
{
  if (!state_initialized_) {
    if (!mag_init_) {
      calc_mag_field_properties(input);
    }
    if (!initialize_heading(input)) {
      return;
    }
    return;
  }

  if (input.imu_new) {
    propagation_step(input);
  }
  measurement_update(input);

  const Eigen::Vector3d euler = state_.eulerAngles();

  output.pn = state_.pn();
  output.pe = state_.pe();
  output.pd = state_.pd();

  const Eigen::Vector3d velocity_body =
    state_.rotationGlobalFromImu().transpose() * state_.velocityGlobal();
  output.vx = velocity_body.x();
  output.vy = velocity_body.y();
  output.vz = velocity_body.z();

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

bool EstimatorEQF::initialize_heading(const Input & input)
{
  if (!mag_init_ || !input.mag_new) {
    return false;
  }

  const std::optional<double> heading = eqf::headingFromMagnetometer(
    Eigen::Vector3d(input.mag_x, input.mag_y, input.mag_z), declination_rad());
  if (!heading.has_value()) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Waiting to initialize EqF heading because the magnetometer sample is invalid.");
    return false;
  }

  reference_state_.setRotationGlobalFromImu(
    eqf::rotationGlobalFromImuAtHeading(*heading));
  observer_state_.reset();
  error_state_.reset();
  state_ = eqf::phi(observer_state_, reference_state_);
  state_initialized_ = true;
  has_last_time_ = false;
  last_imu_stamp_nanoseconds_ = 0;

  RCLCPP_INFO(
    this->get_logger(), "Initialized EqF heading to %.6f rad from magnetometer.", *heading);
  return true;
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

  const int64_t current_stamp_nanoseconds = input.imu_stamp_nanoseconds;
  if (current_stamp_nanoseconds <= 0) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(),
      *this->get_clock(),
      1000,
      "Skipping EqF propagation because the IMU timestamp is invalid.");
    return;
  }

  if (!has_last_time_) {
    last_imu_stamp_nanoseconds_ = current_stamp_nanoseconds;
    has_last_time_ = true;
    state_ = eqf::phi(observer_state_, reference_state_);
    return;
  }

  const double dt =
    static_cast<double>(current_stamp_nanoseconds - last_imu_stamp_nanoseconds_) * 1e-9;

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
  // Equation (10) is evaluated at the beginning of this IMU interval.
  // A_begin_ is retained for covariance propagation once the paper-derived
  // 21x18 noise-input mapping L is available.
  const eqf::State state_begin = eqf::phi(observer_state_, reference_state_);
  A_begin_ = eqf::buildErrorDynamicsA(
    observer_state_, state_begin, imu_input, gravity_in_inertial);
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

  eqf::CovariancePropagationResult covariance_propagation;
  try {
    covariance_propagation =
      eqf::propagateCovarianceSecondOrder(P_, A_begin_, L_, Qc_, dt);
  } catch (const std::exception & exception) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Skipping EqF propagation because covariance propagation failed: %s", exception.what());
    state_ = eqf::phi(observer_state_, reference_state_);
    return;
  }

  last_imu_stamp_nanoseconds_ = current_stamp_nanoseconds;
  observer_state_ = propagation.observer_after;
  state_ = propagation.state_after;
  P_ = covariance_propagation.P_pred;
}

void EstimatorEQF::measurement_update(const Input & input)
{
  if (!params_.get_bool("eqf_measurement_update_enabled")) {
    return;
  }

  // Deterministic sequential policy: each successful update becomes the prior
  // for the next sensor.
  gnss_position_measurement_update_step(input);
  gnss_velocity_measurement_update_step(input);
  magnetometer_measurement_update_step(input);
}

eqf::RobustUpdateConfig EstimatorEQF::robust_update_config(const double alpha)
{
  eqf::RobustUpdateConfig config;
  config.alpha = alpha;
  config.robust_inflation_enabled = params_.get_bool("eqf_robust_inflation_enabled");
  config.innovation_covariance_jitter =
    params_.get_double("eqf_innovation_covariance_jitter");
  config.covariance_jitter = params_.get_double("eqf_covariance_jitter");
  config.maximum_condition_number = params_.get_double("eqf_maximum_condition_number");
  return config;
}

void EstimatorEQF::commit_measurement_update(
  const eqf::RobustUpdateResult & result,
  const char * sensor_name)
{
  const auto & diagnostics = result.diagnostics;
  if (!diagnostics.accepted) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Rejected EqF %s update because its inputs or numerical solve were invalid.", sensor_name);
    return;
  }

  observer_state_ = result.observer_after;
  state_ = result.state_after;
  P_ = result.covariance_after;
  error_state_.reset();

  if (diagnostics.numerical_regularization_used) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Regularized EqF %s measurement covariance.", sensor_name);
  }
  if (diagnostics.post_inflation_normalized_innovation > 1.0 + 1.0e-8) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "EqF %s post-inflation normalized innovation is %.6f (> 1).",
      sensor_name, diagnostics.post_inflation_normalized_innovation);
  }
  RCLCPP_DEBUG_THROTTLE(
    this->get_logger(), *this->get_clock(), 1000,
    "EqF %s update: r=%.6f beta=%.6f alpha=%.3f correction_norm=%.6e",
    sensor_name, diagnostics.normalized_innovation_magnitude, diagnostics.beta,
    diagnostics.alpha, diagnostics.delta_epsilon.norm());
}

void EstimatorEQF::gnss_position_measurement_update_step(const Input & input)
{
  if (!input.gps_new || !gps_init_) {
    return;
  }
  const Eigen::Vector3d measurement(input.gps_n, input.gps_e, -input.gps_h);
  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  const double sigma_n = params_.get_double("sigma_n_gps");
  const double sigma_e = params_.get_double("sigma_e_gps");
  const double sigma_d = params_.get_double("sigma_d_gps");
  if (!std::isfinite(sigma_n) || !std::isfinite(sigma_e) || !std::isfinite(sigma_d) ||
    sigma_n < 0.0 || sigma_e < 0.0 || sigma_d < 0.0)
  {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Skipping EqF GNSS position update because its noise parameters are invalid.");
    return;
  }
  covariance.diagonal() << sigma_n * sigma_n, sigma_e * sigma_e, sigma_d * sigma_d;

  try {
    const eqf::RobustUpdateResult result = eqf::applyRobustEqfCorrection(
      eqf::MeasurementType::kGnssPosition,
      eqf::gnssPositionInnovation(state_, measurement),
      eqf::computeGnssPositionCStar(observer_state_, measurement),
      covariance,
      robust_update_config(params_.get_double("eqf_robust_alpha_gnss_position")),
      reference_state_, observer_state_, P_);
    commit_measurement_update(result, "GNSS position");
  } catch (const std::exception & exception) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Skipping EqF GNSS position update: %s", exception.what());
  }
}

void EstimatorEQF::gnss_velocity_measurement_update_step(const Input & input)
{
  if (!input.gps_new || !gps_init_) {
    return;
  }
  const Eigen::Vector3d measurement(input.gps_vn, input.gps_ve, input.gps_vd);
  const Eigen::Vector3d measured_rate(input.gyro_x, input.gyro_y, input.gyro_z);
  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  const double sigma_n = params_.get_double("sigma_vn_gps");
  const double sigma_e = params_.get_double("sigma_ve_gps");
  const double sigma_d = params_.get_double("sigma_vd_gps");
  if (!std::isfinite(sigma_n) || !std::isfinite(sigma_e) || !std::isfinite(sigma_d) ||
    sigma_n < 0.0 || sigma_e < 0.0 || sigma_d < 0.0)
  {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Skipping EqF GNSS velocity update because its noise parameters are invalid.");
    return;
  }
  covariance.diagonal() << sigma_n * sigma_n, sigma_e * sigma_e, sigma_d * sigma_d;

  try {
    const Eigen::Vector3d corrected_rate = state_.correctedAngularRate(measured_rate);
    const eqf::RobustUpdateResult result = eqf::applyRobustEqfCorrection(
      eqf::MeasurementType::kGnssVelocity,
      eqf::gnssVelocityInnovation(state_, measurement, corrected_rate),
      eqf::computeGnssVelocityCStar(observer_state_, measurement, corrected_rate),
      covariance,
      robust_update_config(params_.get_double("eqf_robust_alpha_gnss_velocity")),
      reference_state_, observer_state_, P_);
    commit_measurement_update(result, "GNSS velocity");
  } catch (const std::exception & exception) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Skipping EqF GNSS velocity update: %s", exception.what());
  }
}

void EstimatorEQF::magnetometer_measurement_update_step(const Input & input)
{
  if (!input.mag_new || !mag_init_) {
    return;
  }
  const Eigen::Vector3d measurement(input.mag_x, input.mag_y, input.mag_z);
  const double minimum_norm = params_.get_double("eqf_minimum_direction_norm");
  if (!measurement.allFinite() || !std::isfinite(minimum_norm) ||
    minimum_norm <= 0.0 || measurement.norm() <= minimum_norm)
  {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Skipping EqF magnetometer update because the direction is invalid.");
    return;
  }

  const double sigma = params_.get_double("sigma_mag");
  if (!std::isfinite(sigma) || sigma < 0.0) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Skipping EqF magnetometer update because sigma_mag is invalid.");
    return;
  }
  const Eigen::Matrix3d covariance = sigma * sigma * Eigen::Matrix3d::Identity();
  try {
    const Eigen::Vector3d prediction =
      eqf::predictedMagnetometerDirection(state_, magnetic_reference_G_);
    const eqf::RobustUpdateResult result = eqf::applyRobustEqfCorrection(
      eqf::MeasurementType::kMagnetometer,
      eqf::magnetometerResidual(observer_state_, measurement, magnetic_reference_G_),
      eqf::computeMagnetometerCStar(
        observer_state_, prediction, magnetic_reference_G_),
      covariance,
      robust_update_config(params_.get_double("eqf_robust_alpha_magnetometer")),
      reference_state_, observer_state_, P_);
    commit_measurement_update(result, "magnetometer");
  } catch (const std::exception & exception) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 1000,
      "Skipping EqF magnetometer update: %s", exception.what());
  }
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

  auto configured_mag_is_valid = [&]() {
      return std::isfinite(configured_declination) &&
             std::isfinite(configured_inclination) &&
             configured_declination >= -180.0 && configured_declination <= 180.0 &&
             configured_inclination >= -90.0 && configured_inclination <= 90.0;
    };

  auto use_configured_mag = [&]() {
      if (!configured_mag_is_valid()) {
        RCLCPP_ERROR(
          this->get_logger(),
          "Configured magnetic inclination or declination is outside its valid range.");
        return false;
      }
      declination_ = configured_declination;
      inclination_ = configured_inclination;
      magnetic_reference_G_ = calculate_magnetic_reference(declination_rad(), inclination_rad());
      mag_init_ = true;
      return true;
    };

  if (force_param) {
    return use_configured_mag();
  }

  if (require_gnss && !has_fix_) {
    return false;
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
