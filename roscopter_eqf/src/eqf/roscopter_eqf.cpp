#include "eqf/roscopter_eqf.hpp"

#include "ekf/geomag.h"

#include <cmath>
#include <limits>
#include <memory>

EstimatorEQF::EstimatorEQF()
: roscopter::EstimatorROS()
{
  declare_parameters();
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
}

void EstimatorEQF::estimate(const Input & input, Output & output)
{
  if (!mag_init_) {
    calc_mag_field_properties(input);
  }

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
