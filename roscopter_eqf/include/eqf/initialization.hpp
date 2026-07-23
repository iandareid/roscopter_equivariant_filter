#ifndef ROSCOPTER_EQF_INITIALIZATION_HPP
#define ROSCOPTER_EQF_INITIALIZATION_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <optional>

namespace eqf
{

inline std::optional<double> headingFromMagnetometer(
  const Eigen::Vector3d & magnetic_field_imu,
  const double declination_rad,
  const double minimum_horizontal_norm = 1.0e-12)
{
  if (!magnetic_field_imu.allFinite() || !std::isfinite(declination_rad) ||
    !std::isfinite(minimum_horizontal_norm) || minimum_horizontal_norm <= 0.0)
  {
    return std::nullopt;
  }

  const double horizontal_norm = magnetic_field_imu.head<2>().norm();
  if (!std::isfinite(horizontal_norm) || horizontal_norm <= minimum_horizontal_norm) {
    return std::nullopt;
  }

  return -std::atan2(magnetic_field_imu.y(), magnetic_field_imu.x()) + declination_rad;
}

inline Eigen::Matrix3d rotationGlobalFromImuAtHeading(const double heading_rad)
{
  return Eigen::AngleAxisd(heading_rad, Eigen::Vector3d::UnitZ()).toRotationMatrix();
}

}  // namespace eqf

#endif  // ROSCOPTER_EQF_INITIALIZATION_HPP
