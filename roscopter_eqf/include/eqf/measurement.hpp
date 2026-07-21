#ifndef ROSCOPTER_EQF_MEASUREMENT_HPP
#define ROSCOPTER_EQF_MEASUREMENT_HPP

#include <Eigen/Core>

#include <cmath>
#include <stdexcept>

#include "eqf/equivariant_state.hpp"
#include "eqf/utils.hpp"

namespace eqf
{

constexpr double kDirectionNormTolerance = 1e-12;

[[nodiscard]] inline EquivariantState groupInverse(const EquivariantState & X)
{
  const Eigen::Matrix3d A_transpose = X.rotationBlock().transpose();

  EquivariantState inverse;
  EquivariantState::SE23Matrix C_inverse = EquivariantState::SE23Matrix::Identity();
  C_inverse.block<3, 3>(0, 0) = A_transpose;
  C_inverse.block<3, 1>(0, 3) = -A_transpose * X.firstTranslationColumn();
  C_inverse.block<3, 1>(0, 4) = -A_transpose * X.secondTranslationColumn();
  inverse.setC(C_inverse);
  inverse.setDelta(-A_transpose * X.delta());
  inverse.setE(X.E().transpose());

  const Vector6d inverse_gamma =
    invertBiasLinearPart(inverse, -biasVector(X));
  inverse.setGammaGyro(inverse_gamma.head<3>());
  inverse.setGammaAccel(inverse_gamma.tail<3>());
  return inverse;
}

[[nodiscard]] inline Eigen::Vector3d rhoMagnetometer(
  const EquivariantState & X,
  const Eigen::Vector3d & measured_direction)
{
  if (!measured_direction.allFinite()) {
    throw std::invalid_argument("Magnetometer direction contains NaN or infinity");
  }
  return X.E().transpose() * measured_direction;
}

[[nodiscard]] inline Eigen::Vector3d rhoGnssPosition(
  const EquivariantState & X,
  const Eigen::Vector3d & configured_position_output)
{
  if (!configured_position_output.allFinite()) {
    throw std::invalid_argument("GNSS position output contains NaN or infinity");
  }
  return X.rotationBlock().transpose() *
         (configured_position_output - X.secondTranslationColumn() + X.delta());
}

[[nodiscard]] inline Eigen::Vector3d rhoGnssVelocity(
  const EquivariantState & X,
  const Eigen::Vector3d & configured_velocity_output,
  const Eigen::Vector3d & angular_rate_imu)
{
  if (!configured_velocity_output.allFinite() || !angular_rate_imu.allFinite()) {
    throw std::invalid_argument("GNSS velocity action input contains NaN or infinity");
  }
  return X.rotationBlock().transpose() *
         (configured_velocity_output - X.firstTranslationColumn() -
         X.delta().cross(angular_rate_imu));
}

[[nodiscard]] inline Eigen::Vector3d directionResidualMap(
  const Eigen::Vector3d & transformed_direction,
  const Eigen::Vector3d & reference_direction)
{
  if (!transformed_direction.allFinite() || !reference_direction.allFinite()) {
    throw std::invalid_argument("Direction residual input contains NaN or infinity");
  }
  const double transformed_norm = transformed_direction.norm();
  const double reference_norm = reference_direction.norm();
  if (transformed_norm <= kDirectionNormTolerance ||
    reference_norm <= kDirectionNormTolerance)
  {
    throw std::invalid_argument("Direction residual requires nonzero directions");
  }
  return (reference_direction / reference_norm).cross(
    transformed_direction / transformed_norm);
}

[[nodiscard]] inline Eigen::Vector3d positionResidualMap(
  const Eigen::Vector3d & transformed_position)
{
  if (!transformed_position.allFinite()) {
    throw std::invalid_argument("Position residual contains NaN or infinity");
  }
  return transformed_position;
}

[[nodiscard]] inline Eigen::Vector3d velocityResidualMap(
  const Eigen::Vector3d & transformed_velocity)
{
  if (!transformed_velocity.allFinite()) {
    throw std::invalid_argument("Velocity residual contains NaN or infinity");
  }
  return transformed_velocity;
}

[[nodiscard]] inline Eigen::Vector3d magnetometerResidual(
  const EquivariantState & observer_state,
  const Eigen::Vector3d & measured_direction,
  const Eigen::Vector3d & reference_direction)
{
  if (!measured_direction.allFinite() || !reference_direction.allFinite() ||
    measured_direction.norm() <= kDirectionNormTolerance ||
    reference_direction.norm() <= kDirectionNormTolerance)
  {
    throw std::invalid_argument("Magnetometer residual requires finite nonzero directions");
  }

  const Eigen::Vector3d normalized_measurement = measured_direction.normalized();
  const Eigen::Vector3d transformed =
    rhoMagnetometer(groupInverse(observer_state), normalized_measurement);
  return directionResidualMap(transformed, reference_direction);
}

[[nodiscard]] inline Eigen::Vector3d gnssPositionResidual(
  const EquivariantState & observer_state,
  const Eigen::Vector3d & configured_position_output)
{
  return positionResidualMap(
    rhoGnssPosition(groupInverse(observer_state), configured_position_output));
}

/** angular_rate_imu is the bias-corrected IMU-frame rate used by the lift. */
[[nodiscard]] inline Eigen::Vector3d gnssVelocityResidual(
  const EquivariantState & observer_state,
  const Eigen::Vector3d & configured_velocity_output,
  const Eigen::Vector3d & angular_rate_imu)
{
  return velocityResidualMap(
    rhoGnssVelocity(
      groupInverse(observer_state), configured_velocity_output, angular_rate_imu));
}

[[nodiscard]] inline Eigen::Matrix3d rotateResidualCovariance(
  const Eigen::Matrix3d & raw_covariance,
  const Eigen::Matrix3d & residual_from_raw_rotation)
{
  if (!raw_covariance.allFinite() || !residual_from_raw_rotation.allFinite()) {
    throw std::invalid_argument("Measurement covariance transform contains NaN or infinity");
  }
  Eigen::Matrix3d transformed =
    residual_from_raw_rotation * raw_covariance * residual_from_raw_rotation.transpose();
  return 0.5 * (transformed + transformed.transpose()).eval();
}

[[nodiscard]] inline Eigen::Matrix3d magnetometerResidualCovariance(
  const Eigen::Matrix3d & transformed_direction_covariance,
  const Eigen::Vector3d & reference_direction)
{
  if (!reference_direction.allFinite() ||
    reference_direction.norm() <= kDirectionNormTolerance)
  {
    throw std::invalid_argument("Magnetometer covariance requires a finite reference direction");
  }
  const Eigen::Matrix3d residual_jacobian =
    skew_matrix(reference_direction.normalized());
  return rotateResidualCovariance(transformed_direction_covariance, residual_jacobian);
}

}  // namespace eqf

#endif  // ROSCOPTER_EQF_MEASUREMENT_HPP
