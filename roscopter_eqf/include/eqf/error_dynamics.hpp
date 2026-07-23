#ifndef ROSCOPTER_EQF_ERROR_DYNAMICS_HPP
#define ROSCOPTER_EQF_ERROR_DYNAMICS_HPP

#include <Eigen/Core>

#include <stdexcept>

#include "eqf/covariance.hpp"
#include "eqf/equivariant_state.hpp"
#include "eqf/state.hpp"
#include "eqf/utils.hpp"

namespace eqf
{

struct ErrorDynamicsQuantities
{
  Eigen::Vector3d omega_q = Eigen::Vector3d::Zero();
  Eigen::Vector3d v_q = Eigen::Vector3d::Zero();
  Eigen::Vector3d omega_t = Eigen::Vector3d::Zero();
};

/**
 * Compute the estimate-dependent vectors in Equation (10) of the ArduPilot
 * EqF paper. The IMU values are the raw inputs used by the lift; bias
 * correction is already represented by the observer's gamma component.
 */
[[nodiscard]] inline ErrorDynamicsQuantities errorDynamicsQuantities(
  const EquivariantState & observer,
  const ImuInput & input,
  const Eigen::Vector3d & gravity_global)
{
  if (!input.angular_velocity.allFinite() || !input.linear_acceleration.allFinite() ||
    !gravity_global.allFinite())
  {
    throw std::invalid_argument("EqF error-dynamics input contains NaN or infinity");
  }

  const Eigen::Matrix3d A_hat = observer.rotationBlock();
  const Eigen::Vector3d rotated_omega = A_hat * input.angular_velocity;

  ErrorDynamicsQuantities result;
  result.omega_q = rotated_omega + observer.gammaGyro();
  result.v_q =
    A_hat * input.linear_acceleration +
    observer.firstTranslationColumn().cross(rotated_omega) +
    gravity_global + observer.gammaAccel();
  result.omega_t = result.omega_q;
  return result;
}

[[nodiscard]] inline Matrix21d buildErrorDynamicsA(
  const Eigen::Vector3d & gravity_global,
  const Eigen::Vector3d & estimated_gyro_bias,
  const Eigen::Vector3d & omega_q,
  const Eigen::Vector3d & v_q,
  const Eigen::Vector3d & omega_t)
{
  if (!gravity_global.allFinite() || !estimated_gyro_bias.allFinite() ||
    !omega_q.allFinite() || !v_q.allFinite() || !omega_t.allFinite())
  {
    throw std::invalid_argument("EqF A matrix input contains NaN or infinity");
  }

  Matrix21d A = Matrix21d::Zero();
  const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

  A.block<3, 3>(ErrorState::kAttitudeOffset, ErrorState::kGammaGyroOffset) = I;
  A.block<3, 3>(ErrorState::kVelocityOffset, ErrorState::kAttitudeOffset) =
    skew_matrix(gravity_global);
  A.block<3, 3>(ErrorState::kVelocityOffset, ErrorState::kGammaAccelOffset) = I;
  A.block<3, 3>(ErrorState::kPositionOffset, ErrorState::kVelocityOffset) = I;
  A.block<3, 3>(ErrorState::kPositionOffset, ErrorState::kGammaGyroOffset) =
    skew_matrix(estimated_gyro_bias);
  A.block<3, 3>(ErrorState::kGammaGyroOffset, ErrorState::kGammaGyroOffset) =
    skew_matrix(omega_q);
  A.block<3, 3>(ErrorState::kGammaAccelOffset, ErrorState::kGammaGyroOffset) =
    skew_matrix(v_q);
  A.block<3, 3>(ErrorState::kGammaAccelOffset, ErrorState::kGammaAccelOffset) =
    skew_matrix(omega_q);
  A.block<3, 3>(ErrorState::kLeverArmOffset, ErrorState::kLeverArmOffset) =
    skew_matrix(omega_t);
  A.block<3, 3>(ErrorState::kMagRotationOffset, ErrorState::kAttitudeOffset) =
    -skew_matrix(omega_t);
  A.block<3, 3>(ErrorState::kMagRotationOffset, ErrorState::kGammaGyroOffset) = I;
  A.block<3, 3>(ErrorState::kMagRotationOffset, ErrorState::kMagRotationOffset) =
    skew_matrix(omega_t);
  return A;
}

/** Build A_t^0 from the beginning-of-step observer and physical estimate. */
[[nodiscard]] inline Matrix21d buildErrorDynamicsA(
  const EquivariantState & observer_begin,
  const State & state_begin,
  const ImuInput & input,
  const Eigen::Vector3d & gravity_global)
{
  const ErrorDynamicsQuantities quantities =
    errorDynamicsQuantities(observer_begin, input, gravity_global);
  return buildErrorDynamicsA(
    gravity_global, state_begin.gyroBiasImu(), quantities.omega_q,
    quantities.v_q, quantities.omega_t);
}

}  // namespace eqf

#endif  // ROSCOPTER_EQF_ERROR_DYNAMICS_HPP
