#ifndef ROSCOPTER_EQF_ROBUST_MEASUREMENT_HPP
#define ROSCOPTER_EQF_ROBUST_MEASUREMENT_HPP

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>

#include "eqf/covariance.hpp"
#include "eqf/error_state.hpp"
#include "eqf/measurement.hpp"
#include "eqf/state.hpp"
#include "eqf/utils.hpp"

namespace eqf
{

using Matrix3x21d = Eigen::Matrix<double, 3, ErrorState::kErrorStateSize>;
using Matrix21x3d = Eigen::Matrix<double, ErrorState::kErrorStateSize, 3>;
using Vector21d = ErrorState::Vector21d;

enum class MeasurementType
{
  kMagnetometer,
  kGnssPosition,
  kGnssVelocity
};

struct RobustUpdateConfig
{
  double alpha = 0.5;
  bool robust_inflation_enabled = true;
  double innovation_covariance_jitter = 1e-9;
  double covariance_jitter = 1e-12;
  double maximum_condition_number = 1e12;
};

struct RobustUpdateDiagnostics
{
  MeasurementType type = MeasurementType::kGnssPosition;
  Eigen::Vector3d innovation = Eigen::Vector3d::Zero();
  Eigen::Matrix3d nominal_innovation_covariance = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d inflated_innovation_covariance = Eigen::Matrix3d::Zero();
  Vector21d delta_epsilon = Vector21d::Zero();
  double normalized_innovation_magnitude = 0.0;
  double post_inflation_normalized_innovation = 0.0;
  double beta = 1.0;
  double alpha = 0.5;
  bool accepted = false;
  bool numerical_regularization_used = false;
};

struct RobustUpdateResult
{
  EquivariantState observer_after;
  State state_after = State::Identity();
  Matrix21d covariance_after = Matrix21d::Identity();
  RobustUpdateDiagnostics diagnostics;
};

[[nodiscard]] inline Eigen::Matrix3d symmetrize3(const Eigen::Matrix3d & matrix)
{
  return 0.5 * (matrix + matrix.transpose()).eval();
}

[[nodiscard]] inline Matrix21d symmetrize21(const Matrix21d & matrix)
{
  return 0.5 * (matrix + matrix.transpose()).eval();
}

[[nodiscard]] inline Eigen::Vector3d predictedMagnetometerDirection(
  const State & state,
  const Eigen::Vector3d & magnetic_reference_global)
{
  if (!magnetic_reference_global.allFinite() ||
    magnetic_reference_global.norm() <= kDirectionNormTolerance)
  {
    throw std::invalid_argument("Magnetic reference must be finite and nonzero");
  }
  return state.rotationImuFromMagnetometer().transpose() *
         state.rotationGlobalFromImu().transpose() *
         magnetic_reference_global.normalized();
}

[[nodiscard]] inline Eigen::Vector3d gnssPositionInnovation(
  const State & state,
  const Eigen::Vector3d & antenna_position_global)
{
  if (!antenna_position_global.allFinite()) {
    throw std::invalid_argument("GNSS position contains NaN or infinity");
  }
  const Eigen::Vector3d predicted = state.positionGlobal() +
    state.rotationGlobalFromImu() * state.gnssLeverArmImu();
  return state.rotationGlobalFromImu().transpose() *
         (antenna_position_global - predicted);
}

[[nodiscard]] inline Eigen::Vector3d gnssVelocityInnovation(
  const State & state,
  const Eigen::Vector3d & antenna_velocity_global,
  const Eigen::Vector3d & angular_rate_imu)
{
  if (!antenna_velocity_global.allFinite() || !angular_rate_imu.allFinite()) {
    throw std::invalid_argument("GNSS velocity input contains NaN or infinity");
  }
  const Eigen::Vector3d predicted = state.velocityGlobal() +
    state.rotationGlobalFromImu() *
    angular_rate_imu.cross(state.gnssLeverArmImu());
  return state.rotationGlobalFromImu().transpose() *
         (antenna_velocity_global - predicted);
}

/** ArduPilot EqF Equation (12), remapped to ErrorState's documented ordering. */
[[nodiscard]] inline Matrix3x21d computeGnssPositionCStar(
  const EquivariantState & observer,
  const Eigen::Vector3d & antenna_position_global)
{
  if (!antenna_position_global.allFinite()) {
    throw std::invalid_argument("GNSS position contains NaN or infinity");
  }
  Matrix3x21d C = Matrix3x21d::Zero();
  const Eigen::Vector3d symmetrized_output = antenna_position_global +
    observer.secondTranslationColumn() - observer.delta();
  C.block<3, 3>(0, ErrorState::kAttitudeOffset) =
    0.5 * skew_matrix(symmetrized_output);
  C.block<3, 3>(0, ErrorState::kPositionOffset) = -Eigen::Matrix3d::Identity();
  C.block<3, 3>(0, ErrorState::kLeverArmOffset) = Eigen::Matrix3d::Identity();
  return C;
}

/** ArduPilot EqF Equation (13), remapped to ErrorState's documented ordering. */
[[nodiscard]] inline Matrix3x21d computeGnssVelocityCStar(
  const EquivariantState & observer,
  const Eigen::Vector3d & antenna_velocity_global,
  const Eigen::Vector3d & angular_rate_imu)
{
  if (!antenna_velocity_global.allFinite() || !angular_rate_imu.allFinite()) {
    throw std::invalid_argument("GNSS velocity input contains NaN or infinity");
  }
  Matrix3x21d C = Matrix3x21d::Zero();
  const Eigen::Vector3d symmetrized_output = antenna_velocity_global +
    observer.firstTranslationColumn() - angular_rate_imu.cross(observer.delta());
  C.block<3, 3>(0, ErrorState::kAttitudeOffset) =
    0.5 * skew_matrix(symmetrized_output);
  C.block<3, 3>(0, ErrorState::kVelocityOffset) = -Eigen::Matrix3d::Identity();
  C.block<3, 3>(0, ErrorState::kLeverArmOffset) = skew_matrix(angular_rate_imu);
  return C;
}

/**
 * ArduPilot EqF Equation (11), expressed in the established tangent residual
 * coordinates m_ref x rho_m(X^-1, y_m).  The paper's calibration block is
 * placed at ErrorState::kMagRotationOffset rather than relying on typeset
 * block position.
 */
[[nodiscard]] inline Matrix3x21d computeMagnetometerCStar(
  const EquivariantState & observer,
  const Eigen::Vector3d & predicted_direction_magnetometer,
  const Eigen::Vector3d & magnetic_reference_global)
{
  if (!predicted_direction_magnetometer.allFinite() ||
    !magnetic_reference_global.allFinite() ||
    predicted_direction_magnetometer.norm() <= kDirectionNormTolerance ||
    magnetic_reference_global.norm() <= kDirectionNormTolerance)
  {
    throw std::invalid_argument("Magnetometer C-star requires finite nonzero directions");
  }
  const Eigen::Vector3d reference = magnetic_reference_global.normalized();
  const Eigen::Vector3d transported_prediction =
    observer.E() * predicted_direction_magnetometer.normalized();
  Matrix3x21d C = Matrix3x21d::Zero();
  C.block<3, 3>(0, ErrorState::kMagRotationOffset) =
    skew_matrix(reference) * 0.5 *
    skew_matrix(reference + transported_prediction);
  return C;
}

/**
 * Implements (D phi|I)^dagger (D vartheta^-1|0) for the existing right-local
 * state coordinates.  This is deliberately not an identity map.
 */
[[nodiscard]] inline Lift mapNormalCorrectionToLieAlgebra(
  const Vector21d & delta_epsilon,
  const State & origin)
{
  if (!delta_epsilon.allFinite()) {
    throw std::invalid_argument("EqF correction contains NaN or infinity");
  }
  const Eigen::Vector3d theta =
    delta_epsilon.segment<3>(ErrorState::kAttitudeOffset);
  const Eigen::Vector3d velocity =
    delta_epsilon.segment<3>(ErrorState::kVelocityOffset);
  const Eigen::Vector3d position =
    delta_epsilon.segment<3>(ErrorState::kPositionOffset);
  const Eigen::Vector3d gyro_bias =
    delta_epsilon.segment<3>(ErrorState::kGammaGyroOffset);
  const Eigen::Vector3d accel_bias =
    delta_epsilon.segment<3>(ErrorState::kGammaAccelOffset);
  const Eigen::Vector3d lever_arm =
    delta_epsilon.segment<3>(ErrorState::kLeverArmOffset);
  const Eigen::Vector3d mag_rotation =
    delta_epsilon.segment<3>(ErrorState::kMagRotationOffset);

  Lift correction;
  correction.lambda1 = se23Wedge(theta, velocity, position);
  correction.lambda2.head<3>() =
    origin.gyroBiasImu().cross(theta) - gyro_bias;
  correction.lambda2.tail<3>() =
    origin.accelBiasImu().cross(theta) +
    origin.gyroBiasImu().cross(velocity) - accel_bias;
  correction.lambda3 = origin.gnssLeverArmImu().cross(theta) + lever_arm;
  correction.lambda4 =
    origin.rotationImuFromMagnetometer().transpose() * theta - mag_rotation;
  return correction;
}

[[nodiscard]] inline double generalizedCovarianceUnionFactor(const double r)
{
  if (!std::isfinite(r) || r < 0.0) {
    throw std::invalid_argument("Normalized innovation magnitude must be finite and nonnegative");
  }
  if (r >= 1.0) {
    return 2.0;
  }
  const double root_r = std::sqrt(r);
  return ((1.0 + root_r) * (1.0 + root_r)) / (1.0 + r);
}

namespace detail
{

struct Factorization3
{
  Eigen::Matrix3d matrix = Eigen::Matrix3d::Zero();
  Eigen::LDLT<Eigen::Matrix3d> ldlt;
  bool valid = false;
  bool regularized = false;
};

[[nodiscard]] inline Factorization3 factorInnovationCovariance(
  const Eigen::Matrix3d & input,
  const RobustUpdateConfig & config)
{
  Factorization3 result;
  result.matrix = symmetrize3(input);
  if (!result.matrix.allFinite()) {
    return result;
  }

  for (int attempt = 0; attempt < 2; ++attempt) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(result.matrix);
    if (eigen_solver.info() == Eigen::Success) {
      const Eigen::Vector3d eigenvalues = eigen_solver.eigenvalues();
      const double min_eigenvalue = eigenvalues.minCoeff();
      const double max_eigenvalue = eigenvalues.maxCoeff();
      const double condition = min_eigenvalue > 0.0 ?
        max_eigenvalue / min_eigenvalue : std::numeric_limits<double>::infinity();
      result.ldlt.compute(result.matrix);
      if (result.ldlt.info() == Eigen::Success && min_eigenvalue > 0.0 &&
        std::isfinite(condition) && condition <= config.maximum_condition_number)
      {
        result.valid = true;
        return result;
      }
    }
    if (attempt == 0) {
      const double scale = std::max(1.0, std::abs(result.matrix.trace()) / 3.0);
      result.matrix += config.innovation_covariance_jitter * scale *
        Eigen::Matrix3d::Identity();
      result.regularized = true;
    }
  }
  return result;
}

[[nodiscard]] inline bool validConfig(const RobustUpdateConfig & config)
{
  return std::isfinite(config.alpha) && config.alpha >= 0.5 && config.alpha <= 1.0 &&
         std::isfinite(config.innovation_covariance_jitter) &&
         config.innovation_covariance_jitter > 0.0 &&
         std::isfinite(config.covariance_jitter) && config.covariance_jitter >= 0.0 &&
         std::isfinite(config.maximum_condition_number) &&
         config.maximum_condition_number > 1.0;
}

}  // namespace detail

/**
 * Robust EqF correction using ArduPilot EqF Section VI and a Joseph-form
 * discrete covariance update.  Failures return the prior atomically.
 */
[[nodiscard]] inline RobustUpdateResult applyRobustEqfCorrection(
  const MeasurementType type,
  const Eigen::Vector3d & innovation,
  const Matrix3x21d & C_star,
  const Eigen::Matrix3d & measurement_covariance,
  const RobustUpdateConfig & config,
  const State & origin,
  const EquivariantState & observer_prior,
  const Matrix21d & covariance_prior)
{
  RobustUpdateResult result;
  result.observer_after = observer_prior;
  result.state_after = phi(observer_prior, origin);
  result.covariance_after = covariance_prior;
  result.diagnostics.type = type;
  result.diagnostics.innovation = innovation;
  result.diagnostics.alpha = config.alpha;

  if (!detail::validConfig(config) || !innovation.allFinite() ||
    !C_star.allFinite() || !measurement_covariance.allFinite() ||
    !covariance_prior.allFinite())
  {
    return result;
  }

  const Eigen::Matrix3d R = symmetrize3(measurement_covariance);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> R_eigen(R);
  Eigen::SelfAdjointEigenSolver<Matrix21d> P_eigen(symmetrize21(covariance_prior));
  if (R_eigen.info() != Eigen::Success || P_eigen.info() != Eigen::Success ||
    R_eigen.eigenvalues().minCoeff() < -config.covariance_jitter ||
    P_eigen.eigenvalues().minCoeff() < -config.covariance_jitter)
  {
    return result;
  }

  const Matrix21d P = symmetrize21(covariance_prior);
  const Eigen::Matrix3d Pyy = symmetrize3(C_star * P * C_star.transpose());
  const detail::Factorization3 nominal =
    detail::factorInnovationCovariance(Pyy + R, config);
  result.diagnostics.nominal_innovation_covariance = nominal.matrix;
  result.diagnostics.numerical_regularization_used = nominal.regularized;
  if (!nominal.valid) {
    return result;
  }

  const Eigen::Vector3d nominal_solution = nominal.ldlt.solve(innovation);
  if (nominal.ldlt.info() != Eigen::Success || !nominal_solution.allFinite()) {
    return result;
  }
  const double r = std::max(0.0, innovation.dot(nominal_solution));
  const double beta = config.robust_inflation_enabled ?
    generalizedCovarianceUnionFactor(r) : 1.0;
  const double alpha_term = config.robust_inflation_enabled ? config.alpha : 0.0;
  result.diagnostics.normalized_innovation_magnitude = r;
  result.diagnostics.beta = beta;

  const Eigen::Matrix3d innovation_outer = innovation * innovation.transpose();
  const detail::Factorization3 inflated = detail::factorInnovationCovariance(
    beta * (Pyy + alpha_term * innovation_outer) + R, config);
  result.diagnostics.inflated_innovation_covariance = inflated.matrix;
  result.diagnostics.numerical_regularization_used =
    result.diagnostics.numerical_regularization_used || inflated.regularized;
  if (!inflated.valid) {
    return result;
  }

  const Eigen::Matrix<double, 3, ErrorState::kErrorStateSize> gain_transpose =
    inflated.ldlt.solve(C_star * P);
  const Eigen::Vector3d inflated_solution = inflated.ldlt.solve(innovation);
  if (inflated.ldlt.info() != Eigen::Success || !gain_transpose.allFinite() ||
    !inflated_solution.allFinite())
  {
    return result;
  }
  const Matrix21x3d gain = gain_transpose.transpose();
  const Vector21d delta_epsilon = gain * innovation;
  result.diagnostics.delta_epsilon = delta_epsilon;
  result.diagnostics.post_inflation_normalized_innovation =
    std::max(0.0, innovation.dot(inflated_solution));

  try {
    const Lift correction = mapNormalCorrectionToLieAlgebra(delta_epsilon, origin);
    // C-star and the concrete innovations describe truth-minus-estimate error,
    // while this implementation stores the observer through the established
    // right state action.  Negating exactly here gives the correction side/sign
    // selected by the small-error tests.
    const EquivariantState observer_after =
      groupMultiply(expGroup(scaleLift(correction, -1.0)), observer_prior);
    const State state_after = phi(observer_after, origin);

    const Eigen::Matrix3d R_effective = symmetrize3(
      (beta - 1.0) * Pyy + beta * alpha_term * innovation_outer + R);
    const Matrix21d joseph_factor =
      Matrix21d::Identity() - gain * C_star;
    Matrix21d covariance_after = symmetrize21(
      joseph_factor * P * joseph_factor.transpose() +
      gain * R_effective * gain.transpose());
    Eigen::SelfAdjointEigenSolver<Matrix21d> covariance_eigen(covariance_after);
    if (covariance_eigen.info() != Eigen::Success) {
      return result;
    }
    const double minimum_eigenvalue = covariance_eigen.eigenvalues().minCoeff();
    if (minimum_eigenvalue < -config.covariance_jitter) {
      return result;
    }
    if (minimum_eigenvalue < 0.0) {
      covariance_after += (config.covariance_jitter - minimum_eigenvalue) *
        Matrix21d::Identity();
      covariance_after = symmetrize21(covariance_after);
      result.diagnostics.numerical_regularization_used = true;
    }

    result.observer_after = observer_after;
    result.state_after = state_after;
    result.covariance_after = covariance_after;
    result.diagnostics.accepted = true;
  } catch (const std::exception &) {
    return result;
  }
  return result;
}

}  // namespace eqf

#endif  // ROSCOPTER_EQF_ROBUST_MEASUREMENT_HPP
