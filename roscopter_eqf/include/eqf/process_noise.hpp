#ifndef ROSCOPTER_EQF_PROCESS_NOISE_HPP
#define ROSCOPTER_EQF_PROCESS_NOISE_HPP

#include <Eigen/Core>

#include <stdexcept>

#include "eqf/error_state.hpp"

namespace eqf
{

using Matrix18d = Eigen::Matrix<double, 18, 18>;
using Matrix21x18d = Eigen::Matrix<double, ErrorState::kErrorStateSize, 18>;

struct ProcessNoiseIndex
{
  static constexpr Eigen::Index kGyroNoise = 0;
  static constexpr Eigen::Index kAccelNoise = 3;
  static constexpr Eigen::Index kGyroBiasRandomWalk = 6;
  static constexpr Eigen::Index kAccelBiasRandomWalk = 9;
  static constexpr Eigen::Index kGnssLeverArmRandomWalk = 12;
  static constexpr Eigen::Index kMagCalibrationRandomWalk = 15;
};

struct ProcessNoise
{
  Eigen::Vector3d gyro_noise_density = Eigen::Vector3d::Zero();
  Eigen::Vector3d accel_noise_density = Eigen::Vector3d::Zero();
  Eigen::Vector3d gyro_bias_random_walk = Eigen::Vector3d::Zero();
  Eigen::Vector3d accel_bias_random_walk = Eigen::Vector3d::Zero();
  Eigen::Vector3d gnss_lever_arm_random_walk = Eigen::Vector3d::Zero();
  Eigen::Vector3d mag_calibration_random_walk = Eigen::Vector3d::Zero();
};

/**
 * Fixed process-noise convention for this EqF implementation. The ArduPilot
 * EqF paper does not print this matrix explicitly. Sensor measurement noise
 * enters the true-minus-estimated navigation error with a negative sign;
 * calibration random walks enter their corresponding error blocks directly.
 */
[[nodiscard]] inline Matrix21x18d makeProcessNoiseMapL()
{
  Matrix21x18d L = Matrix21x18d::Zero();
  const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  L.block<3, 3>(ErrorState::kAttitudeOffset, ProcessNoiseIndex::kGyroNoise) = -I;
  L.block<3, 3>(ErrorState::kVelocityOffset, ProcessNoiseIndex::kAccelNoise) = -I;
  L.block<3, 3>(
    ErrorState::kGammaGyroOffset, ProcessNoiseIndex::kGyroBiasRandomWalk) = I;
  L.block<3, 3>(
    ErrorState::kGammaAccelOffset, ProcessNoiseIndex::kAccelBiasRandomWalk) = I;
  L.block<3, 3>(
    ErrorState::kLeverArmOffset, ProcessNoiseIndex::kGnssLeverArmRandomWalk) = I;
  L.block<3, 3>(
    ErrorState::kMagRotationOffset, ProcessNoiseIndex::kMagCalibrationRandomWalk) = I;
  return L;
}

[[nodiscard]] inline Matrix18d makeContinuousProcessNoiseCovariance(
  const ProcessNoise & noise)
{
  const Eigen::Matrix<double, 18, 1> standard_deviations =
    (Eigen::Matrix<double, 18, 1>() <<
    noise.gyro_noise_density,
    noise.accel_noise_density,
    noise.gyro_bias_random_walk,
    noise.accel_bias_random_walk,
    noise.gnss_lever_arm_random_walk,
    noise.mag_calibration_random_walk).finished();

  if (!standard_deviations.allFinite() ||
    (standard_deviations.array() < 0.0).any())
  {
    throw std::invalid_argument(
            "EqF process-noise densities must be finite and nonnegative");
  }

  return standard_deviations.array().square().matrix().asDiagonal();
}

}  // namespace eqf

#endif  // ROSCOPTER_EQF_PROCESS_NOISE_HPP
