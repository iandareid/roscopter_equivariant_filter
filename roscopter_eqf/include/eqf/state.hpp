#ifndef ROSCOPTER_EQF_STATE_HPP
#define ROSCOPTER_EQF_STATE_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace eqf
{

class EqfState
{
public:
  using SE23Matrix = Eigen::Matrix<double, 5, 5>;

  EqfState()
  {
    reset();
  }

  EqfState(
    const SE23Matrix & T,
    const Eigen::Vector3d & gyro_bias_I,
    const Eigen::Vector3d & accel_bias_I,
    const Eigen::Vector3d & gnss_lever_arm_I,
    const Eigen::Matrix3d & R_IM)
  {
    setNavigationState(T);
    setGyroBiasImu(gyro_bias_I);
    setAccelBiasImu(accel_bias_I);
    setGnssLeverArmImu(gnss_lever_arm_I);
    setRotationImuFromMagnetometer(R_IM);
  }

  void reset()
  {
    T_.setIdentity();
    gyro_bias_I_.setZero();
    accel_bias_I_.setZero();
    gnss_lever_arm_I_.setZero();
    R_IM_.setIdentity();
  }
  
  static EqfState Identity()
  {
    EqfState state;
    state.reset();

    return state;
  }

  const SE23Matrix & navigationState() const
  {
    return T_;
  }

  Eigen::Matrix3d rotationGlobalFromImu() const
  {
    return T_.block<3, 3>(0, 0);
  }

  Eigen::Vector3d velocityGlobal() const
  {
    return T_.block<3, 1>(0, 3);
  }

  Eigen::Vector3d positionGlobal() const
  {
    return T_.block<3, 1>(0, 4);
  }

  Eigen::Quaterniond quaternionGlobalFromImu() const
  {
    return Eigen::Quaterniond(rotationGlobalFromImu()).normalized();
  }

  Eigen::Vector3d eulerAngles() const
  {
    const Eigen::Matrix3d R = rotationGlobalFromImu();
    const double phi = std::atan2(R(2, 1), R(2, 2));
    const double theta = std::asin(std::clamp(-R(2, 0), -1.0, 1.0));
    const double psi = std::atan2(R(1, 0), R(0, 0));
    return Eigen::Vector3d(phi, theta, psi);
  }

  const Eigen::Vector3d & gyroBiasImu() const
  {
    return gyro_bias_I_;
  }

  const Eigen::Vector3d & accelBiasImu() const
  {
    return accel_bias_I_;
  }

  const Eigen::Vector3d & gnssLeverArmImu() const
  {
    return gnss_lever_arm_I_;
  }

  const Eigen::Matrix3d & rotationImuFromMagnetometer() const
  {
    return R_IM_;
  }

  double pn() const {return positionGlobal().x();}
  double pe() const {return positionGlobal().y();}
  double pd() const {return positionGlobal().z();}

  double v_x() const {return velocityGlobal().x();}
  double v_y() const {return velocityGlobal().y();}
  double v_z() const {return velocityGlobal().z();}

  double b_gx() const {return gyro_bias_I_.x();}
  double b_gy() const {return gyro_bias_I_.y();}
  double b_gz() const {return gyro_bias_I_.z();}

  double b_ax() const {return accel_bias_I_.x();}
  double b_ay() const {return accel_bias_I_.y();}
  double b_az() const {return accel_bias_I_.z();}

  void setNavigationState(const SE23Matrix & T)
  {
    validateSE23OrThrow(T);
    T_ = T;
  }

  void setRotationGlobalFromImu(const Eigen::Matrix3d & R_GI)
  {
    validateRotationOrThrow(R_GI);
    T_.block<3, 3>(0, 0) = R_GI;
  }

  void setVelocityGlobal(const Eigen::Vector3d & v_GI)
  {
    validateVectorOrThrow(v_GI, "velocityGlobal");
    T_.block<3, 1>(0, 3) = v_GI;
  }

  void setPositionGlobal(const Eigen::Vector3d & p_GI)
  {
    validateVectorOrThrow(p_GI, "positionGlobal");
    T_.block<3, 1>(0, 4) = p_GI;
  }

  void setGyroBiasImu(const Eigen::Vector3d & gyro_bias_I)
  {
    validateVectorOrThrow(gyro_bias_I, "gyroBiasImu");
    gyro_bias_I_ = gyro_bias_I;
  }

  void setAccelBiasImu(const Eigen::Vector3d & accel_bias_I)
  {
    validateVectorOrThrow(accel_bias_I, "accelBiasImu");
    accel_bias_I_ = accel_bias_I;
  }

  void setGnssLeverArmImu(const Eigen::Vector3d & gnss_lever_arm_I)
  {
    validateVectorOrThrow(gnss_lever_arm_I, "gnssLeverArmImu");
    gnss_lever_arm_I_ = gnss_lever_arm_I;
  }

  void setRotationImuFromMagnetometer(const Eigen::Matrix3d & R_IM)
  {
    validateRotationOrThrow(R_IM);
    R_IM_ = R_IM;
  }

  Eigen::Vector3d correctedAngularRate(
    const Eigen::Vector3d & measured_angular_rate_I) const
  {
    validateVectorOrThrow(measured_angular_rate_I, "measured_angular_rate_I");
    return measured_angular_rate_I - gyro_bias_I_;
  }

  Eigen::Vector3d correctedSpecificForce(
    const Eigen::Vector3d & measured_specific_force_I) const
  {
    validateVectorOrThrow(measured_specific_force_I, "measured_specific_force_I");
    return measured_specific_force_I - accel_bias_I_;
  }

  static bool isValidRotation(
    const Eigen::Matrix3d & R,
    const double tolerance = 1e-6)
  {
    if (!R.allFinite()) {
      return false;
    }

    const Eigen::Matrix3d orthogonality_error =
      R.transpose() * R - Eigen::Matrix3d::Identity();
    return orthogonality_error.norm() <= tolerance &&
           std::abs(R.determinant() - 1.0) <= tolerance;
  }

  static bool isValidSE23(
    const SE23Matrix & T,
    const double tolerance = 1e-6)
  {
    if (!T.allFinite()) {
      return false;
    }

    Eigen::Matrix<double, 1, 5> bottom_row_0;
    bottom_row_0 << 0.0, 0.0, 0.0, 1.0, 0.0;

    Eigen::Matrix<double, 1, 5> bottom_row_1;
    bottom_row_1 << 0.0, 0.0, 0.0, 0.0, 1.0;

    return isValidRotation(T.block<3, 3>(0, 0), tolerance) &&
           (T.block<1, 5>(3, 0) - bottom_row_0).norm() <= tolerance &&
           (T.block<1, 5>(4, 0) - bottom_row_1).norm() <= tolerance;
  }

private:
  static void validateRotationOrThrow(const Eigen::Matrix3d & R)
  {
    if (!isValidRotation(R)) {
      throw std::invalid_argument("Invalid rotation matrix");
    }
  }

  static void validateSE23OrThrow(const SE23Matrix & T)
  {
    if (!isValidSE23(T)) {
      throw std::invalid_argument("Invalid SE_2(3) navigation state matrix");
    }
  }

  static void validateVectorOrThrow(const Eigen::Vector3d & vector, const char * name)
  {
    if (!vector.allFinite()) {
      throw std::invalid_argument(std::string(name) + " contains NaN or infinity");
    }
  }

  // SE_2(3) layout:
  // T_.block<3, 3>(0, 0) = global-from-IMU rotation R_GI
  // T_.block<3, 1>(0, 3) = global-frame velocity v_GI
  // T_.block<3, 1>(0, 4) = global-frame position p_GI
  // Bottom rows are fixed to [0 0 0 1 0] and [0 0 0 0 1].
  SE23Matrix T_;

  Eigen::Vector3d gyro_bias_I_;
  Eigen::Vector3d accel_bias_I_;

  Eigen::Vector3d gnss_lever_arm_I_;
  Eigen::Matrix3d R_IM_;
};

using State = EqfState;

}  // namespace eqf

#endif  // ROSCOPTER_EQF_STATE_HPP
