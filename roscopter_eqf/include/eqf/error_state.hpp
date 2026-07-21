#ifndef ROSCOPTER_EQF_ERROR_STATE_HPP
#define ROSCOPTER_EQF_ERROR_STATE_HPP

#include <Eigen/Core>

#include <stdexcept>
#include <string>

namespace eqf
{

class EqfErrorState
{
public:
  using Vector21d = Eigen::Matrix<double, 21, 1>;
  using Vector9d = Eigen::Matrix<double, 9, 1>;
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  using Se23AlgebraMatrix = Eigen::Matrix<double, 5, 5>;
  using Se3AlgebraMatrix = Eigen::Matrix<double, 4, 4>;

  static constexpr Eigen::Index kAttitudeOffset = 0;
  static constexpr Eigen::Index kVelocityOffset = 3;
  static constexpr Eigen::Index kPositionOffset = 6;
  static constexpr Eigen::Index kGammaGyroOffset = 9;
  static constexpr Eigen::Index kGammaAccelOffset = 12;
  static constexpr Eigen::Index kLeverArmOffset = 15;
  static constexpr Eigen::Index kMagRotationOffset = 18;

  static constexpr Eigen::Index kVector3Size = 3;
  static constexpr Eigen::Index kErrorStateSize = 21;

  EqfErrorState()
  {
    reset();
  }

  explicit EqfErrorState(const Vector21d & epsilon)
  {
    setVector(epsilon);
  }

  EqfErrorState(
    const Eigen::Vector3d & delta_theta,
    const Eigen::Vector3d & delta_velocity,
    const Eigen::Vector3d & delta_position,
    const Eigen::Vector3d & delta_gamma_gyro,
    const Eigen::Vector3d & delta_gamma_accel,
    const Eigen::Vector3d & delta_lever_arm,
    const Eigen::Vector3d & delta_mag_rotation)
  {
    reset();
    setAttitudeError(delta_theta);
    setVelocityError(delta_velocity);
    setPositionError(delta_position);
    setGammaGyroError(delta_gamma_gyro);
    setGammaAccelError(delta_gamma_accel);
    setGnssLeverArmError(delta_lever_arm);
    setMagnetometerRotationError(delta_mag_rotation);
  }

  void reset()
  {
    epsilon_.setZero();
  }

  bool isZero(const double tolerance = 1e-12) const
  {
    return epsilon_.norm() <= tolerance;
  }

  const Vector21d & vector() const
  {
    return epsilon_;
  }

  void setVector(const Vector21d & epsilon)
  {
    validateFiniteOrThrow(epsilon);
    epsilon_ = epsilon;
  }

  Eigen::Vector3d attitudeError() const
  {
    return epsilon_.segment<kVector3Size>(kAttitudeOffset);
  }

  Eigen::Vector3d velocityError() const
  {
    return epsilon_.segment<kVector3Size>(kVelocityOffset);
  }

  Eigen::Vector3d positionError() const
  {
    return epsilon_.segment<kVector3Size>(kPositionOffset);
  }

  Eigen::Vector3d gammaGyroError() const
  {
    return epsilon_.segment<kVector3Size>(kGammaGyroOffset);
  }

  Eigen::Vector3d gammaAccelError() const
  {
    return epsilon_.segment<kVector3Size>(kGammaAccelOffset);
  }

  Eigen::Vector3d gnssLeverArmError() const
  {
    return epsilon_.segment<kVector3Size>(kLeverArmOffset);
  }

  Eigen::Vector3d magnetometerRotationError() const
  {
    return epsilon_.segment<kVector3Size>(kMagRotationOffset);
  }

  void setAttitudeError(const Eigen::Vector3d & delta_theta)
  {
    setSegment(kAttitudeOffset, delta_theta, "delta_theta");
  }

  void setVelocityError(const Eigen::Vector3d & delta_velocity)
  {
    setSegment(kVelocityOffset, delta_velocity, "delta_velocity");
  }

  void setPositionError(const Eigen::Vector3d & delta_position)
  {
    setSegment(kPositionOffset, delta_position, "delta_position");
  }

  void setGammaGyroError(const Eigen::Vector3d & delta_gamma_gyro)
  {
    setSegment(kGammaGyroOffset, delta_gamma_gyro, "delta_gamma_gyro");
  }

  void setGammaAccelError(const Eigen::Vector3d & delta_gamma_accel)
  {
    setSegment(kGammaAccelOffset, delta_gamma_accel, "delta_gamma_accel");
  }

  void setGnssLeverArmError(const Eigen::Vector3d & delta_lever_arm)
  {
    setSegment(kLeverArmOffset, delta_lever_arm, "delta_lever_arm");
  }

  void setMagnetometerRotationError(const Eigen::Vector3d & delta_mag_rotation)
  {
    setSegment(kMagRotationOffset, delta_mag_rotation, "delta_mag_rotation");
  }

  Vector9d navigationError() const
  {
    return epsilon_.segment<9>(kAttitudeOffset);
  }

  Vector6d biasGroupError() const
  {
    return epsilon_.segment<6>(kGammaGyroOffset);
  }

  Eigen::Vector3d calibrationTranslationError() const
  {
    return gnssLeverArmError();
  }

  Eigen::Vector3d calibrationRotationError() const
  {
    return magnetometerRotationError();
  }

  void setNavigationError(const Vector9d & epsilon_C)
  {
    validateFiniteOrThrow(epsilon_C, "epsilon_C");
    epsilon_.segment<9>(kAttitudeOffset) = epsilon_C;
  }

  void setBiasGroupError(const Vector6d & epsilon_gamma)
  {
    validateFiniteOrThrow(epsilon_gamma, "epsilon_gamma");
    epsilon_.segment<6>(kGammaGyroOffset) = epsilon_gamma;
  }

  void setCalibrationTranslationError(const Eigen::Vector3d & epsilon_t)
  {
    setGnssLeverArmError(epsilon_t);
  }

  void setCalibrationRotationError(const Eigen::Vector3d & epsilon_M)
  {
    setMagnetometerRotationError(epsilon_M);
  }

  Se23AlgebraMatrix navigationAlgebraMatrix() const
  {
    Se23AlgebraMatrix result = Se23AlgebraMatrix::Zero();
    result.block<3, 3>(0, 0) = hat(attitudeError());
    result.block<3, 1>(0, 3) = velocityError();
    result.block<3, 1>(0, 4) = positionError();
    return result;
  }

  Se3AlgebraMatrix biasAlgebraMatrix() const
  {
    Se3AlgebraMatrix result = Se3AlgebraMatrix::Zero();
    result.block<3, 3>(0, 0) = hat(gammaGyroError());
    result.block<3, 1>(0, 3) = gammaAccelError();
    return result;
  }

  EqfErrorState & operator+=(const EqfErrorState & other)
  {
    epsilon_ += other.epsilon_;
    return *this;
  }

  EqfErrorState & operator-=(const EqfErrorState & other)
  {
    epsilon_ -= other.epsilon_;
    return *this;
  }

  EqfErrorState & operator*=(const double scalar)
  {
    epsilon_ *= scalar;
    return *this;
  }

  static bool isFinite(const Vector21d & epsilon)
  {
    return epsilon.array().isFinite().all();
  }

  static Eigen::Matrix3d hat(const Eigen::Vector3d & vector)
  {
    Eigen::Matrix3d result;
    result << 0.0, -vector.z(), vector.y(),
      vector.z(), 0.0, -vector.x(),
      -vector.y(), vector.x(), 0.0;
    return result;
  }

private:
  template<typename Derived>
  static void validateFiniteOrThrow(
    const Eigen::MatrixBase<Derived> & vector,
    const char * name = "epsilon")
  {
    if (!vector.array().isFinite().all()) {
      throw std::invalid_argument(std::string(name) + " contains NaN or infinity");
    }
  }

  void setSegment(
    const Eigen::Index offset,
    const Eigen::Vector3d & value,
    const char * name)
  {
    validateFiniteOrThrow(value, name);
    epsilon_.segment<kVector3Size>(offset) = value;
  }

  // Local error-coordinate layout:
  // 0:2   delta_theta
  // 3:5   delta_velocity
  // 6:8   delta_position
  // 9:11  delta_gamma_gyro
  // 12:14 delta_gamma_accel
  // 15:17 delta_lever_arm
  // 18:20 delta_mag_rotation
  Vector21d epsilon_;
};

inline EqfErrorState operator+(EqfErrorState lhs, const EqfErrorState & rhs)
{
  lhs += rhs;
  return lhs;
}

inline EqfErrorState operator-(EqfErrorState lhs, const EqfErrorState & rhs)
{
  lhs -= rhs;
  return lhs;
}

inline EqfErrorState operator*(EqfErrorState state, const double scalar)
{
  state *= scalar;
  return state;
}

inline EqfErrorState operator*(const double scalar, EqfErrorState state)
{
  state *= scalar;
  return state;
}

using ErrorState = EqfErrorState;

}  // namespace eqf

#endif  // ROSCOPTER_EQF_ERROR_STATE_HPP
