#ifndef ROSCOPTER_EQF_EQUIVARIANT_STATE_HPP
#define ROSCOPTER_EQF_EQUIVARIANT_STATE_HPP

#include <Eigen/Core>

#include <cmath>
#include <stdexcept>
#include <string>

namespace eqf
{

class EqfEquivariantState
{
public:
  using SE23Matrix = Eigen::Matrix<double, 5, 5>;
  using Se3AlgebraMatrix = Eigen::Matrix<double, 4, 4>;

  EqfEquivariantState()
  {
    reset();
  }

  EqfEquivariantState(
    const SE23Matrix & C,
    const Se3AlgebraMatrix & gamma,
    const Eigen::Vector3d & delta,
    const Eigen::Matrix3d & E)
  {
    setC(C);
    setGamma(gamma);
    setDelta(delta);
    setE(E);
  }

  void reset()
  {
    C_.setIdentity();
    gamma_.setZero();
    delta_.setZero();
    E_.setIdentity();
  }

  const SE23Matrix & C() const
  {
    return C_;
  }

  const Se3AlgebraMatrix & gamma() const
  {
    return gamma_;
  }

  const Eigen::Vector3d & delta() const
  {
    return delta_;
  }

  const Eigen::Matrix3d & E() const
  {
    return E_;
  }

  void setC(const SE23Matrix & C)
  {
    validateSE23OrThrow(C);
    C_ = C;
  }

  void setGamma(const Se3AlgebraMatrix & gamma)
  {
    validateSe3AlgebraOrThrow(gamma);
    gamma_ = gamma;
  }

  void setDelta(const Eigen::Vector3d & delta)
  {
    validateVectorOrThrow(delta, "delta");
    delta_ = delta;
  }

  void setE(const Eigen::Matrix3d & E)
  {
    validateRotationOrThrow(E);
    E_ = E;
  }

  Eigen::Matrix3d rotationBlock() const
  {
    return C_.block<3, 3>(0, 0);
  }

  Eigen::Vector3d firstTranslationColumn() const
  {
    return C_.block<3, 1>(0, 3);
  }

  Eigen::Vector3d secondTranslationColumn() const
  {
    return C_.block<3, 1>(0, 4);
  }

  void setRotationBlock(const Eigen::Matrix3d & A)
  {
    validateRotationOrThrow(A);
    C_.block<3, 3>(0, 0) = A;
  }

  void setFirstTranslationColumn(const Eigen::Vector3d & a)
  {
    validateVectorOrThrow(a, "firstTranslationColumn");
    C_.block<3, 1>(0, 3) = a;
  }

  void setSecondTranslationColumn(const Eigen::Vector3d & b)
  {
    validateVectorOrThrow(b, "secondTranslationColumn");
    C_.block<3, 1>(0, 4) = b;
  }

  Eigen::Vector3d gammaGyro() const
  {
    return vee(gamma_.block<3, 3>(0, 0));
  }

  Eigen::Vector3d gammaAccel() const
  {
    return gamma_.block<3, 1>(0, 3);
  }

  void setGammaGyro(const Eigen::Vector3d & gamma_gyro)
  {
    validateVectorOrThrow(gamma_gyro, "gamma_gyro");
    gamma_.block<3, 3>(0, 0) = hat(gamma_gyro);
    gamma_.block<1, 4>(3, 0).setZero();
  }

  void setGammaAccel(const Eigen::Vector3d & gamma_accel)
  {
    validateVectorOrThrow(gamma_accel, "gamma_accel");
    gamma_.block<3, 1>(0, 3) = gamma_accel;
    gamma_.block<1, 4>(3, 0).setZero();
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
    const SE23Matrix & C,
    const double tolerance = 1e-6)
  {
    if (!C.allFinite()) {
      return false;
    }

    Eigen::Matrix<double, 1, 5> bottom_row_0;
    bottom_row_0 << 0.0, 0.0, 0.0, 1.0, 0.0;

    Eigen::Matrix<double, 1, 5> bottom_row_1;
    bottom_row_1 << 0.0, 0.0, 0.0, 0.0, 1.0;

    return isValidRotation(C.block<3, 3>(0, 0), tolerance) &&
           (C.block<1, 5>(3, 0) - bottom_row_0).norm() <= tolerance &&
           (C.block<1, 5>(4, 0) - bottom_row_1).norm() <= tolerance;
  }

  static bool isValidSe3Algebra(
    const Se3AlgebraMatrix & gamma,
    const double tolerance = 1e-6)
  {
    if (!gamma.allFinite()) {
      return false;
    }

    const Eigen::Matrix3d omega = gamma.block<3, 3>(0, 0);
    const Eigen::Matrix<double, 1, 4> bottom_row = gamma.block<1, 4>(3, 0);
    return (omega + omega.transpose()).norm() <= tolerance &&
           bottom_row.norm() <= tolerance;
  }

private:
  static Eigen::Matrix3d hat(const Eigen::Vector3d & vector)
  {
    Eigen::Matrix3d result;
    result << 0.0, -vector.z(), vector.y(),
              vector.z(), 0.0, -vector.x(),
              -vector.y(), vector.x(), 0.0;
    return result;
  }

  static Eigen::Vector3d vee(const Eigen::Matrix3d & skew)
  {
    return Eigen::Vector3d(skew(2, 1), skew(0, 2), skew(1, 0));
  }

  static void validateRotationOrThrow(const Eigen::Matrix3d & R)
  {
    if (!isValidRotation(R)) {
      throw std::invalid_argument("Invalid rotation matrix");
    }
  }

  static void validateSE23OrThrow(const SE23Matrix & C)
  {
    if (!isValidSE23(C)) {
      throw std::invalid_argument("Invalid SE_2(3) matrix");
    }
  }

  static void validateSe3AlgebraOrThrow(const Se3AlgebraMatrix & gamma)
  {
    if (!isValidSe3Algebra(gamma)) {
      throw std::invalid_argument("Invalid se(3) algebra matrix");
    }
  }

  static void validateVectorOrThrow(const Eigen::Vector3d & vector, const char * name)
  {
    if (!vector.allFinite()) {
      throw std::invalid_argument(std::string(name) + " contains NaN or infinity");
    }
  }

  // Equivariant observer state X = (C, gamma, delta, E).
  // C is stored as one SE_2(3) matrix with fixed bottom rows.
  // gamma is stored as one se(3) algebra matrix with zero bottom row.
  SE23Matrix C_;
  Se3AlgebraMatrix gamma_;
  Eigen::Vector3d delta_;
  Eigen::Matrix3d E_;
};

using EquivariantState = EqfEquivariantState;

}  // namespace eqf

#endif  // ROSCOPTER_EQF_EQUIVARIANT_STATE_HPP
