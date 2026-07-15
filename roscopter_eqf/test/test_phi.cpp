#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "eqf/equivariant_state.hpp"
#include "eqf/state.hpp"
#include "eqf/utils.hpp"

namespace
{

constexpr double kTolerance = 1e-10;

using SE23Matrix = eqf::State::SE23Matrix;

Eigen::Matrix3d rotationFromVector(const Eigen::Vector3d & vector)
{
  return Eigen::AngleAxisd(vector.norm(), vector.normalized()).toRotationMatrix();
}

SE23Matrix makeSE23(
  const Eigen::Matrix3d & R,
  const Eigen::Vector3d & v,
  const Eigen::Vector3d & p)
{
  SE23Matrix T = SE23Matrix::Identity();
  T.block<3, 3>(0, 0) = R;
  T.block<3, 1>(0, 3) = v;
  T.block<3, 1>(0, 4) = p;
  return T;
}

eqf::State makeState()
{
  return eqf::State(
    makeSE23(
      rotationFromVector(Eigen::Vector3d(0.2, -0.4, 0.3)),
      Eigen::Vector3d(1.2, -0.5, 0.7),
      Eigen::Vector3d(-2.0, 0.8, 1.5)),
    Eigen::Vector3d(0.03, -0.02, 0.05),
    Eigen::Vector3d(-0.4, 0.6, 0.2),
    Eigen::Vector3d(0.7, -0.3, 0.4),
    rotationFromVector(Eigen::Vector3d(-0.3, 0.1, 0.25)));
}

eqf::EquivariantState makeSymmetry()
{
  eqf::EquivariantState X;
  X.setC(
    makeSE23(
      rotationFromVector(Eigen::Vector3d(-0.15, 0.25, 0.35)),
      Eigen::Vector3d(0.4, -0.2, 0.9),
      Eigen::Vector3d(-0.6, 0.1, 0.3)));
  X.setGammaGyro(Eigen::Vector3d(0.01, -0.04, 0.02));
  X.setGammaAccel(Eigen::Vector3d(-0.3, 0.2, 0.5));
  X.setDelta(Eigen::Vector3d(0.8, -0.7, 0.6));
  X.setE(rotationFromVector(Eigen::Vector3d(0.12, -0.22, 0.18)));
  return X;
}

void expectStatesApprox(const eqf::State & lhs, const eqf::State & rhs)
{
  EXPECT_TRUE(lhs.navigationState().isApprox(rhs.navigationState(), kTolerance));
  EXPECT_TRUE(lhs.gyroBiasImu().isApprox(rhs.gyroBiasImu(), kTolerance));
  EXPECT_TRUE(lhs.accelBiasImu().isApprox(rhs.accelBiasImu(), kTolerance));
  EXPECT_TRUE(lhs.gnssLeverArmImu().isApprox(rhs.gnssLeverArmImu(), kTolerance));
  EXPECT_TRUE(
    lhs.rotationImuFromMagnetometer().isApprox(
      rhs.rotationImuFromMagnetometer(), kTolerance));
}

Eigen::Matrix<double, 6, 1> applyBiasLinearPart6(
  const eqf::EquivariantState & X,
  const Eigen::Vector3d & gyro,
  const Eigen::Vector3d & accel)
{
  const Eigen::Matrix3d A = X.rotationBlock();
  const Eigen::Vector3d c_v = X.firstTranslationColumn();
  Eigen::Matrix<double, 6, 1> transformed;
  transformed.head<3>() = A.transpose() * gyro;
  transformed.tail<3>() = A.transpose() * (accel - c_v.cross(gyro));
  return transformed;
}

Eigen::Matrix<double, 6, 1> invertBiasLinearPart6(
  const eqf::EquivariantState & X,
  const Eigen::Matrix<double, 6, 1> & transformed)
{
  const Eigen::Matrix3d A = X.rotationBlock();
  const Eigen::Vector3d c_v = X.firstTranslationColumn();
  Eigen::Matrix<double, 6, 1> original;
  original.head<3>() = A * transformed.head<3>();
  original.tail<3>() =
    A * transformed.tail<3>() + c_v.cross(original.head<3>());
  return original;
}

eqf::EquivariantState composeForRightAction(
  const eqf::EquivariantState & X,
  const eqf::EquivariantState & Y)
{
  eqf::EquivariantState XY;
  XY.setC(X.C() * Y.C());

  const Eigen::Matrix3d A_x = X.rotationBlock();
  XY.setDelta(X.delta() + A_x * Y.delta());
  XY.setE(X.E() * Y.E());

  Eigen::Matrix<double, 6, 1> gamma_x;
  gamma_x.head<3>() = X.gammaGyro();
  gamma_x.tail<3>() = X.gammaAccel();

  const Eigen::Matrix<double, 6, 1> linear_y_gamma =
    applyBiasLinearPart6(Y, Y.gammaGyro(), Y.gammaAccel());
  const Eigen::Matrix<double, 6, 1> gamma_xy =
    gamma_x + invertBiasLinearPart6(XY, linear_y_gamma);

  XY.setGammaGyro(gamma_xy.head<3>());
  XY.setGammaAccel(gamma_xy.tail<3>());
  return XY;
}

}  // namespace

TEST(EqfPhi, IdentityActionLeavesStateUnchanged)
{
  const eqf::State xi = makeState();
  const eqf::EquivariantState identity;

  expectStatesApprox(eqf::phi(identity, xi), xi);
}

TEST(EqfPhi, OriginStateMapsToSymmetryComponents)
{
  const eqf::State origin = eqf::State::Identity();
  const eqf::EquivariantState X = makeSymmetry();
  const eqf::State out = eqf::phi(X, origin);

  const Eigen::Matrix3d A = X.rotationBlock();
  const Eigen::Vector3d c_v = X.firstTranslationColumn();
  const Eigen::Vector3d q_gyro = -X.gammaGyro();
  const Eigen::Vector3d q_accel = -X.gammaAccel();

  EXPECT_TRUE(out.navigationState().isApprox(X.C(), kTolerance));
  EXPECT_TRUE(out.gyroBiasImu().isApprox(A.transpose() * q_gyro, kTolerance));
  EXPECT_TRUE(
    out.accelBiasImu().isApprox(
      A.transpose() * (q_accel - c_v.cross(q_gyro)), kTolerance));
  EXPECT_TRUE(
    out.gnssLeverArmImu().isApprox(-A.transpose() * X.delta(), kTolerance));
  EXPECT_TRUE(out.rotationImuFromMagnetometer().isApprox(A.transpose() * X.E(), kTolerance));
}

TEST(EqfPhi, ExtendedPoseUsesRightMultiplication)
{
  const eqf::State xi = makeState();
  const eqf::EquivariantState X = makeSymmetry();
  const eqf::State out = eqf::phi(X, xi);

  const Eigen::Matrix3d R = xi.rotationGlobalFromImu();
  const Eigen::Vector3d v = xi.velocityGlobal();
  const Eigen::Vector3d p = xi.positionGlobal();
  const Eigen::Matrix3d A = X.rotationBlock();
  const Eigen::Vector3d c_v = X.firstTranslationColumn();
  const Eigen::Vector3d c_p = X.secondTranslationColumn();

  EXPECT_TRUE(out.rotationGlobalFromImu().isApprox(R * A, kTolerance));
  EXPECT_TRUE(out.velocityGlobal().isApprox(v + R * c_v, kTolerance));
  EXPECT_TRUE(out.positionGlobal().isApprox(p + R * c_p, kTolerance));
}

TEST(EqfPhi, AccelBiasUsesVelocityCrossGyroTerm)
{
  eqf::State xi = makeState();
  eqf::EquivariantState X = makeSymmetry();
  X.setC(
    makeSE23(
      Eigen::Matrix3d::Identity(),
      Eigen::Vector3d(0.3, -0.8, 0.4),
      Eigen::Vector3d::Zero()));
  X.setGammaGyro(Eigen::Vector3d(-0.2, 0.1, 0.05));
  X.setGammaAccel(Eigen::Vector3d(0.7, -0.4, 0.2));

  const eqf::State out = eqf::phi(X, xi);
  const Eigen::Vector3d q_gyro = xi.gyroBiasImu() - X.gammaGyro();
  const Eigen::Vector3d q_accel = xi.accelBiasImu() - X.gammaAccel();
  const Eigen::Vector3d expected =
    q_accel - X.firstTranslationColumn().cross(q_gyro);

  EXPECT_TRUE(out.accelBiasImu().isApprox(expected, kTolerance));
}

TEST(EqfPhi, MagnetometerCalibrationOrderIsATransposeSE)
{
  const eqf::State xi = makeState();
  const eqf::EquivariantState X = makeSymmetry();
  const eqf::State out = eqf::phi(X, xi);

  const Eigen::Matrix3d expected =
    X.rotationBlock().transpose() * xi.rotationImuFromMagnetometer() * X.E();

  EXPECT_TRUE(out.rotationImuFromMagnetometer().isApprox(expected, kTolerance));
}

TEST(EqfPhi, RightActionCompositionHolds)
{
  const eqf::State xi = makeState();
  const eqf::EquivariantState X = makeSymmetry();

  eqf::EquivariantState Y;
  Y.setC(
    makeSE23(
      rotationFromVector(Eigen::Vector3d(0.31, 0.17, -0.23)),
      Eigen::Vector3d(-0.5, 0.25, 0.35),
      Eigen::Vector3d(0.45, -0.15, 0.55)));
  Y.setGammaGyro(Eigen::Vector3d(-0.06, 0.03, 0.07));
  Y.setGammaAccel(Eigen::Vector3d(0.25, -0.35, 0.15));
  Y.setDelta(Eigen::Vector3d(-0.2, 0.9, -0.4));
  Y.setE(rotationFromVector(Eigen::Vector3d(-0.19, 0.21, 0.11)));

  const eqf::EquivariantState XY = composeForRightAction(X, Y);
  const eqf::State sequential = eqf::phi(Y, eqf::phi(X, xi));
  const eqf::State direct = eqf::phi(XY, xi);

  expectStatesApprox(sequential, direct);
}
