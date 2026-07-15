#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "eqf/state.hpp"
#include "eqf/utils.hpp"

namespace
{

constexpr double kTolerance = 1e-10;

Eigen::Matrix3d rotationFromVector(const Eigen::Vector3d & vector)
{
  return Eigen::AngleAxisd(vector.norm(), vector.normalized()).toRotationMatrix();
}

eqf::State::SE23Matrix makeSE23(
  const Eigen::Matrix3d & R,
  const Eigen::Vector3d & v,
  const Eigen::Vector3d & p)
{
  eqf::State::SE23Matrix T = eqf::State::SE23Matrix::Identity();
  T.block<3, 3>(0, 0) = R;
  T.block<3, 1>(0, 3) = v;
  T.block<3, 1>(0, 4) = p;
  return T;
}

eqf::State makeState(
  const Eigen::Matrix3d & R,
  const Eigen::Vector3d & v,
  const Eigen::Vector3d & p,
  const Eigen::Vector3d & b_gyro,
  const Eigen::Vector3d & b_accel,
  const Eigen::Vector3d & lever_arm,
  const Eigen::Matrix3d & S)
{
  return eqf::State(makeSE23(R, v, p), b_gyro, b_accel, lever_arm, S);
}

eqf::State makeGenericState(const Eigen::Vector3d & velocity)
{
  return makeState(
    rotationFromVector(Eigen::Vector3d(0.25, -0.35, 0.45)),
    velocity,
    Eigen::Vector3d(-1.0, 0.8, 0.2),
    Eigen::Vector3d(0.03, -0.04, 0.05),
    Eigen::Vector3d(-0.2, 0.3, 0.4),
    Eigen::Vector3d(0.7, -0.1, 0.5),
    rotationFromVector(Eigen::Vector3d(-0.15, 0.2, 0.3)));
}

eqf::State::SE23Matrix lambda1MatrixReference(
  const eqf::State & state,
  const eqf::ImuInput & input,
  const Eigen::Vector3d & gravity)
{
  const eqf::State::SE23Matrix W = eqf::se23Wedge(
    input.angular_velocity,
    input.linear_acceleration,
    Eigen::Vector3d::Zero());
  const eqf::State::SE23Matrix B = eqf::se23Wedge(
    state.gyroBiasImu(),
    state.accelBiasImu(),
    Eigen::Vector3d::Zero());
  const eqf::State::SE23Matrix G = eqf::se23Wedge(
    Eigen::Vector3d::Zero(),
    gravity,
    Eigen::Vector3d::Zero());
  eqf::State::SE23Matrix N = eqf::State::SE23Matrix::Zero();
  N(3, 4) = 1.0;

  return (W - B + N) +
         state.navigationState().inverse() * (G - N) * state.navigationState();
}

void expectLiftZero(const eqf::Lift & lift)
{
  EXPECT_TRUE(lift.lambda1.isApprox(eqf::State::SE23Matrix::Zero(), kTolerance));
  EXPECT_TRUE(lift.lambda2.isApprox(eqf::Vector6d::Zero(), kTolerance));
  EXPECT_TRUE(lift.lambda3.isApprox(Eigen::Vector3d::Zero(), kTolerance));
  EXPECT_TRUE(lift.lambda4.isApprox(Eigen::Vector3d::Zero(), kTolerance));
}

}  // namespace

TEST(EqfLambda, ZeroStateZeroInputKeepsOnlyGravityInLambda1)
{
  const eqf::State state = eqf::State::Identity();
  const eqf::ImuInput input;
  const Eigen::Vector3d gravity(0.0, 0.0, 9.81);

  const eqf::Lift lift = eqf::lambda(state, input, gravity);

  EXPECT_TRUE(
    lift.lambda1.isApprox(
      eqf::se23Wedge(
        Eigen::Vector3d::Zero(),
        gravity,
        Eigen::Vector3d::Zero()),
      kTolerance));
  EXPECT_TRUE(lift.lambda2.isApprox(eqf::Vector6d::Zero(), kTolerance));
  EXPECT_TRUE(lift.lambda3.isApprox(Eigen::Vector3d::Zero(), kTolerance));
  EXPECT_TRUE(lift.lambda4.isApprox(Eigen::Vector3d::Zero(), kTolerance));
}

TEST(EqfLambda, StationarySpecificForceCancelsGravity)
{
  const Eigen::Matrix3d R = rotationFromVector(Eigen::Vector3d(0.2, -0.1, 0.3));
  const Eigen::Vector3d b_gyro(0.04, -0.03, 0.02);
  const Eigen::Vector3d b_accel(-0.2, 0.1, 0.3);
  const Eigen::Vector3d gravity(0.0, 0.0, 9.81);
  const eqf::State state = makeState(
    R,
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d(0.4, -0.5, 0.6),
    b_gyro,
    b_accel,
    Eigen::Vector3d(0.3, 0.2, -0.1),
    rotationFromVector(Eigen::Vector3d(-0.1, 0.2, 0.15)));

  eqf::ImuInput input;
  input.angular_velocity = b_gyro;
  input.linear_acceleration = -R.transpose() * gravity + b_accel;

  expectLiftZero(eqf::lambda(state, input, gravity));
}

TEST(EqfLambda, Lambda1MatchesDirectMatrixExpression)
{
  const eqf::State state = makeGenericState(Eigen::Vector3d(3.0, -2.0, 1.0));
  eqf::ImuInput input;
  input.angular_velocity = Eigen::Vector3d(0.4, -0.6, 0.8);
  input.linear_acceleration = Eigen::Vector3d(-1.2, 0.7, 0.9);
  const Eigen::Vector3d gravity(0.1, -0.2, 9.7);

  EXPECT_TRUE(
    eqf::lambda1(state, input, gravity).isApprox(
      lambda1MatrixReference(state, input, gravity), kTolerance));
}

TEST(EqfLambda, Lambda2UsesFirstTranslationProjectionOnly)
{
  const eqf::State state_a = makeGenericState(Eigen::Vector3d(1.0e6, -2.0e6, 3.0e6));
  const eqf::State state_b = makeGenericState(Eigen::Vector3d(-4.0e6, 5.0e6, -6.0e6));
  eqf::ImuInput input;
  input.angular_velocity = Eigen::Vector3d(0.5, -0.1, 0.2);
  input.linear_acceleration = Eigen::Vector3d(0.9, -0.8, 0.7);
  const Eigen::Vector3d gravity(0.0, 0.0, 9.81);

  EXPECT_TRUE(
    eqf::lambda2(state_a, input, gravity).isApprox(
      eqf::lambda2(state_b, input, gravity), kTolerance));
  const bool second_translation_matches =
    eqf::lambda1(state_a, input, gravity).block<3, 1>(0, 4).isApprox(
    eqf::lambda1(state_b, input, gravity).block<3, 1>(0, 4), kTolerance);
  EXPECT_FALSE(second_translation_matches);
}

TEST(EqfLambda, Lambda2LieBracketSignIsBiasThenProjection)
{
  const eqf::State state = makeState(
    Eigen::Matrix3d::Identity(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::UnitX(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Matrix3d::Identity());
  eqf::ImuInput input;
  input.angular_velocity = Eigen::Vector3d::UnitX() + Eigen::Vector3d::UnitY();
  input.linear_acceleration = Eigen::Vector3d::Zero();

  EXPECT_TRUE(
    eqf::lambda2(state, input, Eigen::Vector3d::Zero()).head<3>().isApprox(
      Eigen::Vector3d::UnitZ(), kTolerance));
}

TEST(EqfLambda, Lambda3UsesLeverArmCrossCorrectedAngularVelocity)
{
  const eqf::State state = makeState(
    Eigen::Matrix3d::Identity(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::UnitX(),
    Eigen::Matrix3d::Identity());
  eqf::ImuInput input;
  input.angular_velocity = Eigen::Vector3d::UnitY();

  EXPECT_TRUE(eqf::lambda3(state, input).isApprox(Eigen::Vector3d::UnitZ(), kTolerance));
}

TEST(EqfLambda, Lambda4UsesMagnetometerTranspose)
{
  const Eigen::Matrix3d S =
    Eigen::AngleAxisd(0.5, Eigen::Vector3d::UnitZ()).toRotationMatrix();
  const eqf::State state = makeState(
    Eigen::Matrix3d::Identity(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    S);
  eqf::ImuInput input;
  input.angular_velocity = Eigen::Vector3d::UnitX();

  EXPECT_TRUE(
    eqf::lambda4(state, input).isApprox(
      S.transpose() * input.angular_velocity, kTolerance));
  EXPECT_FALSE(eqf::lambda4(state, input).isApprox(S * input.angular_velocity, kTolerance));
}

TEST(EqfLambda, BiasFreeStateHasZeroLambda2)
{
  const eqf::State state = makeState(
    rotationFromVector(Eigen::Vector3d(-0.25, 0.15, 0.4)),
    Eigen::Vector3d(2.0, -1.0, 0.5),
    Eigen::Vector3d(0.4, 0.3, -0.2),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d(0.3, -0.4, 0.5),
    rotationFromVector(Eigen::Vector3d(0.1, 0.2, -0.3)));
  eqf::ImuInput input;
  input.angular_velocity = Eigen::Vector3d(0.7, -0.6, 0.5);
  input.linear_acceleration = Eigen::Vector3d(-0.4, 0.3, 1.2);

  EXPECT_TRUE(
    eqf::lambda2(state, input, Eigen::Vector3d(0.0, 0.0, 9.81)).isApprox(
      eqf::Vector6d::Zero(), kTolerance));
}
