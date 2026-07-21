#ifndef ROSCOPTER_EQF_UTILS_HPP
#define ROSCOPTER_EQF_UTILS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <cmath>

#include "eqf/equivariant_state.hpp"
#include "eqf/state.hpp"

namespace eqf
{

using Vector6d = Eigen::Matrix<double, 6, 1>;

struct ImuInput
{
  Eigen::Vector3d angular_velocity = Eigen::Vector3d::Zero();
  Eigen::Vector3d linear_acceleration = Eigen::Vector3d::Zero();
};

struct Lift
{
  State::SE23Matrix lambda1 = State::SE23Matrix::Zero();
  Vector6d lambda2 = Vector6d::Zero();
  Eigen::Vector3d lambda3 = Eigen::Vector3d::Zero();
  Eigen::Vector3d lambda4 = Eigen::Vector3d::Zero();
};

struct ImuPropagationResult
{
  bool propagated = false;
  double dt = 0.0;
  State state_before = State::Identity();
  Lift lift;
  EquivariantState increment;
  EquivariantState observer_after;
  State state_after = State::Identity();
};

inline Eigen::Matrix3d skew_matrix(const Eigen::Vector3d & vec)
{
  Eigen::Matrix3d skew;
  skew << 0.0, -vec.z(), vec.y(),
    vec.z(), 0.0, -vec.x(),
    -vec.y(), vec.x(), 0.0;
  return skew;
}

[[nodiscard]] inline State::SE23Matrix se23Wedge(
  const Eigen::Vector3d & rotation,
  const Eigen::Vector3d & first_translation,
  const Eigen::Vector3d & second_translation)
{
  State::SE23Matrix result = State::SE23Matrix::Zero();
  result.block<3, 3>(0, 0) = skew_matrix(rotation);
  result.block<3, 1>(0, 3) = first_translation;
  result.block<3, 1>(0, 4) = second_translation;
  return result;
}

[[nodiscard]] inline EquivariantState::Se3AlgebraMatrix se3Wedge(
  const Eigen::Vector3d & rotation,
  const Eigen::Vector3d & translation)
{
  EquivariantState::Se3AlgebraMatrix result =
    EquivariantState::Se3AlgebraMatrix::Zero();
  result.block<3, 3>(0, 0) = skew_matrix(rotation);
  result.block<3, 1>(0, 3) = translation;
  return result;
}

[[nodiscard]] inline Vector6d se3LieBracket(
  const Vector6d & x,
  const Vector6d & y)
{
  Vector6d result;
  result.head<3>() = x.head<3>().cross(y.head<3>());
  result.tail<3>() =
    x.tail<3>().cross(y.head<3>()) + x.head<3>().cross(y.tail<3>());
  return result;
}

[[nodiscard]] inline Lift scaleLift(const Lift & lift, const double scale)
{
  Lift out;
  out.lambda1 = scale * lift.lambda1;
  out.lambda2 = scale * lift.lambda2;
  out.lambda3 = scale * lift.lambda3;
  out.lambda4 = scale * lift.lambda4;
  return out;
}

[[nodiscard]] inline Vector6d biasVector(const EquivariantState & X)
{
  Vector6d gamma;
  gamma.head<3>() = X.gammaGyro();
  gamma.tail<3>() = X.gammaAccel();
  return gamma;
}

[[nodiscard]] inline Vector6d applyBiasLinearPart(
  const EquivariantState & X,
  const Vector6d & bias)
{
  const Eigen::Matrix3d A = X.rotationBlock();
  const Eigen::Vector3d c_v = X.firstTranslationColumn();

  Vector6d transformed;
  transformed.head<3>() = A.transpose() * bias.head<3>();
  transformed.tail<3>() =
    A.transpose() * (bias.tail<3>() - c_v.cross(bias.head<3>()));
  return transformed;
}

[[nodiscard]] inline Vector6d invertBiasLinearPart(
  const EquivariantState & X,
  const Vector6d & transformed)
{
  const Eigen::Matrix3d A = X.rotationBlock();
  const Eigen::Vector3d c_v = X.firstTranslationColumn();

  Vector6d bias;
  bias.head<3>() = A * transformed.head<3>();
  bias.tail<3>() = A * transformed.tail<3>() + c_v.cross(bias.head<3>());
  return bias;
}

[[nodiscard]] inline EquivariantState groupMultiply(
  const EquivariantState & X,
  const EquivariantState & Y)
{
  EquivariantState XY;
  XY.setC(X.C() * Y.C());
  XY.setDelta(X.delta() + X.rotationBlock() * Y.delta());
  XY.setE(X.E() * Y.E());

  const Vector6d gamma_xy =
    biasVector(X) + invertBiasLinearPart(XY, applyBiasLinearPart(Y, biasVector(Y)));
  XY.setGammaGyro(gamma_xy.head<3>());
  XY.setGammaAccel(gamma_xy.tail<3>());
  return XY;
}

[[nodiscard]] inline EquivariantState expGroup(const Lift & lift)
{
  EquivariantState out;
  out.setC(lift.lambda1.exp());
  out.setGamma(se3Wedge(lift.lambda2.head<3>(), lift.lambda2.tail<3>()));
  out.setDelta(lift.lambda3);
  out.setE(skew_matrix(lift.lambda4).exp());
  return out;
}

[[nodiscard]] inline State phi(
  const EquivariantState & X,
  const State & xi)
{
  const Eigen::Matrix3d A = X.rotationBlock();
  const Eigen::Vector3d c_v = X.firstTranslationColumn();

  const Eigen::Vector3d q_gyro = xi.gyroBiasImu() - X.gammaGyro();
  const Eigen::Vector3d q_accel = xi.accelBiasImu() - X.gammaAccel();

  State out;
  out.setNavigationState(xi.navigationState() * X.C());
  out.setGyroBiasImu(A.transpose() * q_gyro);
  out.setAccelBiasImu(A.transpose() * (q_accel - c_v.cross(q_gyro)));
  out.setGnssLeverArmImu(A.transpose() * (xi.gnssLeverArmImu() - X.delta()));
  out.setRotationImuFromMagnetometer(
    A.transpose() * xi.rotationImuFromMagnetometer() * X.E());

  return out;
}

[[nodiscard]] inline Eigen::Vector3d correctedAngularVelocity(
  const State & state,
  const ImuInput & input)
{
  return input.angular_velocity - state.gyroBiasImu();
}

[[nodiscard]] inline Eigen::Vector3d lambdaAlpha(
  const State & state,
  const ImuInput & input,
  const Eigen::Vector3d & gravity_in_inertial)
{
  return input.linear_acceleration - state.accelBiasImu() +
         state.rotationGlobalFromImu().transpose() * gravity_in_inertial;
}

[[nodiscard]] inline State::SE23Matrix lambda1(
  const State & state,
  const ImuInput & input,
  const Eigen::Vector3d & gravity_in_inertial)
{
  return se23Wedge(
    correctedAngularVelocity(state, input),
    lambdaAlpha(state, input, gravity_in_inertial),
    state.rotationGlobalFromImu().transpose() * state.velocityGlobal());
}

[[nodiscard]] inline Vector6d lambda2(
  const State & state,
  const ImuInput & input,
  const Eigen::Vector3d & gravity_in_inertial)
{
  Vector6d bias;
  bias.head<3>() = state.gyroBiasImu();
  bias.tail<3>() = state.accelBiasImu();

  Vector6d projected_lambda1;
  projected_lambda1.head<3>() = correctedAngularVelocity(state, input);
  projected_lambda1.tail<3>() = lambdaAlpha(state, input, gravity_in_inertial);

  return se3LieBracket(bias, projected_lambda1);
}

[[nodiscard]] inline Eigen::Vector3d lambda3(
  const State & state,
  const ImuInput & input)
{
  return state.gnssLeverArmImu().cross(correctedAngularVelocity(state, input));
}

[[nodiscard]] inline Eigen::Vector3d lambda4(
  const State & state,
  const ImuInput & input)
{
  return state.rotationImuFromMagnetometer().transpose() *
         correctedAngularVelocity(state, input);
}

[[nodiscard]] inline Lift lambda(
  const State & state,
  const ImuInput & input,
  const Eigen::Vector3d & gravity_in_inertial)
{
  Lift out;
  out.lambda1 = lambda1(state, input, gravity_in_inertial);
  out.lambda2 = lambda2(state, input, gravity_in_inertial);
  out.lambda3 = lambda3(state, input);
  out.lambda4 = lambda4(state, input);
  return out;
}

[[nodiscard]] inline ImuPropagationResult propagateImuBeforeA(
  const EquivariantState & observer_state,
  const State & reference_state,
  const ImuInput & input,
  const Eigen::Vector3d & gravity_in_inertial,
  const double dt)
{
  ImuPropagationResult result;
  result.dt = dt;
  result.state_before = phi(observer_state, reference_state);
  result.observer_after = observer_state;
  result.state_after = result.state_before;

  if (!std::isfinite(dt) || dt <= 0.0 ||
    !input.angular_velocity.allFinite() ||
    !input.linear_acceleration.allFinite())
  {
    return result;
  }

  result.lift = lambda(result.state_before, input, gravity_in_inertial);
  result.increment = expGroup(scaleLift(result.lift, dt));
  result.observer_after = groupMultiply(observer_state, result.increment);
  result.state_after = phi(result.observer_after, reference_state);
  result.propagated = true;
  return result;
}

}  // namespace eqf

#endif  // ROSCOPTER_EQF_UTILS_HPP
