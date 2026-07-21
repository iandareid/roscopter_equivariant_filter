#ifndef ROSCOPTER_EQF_COVARIANCE_HPP
#define ROSCOPTER_EQF_COVARIANCE_HPP

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <stdexcept>

#include "eqf/error_state.hpp"

namespace eqf
{

using Matrix21d = Eigen::Matrix<double, 21, 21>;

struct CovariancePropagationResult
{
  Matrix21d F = Matrix21d::Identity();
  Matrix21d Qd = Matrix21d::Zero();
  Matrix21d P_pred = Matrix21d::Zero();
};

/**
 * Propagate covariance in the EqF normal-coordinate ordering documented by
 * ErrorState: attitude, velocity, position, gyro-bias group, accel-bias
 * group, GNSS lever arm, and magnetometer rotation (three states each).
 *
 * Qc is a continuous-time spectral-density matrix. Its columns must use the
 * same noise-channel ordering as the columns of L; this routine deliberately
 * does not assign physical meanings or units to those channels.
 */
[[nodiscard]] inline CovariancePropagationResult propagateCovarianceSecondOrder(
  const Matrix21d & P,
  const Matrix21d & A,
  const Eigen::Ref<const Eigen::MatrixXd> & L,
  const Eigen::Ref<const Eigen::MatrixXd> & Qc,
  const double dt)
{
  if (!std::isfinite(dt) || dt <= 0.0) {
    throw std::invalid_argument("EqF covariance propagation requires finite dt > 0");
  }
  if (L.rows() != ErrorState::kErrorStateSize || Qc.rows() != Qc.cols() ||
    L.cols() != Qc.rows())
  {
    throw std::invalid_argument("EqF covariance propagation matrix dimensions are incompatible");
  }
  if (!P.allFinite() || !A.allFinite() || !L.allFinite() || !Qc.allFinite()) {
    throw std::invalid_argument("EqF covariance propagation input contains NaN or infinity");
  }

  CovariancePropagationResult result;
  const Matrix21d A_dt = A * dt;
  result.F += A_dt + 0.5 * (A_dt * A_dt);

  result.Qd = (L * Qc * L.transpose()) * dt;
  result.Qd = 0.5 * (result.Qd + result.Qd.transpose()).eval();

  result.P_pred = result.F * P * result.F.transpose() + result.Qd;
  result.P_pred = 0.5 * (result.P_pred + result.P_pred.transpose()).eval();

  if (!result.F.allFinite() || !result.Qd.allFinite() || !result.P_pred.allFinite()) {
    throw std::runtime_error("EqF covariance propagation produced NaN or infinity");
  }

#ifndef NDEBUG
  const Eigen::SelfAdjointEigenSolver<Matrix21d> solver(result.P_pred);
  if (solver.info() != Eigen::Success || solver.eigenvalues().minCoeff() < -1e-10) {
    throw std::runtime_error("EqF covariance propagation produced an indefinite covariance");
  }
#endif

  return result;
}

}  // namespace eqf

#endif  // ROSCOPTER_EQF_COVARIANCE_HPP
