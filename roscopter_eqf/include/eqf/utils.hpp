#ifndef ROSCOPTER_EQF_UTILS_HPP
#define ROSCOPTER_EQF_UTILS_HPP

#include <Eigen/Dense>

namespace eqf
{

inline Eigen::Matrix3d skew_matrix(const Eigen::Vector3d & vec)
{
  Eigen::Matrix3d skew;
  skew << 0.0, -vec.z(), vec.y(),
          vec.z(), 0.0, -vec.x(),
          -vec.y(), vec.x(), 0.0;
  return skew;
}

}  // namespace eqf

#endif  // ROSCOPTER_EQF_UTILS_HPP
