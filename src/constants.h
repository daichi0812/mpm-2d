#ifndef constants_h
#define constants_h

#include <Eigen/Core>

enum RigidBodyType2D
{
  POLYGON,
  CIRCLE
};

inline Eigen::Matrix2d rot90()
{
  Eigen::Matrix2d res;
  res << 0.0, -1.0, 1.0, 0.0;
  return res;
};

inline Eigen::Matrix2d rotationMat( const double in_theta )
{
  Eigen::Matrix2d res;
  res << cos( in_theta ), -sin( in_theta ), sin( in_theta ), cos( in_theta );
  return res;
}

inline double cross2d( const Eigen::Vector2d& in_a, const Eigen::Vector2d& in_b )
{
  return in_a.x() * in_b.y() - in_a.y() * in_b.x();
}

inline Eigen::Vector2d getCurrentPosition( const Eigen::Vector2d& in_x0, const double in_scale, const double in_theta, const Eigen::Vector2d& in_center_current )
{
  const Eigen::Matrix2d mat{ rotationMat( in_theta ) };
  return mat * in_x0 * in_scale + in_center_current;
}

inline Eigen::Vector2d getMaterialPosition( const Eigen::Vector2d& in_x_current, const double in_scale, const double in_theta, const Eigen::Vector2d& in_center_current )
{
  const Eigen::Matrix2d mat{ rotationMat( in_theta ) };
  return mat.inverse() * ( in_x_current - in_center_current ) / in_scale;
}

#endif
