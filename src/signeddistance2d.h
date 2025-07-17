#ifndef signed_distance_2d_h
#define signed_distance_2d_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <Eigen/Core>

class ShapeTemplate2D;

class SignedDistanceFunction2D
{
	SignedDistanceFunction2D();
public:
  SignedDistanceFunction2D( const Eigen::Vector2d& in_MinVertex, const double in_dx, const Eigen::Vector2i& in_Resolution );
  ~SignedDistanceFunction2D();

  double windingNumber( const Eigen::Vector2i& idx ) const;
  double minimumDistance( const Eigen::Vector2i& idx ) const;
  double signedDistance( const Eigen::Vector2i& idx ) const;
  double signedDistance( const Eigen::Vector2d& in_x ) const;
  Eigen::Vector2d normal( const Eigen::Vector2d& in_x ) const;
  Eigen::Vector2d minVertex() const;
  double dx() const;
  Eigen::Vector2i resolution() const;
  void computeSignedDistanceFunction( const ShapeTemplate2D* in_Obj );

private:
  Eigen::Vector2d m_MinVertex;
	double m_dx;
  Eigen::Vector2i m_Resolution;
	double* m_WindingNumber;
	double* m_MinimumDistance;
	double* m_SignedDistance;
};

#endif

