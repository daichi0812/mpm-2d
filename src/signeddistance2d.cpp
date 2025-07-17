//
//  signeddistance2d.cpp
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/03.
//

#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE

#include "signeddistance2d.h"
#include "RigidBody2DData.h"

SignedDistanceFunction2D::SignedDistanceFunction2D( const Eigen::Vector2d& in_MinVertex, const double in_dx, const Eigen::Vector2i& in_Resolution )
  : m_MinVertex( in_MinVertex ), m_dx( in_dx ), m_Resolution( in_Resolution ), m_WindingNumber( nullptr ), m_MinimumDistance( nullptr ), m_SignedDistance( nullptr )
{
  m_WindingNumber = ( double* )malloc( sizeof(double) * ( m_Resolution(0) + 1 ) * ( m_Resolution(1) + 1 ) );
  m_MinimumDistance = ( double* )malloc( sizeof(double) * ( m_Resolution(0) + 1 ) * ( m_Resolution(1) + 1 ) );
  m_SignedDistance = ( double* )malloc( sizeof(double) * ( m_Resolution(0) + 1 ) * ( m_Resolution(1) + 1 ) );
}

SignedDistanceFunction2D::~SignedDistanceFunction2D()
{
  if (m_WindingNumber)
  {
    free(m_WindingNumber);
    m_WindingNumber = nullptr;
  }
  if (m_MinimumDistance)
  {
    free(m_MinimumDistance);
    m_MinimumDistance = nullptr;
  }
  if (m_SignedDistance)
  {
    free(m_SignedDistance);
    m_SignedDistance = nullptr;
  }
}

double SignedDistanceFunction2D::windingNumber( const Eigen::Vector2i& idx ) const
{
  if (idx(0) < 0 || idx(0) > m_Resolution(0) || idx(1) < 0 || idx(1) > m_Resolution(1))
    return 0.0;

  return m_WindingNumber[idx(1) * (m_Resolution(0)) + idx(0)];
}

double SignedDistanceFunction2D::minimumDistance( const Eigen::Vector2i& idx ) const
{
  if (idx(0) < 0 || idx(0) > m_Resolution(0) || idx(1) < 0 || idx(1) > m_Resolution(1))
    return 0.0;

  return m_MinimumDistance[idx(1) * (m_Resolution(0)) + idx(0)];
}

double SignedDistanceFunction2D::signedDistance( const Eigen::Vector2i& idx ) const
{
  if (idx(0) < 0 || idx(0) > m_Resolution(0) || idx(1) < 0 || idx(1) > m_Resolution(1))
    return 1.0e33;

  return m_SignedDistance[idx(1) * (m_Resolution(0)) + idx(0)];
}

double SignedDistanceFunction2D::signedDistance( const Eigen::Vector2d& in_x ) const
{
  const Eigen::Vector2d _x = ( in_x - m_MinVertex ) / m_dx;
  const int i = int( floor( _x(0) ) ); const int j = int( floor( _x(1) ) );
  
  if( i < 0 || i >= m_Resolution(0) || j < 0 || j >= m_Resolution(1) )
    return 1.0e33;
  
  const double s = _x(0) - i; const double t = _x(1) - j;
  const int flat_idx_00 = j * m_Resolution(0) + i;
  const int flat_idx_01 = j * m_Resolution(0) + i + 1;
  const int flat_idx_10 = (j + 1) * m_Resolution(0) + i;
  const int flat_idx_11 = (j + 1) * m_Resolution(0) + i + 1;
  
  return (1-t) * ( (1-s) * m_SignedDistance[flat_idx_00] + s * m_SignedDistance[flat_idx_01] )
    + t * ( (1-s) * m_SignedDistance[flat_idx_10] + s * m_SignedDistance[flat_idx_11] );
}

Eigen::Vector2d SignedDistanceFunction2D::normal( const Eigen::Vector2d& in_x ) const
{
  const Eigen::Vector2d _x = ( in_x - m_MinVertex ) / m_dx;
  const int i = int( floor( _x(0) ) ); const int j = int( floor( _x(1) ) );
  
  if( i < 0 || i >= m_Resolution(0) || j < 0 || j >= m_Resolution(1) )
    return Eigen::Vector2d::Zero();
  
  const double s = _x(0) - i; const double t = _x(1) - j;
  const int flat_idx_00 = j * m_Resolution(0) + i;
  const int flat_idx_01 = j * m_Resolution(0) + i + 1;
  const int flat_idx_10 = (j + 1) * m_Resolution(0) + i;
  const int flat_idx_11 = (j + 1) * m_Resolution(0) + i + 1;
  
  const double sd00 = m_SignedDistance[flat_idx_00];
  const double sd01 = m_SignedDistance[flat_idx_01];
  const double sd10 = m_SignedDistance[flat_idx_10];
  const double sd11 = m_SignedDistance[flat_idx_11];
  
  Eigen::Vector2d _n{ (1-t)*(sd01 - sd00) + t*(sd11-sd10), (1-s)*(sd10-sd00) + s*(sd11-sd01) };
  _n.normalize();
  return _n;
}

Eigen::Vector2d SignedDistanceFunction2D::minVertex() const
{
  return m_MinVertex;
}

double SignedDistanceFunction2D::dx() const
{
  return m_dx;
}

Eigen::Vector2i SignedDistanceFunction2D::resolution() const
{
  return m_Resolution;
}

void SignedDistanceFunction2D::computeSignedDistanceFunction( const ShapeTemplate2D* in_Obj )
{
  for( int j = 0; j <= m_Resolution(1); j++ )
  {
    for( int i = 0; i <= m_Resolution(0); i++ )
    {
      const Eigen::Vector2d x = m_MinVertex + Eigen::Vector2d{ i, j } * m_dx;

      const int flat_idx = j * m_Resolution(0) + i;

      m_WindingNumber[flat_idx] = in_Obj->computeWindingNumber(x);
      m_MinimumDistance[flat_idx] = in_Obj->computeMinimumDistance(x);

      if( m_WindingNumber[flat_idx] > 0.5 )
        m_SignedDistance[flat_idx] = -m_MinimumDistance[flat_idx];
      else
        m_SignedDistance[flat_idx] = m_MinimumDistance[flat_idx];
    }
  }
}
