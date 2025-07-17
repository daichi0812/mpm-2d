//
//  PolygonTemplate2D.cpp
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/08.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE

#include "PolygonTemplate2D.h"
#include "signeddistance2d.h"

PolygonTemplate2D::PolygonTemplate2D( const std::string& in_name, const double& in_density, const Eigen::Matrix2Xd& in_vertices, const double& in_size_std )
: m_Name( in_name )
, m_Density( in_density )
, m_Vertices( in_vertices )
, m_SizeStd( in_size_std )
, m_SDF( nullptr )
{
  setCenterOfMassToOrigin();
}

PolygonTemplate2D::~PolygonTemplate2D()
{
  if( m_SDF )
  {
    delete m_SDF;
    m_SDF = nullptr;
  }
}

void PolygonTemplate2D::computeSignedDistanceFunction( const double dx )
{
  double min_x = 1.0e33;
  double min_y = 1.0e33;
  double max_x = -1.0e33;
  double max_y = -1.0e33;

  for( int i = 0; i < m_Vertices.cols(); i++ )
  {
    min_x = std::min<double>( min_x, m_Vertices.col(i)(0) );
    min_y = std::min<double>( min_y, m_Vertices.col(i)(1) );
    max_x = std::max<double>( max_x, m_Vertices.col(i)(0) );
    max_y = std::max<double>( max_y, m_Vertices.col(i)(1) );
  }

  const Eigen::Vector2d min_C{ min_x, min_y };
  const Eigen::Vector2d max_C{ max_x, max_y };

  const Eigen::Vector2d center = ( max_C + min_C ) * 0.5;
  const Eigen::Vector2d half_width = ( max_C - min_C ) * 0.5 * 2.0;

  const Eigen::Vector2i resolution = ( half_width * 2.0 / dx ).cast<int>();

  const Eigen::Vector2d gridMin = center - resolution.cast<double>() * dx * 0.5;
  SignedDistanceFunction2D* _sdf = new SignedDistanceFunction2D( gridMin, dx, resolution );
  _sdf->computeSignedDistanceFunction( this );
  m_SDF = _sdf;
}

const SignedDistanceFunction2D* PolygonTemplate2D::getSDF() const
{
  return m_SDF;
}

