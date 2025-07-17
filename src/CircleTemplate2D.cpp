//
//  CircleTemplate2D.cpp
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/08.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE

#include "CircleTemplate2D.h"
#include "signeddistance2d.h"

CircleTemplate2D::CircleTemplate2D( const std::string& in_name, const double& in_density, const double& in_radius, const double& in_radius_std )
: m_Name( in_name )
, m_Density( in_density )
, m_RadiusMean( in_radius )
, m_RadiusStd( in_radius_std )
, m_SDF( nullptr )
{}

CircleTemplate2D::~CircleTemplate2D()
{
  if( m_SDF )
  {
    delete m_SDF;
    m_SDF = nullptr;
  }
}

void CircleTemplate2D::computeSignedDistanceFunction( const double dx )
{
  const Eigen::Vector2d min_C{ -m_RadiusMean, -m_RadiusMean };
  const Eigen::Vector2d max_C{ m_RadiusMean, m_RadiusMean };

  const Eigen::Vector2d center = ( max_C + min_C ) * 0.5;
  const Eigen::Vector2d half_width = ( max_C - min_C ) * 0.5 * 2.0;

  const Eigen::Vector2i resolution = ( half_width * 2.0 / dx ).cast<int>();

  const Eigen::Vector2d gridMin = center - resolution.cast<double>() * dx * 0.5;
  SignedDistanceFunction2D* _sdf = new SignedDistanceFunction2D( gridMin, dx, resolution );
  _sdf->computeSignedDistanceFunction( this );
  m_SDF = _sdf;
}

const SignedDistanceFunction2D* CircleTemplate2D::getSDF() const
{
  return m_SDF;
}
