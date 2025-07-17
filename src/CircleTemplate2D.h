//
//  CircleTemplate2D.h
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/08.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#ifndef __CIRCLE_TEMPLATE_2D_H__
#define __CIRCLE_TEMPLATE_2D_H__

#define _USE_MATH_DEFINES
#include<math.h>

#include <Eigen/Core>
#include <vector>
#include <map>
#include <iostream>
#include "constants.h"
#include "collisionsample2d.h"
#include "ShapeTemplate2D.h"
#include "signeddistance2d.h"

class SignedDistanceFunction2D;

class CircleTemplate2D : public ShapeTemplate2D
{
	CircleTemplate2D() = delete;
public:
  CircleTemplate2D( const std::string& in_name, const double& in_density, const double& in_radius, const double& in_radius_std );
  
  ~CircleTemplate2D();

	RigidBodyType2D type() const { return CIRCLE; }

	void getData( ShapeTemplateData& out_data ) const
	{
		out_data.name = m_Name;
		out_data.density = m_Density;
		out_data.size_mean = m_RadiusMean * 2.0;
		out_data.size_std = m_RadiusStd * 2.0;
		out_data.vertex_list.setZero(2, 1);
	}

	const std::string& getName() const
	{
		return m_Name;
	}

	double getDensity() const
	{
		return m_Density;
	}
  
  double computeWindingNumber( const Eigen::Vector2d& in_x0 ) const
  {
    if( in_x0.norm() > m_RadiusMean ) return 0.0;
    else return 1.0;
  }
  
  double computeMinimumDistance( const Eigen::Vector2d& in_x0 ) const
  {
    const double s = atan2( in_x0(1), in_x0(0) );
    const Eigen::Vector2d p0{ m_RadiusMean * cos(s), m_RadiusMean * sin(s) };
    return ( p0 - in_x0 ).norm();
  }
  
  double getSignedDistanceAt( const Eigen::Vector2d& in_x0 ) const
  {
    double winding_number = computeWindingNumber( in_x0 );
    double minimumdistance = computeMinimumDistance( in_x0 );
    return winding_number < 0.5 ? minimumdistance : -minimumdistance;
    
    //return m_SDF->signedDistance( in_x0 );
  }
  
  Eigen::Vector2d getNormalAt( const Eigen::Vector2d& in_x0 ) const
  {
    //const Eigen::Vector2d n0 = m_SDF->normal( in_x0 );
    const Eigen::Vector2d n0 = in_x0.normalized();
    
    return n0;
  }
  
  void computeSignedDistanceFunction( const double dx );
  
  const SignedDistanceFunction2D* getSDF() const;
  
  void getBoundingBox( BoundingBox2D& out_BB, const double in_scale, const double in_theta, const Eigen::Vector2d& in_center_current ) const
  {
    out_BB.bb_min(0) = -( m_RadiusMean * in_scale ) + in_center_current(0);
    out_BB.bb_min(1) = -( m_RadiusMean * in_scale ) + in_center_current(1);
    out_BB.bb_max(0) =  ( m_RadiusMean * in_scale ) + in_center_current(0);
    out_BB.bb_max(1) =  ( m_RadiusMean * in_scale ) + in_center_current(1);
  }
  
  // scales with length^2
  double computeArea() const { return M_PI * m_RadiusMean * m_RadiusMean; }
  // scales with length^4
  double computeSecondMomentOfArea() const { return 0.5 * M_PI * m_RadiusMean * m_RadiusMean * m_RadiusMean * m_RadiusMean; }
  
  void generateCollisionSamples( std::vector< CollisionSample2D >& out_CollisionSamples, const double in_dx_sample_points ) const
  {
    out_CollisionSamples.clear();
    const double length = 2.0 * M_PI * m_RadiusMean;
    const int nSegs = std::max<int>( 2, int( ceil( length / in_dx_sample_points ) ) );
    
    for( int i=0; i<nSegs; i++ )
    {
      const double theta = 2.0 * M_PI * ( i + 0.5 ) / nSegs;
      const Eigen::Vector2d p { m_RadiusMean * cos( theta ), m_RadiusMean * sin( theta ) };
      CollisionSample2D sample(p);
      out_CollisionSamples.push_back( sample );
    }
  }
  
  double getRadius0() const { return m_RadiusMean; }

private:
  
	std::string m_Name;
	double m_Density;
	double m_RadiusMean;
	double m_RadiusStd;
  
  const SignedDistanceFunction2D* m_SDF;
};

#endif
