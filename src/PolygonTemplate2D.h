//
//  PolygonTemplate2D.h
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/08.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#ifndef __POLYGON_TEMPLATE_2D_H__
#define __POLYGON_TEMPLATE_2D_H__

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
class ShapeTemplate2D;

class PolygonTemplate2D : public ShapeTemplate2D
{
  PolygonTemplate2D() = delete;
public:
  PolygonTemplate2D( const std::string& in_name, const double& in_density, const Eigen::Matrix2Xd& in_vertices, const double& in_size_std );
  ~PolygonTemplate2D();

  RigidBodyType2D type() const { return POLYGON; }

  void getData( ShapeTemplateData& out_data ) const
  {
    Eigen::Vector2d min_coord = m_Vertices.col(0);
    Eigen::Vector2d max_coord = m_Vertices.col(0);
    for( int i=1; i<m_Vertices.cols(); i++ )
    {
      min_coord = min_coord.cwiseMin( m_Vertices.col( i ) );
      max_coord = max_coord.cwiseMax( m_Vertices.col( i ) );
    }

    const Eigen::Vector2d size = max_coord - min_coord;

    out_data.name = m_Name;
    out_data.density = m_Density;
    out_data.size_mean = 0.5 * ( size(0) + size(1) );
    out_data.size_std = m_SizeStd;
    out_data.vertex_list = m_Vertices;
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
    double totalRad = 0.0;
    for( int i = 0; i < m_Vertices.cols(); i++ )
    {
      const int ip = ( i + 1 ) % m_Vertices.cols();
      totalRad += computeThetaFromSinCos( in_x0, m_Vertices.col(i), m_Vertices.col(ip) );
    }
    return totalRad / ( 2.0 * M_PI );
  }
  
  double computeMinimumDistance( const Eigen::Vector2d& in_x0 ) const
  {
    double minDistance = 1.0e33;
    for( int i = 0; i < m_Vertices.cols(); i++ )
    {
      const int ip = ( i + 1 ) % m_Vertices.cols();
      double _dist = computeDistanceToSegment( in_x0, m_Vertices.col(i), m_Vertices.col(ip) );
      minDistance = std::min<double>( minDistance, _dist );
    }
    return minDistance;
  }
  
  double getSignedDistanceAt( const Eigen::Vector2d& in_x0 ) const
  {
    double winding_number = computeWindingNumber( in_x0 );
    double minimumdistance = computeMinimumDistance( in_x0 );
    return winding_number < 0.5 ? minimumdistance : -minimumdistance;
    
    // return m_SDF->signedDistance( in_x0 );
  }
  
  Eigen::Vector2d getNormalAt( const Eigen::Vector2d& in_x0 ) const
  {
    const Eigen::Vector2d q0 = computeClosestPoint( in_x0 );
    double winding_number = computeWindingNumber( in_x0 );
    
    const Eigen::Vector2d n0 = winding_number < 0.5 ? ( in_x0 - q0 ).normalized() : ( q0 - in_x0 ).normalized();
    
    // const Eigen::Vector2d n0 = m_SDF->normal( in_x0 );
    
    return n0;
  }
  
  void computeSignedDistanceFunction( const double dx );
  
  const SignedDistanceFunction2D* getSDF() const;
  
  void getBoundingBox( BoundingBox2D& out_BB, const double in_scale, const double in_theta, const Eigen::Vector2d& in_center_current ) const
  {
    out_BB.bb_min(0) = 1.0e33; out_BB.bb_min(1) = 1.0e33;
    out_BB.bb_max(0) = -1.0e33; out_BB.bb_max(1) = -1.0e33;

    for( int i = 0; i < m_Vertices.cols(); i++ )
    {
      const Eigen::Vector2d p = getCurrentPosition( m_Vertices.col(i), in_scale, in_theta, in_center_current );
      out_BB.bb_min(0) = std::min<double>( out_BB.bb_min(0), p(0) );
      out_BB.bb_min(1) = std::min<double>( out_BB.bb_min(1), p(1) );
      out_BB.bb_max(0) = std::max<double>( out_BB.bb_max(0), p(0) );
      out_BB.bb_max(1) = std::max<double>( out_BB.bb_max(1), p(1) );
    }
  }
  
  double computeArea() const
  {
    double area = 0.0;
    for( int i=0; i<m_Vertices.cols(); i++ )
    {
      int ip = ( i+1 ) % m_Vertices.cols();
      int im = ( i-1+m_Vertices.cols() ) % m_Vertices.cols();
      
      area += m_Vertices.col(i)(0) * ( m_Vertices.col(ip)(1) - m_Vertices.col(im)(1) );
    }
    return area * 0.5;
  }
  
  double computeSecondMomentOfArea() const
  {
    double _I = 0.0;
    for( int i=0; i<m_Vertices.cols(); i++ )
    {
      int ip = (i+1) % m_Vertices.cols();
      int im = (i-1+m_Vertices.cols()) % m_Vertices.cols();
      
      const double _Iy_x12 = ( m_Vertices.col(i)(0) * m_Vertices.col(ip)(1) - m_Vertices.col(ip)(0) * m_Vertices.col(i)(1) ) * ( m_Vertices.col(i)(0) * m_Vertices.col(i)(0) + m_Vertices.col(i)(0) * m_Vertices.col(ip)(0) + m_Vertices.col(ip)(0) * m_Vertices.col(ip)(0) );
      
      const double _Ix_x12 = ( m_Vertices.col(i)(0) * m_Vertices.col(ip)(1) - m_Vertices.col(ip)(0) * m_Vertices.col(i)(1) ) * ( m_Vertices.col(i)(1) * m_Vertices.col(i)(1) + m_Vertices.col(i)(1) * m_Vertices.col(ip)(1) + m_Vertices.col(ip)(1) * m_Vertices.col(ip)(1) );
      
      _I += _Ix_x12 + _Iy_x12;
    }
    return _I / 12.0;
  }
  
  void generateCollisionSamples( std::vector< CollisionSample2D >& out_CollisionSamples, const double in_dx_sample_points ) const
  {
    out_CollisionSamples.clear();
    for( int k = 0; k < m_Vertices.cols(); k++ )
    {
      const int kp = ( k + 1 ) % m_Vertices.cols();
      const Eigen::Vector2d xk0 = m_Vertices.col(k);
      const Eigen::Vector2d xkp0 = m_Vertices.col(kp);
      const double length = ( xkp0 - xk0 ).norm();
      
      const int nSegs = std::max<int>( 2, int( ceil( length / in_dx_sample_points ) ) );
      
      for( int i=0; i<nSegs; i++ )
      {
        const Eigen::Vector2d p = xk0 + ( xkp0 - xk0 ) * ( double(i) / double(nSegs) );
        CollisionSample2D sample(p);
        out_CollisionSamples.push_back( sample );
      }
    }
  }
  
  int numVertices() const { return m_Vertices.cols(); }
  
  Eigen::Vector2d getVertexPosition0( int k ) const { return m_Vertices.col( k ); }

private:
  void setCenterOfMassToOrigin()
  {
    Eigen::Vector2d _x0 = Eigen::Vector2d::Zero();
    for( int i = 0; i < m_Vertices.cols(); i++ )
    {
      _x0 += m_Vertices.col(i);
    }
    _x0 = _x0 / m_Vertices.cols();
    
    for( int i = 0; i < m_Vertices.cols(); i++ )
    {
      m_Vertices.col(i) -= _x0;
    }
  }
  
  double computeThetaFromSinCos( const Eigen::Vector2d& in_v, const Eigen::Vector2d& in_x0, const Eigen::Vector2d& in_x1 ) const
  {
    double numeratorSin = 0.0;
    double numeratorCos = 0.0;
    double denominator = 0.0;

    numeratorSin = cross2d( in_x0 - in_v, in_x1 - in_v );
    denominator = ( in_x0 - in_v ).norm() * ( in_x1 - in_v ).norm();
    double sinTheta = numeratorSin / denominator;

    numeratorCos = ( in_x0 - in_v ).dot( in_x1 - in_v );
    double cosTheta = numeratorCos / denominator;

    return atan2( sinTheta, cosTheta );
  }
  
  double computeDistanceToSegment( const Eigen::Vector2d& in_Point, const Eigen::Vector2d& in_SegmentFrom, const Eigen::Vector2d& in_SegmentTo ) const
  {
    const double t = ( in_SegmentTo - in_SegmentFrom ).dot( in_Point - in_SegmentFrom ) / ( in_SegmentTo - in_SegmentFrom ).squaredNorm();
    const Eigen::Vector2d Pt = in_SegmentFrom + ( in_SegmentTo - in_SegmentFrom ) * t;
    const double distance = ( Pt - in_Point ).norm();

    if( t < 0 )
      return ( in_SegmentFrom - in_Point ).norm();
    else if( t > 1 )
      return ( in_SegmentTo - in_Point ).norm();
    else
      return distance;
  }
  
  
  double computeDistanceToSegment( const Eigen::Vector2d& in_Point, const Eigen::Vector2d& in_SegmentFrom, const Eigen::Vector2d& in_SegmentTo, Eigen::Vector2d& out_ClosestPoint ) const
  {
    const double t = ( in_SegmentTo - in_SegmentFrom ).dot( in_Point - in_SegmentFrom ) / ( in_SegmentTo - in_SegmentFrom ).squaredNorm();
    const Eigen::Vector2d Pt = in_SegmentFrom + ( in_SegmentTo - in_SegmentFrom ) * t;
    const double distance = ( Pt - in_Point ).norm();

    if( t < 0 )
    {
      out_ClosestPoint = in_SegmentFrom;
      return ( in_SegmentFrom - in_Point ).norm();
    }
    else if( t > 1 )
    {
      out_ClosestPoint = in_SegmentTo;
      return ( in_SegmentTo - in_Point ).norm();
    }
    else
    {
      out_ClosestPoint = Pt;
      return distance;
    }
  }
  
  Eigen::Vector2d computeClosestPoint( const Eigen::Vector2d& in_x0 ) const
  {
    double closest_dist = 1.0e33;
    Eigen::Vector2d res;
    for( int i=0; i<m_Vertices.cols(); i++ )
    {
      const int ip = ( i + 1 ) % m_Vertices.cols();
      Eigen::Vector2d cp;
      const double _dist = computeDistanceToSegment( in_x0, m_Vertices.col(i), m_Vertices.col(ip), cp );
      if( _dist < closest_dist )
      {
        closest_dist = _dist;
        res = cp;
      }
    }
    
    return res;
  }
  
  std::string m_Name;
  double m_Density;
  Eigen::Matrix2Xd m_Vertices;
  double m_SizeStd;
  
  const SignedDistanceFunction2D* m_SDF;
};

#endif
