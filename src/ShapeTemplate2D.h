//
//  ShapeTemplate2D.h
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/08.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#ifndef __SHAPE_TEMPLATE_2D_H__
#define __SHAPE_TEMPLATE_2D_H__

#include <Eigen/Core>
#include <vector>
#include <map>
#include <iostream>
#include "constants.h"
#include "collisionsample2d.h"

class SignedDistanceFunction2D;

struct ShapeTemplateData
{
  std::string name;
  double density;
  double size_mean;
  double size_std;
  Eigen::Matrix2Xd vertex_list;
};

class ShapeTemplate2D
{
public:
  virtual ~ShapeTemplate2D() {};
  virtual RigidBodyType2D type() const = 0;

  virtual void getData( ShapeTemplateData& out_data ) const = 0;
  virtual const std::string& getName() const = 0;
  virtual double getDensity() const = 0;
  
  virtual double computeWindingNumber( const Eigen::Vector2d& in_x0 ) const = 0;
  virtual double computeMinimumDistance( const Eigen::Vector2d& in_x0 ) const = 0;
  virtual double getSignedDistanceAt( const Eigen::Vector2d& in_x0 ) const = 0;
  virtual Eigen::Vector2d getNormalAt( const Eigen::Vector2d& in_x0 ) const = 0;
  
  virtual void computeSignedDistanceFunction( const double dx ) = 0;
  virtual const SignedDistanceFunction2D* getSDF() const = 0;
  
  virtual void getBoundingBox( BoundingBox2D& out_BB, const double in_scale, const double in_theta, const Eigen::Vector2d& in_center_current ) const = 0;
  
  virtual double computeArea() const = 0;
  virtual double computeSecondMomentOfArea() const = 0;
  
  virtual void generateCollisionSamples( std::vector< CollisionSample2D >& out_CollisionSamples, const double in_dx_sample_points ) const = 0;
};

ShapeTemplate2D* initializeTemplate( const ShapeTemplateData& in_data );

#endif
