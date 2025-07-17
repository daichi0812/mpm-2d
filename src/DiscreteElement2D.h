//
//  DiscreteElement2D.h
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/08.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#ifndef __DISCRETE_ELEMENT_2D_H__
#define __DISCRETE_ELEMENT_2D_H__

#include <Eigen/Core>
#include <vector>
#include <map>
#include <iostream>
#include "constants.h"
#include "collisionsample2d.h"
#include "ShapeTemplate2D.h"

class SignedDistanceFunction2D;

struct Element
{
	Element()
	: template_idx( -1 ), size_ratio( 0.1 ), center_of_mass( Eigen::Vector2d::Zero() )
	, rotation_angle( 0.0 ), velocity( Eigen::Vector2d::Zero() ), angular_velocity( 0.0 )
	, is_static( true )
	{}
  
  void initialize( const int in_index, const ShapeTemplate2D* in_template_ptr)
  {
    index = in_index;
    template_ptr = in_template_ptr;
    clearForceAndTorque();
    
    density = template_ptr-> getDensity();
    mass = template_ptr->computeArea() * density * size_ratio * size_ratio;
    
    /*
    std::cout << "initializing idx [" << in_index << "]:" << std::endl;
    std::cout << "  type: " << ( in_template_ptr->type() == CIRCLE ? "circle" : "polygon" ) << std::endl;
    std::cout << "  size_ratio: " << size_ratio << std::endl;
    std::cout << "  density: " << template_ptr->getDensity() << std::endl;
    std::cout << "  area: " << template_ptr->computeArea() << std::endl;
    std::cout << "  mass: " << mass << std::endl;
    std::cout << "  2nd moment: " << template_ptr->computeSecondMomentOfArea() << std::endl;
    std::cout << "  inertia: " << inertia << std::endl;
    //*/
  }
  
  void getBoundingBox( BoundingBox2D& out_BB ) const
  {
    template_ptr->getBoundingBox( out_BB, size_ratio, rotation_angle, center_of_mass );
  }
  
  // For signed distance function:
  double getSignedDistanceAt( const Eigen::Vector2d& in_p_current ) const
  {
    const Eigen::Vector2d p0 = getMaterialPosition( in_p_current, size_ratio, rotation_angle, center_of_mass );
    return template_ptr->getSignedDistanceAt( p0 ) * size_ratio;
  }
  
  Eigen::Vector2d getNormalAt( const Eigen::Vector2d& in_p_current ) const
  {
    const Eigen::Matrix2d mat{ rotationMat( rotation_angle ) };
    const Eigen::Vector2d p0 = getMaterialPosition( in_p_current, size_ratio, rotation_angle, center_of_mass );
    
    return mat * template_ptr->getNormalAt( p0 );
  }
  
  // For collisions:
  int numCollisionSamples() const { return collision_samples.size(); }
  const CollisionSample2D& collisionSample( int in_Idx ) const { return collision_samples[ in_Idx ]; }
  CollisionSample2D& collisionSample( int in_Idx ) { return collision_samples[ in_Idx ]; }

  // For dynamics:
  int getIndex() const { return index; }
  int getTemplateIndex() const { return template_idx; }
  const ShapeTemplate2D* getTemplatePtr() const { return template_ptr; }
  
  bool getIsStatic() const { return is_static; }
  
  double getMass() const { return mass; }
  double getInertia() const { return inertia; }
  Eigen::Vector2d getCenterOfMass() const { return center_of_mass; }
  Eigen::Vector2d getVelocity() const { return velocity; }
  double getRotationAngle() const { return rotation_angle; }
  double getAngularVelocity() const { return angular_velocity; }

  void clearForceAndTorque()
  {
    force.setZero();
    torque = 0.0;
  }
  
  void accumulateForceAndTorque( const Eigen::Vector2d& in_Force, const double in_Torque )
  {
    force += in_Force;
    torque += in_Torque;
  }
  
  void stepSymplecticEuler( const double in_dt )
  {
    velocity += force * in_dt / mass;
    angular_velocity += torque * in_dt / inertia;
    
    center_of_mass += velocity * in_dt;
    rotation_angle += angular_velocity * in_dt;
  }

  // For mapping reference frame to current frame:
  Eigen::Vector2d getCurrentPosition( const Eigen::Vector2d& in_p0 ) const { return ::getCurrentPosition( in_p0, size_ratio, rotation_angle, center_of_mass ); }
  double scale() const { return size_ratio; }

	int template_idx;
	double size_ratio;
	Eigen::Vector2d center_of_mass;
	double rotation_angle;
	Eigen::Vector2d velocity;
	double angular_velocity;
	bool is_static;
  
  int index;
  const ShapeTemplate2D* template_ptr;
  double mass;
  double inertia;
  Eigen::Vector2d force;
  double torque;
  double density;
  
  std::vector< CollisionSample2D > collision_samples;
};

#endif
