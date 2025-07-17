//
//  RigidBody2DData.h
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/02.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#ifndef __RIGID_BODY_2D_DATA_H__
#define __RIGID_BODY_2D_DATA_H__

#include <Eigen/Core>
#include <vector>
#include <map>
#include <iostream>
#include "constants.h"
#include "collisionsample2d.h"
#include "DiscreteElement2D.h"
#include "ShapeTemplate2D.h"
#include "HDF5io.h"

class SignedDistanceFunction2D;

struct ProgressData
{
  ProgressData()
  : dt( -1.0 ), time( 0.0 ), tick_count( 0 ), serialization_idx( -1 ) {}
  
  double dt;
  double time;
  long long tick_count;
  int serialization_idx;
};

struct HomogenizeData
{
  HomogenizeData()
  : sigma( Eigen::Matrix2d::Zero() ), strain( Eigen::Matrix2d::Zero() ) {}
  
  Eigen::Matrix2d sigma;
  Eigen::Matrix2d strain;
};

class RigidBody2DData
{
private:
  struct RawShapeTemplateData
  {
  	RawShapeTemplateData()
  	: density( 1.0 ), size_mean( 0.1 ), size_std( 0.01 )
  	{ vertex_list.setZero(2, 1); }

  	double density;
  	double size_mean;
  	double size_std;
  	Eigen::Matrix2Xd vertex_list;
  };

public:
  RigidBody2DData();
  ~RigidBody2DData();
  
  void clearData();

	void deserializeData( const std::string& in_template_fn, const std::string& in_element_fn );
	void serializeData( const std::string& in_template_fn, const std::string& in_element_fn ) const;
	void printData() const;

  void deserializeTemplates( const std::string& in_template_fn );
	void deserializeElements( const std::string& in_element_fn, const std::string& in_time_step = "" );
  void deserializeElements( HDF5File& h5, const std::string& in_time_step = "");
  void deserializeInterpolatedData( HDF5File& h5, const std::string& in_time_step = "");
	void serializeTemplates( const std::string& in_template_fn ) const;
	void serializeElements( const std::string& in_element_fn ) const;
 
  void serializeForces( HDF5File& h5, const std::string& in_time_step = "");
  
  inline int numTemplates() const { return m_Templates.size(); }
  inline int numElements() const { return m_Elements.size(); }
  inline int numHomogenizeData() const { return homogenizeDataArray.size(); }

  inline const ShapeTemplate2D* getTemplate( const int idx ) const { return m_Templates[idx]; }
  inline ShapeTemplate2D* getTemplate( const int idx ){ return m_Templates[idx]; }
  inline const Element& getElement( const int idx ) const { return m_Elements[idx]; }
  inline Element& getElement( const int idx ){ return m_Elements[idx]; }
  
  inline const ProgressData& getProgressData() const { return m_ProgressData; }
  inline ProgressData& getProgressData() { return m_ProgressData; }
  
  inline const std::vector<Element>& getElements() const { return m_Elements; }
  inline std::vector<Element>& getElements() { return m_Elements; }
  
  inline const HomogenizeData& getHomogenizeData( const int idx ) const { return homogenizeDataArray[idx]; }
  inline HomogenizeData& getHomogenizeData( const int idx ){ return homogenizeDataArray[idx]; }
  
  inline const Eigen::Vector2i& getDemResolution() const { return DEM_resolution; }
  inline Eigen::Vector2i& getDemResolution() { return DEM_resolution; }
  
  inline const Eigen::Vector2d& getGridStart() const { return grid_start; }
  inline Eigen::Vector2d& getGridStart() { return grid_start; }

private:
  void setupTemplates( const std::vector<std::string>& in_template_names );
  void _serializeElements( HDF5File& h5, const std::string& in_time_step = "" ) const;

	std::map<std::string, RawShapeTemplateData> m_RawTemplates;
  std::map<std::string, int> m_TemplateNameDict;
  std::vector<ShapeTemplate2D*> m_Templates;
  std::vector<Element> m_Elements;
  
  Eigen::Vector2d grid_start;
  Eigen::Vector2i DEM_resolution;
  std::vector<HomogenizeData> homogenizeDataArray;
  
  ProgressData m_ProgressData;
};

#endif
