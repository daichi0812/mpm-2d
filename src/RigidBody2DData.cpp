//
//  RigidBody2DData.cpp
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/02.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE

#include "RigidBody2DData.h"
#include "HDF5io.h"

RigidBody2DData::RigidBody2DData()
{

}

RigidBody2DData::~RigidBody2DData()
{

}

void RigidBody2DData::clearData()
{
  for( int i=0; i<m_Templates.size(); i++ )
  {
    delete m_Templates[i];
    m_Templates[i] = nullptr;
  }
  m_Templates.clear();

  m_RawTemplates.clear();
  m_TemplateNameDict.clear();
  m_Elements.clear();
}

void RigidBody2DData::deserializeData( const std::string& in_template_fn, const std::string& in_element_fn )
{
  deserializeTemplates( in_template_fn );
  deserializeElements( in_element_fn );
}

void RigidBody2DData::serializeData( const std::string& in_template_fn, const std::string& in_element_fn ) const
{
  serializeTemplates( in_template_fn );
  serializeElements( in_element_fn );
}

void RigidBody2DData::deserializeTemplates( const std::string& in_template_fn )
{
  HDF5File h5( in_template_fn, HDF5AccessType::READ_ONLY );
  m_RawTemplates.clear();
  m_TemplateNameDict.clear();
  m_Templates.clear();

  std::vector<std::string> group_names;
  h5.enumerateGroups( "shape_templates_2d", group_names );

  for( auto s: group_names )
  {
    std::string shape_template_name = std::string( "shape_templates_2d" ) + "/" + s;
    RawShapeTemplateData st;
    h5.readScalar( shape_template_name, "density", st.density );
    h5.readScalar( shape_template_name, "size_mean", st.size_mean );
    h5.readScalar( shape_template_name, "size_std", st.size_std );
    h5.readMatrix( shape_template_name, "vertex_list", st.vertex_list );

    m_RawTemplates.insert( std::pair<std::string, RawShapeTemplateData>( s, st ) );
  }
}

void RigidBody2DData::deserializeElements( const std::string& in_element_fn, const std::string& in_time_step )
{
  HDF5File h5( in_element_fn, HDF5AccessType::READ_ONLY );
  deserializeElements( h5, in_time_step );
}

void RigidBody2DData::deserializeElements( HDF5File& h5, const std::string& in_time_step ) {
  m_Elements.clear();
  
  std::string group_name_progress = in_time_step + "/progress_data";
  if( h5.doesGroupExist( group_name_progress ) )
  {
    h5.readScalar( group_name_progress, "dt", m_ProgressData.dt );
    h5.readScalar( group_name_progress, "time", m_ProgressData.time );
    h5.readScalar( group_name_progress, "tick_count", m_ProgressData.tick_count );
    h5.readScalar( group_name_progress, "serialization_idx", m_ProgressData.serialization_idx );
  }
  else
  {
    m_ProgressData.dt = -1.0;
    m_ProgressData.time = 0.0;
    m_ProgressData.tick_count = 0;
    m_ProgressData.serialization_idx = -1;
  }

  std::vector<std::string> template_names;
  std::string group_name_template = in_time_step + "/template_name_dict";
  for( int i=0; ; i++ )
  {
    try
    {
      std::string dataset_name = std::to_string(i);
      if( !h5.doesDatasetExist( group_name_template, dataset_name ) )
        break;
      std::string data_str;
      h5.readString( group_name_template, dataset_name, data_str );
      template_names.push_back( data_str );
    }
    catch( const std::string& e )
    {
      std::cout << "error: " << e << std::endl;
      break;
    }
  }

  setupTemplates( template_names );

  Eigen::VectorXd angular_velocity_array;
  Eigen::Matrix2Xd center_of_mass_array;
  Eigen::VectorXd rotation_angle_array;
  Eigen::VectorXd size_ratio_array;
  Eigen::Matrix<char, Eigen::Dynamic, 1> isStatic_array;
  Eigen::VectorXi template_idx;
  Eigen::Matrix2Xd velocity_array;

  std::string group_name_elements = in_time_step + "/elements_2d";
  try
  {
    h5.readMatrix( group_name_elements, "angular_velocity", angular_velocity_array );
    h5.readMatrix( group_name_elements, "center_of_mass", center_of_mass_array );
    h5.readMatrix( group_name_elements, "rotation_angle", rotation_angle_array );
    h5.readMatrix( group_name_elements, "size_ratio", size_ratio_array );
    h5.readMatrix( group_name_elements, "static", isStatic_array );
    h5.readMatrix( group_name_elements, "template_idx", template_idx );
    h5.readMatrix( group_name_elements, "velocity", velocity_array );
  }
  catch( const std::string& e )
  {
    std::cout << "error: " << e << std::endl;
    return;
  }

  const int num_elems = template_idx.size();
  m_Elements.resize( num_elems );
  
  for( int i=0; i<num_elems; i++ )
  {
    m_Elements[i].template_idx = template_idx(i);
    m_Elements[i].size_ratio = size_ratio_array(i);
    m_Elements[i].center_of_mass = center_of_mass_array.col(i);
    m_Elements[i].rotation_angle = rotation_angle_array(i);
    m_Elements[i].velocity = velocity_array.col(i);
    m_Elements[i].angular_velocity = angular_velocity_array(i);
    m_Elements[i].is_static = isStatic_array(i) == false ? 0 : 1;
  }
}

void RigidBody2DData::deserializeInterpolatedData( HDF5File& h5, const std::string& in_time_step ) {
  homogenizeDataArray.clear();
  
  Eigen::Matrix2Xd sigma_array;
  Eigen::Matrix2Xd grid_p_array;
  Eigen::Matrix2Xi resolution_array;

  std::string group_name_elements = in_time_step + "/homogenization";
  try
  {
    h5.readSquareMatrixArray( group_name_elements, "sigma", sigma_array );
    h5.readMatrix( group_name_elements, "grid_p", grid_p_array );
    h5.readMatrix( group_name_elements, "resolution", resolution_array );
  }
  catch( const std::string& e )
  {
    std::cout << "error: " << e << std::endl;
    return;
  }
  
  const int num_elems = resolution_array.cols();
  homogenizeDataArray.resize( num_elems );
  
  for( int i=0; i<num_elems; i++ )
  {
    homogenizeDataArray[i].sigma.col(0) = sigma_array.col(i);
    homogenizeDataArray[i].sigma.col(1) = sigma_array.col(num_elems + i);
  }
  
  DEM_resolution = resolution_array.col(0);
  grid_start = grid_p_array.col(0);
}

void RigidBody2DData::serializeTemplates( const std::string& in_template_fn ) const
{
  HDF5File h5( in_template_fn, HDF5AccessType::READ_WRITE );

  for( auto t: m_Templates )
  {
    ShapeTemplateData data;
    t->getData( data );

    std::string shape_template_name = std::string( "shape_templates_2d" ) + "/" + data.name;
    h5.writeScalar( shape_template_name, "density", data.density );
    h5.writeScalar( shape_template_name, "size_mean", data.size_mean );
    h5.writeScalar( shape_template_name, "size_std", data.size_std );
    h5.writeMatrix( shape_template_name, "vertex_list", data.vertex_list );
  }
}

void RigidBody2DData::_serializeElements( HDF5File& h5, const std::string& in_time_step ) const
{
  const std::string group_name_progress = in_time_step + "/progress_data";
  h5.writeScalar( group_name_progress, "dt", m_ProgressData.dt );
  h5.writeScalar( group_name_progress, "time", m_ProgressData.time );
  h5.writeScalar( group_name_progress, "tick_count", m_ProgressData.tick_count );
  h5.writeScalar( group_name_progress, "serialization_idx", m_ProgressData.serialization_idx );

  for( auto d: m_TemplateNameDict )
  {
    std::string dataset_name = std::to_string( d.second );
    h5.writeString( in_time_step + "/template_name_dict", dataset_name, d.first );
  }

  Eigen::VectorXd angular_velocity_array; angular_velocity_array.resize( m_Elements.size() );
  Eigen::Matrix2Xd center_of_mass_array; center_of_mass_array.resize( 2, m_Elements.size() );
  Eigen::VectorXd rotation_angle_array; rotation_angle_array.resize( m_Elements.size() );
  Eigen::VectorXd size_ratio_array; size_ratio_array.resize( m_Elements.size() );
  Eigen::Matrix<char, Eigen::Dynamic, 1> isStatic_array; isStatic_array.resize( m_Elements.size(), 1 );
  Eigen::VectorXi template_idx; template_idx.resize( m_Elements.size() );
  Eigen::Matrix2Xd velocity_array; velocity_array.resize( 2, m_Elements.size() );

  const int num_elems = m_Elements.size();
  for( int i=0; i<num_elems; i++ )
  {
    angular_velocity_array(i) = m_Elements[i].angular_velocity;
    center_of_mass_array.col(i) = m_Elements[i].center_of_mass;
    rotation_angle_array(i) = m_Elements[i].rotation_angle;
    size_ratio_array(i) = m_Elements[i].size_ratio;
    isStatic_array(i) = m_Elements[i].is_static ? 1 : 0;
    template_idx(i) = m_Elements[i].template_idx;
    velocity_array.col(i) = m_Elements[i].velocity;
  }

  const std::string group_name_elements = in_time_step + "/elements_2d";
  h5.writeVector( group_name_elements, "angular_velocity", angular_velocity_array );
  h5.writeMatrix( group_name_elements, "center_of_mass", center_of_mass_array );
  h5.writeVector( group_name_elements, "rotation_angle", rotation_angle_array );
  h5.writeVector( group_name_elements, "size_ratio", size_ratio_array );
  h5.writeVector( group_name_elements, "static", isStatic_array );
  h5.writeVector( group_name_elements, "template_idx", template_idx );
  h5.writeMatrix( group_name_elements, "velocity", velocity_array );
}

void RigidBody2DData::serializeElements( const std::string& in_element_fn ) const
{
  HDF5File h5( in_element_fn, HDF5AccessType::READ_WRITE );
  _serializeElements( h5 );
}

void RigidBody2DData::setupTemplates( const std::vector<std::string>& in_template_names )
{
  m_Templates.clear();
  m_TemplateNameDict.clear();

  for( int i=0; i<in_template_names.size(); i++ )
  {
    auto p = m_RawTemplates.find( in_template_names[i] );
    if( p == m_RawTemplates.end() )
    {
      std::cout << "Error: cannot find template name " << in_template_names[i] << " in setupTemplates()" << std::endl;
      exit(-1);
    }

    ShapeTemplateData st_data;
    st_data.name = in_template_names[i];
    st_data.density = p->second.density;
  	st_data.size_mean = p->second.size_mean;
  	st_data.size_std = p->second.size_std;
  	st_data.vertex_list = p->second.vertex_list;
    ShapeTemplate2D* st = initializeTemplate( st_data );

    m_Templates.push_back( st );

    m_TemplateNameDict.insert( std::pair<std::string, int>( in_template_names[i], i ) );
  }
}

void RigidBody2DData::printData() const
{
  std::cout << "Templates: " << std::endl;
  for( auto t: m_Templates )
  {
    ShapeTemplateData data;
    t->getData( data );

    std::cout << data.name << std::endl;
    std::cout << "  density: " << data.density << std::endl;
    std::cout << "  size mean: " << data.size_mean << std::endl;
    std::cout << "  size std: " << data.size_std << std::endl;
    std::cout << "  vertex list: " << data.vertex_list << std::endl;
  }

  std::cout << "Elements: " << std::endl;
  for( auto e: m_Elements )
  {
    std::cout << "<E>" << std::endl;
    std::cout << "  template_name: " << m_Templates[e.template_idx]->getName() << std::endl;
    std::cout << "  size ratio: " << e.size_ratio << std::endl;
    std::cout << "  center of mass: " << e.center_of_mass << std::endl;
    std::cout << "  rotation angle: " << e.rotation_angle << std::endl;
    std::cout << "  velocity: " << e.velocity << std::endl;
    std::cout << "  angular velocity: " << e.angular_velocity << std::endl;
    std::cout << "  is static: " << e.is_static << std::endl;
  }
}
