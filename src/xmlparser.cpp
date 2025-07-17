//
//  xmlparser.cpp
//  SimpleMPM2D
//
//  Created by Yonghao Yue on 2020/09/04.
//

#include "xmlparser.h"
#include <vector>
#include <iostream>
#include <fstream>

void setDefaultSetting( Setting& setting )
{
  setting.gridRegion << 0.0, 0.0, 1.0, 1.0 ; //minGridX, minGridY, maxGridX, maxGridY
  setting.gridRes << 128, 128; //gridResX, gridResY
  setting.particleCountPerCellEdge = 2; //particle per cell *edge* ( = sqrt particle per cell in 2D )
  setting.dt = 0.00001;
  setting.valpha = 0.95;
  setting.maxSimulationTime = 10.0;
  
  setting.capture = false;
  setting.captureInterval = 0.0333333;
  
  setting.obstacleMapFileName = "";
  
  setting.initialParticlesFileName = "";
  setting.materialDensity = 400.0;
  
  setting.xi = 10.0;
  setting.theta_c = 0.025;
  setting.theta_s = 0.0075;
  setting.YoungModulus = 140000.0;
  setting.PoissonRatio = 0.2;
  
  setting.fric = 0.1;
  setting.g = 9.8;
  
  setting.templates_file_name_to_resume = "";
  setting.forces_file_name_to_resume = "";
  setting.homogenization_file_name = "";
  setting.serialization_folder = "Save";
  setting.forces_file_name_for_serialization = "";
  
  setting.cell_height = 0.04;
  setting.compute_sigma_only = false;
}

void showSettings( const Setting& setting )
{
  std::cout << "[MPM2D] ===== Integrator =====" << std::endl;
  std::cout << "[MPM2D]   grid_min: " << setting.gridRegion.segment<2>(0).transpose() << std::endl;
  std::cout << "[MPM2D]   grid_max: " << setting.gridRegion.segment<2>(2).transpose() << std::endl;
  std::cout << "[MPM2D]   grid_res: " << setting.gridRes.transpose() << std::endl;
  std::cout << "[MPM2D]   particle_count_per_cell_edge: " << setting.particleCountPerCellEdge << std::endl;
  std::cout << "[MPM2D]   dt: " << setting.dt << std::endl;
  std::cout << "[MPM2D]   valpha: " << setting.valpha << std::endl;
  std::cout << "[MPM2D]   max_simulation_time: " << setting.maxSimulationTime << std::endl;
  
  std::cout << "[MPM2D] ===== Capture =====" << std::endl;
  std::cout << "[MPM2D]   capture: " << setting.capture << std::endl;
  std::cout << "[MPM2D]   interval: " << setting.captureInterval << std::endl;
  
  std::cout << "[MPM2D] ===== Obstacle =====" << std::endl;
  std::cout << "[MPM2D]   filename: " << setting.obstacleMapFileName << std::endl;
  
  std::cout << "[MPM2D] ===== Particles =====" << std::endl;
  std::cout << "[MPM2D]   filename: " << setting.initialParticlesFileName << std::endl;
  std::cout << "[MPM2D]   density: " << setting.materialDensity << std::endl;
  
  std::cout << "[MPM2D] ===== Constitutive model =====" << std::endl;
  std::cout << "[MPM2D]   xi: " << setting.xi << std::endl;
  std::cout << "[MPM2D]   theta_c: " << setting.theta_c << std::endl;
  std::cout << "[MPM2D]   theta_s: " << setting.theta_s << std::endl;
  std::cout << "[MPM2D]   youngs_modulus: " << setting.YoungModulus << std::endl;
  std::cout << "[MPM2D]   poissons_ratio: " << setting.PoissonRatio << std::endl;
  
  std::cout << "[MPM2D] ===== External force =====" << std::endl;
  std::cout << "[MPM2D]   friction: " << setting.fric << std::endl;
  std::cout << "[MPM2D]   g: " << setting.g << std::endl;
  
  std::cout << "[MPM2D] ===== input data =====" << std::endl;
  std::cout << "[MPM2D]   template: " << setting.templates_file_name_to_resume << std::endl;
  std::cout << "[MPM2D]   forces: " << setting.forces_file_name_to_resume << std::endl;
  std::cout << "[MPM2D]   homogenization: " << setting.homogenization_file_name << std::endl;
  
  std::cout << "[MPM2D] ===== Serialize data =====" << std::endl;
  std::cout << "[MPM2D]   folder: " << setting.serialization_folder << std::endl;
  std::cout << "[MPM2D]   forces: " << setting.forces_file_name_for_serialization << std::endl;
  std::cout << "[MPM2D]   h: " << setting.cell_height << std::endl;
  std::cout << "[MPM2D]   compute_sigma_only: " << setting.compute_sigma_only << std::endl;
}

template<class T>
bool extractFromString( const std::string& str, T& res )
{
  std::stringstream input_strm( str );
  input_strm >> res;
  return !input_strm.fail();
}

template<class T, int N>
bool extractFromString( const std::string& str, Eigen::Matrix<T, N, 1>& vec )
{
  std::stringstream input_strm( str );
  for( int i=0; i<vec.size(); i++ )
    input_strm >> vec(i);
  return !input_strm.fail();
}

bool loadIntegrator( const rapidxml::xml_node<>& node, Setting& setting )
{
  const rapidxml::xml_node<>* const integrator_node{ node.first_node( "integrator" ) };
  if( integrator_node == nullptr )
  {
    std::cerr << "Failed to locate integrator node." << std::endl;
    return false;
  }
  
  Eigen::Vector2d grid_region_min;
  const rapidxml::xml_attribute<>* const gridMin_attrib{ integrator_node->first_attribute( "grid_min" ) };
  if( !gridMin_attrib )
  {
    std::cerr << "Failed to locate grid_min attribute for integrator node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ gridMin_attrib->value() }, grid_region_min ) )
  {
    std::cerr << "Failed to load grid_min attribute for integrator node. Must provide a two-vector." << std::endl;
    return false;
  }
  setting.gridRegion.segment<2>(0) = grid_region_min;
  
  Eigen::Vector2d grid_region_max;
  const rapidxml::xml_attribute<>* const gridMax_attrib{ integrator_node->first_attribute( "grid_max" ) };
  if( !gridMax_attrib )
  {
    std::cerr << "Failed to locate grid_max attribute for integrator node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ gridMax_attrib->value() }, grid_region_max ) )
  {
    std::cerr << "Failed to load grid_max attribute for integrator node. Must provide a two-vector." << std::endl;
    return false;
  }
  setting.gridRegion.segment<2>(2) = grid_region_max;
  
  const rapidxml::xml_attribute<>* const gridRes_attrib{ integrator_node->first_attribute( "grid_res" ) };
  if( !gridRes_attrib )
  {
    std::cerr << "Failed to locate grid_res attribute for integrator node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ gridRes_attrib->value() }, setting.gridRes ) )
  {
    std::cerr << "Failed to load grid_res attribute for integrator node. Must provide an integer two-vector." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const particleCountPerCellEdge_attrib{ integrator_node->first_attribute( "particle_count_per_cell_edge" ) };
  if( !particleCountPerCellEdge_attrib )
  {
    std::cerr << "Failed to locate particle_count_per_cell_edge attribute for integrator node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ particleCountPerCellEdge_attrib->value() }, setting.particleCountPerCellEdge ) )
  {
    std::cerr << "Failed to load particle_count_per_cell_edge attribute for integrator node. Must provide an integer." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const dt_attrib{ integrator_node->first_attribute( "dt" ) };
  if( !dt_attrib )
  {
    std::cerr << "Failed to locate dt attribute for integrator node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ dt_attrib->value() }, setting.dt ) )
  {
    std::cerr << "Failed to load dt attribute for integrator node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const valpha_attrib{ integrator_node->first_attribute( "valpha" ) };
  if( !valpha_attrib )
  {
    std::cerr << "Failed to locate valpha attribute for integrator node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ valpha_attrib->value() }, setting.valpha ) )
  {
    std::cerr << "Failed to load valpha attribute for integrator node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const maxSimulationTime_attrib{ integrator_node->first_attribute( "max_simulation_time" ) };
  if( !maxSimulationTime_attrib )
  {
    std::cerr << "Failed to locate max_simulation_time attribute for integrator node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ maxSimulationTime_attrib->value() }, setting.maxSimulationTime ) )
  {
    std::cerr << "Failed to load max_simulation_time attribute for integrator node. Must provide a real number." << std::endl;
    return false;
  }
  
  return true;
}

bool loadCapture( const rapidxml::xml_node<>& node, Setting& setting )
{
  const rapidxml::xml_node<>* const capture_node{ node.first_node( "capture" ) };
  if( capture_node == nullptr )
  {
    setting.capture = false;
    return true;
  }
    
  const rapidxml::xml_attribute<>* const captureInterval_attrib{ capture_node->first_attribute( "interval" ) };
  if( !captureInterval_attrib )
  {
    std::cerr << "Failed to locate interval attribute for capture node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ captureInterval_attrib->value() }, setting.captureInterval ) )
  {
    std::cerr << "Failed to load interval attribute for capture node. Must provide a real number." << std::endl;
    return false;
  }
  else
  {
    setting.capture = true;
  }
  
  return true;
}

bool loadObstacle( const rapidxml::xml_node<>& node, Setting& setting )
{
  const rapidxml::xml_node<>* const obstacle_node{ node.first_node( "obstacle" ) };
  if( obstacle_node == nullptr )
  {
    std::cerr << "Failed to locate obstacle node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const filename_attrib{ obstacle_node->first_attribute( "filename" ) };
  if( !filename_attrib )
  {
    std::cerr << "Failed to locate filename attribute for obstacle node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ filename_attrib->value() }, setting.obstacleMapFileName ) )
  {
    std::cerr << "Failed to load filename attribute for obstacle node." << std::endl;
    return false;
  }
  
  return true;
}

bool loadParticles( const rapidxml::xml_node<>& node, Setting& setting )
{
  const rapidxml::xml_node<>* const particles_node{ node.first_node( "particles" ) };
  if( particles_node == nullptr )
  {
    std::cerr << "Failed to locate particles node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const filename_attrib{ particles_node->first_attribute( "filename" ) };
  if( !filename_attrib )
  {
    std::cerr << "Failed to locate filename attribute for particles node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ filename_attrib->value() }, setting.initialParticlesFileName ) )
  {
    std::cerr << "Failed to load filename attribute for particles node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const density_attrib{ particles_node->first_attribute( "density" ) };
  if( !density_attrib )
  {
    std::cerr << "Failed to locate density attribute for particles node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ density_attrib->value() }, setting.materialDensity ) )
  {
    std::cerr << "Failed to load density attribute for particles node. Must provide a real number." << std::endl;
    return false;
  }
  
  return true;
}

bool loadConstitutiveModelSnow( const rapidxml::xml_node<>& constitutiveModel_node, Setting& setting )
{
  const rapidxml::xml_attribute<>* const xi_attrib{ constitutiveModel_node.first_attribute( "xi" ) };
  if( !xi_attrib )
  {
    std::cerr << "Failed to locate xi attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ xi_attrib->value() }, setting.xi ) )
  {
    std::cerr << "Failed to load xi attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const theta_c_attrib{ constitutiveModel_node.first_attribute( "theta_c" ) };
  if( !theta_c_attrib )
  {
    std::cerr << "Failed to locate theta_c attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ theta_c_attrib->value() }, setting.theta_c ) )
  {
    std::cerr << "Failed to load theta_c attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const theta_s_attrib{ constitutiveModel_node.first_attribute( "theta_s" ) };
  if( !theta_s_attrib )
  {
    std::cerr << "Failed to locate theta_s attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ theta_s_attrib->value() }, setting.theta_s ) )
  {
    std::cerr << "Failed to load theta_s attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const youngs_modulus_attrib{ constitutiveModel_node.first_attribute( "youngs_modulus" ) };
  if( !youngs_modulus_attrib )
  {
    std::cerr << "Failed to locate youngs_modulus attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ youngs_modulus_attrib->value() }, setting.YoungModulus ) )
  {
    std::cerr << "Failed to load youngs_modulus attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const poissons_ratio_attrib{ constitutiveModel_node.first_attribute( "poissons_ratio" ) };
  if( !poissons_ratio_attrib )
  {
    std::cerr << "Failed to locate poissons_ratio attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ poissons_ratio_attrib->value() }, setting.PoissonRatio ) )
  {
    std::cerr << "Failed to load poissons_ratio attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  return true;
}

bool loadConstitutiveModelHerschelBulkley( const rapidxml::xml_node<>& constitutiveModel_node, Setting& setting )
{
  const rapidxml::xml_attribute<>* const kappa_attrib{ constitutiveModel_node.first_attribute( "kappa" ) };
  if( !kappa_attrib )
  {
    std::cerr << "Failed to locate kappa attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ kappa_attrib->value() }, setting.kappa ) )
  {
    std::cerr << "Failed to load kappa attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const mu_attrib{ constitutiveModel_node.first_attribute( "mu" ) };
  if( !mu_attrib )
  {
    std::cerr << "Failed to locate mu attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ mu_attrib->value() }, setting.mu ) )
  {
    std::cerr << "Failed to load mu attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const n_attrib{ constitutiveModel_node.first_attribute( "n" ) };
  if( !n_attrib )
  {
    std::cerr << "Failed to locate n attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ n_attrib->value() }, setting.n ) )
  {
    std::cerr << "Failed to load n attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const sigma_y_attrib{ constitutiveModel_node.first_attribute( "sigma_y" ) };
  if( !sigma_y_attrib )
  {
    std::cerr << "Failed to locate sigma_y attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ sigma_y_attrib->value() }, setting.sigma_y ) )
  {
    std::cerr << "Failed to load sigma_y attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const eta_attrib{ constitutiveModel_node.first_attribute( "eta" ) };
  if( !eta_attrib )
  {
    std::cerr << "Failed to locate eta attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ eta_attrib->value() }, setting.eta ) )
  {
    std::cerr << "Failed to load eta attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  return true;
}

bool loadConstitutiveModelThixotropicHerschelBulkley( const rapidxml::xml_node<>& constitutiveModel_node, Setting& setting )
{
  const rapidxml::xml_attribute<>* const kappa_attrib{ constitutiveModel_node.first_attribute( "kappa" ) };
  if( !kappa_attrib )
  {
    std::cerr << "Failed to locate kappa attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ kappa_attrib->value() }, setting.kappa ) )
  {
    std::cerr << "Failed to load kappa attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const mu_attrib{ constitutiveModel_node.first_attribute( "mu" ) };
  if( !mu_attrib )
  {
    std::cerr << "Failed to locate mu attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ mu_attrib->value() }, setting.mu ) )
  {
    std::cerr << "Failed to load mu attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const n_attrib{ constitutiveModel_node.first_attribute( "n" ) };
  if( !n_attrib )
  {
    std::cerr << "Failed to locate n attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ n_attrib->value() }, setting.n ) )
  {
    std::cerr << "Failed to load n attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const sigma_y_attrib{ constitutiveModel_node.first_attribute( "sigma_y" ) };
  if( !sigma_y_attrib )
  {
    std::cerr << "Failed to locate sigma_y attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ sigma_y_attrib->value() }, setting.sigma_y ) )
  {
    std::cerr << "Failed to load sigma_y attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const eta_attrib{ constitutiveModel_node.first_attribute( "eta" ) };
  if( !eta_attrib )
  {
    std::cerr << "Failed to locate eta attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ eta_attrib->value() }, setting.eta ) )
  {
    std::cerr << "Failed to load eta attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const k1_attrib{ constitutiveModel_node.first_attribute( "k1" ) };
  if( !k1_attrib )
  {
    std::cerr << "Failed to locate k1 attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ k1_attrib->value() }, setting.k1 ) )
  {
    std::cerr << "Failed to load k1 attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const k2_attrib{ constitutiveModel_node.first_attribute( "k2" ) };
  if( !k2_attrib )
  {
    std::cerr << "Failed to locate k2 attribute for constitutive_model node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ k2_attrib->value() }, setting.k2 ) )
  {
    std::cerr << "Failed to load k2 attribute for constitutive_model node. Must provide a real number." << std::endl;
    return false;
  }
  
  return true;
}
  
bool loadConstitutiveModel( const rapidxml::xml_node<>& node, Setting& setting )
{
  const rapidxml::xml_node<>* const constitutiveModel_node{ node.first_node( "constitutive_model" ) };
  if( constitutiveModel_node == nullptr )
  {
    std::cerr << "Failed to locate constitutive_model node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const type_attrib{ constitutiveModel_node->first_attribute( "type" ) };
  if( !type_attrib )
  {
    std::cerr << "Failed to locate type attribute for constitutive_model node." << std::endl;
    return false;
  }
  
  if( std::string{ type_attrib->value() }.compare( "Snow" ) == 0 )
  {
    if( !loadConstitutiveModelSnow( *constitutiveModel_node, setting ) )
    {
      std::cerr << "Failed to load snow constitutive model" << std::endl;
      return false;
    }
    setting.constitutiveModelType = CMT_SNOW;
  }
  else if( std::string{ type_attrib->value() }.compare( "Herschel_Bulkley" ) == 0 )
  {
    if( !loadConstitutiveModelHerschelBulkley( *constitutiveModel_node, setting ) )
    {
      std::cerr << "Failed to load Herschel_Bulkley constitutive model" << std::endl;
      return false;
    }
    setting.constitutiveModelType = CMT_HERSCHEL_BULKLEY;
  }
  else if( std::string{ type_attrib->value() }.compare( "Thixotropic_Herschel_Bulkley" ) == 0 )
  {
    if( !loadConstitutiveModelThixotropicHerschelBulkley( *constitutiveModel_node, setting ) )
    {
      std::cerr << "Failed to load Thixotropic_Herschel_Bulkley constitutive model" << std::endl;
      return false;
    }
    setting.constitutiveModelType = CMT_THIXOTROPIC_HERSCHEL_BULKLEY;
  }
  else
  {
    std::cerr << "Unknown constitutive model type" << std::endl;
    return false;
  }
    
  return true;
}

bool loadExternalForce( const rapidxml::xml_node<>& node, Setting& setting )
{
  const rapidxml::xml_node<>* const externalForce_node{ node.first_node( "external_force" ) };
  if( externalForce_node == nullptr )
  {
    std::cerr << "Failed to locate external_force node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const friction_attrib{ externalForce_node->first_attribute( "friction" ) };
  if( !friction_attrib )
  {
    std::cerr << "Failed to locate friction attribute for external_force node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ friction_attrib->value() }, setting.fric ) )
  {
    std::cerr << "Failed to load friction attribute for external_force node. Must provide a real number." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const g_attrib{ externalForce_node->first_attribute( "g" ) };
  if( !g_attrib )
  {
    std::cerr << "Failed to locate g attribute for external_force node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ g_attrib->value() }, setting.g ) )
  {
    std::cerr << "Failed to load g attribute for external_force node. Must provide a real number." << std::endl;
    return false;
  }
  
  return true;
}

bool loadResume( const rapidxml::xml_node<>& node, Setting& setting )
{
  const rapidxml::xml_node<>* const resume_node{ node.first_node( "resume" ) };
  if( resume_node == nullptr )
  {
    std::cerr << "Failed to locate resume node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const templates_attrib{ resume_node->first_attribute( "templates" ) };
  if( !templates_attrib )
  {
    std::cerr << "Failed to locate templates attribute for resume node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ templates_attrib->value() }, setting.templates_file_name_to_resume ) )
  {
    std::cerr << "Failed to load templates attribute for resume node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const forces_attrib{ resume_node->first_attribute( "forces" ) };
  if( !forces_attrib )
  {
    std::cerr << "Failed to locate forces attribute for resume node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ forces_attrib->value() }, setting.forces_file_name_to_resume ) )
  {
    std::cerr << "Failed to load forces attribute for resume node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const homogenization_attrib{ resume_node->first_attribute( "homogenization" ) };
  if( !homogenization_attrib )
  {
    std::cerr << "Failed to locate homogenization attribute for resume node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ homogenization_attrib->value() }, setting.homogenization_file_name ) )
  {
    std::cerr << "Failed to load homogenization attribute for resume node." << std::endl;
    return false;
  }
  
  return true;
}

bool loadSerialization( const rapidxml::xml_node<>& node, Setting& setting )
{
  const rapidxml::xml_node<>* const serialization_node{ node.first_node( "serialization" ) };
  if( serialization_node == nullptr )
  {
    std::cerr << "Failed to locate resume node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const folder_attrib{ serialization_node->first_attribute( "folder" ) };
  if( !folder_attrib )
  {
    std::cerr << "Failed to locate folder attribute for serialization node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ folder_attrib->value() }, setting.serialization_folder ) )
  {
    std::cerr << "Failed to load folder attribute for serialization node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const forces_attrib{ serialization_node->first_attribute( "forces" ) };
  if( !forces_attrib )
  {
    std::cerr << "Failed to locate forces attribute for serialization node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ forces_attrib->value() }, setting.forces_file_name_for_serialization ) )
  {
    std::cerr << "Failed to load forces attribute for serialization node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const cellHeight_attrib{ serialization_node->first_attribute( "cell_height" ) };
  if( !cellHeight_attrib )
  {
    std::cerr << "Failed to locate cell_height attribute for serialization node." << std::endl;
    return false;
  }
  if( !extractFromString( std::string{ cellHeight_attrib->value() }, setting.cell_height ) )
  {
    std::cerr << "Failed to load cell_height attribute for serialization node." << std::endl;
    return false;
  }
  
  const rapidxml::xml_attribute<>* const mode_attrib{ serialization_node->first_attribute( "mode" ) };
  if( !mode_attrib )
  {
    std::cerr << "Failed to locate mode attribute for serialization node." << std::endl;
    return false;
  }
  
  int num;
  if( !extractFromString( std::string{ mode_attrib->value() }, num ) )
  {
    std::cerr << "Failed to load mode attribute for serialization node." << std::endl;
    return false;
  }
  
  setting.compute_sigma_only = num == 0.0 ? false : true;
  
  return true;
}

bool openXMLFile( const std::string& filename, Setting& setting )
{
  setDefaultSetting( setting );
  
  rapidxml::xml_document<> doc;
  doc.clear();
  // Read the xml file into a vector
  std::ifstream theFile( filename );
  if( !theFile.is_open() )
  {
    return false;
  }
  
  std::vector<char> buffer( ( std::istreambuf_iterator<char>( theFile ) ), std::istreambuf_iterator<char>() );
  buffer.push_back('\0');
  // Parse the buffer using the xml file parsing library into doc.
  try
  {
    doc.parse<0>( &buffer[0] );
  }
  catch( const rapidxml::parse_error& e )
  {
    std::cerr << "Failed to open xml file: " << filename << std::endl;
    std::cerr << "Error message: " << e.what() << std::endl;
    return false;
  }
  
  if( doc.first_node( "mpm2d" ) == nullptr )
  {
    std::cerr << "Failed to locate root node 'mpm2d' in xml file: " << filename << std::endl;
    return false;
  }
  const rapidxml::xml_node<>& root_node{ *doc.first_node( "mpm2d" ) };
  
  if( !loadIntegrator( root_node, setting ) )
  {
    return false;
  }
  
  if( !loadCapture( root_node, setting ) )
  {
    return false;
  }
  
  if( !loadObstacle( root_node, setting ) )
  {
    return false;
  }
  
  if( !loadParticles( root_node, setting ) )
  {
    return false;
  }
  
  if( !loadConstitutiveModel( root_node, setting ) )
  {
    return false;
  }
  
  if( !loadExternalForce( root_node, setting ) )
  {
    return false;
  }
  
  if( !loadResume( root_node, setting ) )
  {
    return false;
  }
  
  if( !loadSerialization( root_node, setting ) )
  {
    return false;
  }
  
  showSettings( setting );
  
  return true;
}
