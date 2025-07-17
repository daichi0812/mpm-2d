//
//  xmlparser.h
//  SimpleMPM2D
//
//  Created by Yonghao Yue on 2020/09/04.
//

#ifndef xmlparser_h
#define xmlparser_h

#include <Eigen/Core>
#include <rapidxml.hpp>

enum ConStitutiveModelType
{
  CMT_SNOW,
  CMT_HERSCHEL_BULKLEY,
  CMT_THIXOTROPIC_HERSCHEL_BULKLEY,
};

struct Setting
{
  // integrator node
  Eigen::Vector4d gridRegion; //minGridX, minGridY, maxGridX, maxGridY
  Eigen::Vector2i gridRes; //gridResX, gridResY
  int particleCountPerCellEdge; //particle per cell *edge* ( = sqrt particle per cell in 2D )
  double dt;
  double valpha;
  double maxSimulationTime;
  double cell_height;
  
  // capture node (optional)
  bool capture;
  double captureInterval;
  
  // obstacle node
  std::string obstacleMapFileName;
  
  // particle node
  std::string initialParticlesFileName;
  double materialDensity;
  
  // constitutive_model node
  ConStitutiveModelType constitutiveModelType;
  // snow
  double xi;
  double theta_c;
  double theta_s;
  double YoungModulus;
  double PoissonRatio;
  // Herschel bulkley
  double kappa;
  double mu;
  double n;
  double sigma_y;
  double eta;
  // Thixotropic Herschel bulkley
  double k1;
  double k2;
  
  // external_force node
  double fric;
  double g;
  
  // input data file
  std::string templates_file_name_to_resume;
  std::string forces_file_name_to_resume;
  std::string homogenization_file_name;
  
  // output data file
  std::string serialization_folder;
  std::string forces_file_name_for_serialization;
  
  bool compute_sigma_only;
};

bool openXMLFile( const std::string& filename, Setting& setting );

#endif /* xmlparser_h */
