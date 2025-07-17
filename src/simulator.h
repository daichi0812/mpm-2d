#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

#include "image.h"

#include <map>
#include <vector>
using namespace std;

#include "materialpoints.h"
#include "grid.h"
#include "xmlparser.h"
#include "RigidBody2DData.h"
#include "ForceData2D.h"
#include <Eigen/Core>

class LinearBasisFunction;
class ConstitutiveModel;

struct Simulator
{
	MaterialPoints points;
	Grid grid;
	Image<double, 1> obstacleSDF;
	
	Image<double, 1> fixedObjectBinary;
	Image<double, 1> objectBinary;
	Image<double, 1> fixedAdhesionMap;
	Image<double, 1> adhesionMap;
	Image<double, 2> velocityMap;
	
	int timestep;
	double deltaT;

  double valpha;
  
  double fric;
  Eigen::Vector2d g;
	  
  LinearBasisFunction* basisFunction;
  ConstitutiveModel* constitutiveModel;
  
  RigidBody2DData data;
  ForceData2D force_data;
};

void initSimulator( Simulator& io_Simulator, const Setting& in_Setting );
void initSimulatorWithData( Simulator& io_Simulator, const Setting& in_Setting );
void testStrainValue( Simulator& io_Simulator, const Setting& in_Setting );
void advanceStepWithData( Simulator& io_Simulator );
void advanceSingleStep( Simulator& io_Simulator );

#endif
