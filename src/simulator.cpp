#include "simulator.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

#include "imageiohdr.h"

#include "signeddistancefield.h"
#include "timer.h"
#include "mpm.h"
#include <Eigen/SVD>
#include <Eigen/Dense>

#include "basisfunction.h"
#include "constitutivemodel.h"

//inside(sdf<0.0): true
bool getObstacleInfo( const Simulator& in_Simulator, const Eigen::Vector2d& x, Eigen::Vector2d& n, Eigen::Vector2d& vel, double* adhesion )
{
  Eigen::Vector2d _i = ( x - in_Simulator.grid.minGrid ).array() / ( in_Simulator.grid.maxGrid - in_Simulator.grid.minGrid ).array();
	
  int i = int( floor( _i.x() * in_Simulator.obstacleSDF.getWidth() ) );
  int j = int( floor( _i.y() * in_Simulator.obstacleSDF.getHeight() ) );
  in_Simulator.obstacleSDF.clampPixelIndex( i, j );

	bool ret = in_Simulator.obstacleSDF( i, j, 0 ) < 0.0;

  int im = i-1, jm = j-1;
  in_Simulator.obstacleSDF.clampPixelIndex( im, jm );
  int ip = i+1, jp = j+1;
  in_Simulator.obstacleSDF.clampPixelIndex( ip, jp );

	double hx = ( ( im == i ) || ( ip == i ) ) ? 1.0 / in_Simulator.obstacleSDF.getWidth() : 2.0 / in_Simulator.obstacleSDF.getWidth();
	double hy = ( ( jm == j ) || ( jp == j ) ) ? 1.0 / in_Simulator.obstacleSDF.getHeight() : 2.0 / in_Simulator.obstacleSDF.getHeight();
	
	double dx = ( in_Simulator.obstacleSDF( ip, j, 0 ) - in_Simulator.obstacleSDF( im, j, 0 ) ) / hx;
	double dy = ( in_Simulator.obstacleSDF( i, jp, 0 ) - in_Simulator.obstacleSDF( i, jm, 0 ) ) / hy;
	
	vel[0] = in_Simulator.velocityMap( i, j, 0 );
	vel[1] = in_Simulator.velocityMap( i, j, 1 );
	
	*adhesion = in_Simulator.adhesionMap( i, j, 0 );
	
	double len_d = sqrt(dx*dx+dy*dy);
	n[0] = dx/len_d;
	n[1] = dy/len_d;
	
	return ret;
}

void setUpObjects( Simulator& io_Simulator )
{
	if( io_Simulator.timestep == 0 )
		generateSDF( io_Simulator.objectBinary, io_Simulator.obstacleSDF );
}

void updateObjects( Simulator& io_Simulator )
{
}

void transferMassFromPointsToGrid( const LinearBasisFunction* in_BasisFunction, const MaterialPoints& in_Points, Grid& io_Grid )
{
	for( int i=0; i<in_Points.nPoints; i++ )
	{
    auto stencil = in_BasisFunction->getStencil( in_Points.x[i], io_Grid );
    for( auto p = stencil.begin(); p!= stencil.end(); ++p )
		{
      const Eigen::Vector2i& gridIdx = *p;
      io_Grid.mass[ io_Grid.flatIdx( gridIdx ) ] += in_Points.mass0[i] * in_BasisFunction->weight( in_Points.x[i], gridIdx, io_Grid );
		}
	}
}

void transferVelocityFromPointsToGrid( const LinearBasisFunction* in_BasisFunction, const MaterialPoints& in_Points, Grid& io_Grid )
{
	for( int i=0; i<in_Points.nPoints; i++ )
	{
		auto stencil = in_BasisFunction->getStencil( in_Points.x[i], io_Grid );
    for( auto p = stencil.begin(); p!= stencil.end(); ++p )
    {
      const Eigen::Vector2i& gridIdx = *p;
      io_Grid.v[ io_Grid.flatIdx( gridIdx ) ] += in_Points.v[i] * in_Points.mass0[i] * in_BasisFunction->weight( in_Points.x[i], gridIdx, io_Grid ) / io_Grid.mass[ io_Grid.flatIdx( gridIdx ) ];
		}
	}
}

void computeParticleDensityUsingGrid( const LinearBasisFunction* in_BasisFunction, const Grid& in_Grid, MaterialPoints& io_Points )
{
	for( int i=0; i<io_Points.nPoints; i++ )
	{
		io_Points.density0[i] = 0.0;
		
    auto stencil = in_BasisFunction->getStencil( io_Points.x[i], in_Grid );
    for( auto p = stencil.begin(); p!= stencil.end(); ++p )
    {
      const Eigen::Vector2i& gridIdx = *p;
      io_Points.density0[i] += in_Grid.mass[ in_Grid.flatIdx( gridIdx ) ] * in_BasisFunction->weight( io_Points.x[i], gridIdx, in_Grid ) / ( in_Grid.h * in_Grid.h );
		}
	}
}

void compareParticleDensityUsingGrid( const LinearBasisFunction* in_BasisFunction, const Grid& in_Grid, MaterialPoints& io_Points )
{
  double total = 0.0;
  
  for( int i=0; i<io_Points.nPoints; i++ )
  {
    double density = 0.0;
    auto stencil = in_BasisFunction->getStencil( io_Points.x[i], in_Grid );
    for( auto p = stencil.begin(); p!= stencil.end(); ++p )
    {
      const Eigen::Vector2i& gridIdx = *p;
      density += in_Grid.mass[ in_Grid.flatIdx( gridIdx ) ] * in_BasisFunction->weight( io_Points.x[i], gridIdx, in_Grid ) / ( in_Grid.h * in_Grid.h );
    }
    total += io_Points.density0[i] - density ;
  }
  
  std::cout << "average diff: " << total / (double)io_Points.nPoints << std::endl;
}

void computeParticleVolumeUsingGrid( const Grid& in_Grid, MaterialPoints& io_Points )
{
	for( int i=0; i<io_Points.nPoints; i++ )
	{
		io_Points.volume0[i] = io_Points.mass0[i] / io_Points.density0[i];
	}
}

void computeParticleKirchhoffStress( const ConstitutiveModel* in_ConstitutiveModel, MaterialPoints& io_Points )
{
	for( int i=0; i<io_Points.nPoints; i++ )
	{
    io_Points.tau[i] = in_ConstitutiveModel->computeKirchhoffStress( io_Points, i );
	}
}

void computeParticleKirchhoffStressWithData( const ConstitutiveModel* in_ConstitutiveModel, MaterialPoints& io_Points, ForceData2D& io_ForceData )
{
  for( int i=0; i<io_Points.nPoints; i++ )
  {
    const double J = sqrt( io_Points.be[i].determinant() );
    io_Points.tau[i] = in_ConstitutiveModel->computeKirchhoffStress( io_Points, i );
    io_ForceData.storeForce( io_Points.x[i], io_Points.tau[i] / J, io_Points.density0[i], io_Points.volume0[i], J );
  }
}

void computeForceOnGrid( const LinearBasisFunction* in_BasisFunction, const MaterialPoints& in_Points, Grid& io_Grid )
{
	for( int i=0; i<in_Points.nPoints; i++ )
	{
    auto stencil = in_BasisFunction->getStencil( in_Points.x[i], io_Grid );
    for( auto p = stencil.begin(); p!= stencil.end(); ++p )
    {
      const Eigen::Vector2i& gridIdx = *p;
      Eigen::Vector2d wg = in_BasisFunction->weightGrad( in_Points.x[i], gridIdx, io_Grid );
      io_Grid.f[ io_Grid.flatIdx( gridIdx ) ] -= in_Points.volume0[i] * in_Points.tau[i] * wg;
		}
	}
}

void computeVstarOnGrid( Grid& io_Grid, double in_dt, const Eigen::Vector2d& g )
{
  int n = io_Grid.numGridPoints();
	for( int i=0; i<n; i++ )
	{
		if( io_Grid.mass[i] > 0.0 )
		{
			io_Grid.vstar[i] = io_Grid.v[i] + in_dt * io_Grid.f[i] / io_Grid.mass[i] + in_dt * g;
		}
	}
}

void applyGridBasedBodyCollision( Simulator& io_Simulator, double fric )
{
  int N = io_Simulator.grid.numGridPoints();
	for( int i=0; i<N; i++ )
	{
		if( io_Simulator.grid.mass[i] <= 0.0 ) continue;
		
    Eigen::Vector2i _i = io_Simulator.grid.gridIdx(i);
		
    Eigen::Vector2d x = io_Simulator.grid.minGrid + io_Simulator.grid.h * _i.cast<double>();
		
		Eigen::Vector2d n; double adhesion;
		Eigen::Vector2d velObj;
		bool insideObstacle = getObstacleInfo( io_Simulator, x, n, velObj, &adhesion );
		
		if(!insideObstacle) continue;
		
		if(adhesion > 0.5)
		{
			io_Simulator.grid.vstar[i] = velObj;
			continue;
		}

    Eigen::Vector2d velRel = io_Simulator.grid.vstar[i] - velObj;

    double vn = n.dot( velRel );
		if( vn > 0.0 ) continue;
		
		Eigen::Vector2d vt = velRel - vn * n;
    double len_vt = vt.norm();
		
		if( len_vt <= -fric * vn )
		{
			io_Simulator.grid.vstar[i] = velObj;
		}
		else
		{
			io_Simulator.grid.vstar[i] = velObj + vt + fric * vn * vt / len_vt;
		}
	}
}

void updateVelocityOnGrid( Simulator& io_Simulator )
{
  int n = io_Simulator.grid.numGridPoints();
  for( int i=0; i<n; i++ )
  {
    io_Simulator.grid.vnew[i] = io_Simulator.grid.vstar[i];
  }
}

void updateDeformationStatus( const ConstitutiveModel* in_ConstitutiveModel, const LinearBasisFunction* in_BasisFunction, const Grid& in_Grid, MaterialPoints& io_Points, double in_dt )
{
	for( int i=0; i<io_Points.nPoints; i++ )
	{
		Eigen::Matrix2d velocity_gradient = Eigen::Matrix2d::Zero();
		
		auto stencil = in_BasisFunction->getStencil( io_Points.x[i], in_Grid );
    for( auto p = stencil.begin(); p!= stencil.end(); ++p )
    {
      const Eigen::Vector2i& gridIdx = *p;
      Eigen::Vector2d wg = in_BasisFunction->weightGrad( io_Points.x[i], gridIdx, in_Grid );
      velocity_gradient += in_Grid.vnew[ in_Grid.flatIdx( gridIdx ) ] * wg.transpose();
		}
		
    in_ConstitutiveModel->updateDeformationStatus( io_Points, i, in_dt, velocity_gradient );
	}
}

void updateDeformationStatusWithData( const ConstitutiveModel* in_ConstitutiveModel, const LinearBasisFunction* in_BasisFunction, const Grid& in_Grid, MaterialPoints& io_Points, double in_dt )
{
  for( int i=0; i<io_Points.nPoints; i++ )
  {
    Eigen::Matrix2d velocity_gradient = Eigen::Matrix2d::Zero();
    
    auto stencil = in_BasisFunction->getStencil( io_Points.x[i], in_Grid );
    for( auto p = stencil.begin(); p!= stencil.end(); ++p )
    {
      const Eigen::Vector2i& gridIdx = *p;
      Eigen::Vector2d wg = in_BasisFunction->weightGrad( io_Points.x[i], gridIdx, in_Grid );
      velocity_gradient += in_Grid.v[ in_Grid.flatIdx( gridIdx ) ] * wg.transpose();
    }
    
    io_Points.be[i] = io_Points.be[i] + in_dt * ( velocity_gradient * io_Points.be[i] + io_Points.be[i] * velocity_gradient.transpose() );
  }
}

void updateParticleVelocity( const LinearBasisFunction* in_BasisFunction, const Grid& in_Grid, MaterialPoints& io_Points, double alpha )
{
	for( int i=0; i<io_Points.nPoints; i++ )
	{
		Eigen::Vector2d vpic = Eigen::Vector2d::Zero();
		Eigen::Vector2d vflip = io_Points.v[i];
		
		auto stencil = in_BasisFunction->getStencil( io_Points.x[i], in_Grid );
    for( auto p = stencil.begin(); p!= stencil.end(); ++p )
    {
      const Eigen::Vector2i& gridIdx = *p;
      double w = in_BasisFunction->weight( io_Points.x[i], gridIdx, in_Grid );
				
      vpic += in_Grid.vnew[ in_Grid.flatIdx( gridIdx ) ] * w;
      vflip += w * ( in_Grid.vnew[ in_Grid.flatIdx( gridIdx ) ] - in_Grid.v[ in_Grid.flatIdx( gridIdx ) ] );
		}
		
		io_Points.v[i] = ( 1.0 - alpha ) * vpic + alpha * vflip;
	}
}

void applyParticleBasedBodyCollision( Simulator& io_Simulator, double fric )
{
	for( int i=0; i<io_Simulator.points.nPoints; i++ )
	{
		Eigen::Vector2d n; double adhesion;
		Eigen::Vector2d velObj;
		bool insideObstacle = getObstacleInfo( io_Simulator, io_Simulator.points.x[i], n, velObj, &adhesion);
		
		if( !insideObstacle ) continue;
		
		if( adhesion > 0.5 )
		{
			io_Simulator.points.v[i] = velObj;
			continue;
		}
		
		Eigen::Vector2d velRel = io_Simulator.points.v[i] - velObj;
		
    double vn = n.dot( velRel );
		if( vn > 0.0 ) continue;
		
		Eigen::Vector2d vt = velRel - vn * n;
    double len_vt = vt.norm();
		
		if( len_vt <= -fric * vn )
		{
			io_Simulator.points.v[i] = velObj;
		}
		else
		{
			io_Simulator.points.v[i] = velObj + vt + fric * vn * vt / len_vt;
		}
	}
}

void updateParticlePosition( MaterialPoints& io_Points, double in_dt )
{
	for( int i=0; i<io_Points.nPoints; i++ )
	{
		io_Points.x[i] += in_dt * io_Points.v[i];
	}
}

void initParticlesWithData( Simulator& io_Simulator )
{
  int nParticles = 0;
  
  for(int i = 0; i < io_Simulator.data.numElements(); i++)
  {
    const Element element = io_Simulator.data.getElement(i);
    if( !element.is_static && io_Simulator.grid.isInsideGrid( element.center_of_mass ) )
      nParticles++;
  }
  
  io_Simulator.points.nPoints = nParticles;
  io_Simulator.points.nPointsAlloc = nParticles;
  
  io_Simulator.points.mass0.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.density0.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.volume0.resize( io_Simulator.points.nPointsAlloc );
  
  io_Simulator.points.x.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.v.resize( io_Simulator.points.nPointsAlloc );
  
  io_Simulator.points.F.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.Fe.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.Fp.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.tau.resize( io_Simulator.points.nPointsAlloc );
  
  io_Simulator.points.be.resize( io_Simulator.points.nPointsAlloc );
  
  io_Simulator.points.xi.resize( io_Simulator.points.nPointsAlloc );
  
  //set up particles
  
  int idx = 0;
  
  for(int i = 0; i < io_Simulator.data.numElements(); i++)
  {
    const Element element = io_Simulator.data.getElement(i);
    
    if( !element.is_static && io_Simulator.grid.isInsideGrid( element.center_of_mass ) )
    {
      io_Simulator.points.x[idx] = element.center_of_mass;
      
      io_Simulator.points.mass0[idx] = element.mass;
      
      io_Simulator.points.density0[idx] = element.density;
      
      io_Simulator.points.v[idx] = element.velocity;

      io_Simulator.points.F[idx].setIdentity();

      io_Simulator.points.Fe[idx].setIdentity();

      io_Simulator.points.Fp[idx].setIdentity();

      io_Simulator.points.tau[idx].setZero();
      
      io_Simulator.points.xi[idx] = 1.0;

      io_Simulator.points.be[idx].setZero();
      
//      auto stencil = io_Simulator.basisFunction->getCenterStencil( element.center_of_mass, io_Simulator.grid );
//      for( auto p = stencil.begin(); p != stencil.end(); ++p )
//      {
//        const Eigen::Vector2i& gridIdx = *p;
//        const Eigen::Vector2d gridCenterPos = ( io_Simulator.grid.minGrid + gridIdx.cast<double>() * io_Simulator.grid.h ).array() + 0.5 * io_Simulator.grid.h;
//        const double w = io_Simulator.basisFunction->weight( element.center_of_mass, gridCenterPos, io_Simulator.grid );
//        const int flatIdx = io_Simulator.grid.flatCenterIdx( gridIdx );
//        io_Simulator.points.be[idx] += io_Simulator.data.getHomogenizeData( flatIdx ).strain * w;
//      }

      const Eigen::Vector2i gridIdx = io_Simulator.grid.gridCenterIdx( element.center_of_mass );
      const int flatIdx = io_Simulator.grid.flatCenterIdx( gridIdx );
      io_Simulator.points.be[idx] = io_Simulator.data.getHomogenizeData( flatIdx ).strain;
      
      idx++;
    }
  }
};

void initParticles( Simulator& io_Simulator, const Setting& in_Setting )
{
	Image<double, 1> the_Image;
  loadImageHdr( in_Setting.initialParticlesFileName.c_str(), the_Image );
	
	int minImgRegionX = the_Image.getWidth();
	int minImgRegionY = the_Image.getHeight();
	int maxImgRegionX = 0;
	int maxImgRegionY = 0;
	
	for( int j=0; j<the_Image.getHeight(); j++ )
	{
		for( int i=0; i<the_Image.getWidth(); i++ )
		{
			double val = the_Image( i, j, 0 );
			if(val < 0.5)
			{
				minImgRegionX = min( minImgRegionX, i ); minImgRegionY = min( minImgRegionY, j );
				maxImgRegionX = max( maxImgRegionX, i ); maxImgRegionY = max( maxImgRegionY, j );
			}
		}
	}
	
  int minGridX = max( 0, int( floor( ( double( minImgRegionX ) / double( the_Image.getWidth() ) ) * io_Simulator.grid.gridRes.x() ) ) );
  int minGridY = max( 0, int( floor( ( double( minImgRegionY ) / double( the_Image.getHeight() ) ) * io_Simulator.grid.gridRes.y() ) ) );
  int maxGridX = min( io_Simulator.grid.gridRes.x(), int( ceil( ( double( maxImgRegionX+1 ) / double( the_Image.getWidth() ) ) * io_Simulator.grid.gridRes.x() ) ) );
  int maxGridY = min( io_Simulator.grid.gridRes.y(), int( ceil( ( double( maxImgRegionY+1 ) / double( the_Image.getHeight() ) ) * io_Simulator.grid.gridRes.y() ) ) );
	
	//count particles
	
	int nParticles = 0;
	
	for( int j=minGridY; j<maxGridY; j++ )
	{
		for( int i=minGridX; i<maxGridX; i++ )
		{
      for( int v=0; v<in_Setting.particleCountPerCellEdge; v++ )
			{
				for( int u=0; u<in_Setting.particleCountPerCellEdge; u++ )
				{
          double s = ( i + (u+0.5)/in_Setting.particleCountPerCellEdge ) / io_Simulator.grid.gridRes.x();
          double t = ( j + (v+0.5)/in_Setting.particleCountPerCellEdge ) / io_Simulator.grid.gridRes.y();
					
          double val = the_Image.atUVNearest( s, t, 0 );
					
					if( val < 0.5 ) nParticles++;
				}
			}
		}
	}
	
	//alloc data
	
	io_Simulator.points.nPoints = nParticles;
	io_Simulator.points.nPointsAlloc = nParticles;
	
  io_Simulator.points.mass0.resize( io_Simulator.points.nPointsAlloc );
	io_Simulator.points.density0.resize( io_Simulator.points.nPointsAlloc );
	io_Simulator.points.volume0.resize( io_Simulator.points.nPointsAlloc );
	
	
  io_Simulator.points.x.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.v.resize( io_Simulator.points.nPointsAlloc );
	
	io_Simulator.points.F.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.Fe.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.Fp.resize( io_Simulator.points.nPointsAlloc );
  io_Simulator.points.tau.resize( io_Simulator.points.nPointsAlloc );
  
  io_Simulator.points.be.resize( io_Simulator.points.nPointsAlloc );
  
  io_Simulator.points.xi.resize( io_Simulator.points.nPointsAlloc );
	
	double dens = in_Setting.materialDensity; //1000.0 * 0.02;
  double particleMass = dens * io_Simulator.grid.singleCellArea() / ( in_Setting.particleCountPerCellEdge * in_Setting.particleCountPerCellEdge );
	
	//set up particles
	
	int idx = 0;
	
	for( int j=minGridY; j<maxGridY; j++ )
	{
		for( int i=minGridX; i<maxGridX; i++ )
		{
			for( int v=0; v<in_Setting.particleCountPerCellEdge; v++ )
			{
				for( int u=0; u<in_Setting.particleCountPerCellEdge; u++ )
				{
          double s = ( i + (u+0.5)/in_Setting.particleCountPerCellEdge ) / io_Simulator.grid.gridRes.x();
          double t = ( j + (v+0.5)/in_Setting.particleCountPerCellEdge ) / io_Simulator.grid.gridRes.y();
					
					double val = the_Image.atUVNearest( s, t, 0 );
					
					if( val < 0.5 )
					{
            io_Simulator.points.x[idx] = Eigen::Vector2d{ s, 1.0 - t }.array() * ( io_Simulator.grid.maxGrid - io_Simulator.grid.minGrid ).array() + io_Simulator.grid.minGrid.array();
            
						io_Simulator.points.mass0[idx] = particleMass;
            io_Simulator.points.v[idx].setZero();
		
            io_Simulator.points.F[idx].setIdentity();
		
						io_Simulator.points.Fe[idx].setIdentity();
		
						io_Simulator.points.Fp[idx].setIdentity();
		
            io_Simulator.points.tau[idx].setZero();
            
            io_Simulator.points.be[idx].setIdentity();
            
            io_Simulator.points.xi[idx] = 1.0;
					
						idx++;
					}
				}
			}
		}
	}
}

void separateObstacleAdhesionMaps( const Image<double, 2>& in_Image, Image<double, 1>& out_ObstacleMap, Image<double, 1>& out_AdhesionMap)
{
  out_ObstacleMap.resizeToFitTargetImage( in_Image );
  out_AdhesionMap.resizeToFitTargetImage( in_Image );
	
  for( int j=0; j<in_Image.getHeight(); j++ )
  {
    for( int i=0; i<in_Image.getWidth(); i++ )
    {
      double r = in_Image( i, j, 0 );
      double g = in_Image( i, j, 1 );
      
      if( ( r < 0.5 ) || ( g < 0.5 ) ) out_ObstacleMap( i, j, 0 ) = 0.0; else out_ObstacleMap( i, j, 0 ) = 1.0;
      if( ( r > 0.5 ) && ( g < 0.5 ) ) out_AdhesionMap( i, j, 0 ) = 1.0; else out_AdhesionMap( i, j, 0 ) = 0.0;
    }
  }
}

void setGridFromData( Simulator& io_Simulator, const Setting& in_Setting )
{
  io_Simulator.grid.setGrid( io_Simulator.data.getGridStart(), in_Setting.cell_height, io_Simulator.data.getDemResolution() );
  io_Simulator.force_data.setGridInfo( io_Simulator.grid );
}

void initSimulator( Simulator& io_Simulator, const Setting& in_Setting )
{
	io_Simulator.timestep = 0;
  io_Simulator.deltaT = in_Setting.dt;

  io_Simulator.fric = in_Setting.fric;
  io_Simulator.g = Eigen::Vector2d{ 0.0, -in_Setting.g };
    
	io_Simulator.valpha = in_Setting.valpha;
  		
	const double lambda = in_Setting.YoungModulus * in_Setting.PoissonRatio / ( (1.0 + in_Setting.PoissonRatio ) * ( 1.0 - 2.0 * in_Setting.PoissonRatio ) );
	double bulkModulus = in_Setting.YoungModulus / ( 3.0 * ( 1.0 - 2.0 * in_Setting.PoissonRatio ) );
	const double mu = bulkModulus - lambda; //2D
	//io_Simulator.mu = 1.5 * ( bulkModulus - lambda ); //3D
  
	printf( "lambda: %e, mu: %e\n", lambda, mu );
  
  if( in_Setting.constitutiveModelType == CMT_SNOW )
  {
    io_Simulator.constitutiveModel = new SnowModel( lambda, mu, in_Setting.theta_c, in_Setting.theta_s, in_Setting.xi );
  }
  else if( in_Setting.constitutiveModelType == CMT_HERSCHEL_BULKLEY )
  {
    io_Simulator.constitutiveModel = new HerschelBulkley( in_Setting.kappa, in_Setting.mu, in_Setting.n, in_Setting.sigma_y, in_Setting.eta );
  }
  else if( in_Setting.constitutiveModelType == CMT_THIXOTROPIC_HERSCHEL_BULKLEY )
  {
    io_Simulator.constitutiveModel = new ThixotropicHerschelBulkley( in_Setting.kappa, in_Setting.mu, in_Setting.n, in_Setting.sigma_y, in_Setting.eta, in_Setting.k1, in_Setting.k2 );
  }
  
  io_Simulator.grid.setGrid( in_Setting.gridRegion.segment<2>(0), in_Setting.gridRegion.segment<2>(2), in_Setting.gridRes );

  io_Simulator.grid.clear();

  initParticles( io_Simulator, in_Setting );

	Image<double, 2> the_Image;

  loadImageHdr( in_Setting.obstacleMapFileName.c_str(), the_Image );
  the_Image.invertY();

	separateObstacleAdhesionMaps( the_Image, io_Simulator.fixedObjectBinary, io_Simulator.fixedAdhesionMap );

  io_Simulator.objectBinary.resizeToFitTargetImage( the_Image );
  io_Simulator.adhesionMap.resizeToFitTargetImage( the_Image );
  io_Simulator.obstacleSDF.resizeToFitTargetImage( the_Image );
  io_Simulator.velocityMap.resizeToFitTargetImage( the_Image );

  io_Simulator.velocityMap.zeroClear();
  io_Simulator.objectBinary.copyDataFrom( io_Simulator.fixedObjectBinary );
  io_Simulator.adhesionMap.copyDataFrom( io_Simulator.fixedAdhesionMap );

	setUpObjects(io_Simulator);
  
  io_Simulator.basisFunction = new LinearBasisFunction();
}

Eigen::Matrix2d dev2( const Eigen::Matrix2d& in_x )
{
  return in_x - in_x.trace() * 0.5 * Eigen::Matrix2d::Identity(); //for 2D
}

void computeStrain( Simulator& io_Simulator, const Setting& in_Setting )
{
  for( int i = 0; i < io_Simulator.data.numHomogenizeData(); i++ )
  {
    const double b = io_Simulator.data.getHomogenizeData(i).sigma.trace() / in_Setting.kappa;
    const double J = ( b + sqrt( b * b + 4.0 ) ) / 2.0;
    
    const Eigen::Matrix2d dev_be_bar = J * dev2( io_Simulator.data.getHomogenizeData(i).sigma ) / in_Setting.mu;
    const double t = - ( ( dev_be_bar(0, 0) + dev_be_bar(1, 1) ) + sqrt( ( dev_be_bar(0, 0) - dev_be_bar(1, 1) ) * ( dev_be_bar(0, 0) - dev_be_bar(1, 1) ) + 4 * ( dev_be_bar(0, 1) * dev_be_bar(1, 0) + 1 ) ) ) / 2.0;
    
    const Eigen::Matrix2d be_bar = dev_be_bar + t * Eigen::Matrix2d::Identity();
    
    const Eigen::Matrix2d be = J * be_bar;
    
    io_Simulator.data.getHomogenizeData(i).strain = be;
  }
}

void testStrainValue( Simulator& io_Simulator, const Setting& in_Setting )
{
  for( int i = 0; i < io_Simulator.data.numHomogenizeData(); i++ )
  {
    const double J = sqrt( io_Simulator.data.getHomogenizeData(i).strain.determinant() );
    const Eigen::Matrix2d be_bar = io_Simulator.data.getHomogenizeData(i).strain / J; // for 2D
    const Eigen::Matrix2d tau = 0.5 * in_Setting.kappa * ( J*J - 1.0 ) * Eigen::Matrix2d::Identity() + in_Setting.mu * dev2( be_bar );
    const Eigen::Matrix2d sigma = tau / J;

    if( ( ( io_Simulator.data.getHomogenizeData(i).sigma - sigma ).array().abs() > 1.0e-10 ).any() )
    {
      std::cout << "---------- " << J << std::endl << std::endl;
      std::cout << io_Simulator.data.getHomogenizeData(i).sigma << std::endl << std::endl << sigma << std::endl << std::endl;
      std::cout << io_Simulator.data.getHomogenizeData(i).sigma - sigma << std::endl << std::endl;
    }
  }
}

void testStrainValue2( Simulator& io_Simulator )
{
//  vector<Eigen::Matrix2d> sigmaArray; sigmaArray.resize( io_Simulator.data.numHomogenizeData() );
//
//  for( int i = 0; i < sigmaArray.size(); i++ )
//    sigmaArray[i].setZero();
//
//  for( int i = 0; i < io_Simulator.points.nPoints; i++ )
//  {
//    const double J = sqrt( io_Simulator.points.be[i].determinant() );
//
//    auto stencil = io_Simulator.basisFunction->getCenterStencil( io_Simulator.points.x[i], io_Simulator.grid );
//    for( auto p = stencil.begin(); p != stencil.end(); ++p )
//    {
//      const Eigen::Vector2i& gridIdx = *p;
//      const Eigen::Vector2d gridCenterPos = ( io_Simulator.grid.minGrid + gridIdx.cast<double>() * io_Simulator.grid.h ).array() + 0.5 * io_Simulator.grid.h;
//      const double w = io_Simulator.basisFunction->weight( io_Simulator.points.x[i], gridCenterPos, io_Simulator.grid );
//      const int flatIdx = io_Simulator.grid.flatCenterIdx( gridIdx );
//      if( !io_Simulator.data.getHomogenizeData(flatIdx).strain.isIdentity() )
//        sigmaArray[flatIdx] += ( io_Simulator.points.tau[i] / J ) * w;
//    }
//  }
//
//  for( int i = 0; i < io_Simulator.data.numHomogenizeData(); i++ )
//  {
//    if( ( ( io_Simulator.data.getHomogenizeData(i).sigma - sigmaArray[i] ).array().abs() > 1.0e-10 ).any() )
//      std::cout << i << std::endl << io_Simulator.data.getHomogenizeData(i).sigma - sigmaArray[i] << std::endl << std::endl;
//  }
  
  for( int i = 0; i < io_Simulator.points.nPoints; i++ )
  {
    const double J = sqrt( io_Simulator.points.be[i].determinant() );
    const auto gridIdx = io_Simulator.grid.gridCenterIdx( io_Simulator.points.x[i] );
    const int flatIdx = io_Simulator.grid.flatCenterIdx( gridIdx );
    if( ( ( io_Simulator.data.getHomogenizeData(flatIdx).sigma - io_Simulator.points.tau[i] / J ).array().abs() > 1.0e-10 ).any() )
      std::cout << i << std::endl << io_Simulator.data.getHomogenizeData(flatIdx).sigma - io_Simulator.points.tau[i] / J << std::endl << std::endl;
  }
}

void initSimulatorWithData( Simulator& io_Simulator, const Setting& in_Setting )
{
  io_Simulator.timestep = 0;
  io_Simulator.deltaT = io_Simulator.data.getProgressData().dt;

  io_Simulator.fric = in_Setting.fric;
  io_Simulator.g = Eigen::Vector2d{ 0.0, -in_Setting.g };
    
  io_Simulator.valpha = in_Setting.valpha;
  
  computeStrain( io_Simulator, in_Setting );
  
  if( io_Simulator.constitutiveModel == nullptr )
    io_Simulator.constitutiveModel = new HerschelBulkley( in_Setting.kappa, in_Setting.mu, in_Setting.n, in_Setting.sigma_y, in_Setting.eta );
    
  setGridFromData( io_Simulator, in_Setting );
  io_Simulator.grid.clear();
  
  initParticlesWithData( io_Simulator );

  if( io_Simulator.basisFunction == nullptr )
    io_Simulator.basisFunction = new LinearBasisFunction();
}

void advanceSingleStep( Simulator& io_Simulator )
{
	Timer the_Timer;
	
	startTimer( the_Timer );
	setUpObjects( io_Simulator );
	stopTimer( the_Timer );
	double eTSetUpObjects = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
  io_Simulator.grid.clear();
	stopTimer( the_Timer );
	double eTClear = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
  transferMassFromPointsToGrid( io_Simulator.basisFunction, io_Simulator.points, io_Simulator.grid );
	transferVelocityFromPointsToGrid( io_Simulator.basisFunction, io_Simulator.points, io_Simulator.grid );
	stopTimer( the_Timer );
	double eTTransfer = getElapsedMs( the_Timer );
	
	if( io_Simulator.timestep == 0 )
	{
		computeParticleDensityUsingGrid( io_Simulator.basisFunction, io_Simulator.grid, io_Simulator.points );
		computeParticleVolumeUsingGrid( io_Simulator.grid, io_Simulator.points );
	}
	
	startTimer( the_Timer );
  computeParticleKirchhoffStress( io_Simulator.constitutiveModel, io_Simulator.points );
	stopTimer( the_Timer );
	double eTStress = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
	computeForceOnGrid( io_Simulator.basisFunction, io_Simulator.points, io_Simulator.grid );
	stopTimer( the_Timer );
	double eTGForce = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
	computeVstarOnGrid( io_Simulator.grid, io_Simulator.deltaT, io_Simulator.g );
	stopTimer( the_Timer );
	double eTGVstar = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
	applyGridBasedBodyCollision( io_Simulator, io_Simulator.fric );
	stopTimer( the_Timer );
	double eTGCollision = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
	updateVelocityOnGrid( io_Simulator );
	stopTimer( the_Timer );
	double eTGVel = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
  updateDeformationStatus( io_Simulator.constitutiveModel, io_Simulator.basisFunction, io_Simulator.grid, io_Simulator.points, io_Simulator.deltaT );
	stopTimer( the_Timer );
	double eTDefGrad = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
	updateParticleVelocity( io_Simulator.basisFunction, io_Simulator.grid, io_Simulator.points, io_Simulator.valpha );
	stopTimer( the_Timer );
	double eTPVel = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
	applyParticleBasedBodyCollision( io_Simulator, io_Simulator.fric );
	stopTimer( the_Timer );
	double eTPCollision = getElapsedMs( the_Timer );
	
	startTimer( the_Timer );
	updateParticlePosition( io_Simulator.points, io_Simulator.deltaT );
	stopTimer( the_Timer );
	double eTPPos = getElapsedMs( the_Timer );
	
	updateObjects( io_Simulator );
		
	printf( "\rtimestep: %d (SU:%3.1f, C:%3.1f, T:%3.1f, S:%3.1f, GF:%3.1f, GVs:%3.1f, GC:%3.1f, GV:%3.1f, DG:%3.1f, PV:%3.1f, PC:%3.1f, PP:%3.1f)   ", io_Simulator.timestep,
		eTSetUpObjects, eTClear, eTTransfer, eTStress, eTGForce, eTGVstar, eTGCollision, eTGVel, eTDefGrad, eTPVel, eTPCollision, eTPPos );
	
	io_Simulator.timestep++;
}

void advanceStepWithData( Simulator& io_Simulator )
{
  Timer the_Timer;
  
  startTimer( the_Timer );
  io_Simulator.grid.clear();
  stopTimer( the_Timer );
  double eTClear = getElapsedMs( the_Timer );
  
  startTimer( the_Timer );
  transferMassFromPointsToGrid( io_Simulator.basisFunction, io_Simulator.points, io_Simulator.grid );
  transferVelocityFromPointsToGrid( io_Simulator.basisFunction, io_Simulator.points, io_Simulator.grid );
  stopTimer( the_Timer );
  double eTTransfer = getElapsedMs( the_Timer );

  computeParticleDensityUsingGrid( io_Simulator.basisFunction, io_Simulator.grid, io_Simulator.points );
  computeParticleVolumeUsingGrid( io_Simulator.grid, io_Simulator.points );
  
  startTimer( the_Timer );
  updateDeformationStatusWithData( io_Simulator.constitutiveModel, io_Simulator.basisFunction, io_Simulator.grid, io_Simulator.points, io_Simulator.deltaT );
  stopTimer( the_Timer );
  double eTDefGrad = getElapsedMs( the_Timer );
  
  startTimer( the_Timer );
  computeParticleKirchhoffStressWithData( io_Simulator.constitutiveModel, io_Simulator.points, io_Simulator.force_data );
  stopTimer( the_Timer );
  double eTStress = getElapsedMs( the_Timer );
  
  /*  Cauchy stress への変換は書き出し用変数への保存時に行っている ( sigma は ForceData2D 側で保持している ） */
  
//  testStrainValue2( io_Simulator );
//  printf( "\rtimestep: %d (C:%3.1f, T:%3.1f, S:%3.1f, DG:%3.1f)\n", io_Simulator.timestep, eTClear, eTTransfer, eTStress, eTDefGrad);
}
