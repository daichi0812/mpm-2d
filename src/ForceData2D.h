//
//  ForceData2D.h
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/11.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#ifndef __FORCE_DATA_2D_H__
#define __FORCE_DATA_2D_H__

#include <Eigen/Core>
#include <vector>
#include "HDF5io.h"

struct SingleForceData2D
{
  SingleForceData2D( const Eigen::Vector2d& in_center_of_mass, const Eigen::Matrix2d& in_sigma, const double& in_density, const double& in_volume, const double& in_J )
  : center_of_mass( in_center_of_mass ), sigma( in_sigma ), density( in_density ), volume( in_volume ), J( in_J ) {}
  
  Eigen::Vector2d center_of_mass;
  Eigen::Matrix2d sigma;
  double density;
  double volume;
  double J;
};

class ForceData2D
{
public:
  ForceData2D() {}
  ~ForceData2D() {}
  
  void clear() { m_Forces.clear(); }
  
  void storeForce( const Eigen::Vector2d& in_center_of_mass, const Eigen::Matrix2d& in_sigma, const double& in_density, const double& in_volume, const double& in_J )
  {
    m_Forces.emplace_back( in_center_of_mass, in_sigma, in_density, in_volume, in_J );
  }
  
  void setGridInfo( Grid& grid )
  {
    gridRes = grid.gridRes;
    minGrid = grid.minGrid;
    maxGrid = grid.maxGrid;
    h = grid.h;
  };
  
  void serializeForces( HDF5File& io_HDF5, const std::string& in_time_step = "" ) const
  {
    Eigen::Matrix2Xd center_of_mass_array; center_of_mass_array.resize( 2, m_Forces.size() );
    Eigen::Matrix2Xd sigma_array; sigma_array.resize( 2, m_Forces.size() * 2 );
    Eigen::VectorXd density_array; density_array.resize( m_Forces.size() );
    Eigen::VectorXd volume_array; volume_array.resize( m_Forces.size() );
    Eigen::VectorXd J_array; J_array.resize( m_Forces.size() );
    
    const int numGridPoints = ( gridRes.x() + 1 ) * ( gridRes.y() + 1 );
    Eigen::Matrix2Xd grid_pos_array; grid_pos_array.resize( 2, numGridPoints );
    
    int idx = 0;
    for( int i = 0; i <= gridRes.y(); i++ )
    {
      for( int j = 0; j <= gridRes.x(); j++ )
      {
        grid_pos_array.col(idx) = minGrid + h * Eigen::Vector2d(j, i);
        idx++;
      }
    }
    
    const int num_forces = m_Forces.size();
    for( int i=0; i<num_forces; i++ )
    {
      center_of_mass_array.col(i) = m_Forces[i].center_of_mass;
      density_array(i) = m_Forces[i].density;
      volume_array(i) = m_Forces[i].volume;
      J_array(i) = m_Forces[i].J;
      /*
       * python対応用の出力形式
       * 2*2行列の場合： 第１列と第２列を行列の総数分だけ開けて追加
       */
      sigma_array.col(i) = m_Forces[i].sigma.col(0);
      sigma_array.col(num_forces + i) = m_Forces[i].sigma.col(1);
    }

    const std::string group_name_homogenization = in_time_step + "/homogenization";
    io_HDF5.writeMatrix( group_name_homogenization, "center_of_mass", center_of_mass_array );
    io_HDF5.writeSquareMatrixArray( group_name_homogenization, "sigma", sigma_array );
    io_HDF5.writeVector( group_name_homogenization, "density", density_array );
    io_HDF5.writeVector( group_name_homogenization, "volume", volume_array );
    io_HDF5.writeVector( group_name_homogenization, "J", J_array );
    io_HDF5.writeVector( group_name_homogenization, "resolution", gridRes );
    io_HDF5.writeMatrix( group_name_homogenization, "grid_p", grid_pos_array );
    io_HDF5.writeScalar( group_name_homogenization, "h", h );
  }
  
  void printData() {
    for( int i = 0; i < m_Forces.size(); i++ ) {
      std::cout << i << std::endl;
      std::cout << "center_of_mass: " << m_Forces[i].center_of_mass.transpose() << std::endl;
      std::cout << "density: " << m_Forces[i].density << std::endl << std::endl;
      std::cout << "volume: " << m_Forces[i].volume << std::endl << std::endl;
      std::cout << "J     : " << m_Forces[i].J << std::endl << std::endl;
      std::cout << m_Forces[i].sigma << std::endl << std::endl;
    }
    
    std::cout << "gridRes: " << gridRes.transpose() << std::endl;
    std::cout << "minGrid: " << minGrid.transpose() << std::endl;
    std::cout << "maxGrid: " << maxGrid.transpose() << std::endl;
    std::cout << "h: " << h << std::endl;
  }
  
private:
  std::vector<SingleForceData2D> m_Forces;
  
  Eigen::Vector2i gridRes;
  Eigen::Vector2d minGrid;
  Eigen::Vector2d maxGrid;
  double h;
};


#endif
