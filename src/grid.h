//
//  grid.h
//  SimpleMPM2D
//
//  Created by Yonghao Yue on 2020/09/04.
//  Copyright © 2020 YY. All rights reserved.
//

#ifndef grid_h
#define grid_h

#include <Eigen/Core>

struct Grid
{
  std::vector<double> mass;
  std::vector<Eigen::Vector2d> v;
  std::vector<Eigen::Vector2d> f;
  std::vector<Eigen::Vector2d> vstar;
  std::vector<Eigen::Vector2d> vnew;

  Eigen::Vector2i gridRes;
  Eigen::Vector2d minGrid;
  Eigen::Vector2d maxGrid;
  double h;
  
  int flatIdx( const Eigen::Vector2i& in_GridIdx ) const
  {
    return in_GridIdx.y() * ( gridRes.x() + 1 ) + in_GridIdx.x();
  }
  
  int flatCenterIdx( const Eigen::Vector2i& in_GridIdx ) const
  {
    return in_GridIdx.y() * gridRes.x() + in_GridIdx.x();
  }
  
  Eigen::Vector2i gridIdx( const int in_FlatIdx ) const
  {
    return Eigen::Vector2i{ in_FlatIdx % ( gridRes.x() + 1 ), in_FlatIdx / ( gridRes.x() + 1 ) };
  }
  
  Eigen::Vector2i gridCenterIdx( const Eigen::Vector2d in_x ) const
  {
    return Eigen::Vector2i{ in_x.x() / h, in_x.y() / h };
  }
  
  bool isInsideGrid( const Eigen::Vector2i& in_GridIdx ) const
  {
    return in_GridIdx.x() < gridRes.x() && in_GridIdx.y() < gridRes.y();
  }
  
  bool isInsideGrid( const Eigen::Vector2d& in_x ) const
  {
    const Eigen::Vector2i gridIdx = gridCenterIdx( in_x );
    return gridIdx.x() < gridRes.x() && gridIdx.y() < gridRes.y();
  }
  
  int numGridPoints() const
  {
    return ( gridRes.x() + 1 ) * ( gridRes.y() + 1 );
  }
  
  void setGrid( const Eigen::Vector2d& in_Min, const Eigen::Vector2d& in_Max, const Eigen::Vector2i& in_Res )
  {
    const Eigen::Vector2d _h = ( in_Max - in_Min ).array() / in_Res.cast<double>().array();
    h = _h.minCoeff();
    const Eigen::Vector2d center = ( in_Max + in_Min ) * 0.5;
    
    gridRes = in_Res;
    minGrid = center - 0.5 * h * in_Res.cast<double>();
    maxGrid = center + 0.5 * h * in_Res.cast<double>();
    
    mass.resize( numGridPoints() );
    v.resize( numGridPoints() );
    f.resize( numGridPoints() );
    vstar.resize( numGridPoints() );
    vnew.resize( numGridPoints() );
  }
  
  void setGrid( const Eigen::Vector2d& in_gridStart, const double in_h, const Eigen::Vector2i& in_Res )
  {
    h = in_h;
    gridRes = in_Res;

    minGrid = in_gridStart;
    maxGrid = in_gridStart + h * gridRes.cast<double>();
    
    mass.resize( numGridPoints() );
    v.resize( numGridPoints() );
    f.resize( numGridPoints() );
    vstar.resize( numGridPoints() );
    vnew.resize( numGridPoints() );
  }
  
  void clear()
  {
    const int n = numGridPoints();
    for( int i=0; i<n; i++ )
    {
      mass[i] = 0.0;
      v[i].setZero();
      f[i].setZero();
      vstar[i].setZero();
      vnew[i].setZero();
    }
  }
  
  double singleCellArea() const
  {
    return ( ( maxGrid.x() - minGrid.x() ) / gridRes.x() ) * ( ( maxGrid.y() - minGrid.y() ) / gridRes.y() );
  }
};

#endif /* grid_h */
