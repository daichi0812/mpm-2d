//
//  basisfunction.h
//  SimpleMPM2D
//
//  Created by Yonghao Yue on 2020/09/04.
//  Copyright © 2020 YY. All rights reserved.
//

#ifndef basisfunction_h
#define basisfunction_h

#include <Eigen/Core>
#include "simulator.h"

class StencilIterator
{
  StencilIterator();
public:
  StencilIterator( const Eigen::Vector2i& in_Min, const Eigen::Vector2i& in_Max, const Eigen::Vector2i& in_Current )
   : m_Min( in_Min ), m_Max( in_Max ), m_Current( in_Current )
  {}
  
  StencilIterator( const StencilIterator& in_Other )
  {
    m_Min = in_Other.m_Min;
    m_Max = in_Other.m_Max;
    m_Current = in_Other.m_Current;
  }
  
  ~StencilIterator() {}
  
  StencilIterator& operator=( const StencilIterator& in_Other )
  {
    m_Min = in_Other.m_Min;
    m_Max = in_Other.m_Max;
    m_Current = in_Other.m_Current;
    return *this;
  }
  
  StencilIterator& operator++()
  {
    m_Current.x()++;
    if( m_Current.x() > m_Max.x() )
    {
      m_Current.x() = m_Min.x();
      m_Current.y()++;
    }
    
    if( m_Current.y() > m_Max.y() )
    {
      m_Current.y() = m_Max.y();
    }
    
    return *this;
  }
  
  StencilIterator operator++( int )
  {
    StencilIterator tmp(*this);
    operator++();
    return tmp;
  }
  
  const Eigen::Vector2i& operator*() const
  {
    return m_Current;
  }
  
  bool operator==( const StencilIterator& in_Other )
  {
    return m_Current.x() == in_Other.m_Current.x() && m_Current.y() == in_Other.m_Current.y();
  }
  
  bool operator!=( const StencilIterator& in_Other )
  {
    return m_Current.x() != in_Other.m_Current.x() || m_Current.y() != in_Other.m_Current.y();
  }
  
protected:
  Eigen::Vector2i m_Min;
  Eigen::Vector2i m_Max;
  Eigen::Vector2i m_Current;
};

class Stencil
{
  Stencil();
public:
  Stencil( const Eigen::Vector2i& in_Min, const Eigen::Vector2i& in_Max )
   : m_Min( in_Min ), m_Max( in_Max )
  {}
  
  StencilIterator begin() const
  {
    return StencilIterator( m_Min, m_Max, m_Min );
  }
  
  StencilIterator end() const
  {
    return StencilIterator( m_Min, m_Max, m_Max );
  }
  
protected:
  Eigen::Vector2i m_Min;
  Eigen::Vector2i m_Max;
};

class LinearBasisFunction
{
public:
  LinearBasisFunction()
  {
    
  }
  
  double weight( const Eigen::Vector2d& in_xp, const Eigen::Vector2i& in_GridIdx, const Grid& in_Grid ) const
  {
    Eigen::Vector2d dx = in_xp - in_GridIdx.cast<double>() * in_Grid.h - in_Grid.minGrid;
    return N( dx.x() / in_Grid.h ) * N( dx.y() / in_Grid.h );
  }
  
  double weight( const Eigen::Vector2d& in_xp, const Eigen::Vector2d& in_xg, const Grid& in_Grid ) const
  {
    Eigen::Vector2d dx = in_xp - in_xg;
    return N( dx.x() / in_Grid.h ) * N( dx.y() / in_Grid.h );
  }

  Eigen::Vector2d weightGrad( const Eigen::Vector2d& in_xp, const Eigen::Vector2i& in_GridIdx, const Grid& in_Grid ) const
  {
    Eigen::Vector2d dx = in_xp - in_GridIdx.cast<double>() * in_Grid.h - in_Grid.minGrid;
    
    double wx = N( dx.x() / in_Grid.h );
    double wy = N( dx.y() / in_Grid.h );
    
    return Eigen::Vector2d{ wy * dN( dx.x() / in_Grid.h ) / in_Grid.h, wx * dN( dx.y() / in_Grid.h ) / in_Grid.h };
  }
  
  Stencil getStencil( const Eigen::Vector2d& in_xp, const Grid& in_Grid ) const
  {
    Eigen::Vector2i im = floor( ( ( in_xp - in_Grid.minGrid ) / in_Grid.h ).array() ).cast<int>() - 1;
    Eigen::Vector2i iM = ceil( ( ( in_xp - in_Grid.minGrid ) / in_Grid.h ).array() ).cast<int>() + 1;
    const Eigen::Vector2i grid_max = in_Grid.gridRes.array() - 1;
    return Stencil{ im.cwiseMax(0).cwiseMin( grid_max ), iM.cwiseMax(0).cwiseMin( grid_max ) };
  }
  
  Stencil getCenterStencil( const Eigen::Vector2d& in_xp, const Grid& in_Grid ) const
  {
    Eigen::Vector2i im = floor( ( ( in_xp - in_Grid.minGrid ) / in_Grid.h ).array() ).cast<int>() - 2;
    Eigen::Vector2i iM = floor( ( ( in_xp - in_Grid.minGrid ) / in_Grid.h ).array() ).cast<int>() + 2;
    const Eigen::Vector2i grid_max = in_Grid.gridRes.array() - 1;
    return Stencil{ im.cwiseMax(0).cwiseMin( grid_max ), iM.cwiseMax(0).cwiseMin( grid_max ) };
  }
  
protected:
  double N( double x ) const
  {
    double absx = fabs(x);
    if( absx < 1.0 ) return absx*absx*absx*0.5 - absx*absx + 2.0/3.0;
    else if( absx < 2.0 ) return -absx*absx*absx/6.0 + absx*absx - 2.0*absx + 4.0/3.0;
    else return 0.0;
  }

  double dN( double x ) const
  {
    double absx = fabs(x);
    double sgnx = (x>=0.0) ? 1.0 : -1.0;
    
    if( absx < 1.0 ) return ( absx*absx*1.5 - 2.0*absx ) * sgnx;
    else if( absx < 2.0 ) return ( -absx*absx*0.5 + 2.0*absx - 2.0 ) * sgnx;
    else return 0.0;
  }
};


#endif /* basisfunction_h */
