//
//  constitutivemodel.h
//  SimpleMPM2D
//
//  Created by Yonghao Yue on 2020/09/04.
//  Copyright © 2020 YY. All rights reserved.
//

#ifndef constitutivemodel_h
#define constitutivemodel_h

#include <Eigen/Core>
#include <Eigen/SVD>

struct MaterialPoints;

class ConstitutiveModel
{
public:
  virtual Eigen::Matrix2d computeKirchhoffStress( const MaterialPoints& in_Points, const int in_Point_idx ) const = 0;
  
  virtual void updateDeformationStatus( MaterialPoints& io_Points, const int in_Point_idx, const double in_dt, const Eigen::Matrix2d& in_VelocityGradient ) const = 0;
};

class SnowModel final : public ConstitutiveModel
{
  SnowModel();
public:
  SnowModel( const double in_Lambda, const double in_Mu, const double in_Theta_c, const double in_Theta_s, const double in_Xi );
  
  Eigen::Matrix2d computeKirchhoffStress( const MaterialPoints& in_Points, const int in_Point_idx ) const override;
  
  void updateDeformationStatus( MaterialPoints& io_Points, const int in_Point_idx, const double in_dt, const Eigen::Matrix2d& in_VelocityGradient ) const override;

protected:
  double m_Lambda;
  double m_Mu;
  double m_Theta_c;
  double m_Theta_s;
  double m_Xi;
};

class HerschelBulkley final : public ConstitutiveModel
{
  HerschelBulkley();
public:
  HerschelBulkley( const double in_kappa, const double in_mu, const double in_n, const double in_sigma_y, const double in_eta );
  
  Eigen::Matrix2d computeKirchhoffStress( const MaterialPoints& in_Points, const int in_Point_idx ) const override;
  
  void updateDeformationStatus( MaterialPoints& io_Points, const int in_Point_idx, const double in_dt, const Eigen::Matrix2d& in_VelocityGradient ) const override;
  
protected:
  inline double Phi( const double s ) const
  {
    return ( s - m_sigma_y ) / m_eta; // for 2D
  }
  
  inline double gamma_dot( const double s ) const
  {
    return pow( std::max( 0.0, Phi(s) ), 1.0 / m_n );
  }
  
  inline double d_gamma_dot_d_s( const double s ) const
  {
    return pow( gamma_dot(s), 1.0 - m_n ) / ( m_n * m_eta );
  }
  
  inline double f( const double s, const double s0, const double dt, const double mu_hat ) const
  {
    return s - s0 + 2.0 * dt * mu_hat * gamma_dot(s);
  }
  
  inline double df( const double s, const double s0, const double dt, const double mu_hat ) const
  {
    return 1.0 + 2.0 * dt * mu_hat * d_gamma_dot_d_s(s);
  }
  
  inline double solve_for_s( const double s, const double dt, const double mu_hat ) const
  {
    //solve for f(s) = 0
    double s_n = s;
    
    while(1)
    {
      const double _df = df( s_n, s, dt, mu_hat );
      //std::cout << "  _df: " << _df << std::endl;
      
      const double s_n1 = s_n - ( 1.0 / _df ) * f( s_n, s, dt, mu_hat );
      
      if( fabs( s_n - s_n1 ) < 1.0e-8 )
        return s_n1;
      
      s_n = s_n1;
    }
    
    return s_n;
  }
  
  double m_kappa;
  double m_mu;
  double m_n;
  double m_sigma_y;
  double m_eta;
};

class ThixotropicHerschelBulkley final : public ConstitutiveModel
{
  ThixotropicHerschelBulkley();
public:
  ThixotropicHerschelBulkley( const double in_kappa, const double in_mu, const double in_n, const double in_sigma_y, const double in_eta, const double in_k1, const double in_k2 );
  
  Eigen::Matrix2d computeKirchhoffStress( const MaterialPoints& in_Points, const int in_Point_idx ) const override;
  
  void updateDeformationStatus( MaterialPoints& io_Points, const int in_Point_idx, const double in_dt, const Eigen::Matrix2d& in_VelocityGradient ) const override;
  
protected:
  inline double s( const Eigen::Vector2d& v ) const
  {
    return v.x();
  }
  
  inline double xi( const Eigen::Vector2d& v ) const
  {
    return v.y();
  }
  
  inline double Phi( const Eigen::Vector2d& v ) const
  {
    return ( s(v) - ( 1.0 + xi(v) ) * m_sigma_y )/( ( 1.0 + xi(v) ) * m_eta );
  }
  
  inline double gamma_dot( const Eigen::Vector2d& v ) const
  {
    return pow( std::max( 0.0, Phi(v) ), 1.0 / m_n );
  }
  
  inline double d_gamma_dot_d_s( const Eigen::Vector2d& v ) const
  {
    return pow( gamma_dot(v), 1.0 - m_n ) / ( m_n * m_eta * ( 1.0 + xi(v) ) );
  }

  inline double d_gamma_dot_d_xi( const Eigen::Vector2d& v ) const
  {
    return -pow( gamma_dot(v), 1.0 - m_n ) * s(v) / ( m_n * m_eta * ( 1.0 + xi(v) ) * ( 1.0 + xi(v) ) );
  }
  
  inline Eigen::Vector2d f( const Eigen::Vector2d& v, const Eigen::Vector2d& v0, const double h, const double mu_hat ) const
  {
    Eigen::Vector2d res;
    res.x() = s(v) - s(v0) + 2.0 * h * mu_hat * gamma_dot(v);
    res.y() = xi(v) - xi(v0) + h * m_k1 * xi(v) * gamma_dot(v) - h * m_k2 * ( 1.0 - xi(v) );
    return res;
  }
  
  inline Eigen::Matrix2d df( const Eigen::Vector2d& v, const Eigen::Vector2d& v0, const double h, const double mu_hat ) const
  {
    Eigen::Matrix2d res;
    res( 0, 0 ) = 1.0 + 2.0 * h * mu_hat * d_gamma_dot_d_s(v);
    res( 0, 1 ) = 2.0 * h * mu_hat * d_gamma_dot_d_xi(v);
    res( 1, 0 ) = h * m_k1 * xi(v) * d_gamma_dot_d_s(v);
    res( 1, 1 ) = 1.0 + h * m_k1 * ( gamma_dot(v) + xi(v) * d_gamma_dot_d_xi(v) ) + h * m_k2;
    
    return res;
  }
  
  inline Eigen::Vector2d update( const Eigen::Vector2d& v, const double h, const double mu_hat ) const
  {
    //solve for f(v) = 0
    Eigen::Vector2d v_n = v;
    
    while(1)
    {
      const Eigen::Matrix2d _df = df( v_n, v, h, mu_hat );
      //std::cout << "  _df: " << _df << std::endl;
      
      const Eigen::Vector2d v_n1 = v_n - _df.inverse() * f( v_n, v, h, mu_hat );
      
      if( ( v_n - v_n1 ).norm() < 1.0e-8 )
        return v_n1;
      
      v_n = v_n1;
    }
    
    return v_n;
  }
  
  double m_kappa;
  double m_mu;
  double m_n;
  double m_sigma_y;
  double m_eta;
  double m_k1;
  double m_k2;
};

#endif /* constitutivemodel_h */
