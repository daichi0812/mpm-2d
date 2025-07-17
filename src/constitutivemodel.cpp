//
//  constitutivemodel.cpp
//  SimpleMPM2D
//
//  Created by Yonghao Yue on 2020/09/04.
//  Copyright © 2020 YY. All rights reserved.
//

#include "constitutivemodel.h"
#include "materialpoints.h"

#include <Eigen/Dense>

inline Eigen::Matrix2d dev( const Eigen::Matrix2d& in_x )
{
  return in_x - in_x.trace() * 0.5 * Eigen::Matrix2d::Identity(); //for 2D
}

SnowModel::SnowModel( const double in_Lambda, const double in_Mu, const double in_Theta_c, const double in_Theta_s, const double in_Xi )
 : m_Lambda( in_Lambda ), m_Mu( in_Mu ), m_Theta_c( in_Theta_c ), m_Theta_s( in_Theta_s ), m_Xi( in_Xi )
{}

Eigen::Matrix2d SnowModel::computeKirchhoffStress( const MaterialPoints& in_Points, const int in_Point_idx ) const
{
  double Je = in_Points.Fe[in_Point_idx].determinant();
  double Jp = in_Points.Fp[in_Point_idx].determinant();
  
  Eigen::JacobiSVD< Eigen::Matrix2d > svd( in_Points.Fe[in_Point_idx], Eigen::ComputeFullU  | Eigen::ComputeFullV );
  // compute polar decomposition F = RS = (UP);
  Eigen::Matrix2d Re = svd.matrixU() * svd.matrixV().transpose();
  //Eigen::Matrix2d Se = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
  
  double _mu = m_Mu * exp( m_Xi * ( 1.0 - Jp ) );
  double _lambda = m_Lambda * exp( m_Xi * ( 1.0 - Jp ) );
  
  return 2.0 * _mu * ( in_Points.Fe[in_Point_idx] - Re ) * in_Points.Fe[in_Point_idx].transpose() + _lambda * ( Je - 1.0 ) * Je * Eigen::Matrix2d::Identity();
}

void SnowModel::updateDeformationStatus( MaterialPoints& io_Points, const int in_Point_idx, const double in_dt, const Eigen::Matrix2d& in_VelocityGradient ) const
{
  Eigen::Matrix2d Fe = io_Points.Fe[in_Point_idx] + in_dt * in_VelocityGradient * io_Points.Fe[in_Point_idx];
  //Eigen::Matrix2d F = Fe * io_Fp;

  Eigen::JacobiSVD< Eigen::Matrix2d > svd( Fe, Eigen::ComputeFullU | Eigen::ComputeFullV );
  Eigen::Matrix2d S = svd.singularValues().asDiagonal();
  S( 0, 0 ) = std::max<double>( 1.0 - m_Theta_c, std::min<double>( 1.0 + m_Theta_s, S( 0, 0 ) ) );
  S( 1, 1 ) = std::max<double>( 1.0 - m_Theta_c, std::min<double>( 1.0 + m_Theta_s, S( 1, 1 ) ) );

  io_Points.Fe[in_Point_idx] = svd.matrixU() * S * svd.matrixV().transpose();
  io_Points.Fp[in_Point_idx] = svd.matrixV() * S.inverse() * svd.matrixU().transpose() * Fe * io_Points.Fp[in_Point_idx];
  
  io_Points.F[in_Point_idx] = io_Points.Fe[in_Point_idx] * io_Points.Fp[in_Point_idx];
}



HerschelBulkley::HerschelBulkley( const double in_kappa, const double in_mu, const double in_n, const double in_sigma_y, const double in_eta )
 : m_kappa( in_kappa ), m_mu( in_mu ), m_n( in_n ), m_sigma_y( in_sigma_y ), m_eta( in_eta )
{}

Eigen::Matrix2d HerschelBulkley::computeKirchhoffStress( const MaterialPoints& in_Points, const int in_Point_idx ) const
{
  const double J = sqrt( in_Points.be[in_Point_idx].determinant() );
  const Eigen::Matrix2d be_bar = in_Points.be[in_Point_idx] / J; // for 2D
  return 0.5 * m_kappa * ( J*J - 1.0 ) * Eigen::Matrix2d::Identity() + m_mu * dev( be_bar );
}

void HerschelBulkley::updateDeformationStatus( MaterialPoints& io_Points, const int in_Point_idx, const double in_dt, const Eigen::Matrix2d& in_VelocityGradient ) const
{
  const Eigen::Matrix2d f = Eigen::Matrix2d::Identity() + in_dt * in_VelocityGradient;
  const Eigen::Matrix2d be_star = f * io_Points.be[in_Point_idx] * f.transpose();
  const double J_star = sqrt( be_star.determinant() );
  const Eigen::Matrix2d be_star_bar = be_star / J_star; // for 2D
  
  const double mu_hat = 0.5 * be_star_bar.trace() * m_mu; // for 2D
  const Eigen::Matrix2d s_star = m_mu * dev( be_star_bar );
  const double scalar_s_star = s_star.norm();
  const Eigen::Matrix2d s_hat_star = s_star / scalar_s_star;
  
  const double s_new = solve_for_s( scalar_s_star, in_dt, mu_hat );
  
  Eigen::Matrix2d be_new_bar = s_new * s_hat_star / m_mu + 0.5 * be_star_bar.trace() * Eigen::Matrix2d::Identity(); // for 2D
  const double J_new_bar = sqrt( be_new_bar.determinant() );
  be_new_bar = be_new_bar / J_new_bar; // for 2D
  
  io_Points.be[in_Point_idx] = J_star * be_new_bar; // for 2D
}


ThixotropicHerschelBulkley::ThixotropicHerschelBulkley( const double in_kappa, const double in_mu, const double in_n, const double in_sigma_y, const double in_eta, const double in_k1, const double in_k2 )
: m_kappa( in_kappa ), m_mu( in_mu ), m_n( in_n ), m_sigma_y( in_sigma_y ), m_eta( in_eta ), m_k1( in_k1 ), m_k2( in_k2 )
{}

Eigen::Matrix2d ThixotropicHerschelBulkley::computeKirchhoffStress( const MaterialPoints& in_Points, const int in_Point_idx ) const
{
  const double J = sqrt( in_Points.be[in_Point_idx].determinant() );
  const double xi = in_Points.xi[in_Point_idx];
  const Eigen::Matrix2d be_bar = in_Points.be[in_Point_idx] / J; // for 2D
  return 0.5 * m_kappa * ( J*J - 1.0 ) * Eigen::Matrix2d::Identity() + ( 1.0 + xi ) * m_mu * dev( be_bar );
}

void ThixotropicHerschelBulkley::updateDeformationStatus( MaterialPoints& io_Points, const int in_Point_idx, const double in_dt, const Eigen::Matrix2d& in_VelocityGradient ) const
{
  const Eigen::Matrix2d f = Eigen::Matrix2d::Identity() + in_dt * in_VelocityGradient;
  const Eigen::Matrix2d be_star = f * io_Points.be[in_Point_idx] * f.transpose();
  const double J_star = sqrt( be_star.determinant() );
  const Eigen::Matrix2d be_star_bar = be_star / J_star; // for 2D
  
  const double mu_hat = 0.5 * be_star_bar.trace() * m_mu; // for 2D
  const Eigen::Matrix2d s_star = m_mu * dev( be_star_bar );
  const double scalar_s_star = s_star.norm();
  const Eigen::Matrix2d s_hat_star = s_star / scalar_s_star;
  
  const double xi = io_Points.xi[in_Point_idx];
  const Eigen::Vector2d _v{ scalar_s_star, xi };
  const Eigen::Vector2d _v_new = update( _v, in_dt, mu_hat );
  
  const double s_new = _v_new(0);
  
  Eigen::Matrix2d be_new_bar = s_new * s_hat_star / m_mu + 0.5 * be_star_bar.trace() * Eigen::Matrix2d::Identity(); // for 2D
  const double J_new_bar = sqrt( be_new_bar.determinant() );
  be_new_bar = be_new_bar / J_new_bar; // for 2D
  
  io_Points.be[in_Point_idx] = J_star * be_new_bar; // for 2D
  io_Points.xi[in_Point_idx] = _v_new(1);
}
