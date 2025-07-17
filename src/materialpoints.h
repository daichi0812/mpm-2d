//
//  materialpoints.h
//  SimpleMPM2D
//
//  Created by Yonghao Yue on 2020/09/04.
//  Copyright © 2020 YY. All rights reserved.
//

#ifndef materialpoints_h
#define materialpoints_h

#include <Eigen/Core>
#include <vector>

struct MaterialPoints
{
  std::vector<double> mass0;
  std::vector<double> density0;
  std::vector<double> volume0;
  std::vector<Eigen::Vector2d> x;
  std::vector<Eigen::Vector2d> v;
  std::vector<Eigen::Matrix2d> F;
  std::vector<Eigen::Matrix2d> Fe;
  std::vector<Eigen::Matrix2d> Fp;
  std::vector<Eigen::Matrix2d> tau;
  std::vector<Eigen::Matrix2d> be;
  
  std::vector<double> xi; // Thixotropy;
  
  int nPoints;
  int nPointsAlloc;
};

#endif /* materialpoints_h */
