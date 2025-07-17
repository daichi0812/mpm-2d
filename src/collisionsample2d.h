#ifndef collision_sample_h
#define collision_sample_h

#include <Eigen/Core>
#include <map>

struct BoundingBox2D
{
  BoundingBox2D() {}
  BoundingBox2D( const Eigen::Vector2d& _min, const Eigen::Vector2d& _max )
    : bb_min( _min ), bb_max( _max ) {}
  
  Eigen::Vector2d bb_min;
  Eigen::Vector2d bb_max;
};

struct CollisionInfo2D
{
  int frame_idx;
  double penetrationDepth;
  Eigen::Vector2d normal;
  Eigen::Vector2d tangent;
  Eigen::Vector2d delta;
};


struct CollisionSample2D
{
  CollisionSample2D( const Eigen::Vector2d& _x0 )
   : x0(_x0)
  {}
  
  Eigen::Vector2d x0;
  std::map< int, CollisionInfo2D > collision_cache;
};

#endif
