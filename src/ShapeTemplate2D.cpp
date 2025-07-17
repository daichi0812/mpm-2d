//
//  ShapeTemplate2D.cpp
//  rigidbody2dsim
//
//  Created by Yonghao Yue on 2022/04/02.
//  Copyright © 2022 Yonghao Yue. All rights reserved.
//

#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE

#include "ShapeTemplate2D.h"
#include "CircleTemplate2D.h"
#include "PolygonTemplate2D.h"

ShapeTemplate2D* initializeTemplate( const ShapeTemplateData& in_data )
{
  if( in_data.vertex_list.cols() < 3 )
    return new CircleTemplate2D( in_data.name, in_data.density, in_data.size_mean * 0.5, in_data.size_std * 0.5 );
  else
    return new PolygonTemplate2D( in_data.name, in_data.density, in_data.vertex_list, in_data.size_std );
}

