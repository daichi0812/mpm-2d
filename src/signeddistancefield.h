#ifndef __SIGNED_DISTANCE_FIELD_H__
#define __SIGNED_DISTANCE_FIELD_H__

#include "image.h"

void generateSDF( const Image<double, 1>& in_BinaryImage, Image<double, 1>& out_SDF );

#endif
