#ifndef __IMAGE_IO_UTIL_H__
#define __IMAGE_IO_UTIL_H__

#include <stdint.h>
#include <cmath>
#include <iostream>
using namespace std;

template<typename T, int N>
struct color2RGB
{
  void operator()( T c[N], uint8_t rgb[3] ) {}
};

template<typename T>
struct color2RGB<T,1>
{
  void operator()( T c[1], uint8_t rgb[3] )
  {
    int v = max( 0, min( 255, int( c[0] * 256 ) ) );
    rgb[0] = rgb[1] = rgb[2] = v;
  }
};

template<typename T>
struct color2RGB<T,2>
{
  void operator()( T c[2], uint8_t rgb[3] )
  {
    int v1 = max( 0, min( 255, int( c[0] * 256 ) ) );
    int v2 = max( 0, min( 255, int( c[1] * 256 ) ) );
    rgb[0] = v1;
    rgb[1] = v2;
    rgb[2] = 0;
  }
};

template<typename T>
struct color2RGB<T,3>
{
  void operator()( T c[3], uint8_t rgb[3] )
  {
    int v1 = max( 0, min( 255, int( c[0] * 256 ) ) );
    int v2 = max( 0, min( 255, int( c[1] * 256 ) ) );
    int v3 = max( 0, min( 255, int( c[2] * 256 ) ) );
    rgb[0] = v1;
    rgb[1] = v2;
    rgb[2] = v3;
  }
};


#endif
