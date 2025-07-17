#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <stdint.h>
#include <cstdlib>
#include <algorithm>

#include <iostream>
using namespace std;

template<typename T, int N>
class Image
{
public:
  Image()
  {
    m_Width = 0;
    m_Height = 0;
    m_Data = nullptr;
  }

  Image( const int32_t in_Width, const int32_t in_Height )
  {
    if( ( in_Width <= 0 ) || ( in_Height <= 0 ) )
    {
      m_Width = 0;
      m_Height = 0;
      m_Data = nullptr;
    }
    else
    {
      m_Width = in_Width;
      m_Height = in_Height;
      m_Data = (T*)malloc( sizeof(T) * in_Width * in_Height * N );
    }
  }
  
  ~Image()
  {
    m_Width = 0;
    m_Height = 0;
    if( m_Data != nullptr )
    {
      free( m_Data );
      m_Data = nullptr;
    }
  }
  
  void zeroClear()
  {
    memset( m_Data, 0, sizeof(T) * m_Width * m_Height * N );
  }
  
  void resize( const int32_t in_Width, const int32_t in_Height )
  {
    if( ( in_Width <= 0 ) || ( in_Height <= 0 ) )
    {
      std::cout << "Invalid parameters for image resizing. Given parameters: width = " << in_Width << ", height = " << in_Height << ". Both width and height should be positive integers. Image resizing will leave the current image unchanged." << std::endl;
    }
    else
    {
      m_Width = in_Width;
      m_Height = in_Height;
      m_Data = (T*)realloc( m_Data, sizeof(T) * in_Width * in_Height * N );
    }
  }
  
  template<int M>
  void resizeToFitTargetImage( const Image<T, M>& in_Target )
  {
    if( m_Width != in_Target.getWidth() || m_Height != in_Target.getHeight() )
      resize( in_Target.getWidth(), in_Target.getHeight() );
  }
  
  T operator()( const int i, const int j, const int c ) const
  {
    return m_Data[ j * m_Width * N + i * N + c ];
  }
  
  T& operator()( const int i, const int j, const int c )
  {
    return m_Data[ j * m_Width * N + i * N + c ];
  }
  
  T atUVNearest( const double s, const double t, const int c )
  {
    int32_t i = s * m_Width; int32_t j = t * m_Height;
    clampPixelIndex( i, j );
    return m_Data[ j * m_Width * N + i * N + c ];
  }
  
  void invertY()
  {
    for( int j=0; j<m_Height/2; j++ )
    {
      for( int i=0; i<m_Width*N; i++ )
      {
        double temp = m_Data[ j * m_Width * N + i ];
        m_Data[ j * m_Width * N + i ] = m_Data[ ( m_Height - 1 - j ) * m_Width * N + i ];
        m_Data[ ( m_Height - 1 - j ) * m_Width * N + i ] = temp;
      }
    }
  }
  
  void copyDataFrom( const Image<T, N>& in_Target )
  {
    if( m_Width != in_Target.m_Width || m_Height != in_Target.m_Height )
    {
      std::cout << "Image::copyDataFrom() requires the image sizes to be identical." << std::endl;
    }
    else
    {
      memcpy( m_Data, in_Target.m_Data, sizeof(T) * m_Width * m_Height * N );
    }
  }
  
  void clampPixelIndex( int32_t& io_i, int32_t& io_j ) const
  {
    io_i = max( 0, min( m_Width - 1, io_i ) );
    io_j = max( 0, min( m_Height - 1, io_j ) );
  }
  
  int32_t getWidth() const
  {
    return m_Width;
  }
  
  int32_t getHeight() const
  {
    return m_Height;
  }
  
  T* getDataPointer()
  {
    return m_Data;
  }
  
  const T* getConstDataPointer() const
  {
    return m_Data;
  }
  
protected:
	int32_t m_Width;
	int32_t m_Height;
	T* m_Data;

	enum {NC = N};
};

#endif
