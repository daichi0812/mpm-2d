#ifndef __IMAGE_IO_BMP_H__
#define __IMAGE_IO_BMP_H__

#include "image.h"

#include "endianutil.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

#include "imageioutil.h"

template<typename T, int N>
void saveImageBmp( const char* in_Filename, const Image<T, N>& in_Image )
{
  ofstream ofs;
  ofs.open( in_Filename, ios::out | ios::trunc | ios::binary );
  
  uint8_t uc;
  uint16_t us;
  uint32_t ul;
  uint32_t file_size;
  uint32_t data_size;
  
  uint32_t line_offset = ( 4 - ( in_Image.getWidth() * 3 ) % 4 ) % 4;
  data_size = ( in_Image.getWidth() + line_offset ) * in_Image.getHeight() * 3;
  file_size = data_size + 54;
  
  uc = 'B';			ofs.write( (char*)&uc, sizeof(uc) );
  uc = 'M';			ofs.write( (char*)&uc, sizeof(uc) );
  
  file_size = LSB4( file_size );
  ofs.write( (char*)&file_size, sizeof(file_size) );
  us = 0;
  us = LSB2(us);				ofs.write( (char*)&us, sizeof(us) );
  us = 0;
  us = LSB2(us);				ofs.write( (char*)&us, sizeof(us) );
  ul = 54;
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  ul = 40;
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  ul = in_Image.getWidth();
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  ul = in_Image.getHeight();
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  us = 1;
  us = LSB2(us);				ofs.write( (char*)&us, sizeof(us) );
  us = 24;
  us = LSB2(us);				ofs.write( (char*)&us, sizeof(us) );
  ul = 0;
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  data_size = LSB4( data_size );
  ofs.write( (char*)&data_size, sizeof(data_size) );
  ul = 0;
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  ul = 0;
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  ul = 0;
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  ul = 0;
  ul = LSB4(ul);				ofs.write( (char*)&ul, sizeof(ul) );
  
  for( int j=0; j<in_Image.getHeight(); j++ )
  {
    for( int i=0; i<in_Image.getWidth(); i++ )
    {
      T color[3] = { in_Image( i, in_Image.getHeight()-1-j, 0 ), in_Image( i, in_Image.getHeight()-1-j, 1 ), in_Image( i, in_Image.getHeight()-1-j, 2 ) };
      uint8_t rgb[3];
      color2RGB<T, N>()( color, rgb );
      ofs.write( (char*)&rgb[2], sizeof(uint8_t) );
      ofs.write( (char*)&rgb[1], sizeof(uint8_t) );
      ofs.write( (char*)&rgb[0], sizeof(uint8_t) );
    }
    for( unsigned int i=0; i<line_offset; i++ )
    {
      uc = 0;
      ofs.write( (char*)&uc, sizeof(uc) );
    }
  }
  
  ofs.close();
}

#endif
