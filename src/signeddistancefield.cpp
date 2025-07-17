#include "signeddistancefield.h"
#include <cmath>
using namespace std;

//http://www.codersnotes.com/notes/signed-distance-fields

struct SDFPoint
{
	double dx, dy;
	double distSq() const { return dx * dx + dy * dy; }
};

struct SDFGrid
{
	SDFPoint* p;
	int width;
	int height;
};

const SDFPoint INSIDE = { 0, 0 };
const SDFPoint EMPTY = { 1<<14, 1<<14 };

SDFPoint getPoint( SDFGrid& g, int x, int y )
{
	if( ( x>=0 ) && ( y>=0 ) && ( x<g.width ) && ( y<g.height ) )
		return g.p[ y * g.width + x ];
	else
		return EMPTY;
}

void putPoint( SDFGrid& g, int x, int y, const SDFPoint& p )
{
	g.p[ y * g.width + x ] = p;
}

void compareValue( SDFGrid& g, SDFPoint& p, int x, int y, int offsetx, int offsety )
{
	SDFPoint other = getPoint( g, x + offsetx, y + offsety );
	other.dx += offsetx;
	other.dy += offsety;
	
	if( other.distSq() < p.distSq() )
		p = other;
}

void generateSDFGrid( SDFGrid& g )
{
	for( int y=0; y<g.height; y++ )
	{
		for( int x=0; x<g.width; x++ )
		{
			SDFPoint p = getPoint( g, x, y );
			compareValue( g, p, x, y, -1,  0 );
			compareValue( g, p, x, y,  0, -1 );
			compareValue( g, p, x, y, -1, -1 );
			compareValue( g, p, x, y,  1, -1 );
			putPoint( g, x, y, p );
		}
		
		for( int x=g.width-1; x>=0; x-- )
		{
			SDFPoint p = getPoint( g, x, y );
			compareValue( g, p, x, y, 1, 0 );
			putPoint( g, x, y, p );
		}
	}
	
	for( int y=g.height-1; y>=0; y-- )
	{
		for( int x=g.width-1; x>=0; x-- )
		{
			SDFPoint p = getPoint( g, x, y );
			compareValue( g, p, x, y,  1,  0 );
			compareValue( g, p, x, y,  0,  1 );
			compareValue( g, p, x, y, -1,  1 );
			compareValue( g, p, x, y,  1,  1 );
			putPoint( g, x, y, p );
		}
		
		for( int x=0; x<g.width; x++ )
		{
			SDFPoint p = getPoint( g, x, y );
			compareValue( g, p, x, y, -1, 0 );
			putPoint( g, x, y, p );
		}
	}
}

void generateSDF( const Image<double, 1>& in_BinaryImage, Image<double, 1>& out_SDF )
{
  out_SDF.resizeToFitTargetImage( in_BinaryImage );
		
	SDFGrid g1, g2;
	g1.width = in_BinaryImage.getWidth();
	g1.height = in_BinaryImage.getHeight();
	g1.p = (SDFPoint*)malloc( sizeof(SDFPoint) * g1.width * g1.height );
	
	g2.width = in_BinaryImage.getWidth();
	g2.height = in_BinaryImage.getHeight();
	g2.p = (SDFPoint*)malloc( sizeof(SDFPoint) * g2.width * g2.height );
	
	for( int y=0; y<in_BinaryImage.getHeight(); y++ )
	{
		for( int x=0; x<in_BinaryImage.getWidth(); x++ )
		{
			if( in_BinaryImage( x, y, 0 ) < 0.5 )
			{
				putPoint( g1, x, y, INSIDE );
				putPoint( g2, x, y, EMPTY );
			}
			else
			{
				putPoint( g2, x, y, INSIDE );
				putPoint( g1, x, y, EMPTY );
			}
		}
	}
	
	generateSDFGrid( g1 );
	generateSDFGrid( g2 );
	
	for( int y=0; y<in_BinaryImage.getHeight(); y++ )
	{
		for( int x=0; x<in_BinaryImage.getWidth(); x++ )
		{
			double dist1 = sqrt( getPoint( g1, x, y ).distSq() );
			double dist2 = sqrt( getPoint( g2, x, y ).distSq() );
			double dist = dist1 - dist2;
			out_SDF( x, y, 0 ) = dist;
		}
	}
	
	free( g1.p );
	free( g2.p );
}
