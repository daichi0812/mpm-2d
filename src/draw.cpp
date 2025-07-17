#define _CRT_SECURE_NO_WARNINGS

#include "draw.h"
#include "mpm.h"
#include "simulator.h"
#include "window.h"

#include "imageiobmp.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include <unistd.h>
#else
#include <GL/freeglut.h>
#endif

DisplayParam g_DisplayParam;

void drawPoints()
{
	glColor3f( 0.0, 0.0, 0.0 );
	glBegin( GL_POINTS );
	for( int i=0; i<g_MPM.simulator.points.nPoints; i++ )
	{
    glVertex2f( g_MPM.simulator.points.x[i].x(), g_MPM.simulator.points.x[i].y() );
	}
	glEnd();
}

void drawGrids()
{
	glColor3f(0.9, 0.9, 0.9);
	glBegin(GL_LINES);
	
  for( int i=0; i<=g_MPM.simulator.grid.gridRes.x(); i++ )
	{
    double x = g_MPM.simulator.grid.minGrid.x() + ( g_MPM.simulator.grid.maxGrid.x() - g_MPM.simulator.grid.minGrid.x() ) * double(i) / double( g_MPM.simulator.grid.gridRes.x() );
    glVertex2f(x, g_MPM.simulator.grid.minGrid.y() );
    glVertex2f(x, g_MPM.simulator.grid.maxGrid.y() );
	}
	
  for( int i=0; i<=g_MPM.simulator.grid.gridRes.y(); i++ )
	{
    double y = g_MPM.simulator.grid.minGrid.y() + (g_MPM.simulator.grid.maxGrid.y() - g_MPM.simulator.grid.minGrid.y() ) * double(i) / double( g_MPM.simulator.grid.gridRes.y() );
    glVertex2f( g_MPM.simulator.grid.minGrid.x(), y );
    glVertex2f( g_MPM.simulator.grid.maxGrid.x(), y );
	}
	
	glEnd();
}

void drawObstacles()
{
	glBegin( GL_QUADS );
	
	for( int j=0; j<g_MPM.simulator.obstacleSDF.getHeight(); j++ )
	{
		for( int i=0; i<g_MPM.simulator.obstacleSDF.getWidth(); i++ )
		{
			double v = g_MPM.simulator.obstacleSDF( i, j, 0 );
			double adh = g_MPM.simulator.adhesionMap( i, j, 0 );
      
			if( v > 0.0 ) continue;
			else
			{
				if( adh > 0.5 ) glColor3f( 0.7, 0.3, 0.3 );
				else glColor3f( 0.7, 0.7, 0.7 );
			}
			
			double s0 = double(i) / g_MPM.simulator.obstacleSDF.getWidth();
			double s1 = double(i+1) / g_MPM.simulator.obstacleSDF.getWidth();
			double t0 = double(j) / g_MPM.simulator.obstacleSDF.getHeight();
			double t1 = double(j+1) / g_MPM.simulator.obstacleSDF.getHeight();
			
      double x0 = s0 * ( g_MPM.simulator.grid.maxGrid.x() - g_MPM.simulator.grid.minGrid.x() ) + g_MPM.simulator.grid.minGrid.x();
      double x1 = s1 * ( g_MPM.simulator.grid.maxGrid.x() - g_MPM.simulator.grid.minGrid.x() ) + g_MPM.simulator.grid.minGrid.x();
      double y0 = t0 * ( g_MPM.simulator.grid.maxGrid.y() - g_MPM.simulator.grid.minGrid.y() ) + g_MPM.simulator.grid.minGrid.y();
      double y1 = t1 * ( g_MPM.simulator.grid.maxGrid.y() - g_MPM.simulator.grid.minGrid.y() ) + g_MPM.simulator.grid.minGrid.y();
			
			glVertex2f( x0, y0 );
			glVertex2f( x1, y0 );
			glVertex2f( x1, y1 );
			glVertex2f( x0, y1 );
		}
	}
	
	glEnd();
}

void draw()
{
  double window_width = g_Window.width * g_Window.framesize_windowsize_scale_x;
  double window_height = g_Window.height * g_Window.framesize_windowsize_scale_y;
  
  glClearColor( 1.0, 1.0, 1.0, 1.0 );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  
	glPushMatrix();

  glViewport( 0, 0, window_width, window_height );
  
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	int mins = min( window_width, window_height );
  
  double xc = ( g_DisplayParam.right + g_DisplayParam.left ) * 0.5; double yc = ( g_DisplayParam.top + g_DisplayParam.bottom ) * 0.5;
  double hw = ( g_DisplayParam.right - g_DisplayParam.left ) * 0.5; double hh = ( g_DisplayParam.top - g_DisplayParam.bottom ) * 0.5;
	double w = max(hw, hh);
	
	glOrtho( ( xc - 1.05 * w ) * double(window_width) / mins, ( xc + 1.05 * w ) * double(window_width) / mins,
		( yc - 1.05 * w ) * double(window_height) / mins, ( yc + 1.05 * w ) * double(window_height) / mins, -1.0, 2.0 );
  
	drawGrids();
  drawObstacles();
	drawPoints();

	glPopMatrix();	
	
	if( g_MPM.capture )
	{
		if( g_MPM.simulator.timestep * g_MPM.simulator.deltaT - g_MPM.lastCaptureTime < g_MPM.captureTimeStep ) return;
		g_MPM.lastCaptureTime = g_MPM.simulator.timestep * g_MPM.simulator.deltaT;
		
		if( ( window_width != g_MPM.captureImage.getWidth() ) || ( window_height != g_MPM.captureImage.getHeight() ) )
      g_MPM.captureImage.resize( window_width, window_height );
			
		glReadPixels( 0, 0, window_width, window_height, GL_RGB, GL_FLOAT, g_MPM.captureImage.getDataPointer() );
		
		char buf[256];
		sprintf( buf, "./capture_%d.bmp", g_MPM.captureID++ );
    g_MPM.captureImage.invertY();
		saveImageBmp( buf, g_MPM.captureImage );
		//sprintf(buf, "./capture_%d.png", g_MPM.captureID++);
		//saveImagePng(buf, g_MPM.captureImage);
	}
    
  glutSwapBuffers();
}
