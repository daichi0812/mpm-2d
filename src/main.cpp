#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include <unistd.h>
#else
#include <GL/freeglut.h>
#endif

#define _USE_MATH_DEFINES
#include<math.h>

#include <iostream>
#include <vector>
using namespace std;

#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "mpm.h"
#include "draw.h"
#include "usermanipulation.h"
#include "xmlparser.h"
#include "window.h"
#include "HDF5io.h"

WindowInfo g_Window;

void idle()
{
	if( g_MPM.animation )
	{
		for( int i=0; i<20; i++ ) advanceSingleStep( g_MPM.simulator );
		
    glutPostRedisplay();
		fflush( stdout );
		
		if( ( g_MPM.maxSimulationTime > 0.0 ) && ( g_MPM.maxSimulationTime < g_MPM.simulator.timestep * g_MPM.simulator.deltaT ) )
		{
			printf( "Max Simulation Time %e has passed. Terminating...\n", g_MPM.maxSimulationTime );
			exit(0);
		}
	}
}

void resize( int w, int h )
{
  g_Window.width = w;
  g_Window.height = h;
}

void computeSigmaWithData( const Setting& in_Setting )
{
  g_MPM.simulator.data.deserializeTemplates( in_Setting.templates_file_name_to_resume );
  
  const std::string fn = in_Setting.serialization_folder + "/" + in_Setting.forces_file_name_for_serialization;
  HDF5File in_h5( in_Setting.forces_file_name_to_resume, HDF5AccessType::READ_ONLY );
  HDF5File in_h5_homogenization( in_Setting.homogenization_file_name, HDF5AccessType::READ_ONLY );
  HDF5File out_h5( fn, HDF5AccessType::READ_WRITE );
  
  int start, end;
  in_h5.readScalar("force_num", "start", start);
  in_h5.readScalar("force_num", "end", end);
  
  std::cout << "==== start simulation (start: " << start << ", end: " << end << ") ====" << std::endl;
  
  for( int n = start; n < end; n++ )
  {
    g_MPM.simulator.force_data.clear();
    
    g_MPM.simulator.data.deserializeElements( in_h5, std::to_string( n + 1 ) );
    g_MPM.simulator.data.deserializeInterpolatedData( in_h5_homogenization, std::to_string( n ) );
    
    const int numElements = g_MPM.simulator.data.numElements();
    
    for( int i = 0; i < numElements; i++ )
    {
      const int template_idx = g_MPM.simulator.data.getElement( i ).template_idx;
      g_MPM.simulator.data.getElement( i ).initialize( i, g_MPM.simulator.data.getTemplate( template_idx ) );
    }
    
    initMPM( in_Setting );
    
//    testStrainValue( g_MPM.simulator, in_Setting );
    
    advanceStepWithData( g_MPM.simulator );
    g_MPM.simulator.force_data.serializeForces( out_h5, std::to_string( n + 1 ) );
  }
}

int main( int argc, char * argv[] )
{
  if( argc != 2 )
  {
    printf( "Usage: <MPM> <setting file>\n" );
    exit(0);
  }

  Setting the_Setting;
  if( !openXMLFile( argv[1], the_Setting ) ) {
    printf( "Failed to read XML\n" );
    exit(1);
  }
  
  if( the_Setting.compute_sigma_only )
  {
    computeSigmaWithData( the_Setting );
    printf( "Simulation done\n" );
    return 0;
  }
  
  initMPM( the_Setting );
  
  g_Window.width = 512;
  g_Window.height = 512;
  g_MPM.captureImage.resize( g_Window.width, g_Window.height );
  
  initDisplayParam();
  
  glutInit( &argc, argv );
  glutInitWindowSize( g_Window.width, g_Window.height );
  glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
    
  glutCreateWindow( "Material Point Method" );
  
  GLint dims[4] = { 0 };
  glGetIntegerv(GL_VIEWPORT, dims);
  g_Window.framesize_windowsize_scale_x = double( dims[2] ) / double( g_Window.width );
  g_Window.framesize_windowsize_scale_y = double( dims[3] ) / double( g_Window.height );
  
  glutDisplayFunc( draw );
  glutReshapeFunc( resize );
  glutKeyboardFunc( keyboard );
  glutIdleFunc( idle );
  
  glutMainLoop();
  return 0;
}


