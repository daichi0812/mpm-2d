#include "usermanipulation.h"
#include "mpm.h"
#include "simulator.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include <unistd.h>
#else
#include <GL/freeglut.h>
#endif

void keyboard( unsigned char key, int x, int y )
{
	if( key == ' ' )
	{
		g_MPM.animation = !g_MPM.animation;
    glutPostRedisplay();
	}
}
