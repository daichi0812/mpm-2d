#ifndef __MATERIAL_POINT_METHOD_H__
#define __MATERIAL_POINT_METHOD_H__

#include "draw.h"
#include "simulator.h"
#include "xmlparser.h"
#include "image.h"

struct MPM
{
	Simulator simulator;
	bool animation;
	
	Image<float, 3> captureImage;
	bool capture;
	double lastCaptureTime;
	int captureID;
	double captureTimeStep;
	double maxSimulationTime;
};

void initDisplayParam();
void initMPM( const Setting& in_Setting );

extern DisplayParam g_DisplayParam;
extern MPM g_MPM;

#endif
