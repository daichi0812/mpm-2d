#include "mpm.h"

MPM g_MPM;

void initDisplayParam()
{
	g_DisplayParam.left = 0.0;
	g_DisplayParam.right = 1.0;
	g_DisplayParam.bottom = 0.0;
	g_DisplayParam.top = 1.0;
}

void initMPM( const Setting& in_Setting )
{
  if( in_Setting.compute_sigma_only )
    initSimulatorWithData( g_MPM.simulator, in_Setting );
  else
    initSimulator( g_MPM.simulator, in_Setting );
  
  g_MPM.animation = false;
	g_MPM.capture = in_Setting.capture;
	g_MPM.lastCaptureTime = -100000.0;
	g_MPM.captureID = 0;
	g_MPM.captureTimeStep = in_Setting.captureInterval;
	g_MPM.maxSimulationTime = in_Setting.maxSimulationTime;
}
