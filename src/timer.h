#ifndef __TIMER_H__
#define __TIMER_H__

#ifdef __APPLE__
#include <sys/time.h>
#else
#include <time.h>
#endif


struct Timer
{
#ifdef __APPLE__
	timeval s;
	timeval t;
#else
  timespec s;
  timespec t;
#endif
};

inline void startTimer( Timer& io_Timer )
{
#ifdef __APPLE__
	gettimeofday( &io_Timer.s, NULL );
#else
  timespec_get( &io_Timer.s, TIME_UTC );
#endif
}

inline void stopTimer( Timer& io_Timer )
{
#ifdef __APPLE__
	gettimeofday( &io_Timer.t, NULL );
#else
  timespec_get( &io_Timer.t, TIME_UTC );
#endif
}

inline double getElapsedMs( Timer& io_Timer )
{
#ifdef __APPLE__
	double elapsedTime = ( io_Timer.t.tv_sec - io_Timer.s.tv_sec ) * 1000.0; //sec to ms
	elapsedTime += ( io_Timer.t.tv_usec - io_Timer.s.tv_usec ) / 1000.0; //us to ms
#else
  double elapsedTime = ( io_Timer.t.tv_sec - io_Timer.s.tv_sec ) * 1000.0; //sec to ms
  elapsedTime += ( io_Timer.t.tv_nsec - io_Timer.s.tv_nsec ) / 1000000.0; //ns to ms
#endif
  return elapsedTime;
}

#endif
