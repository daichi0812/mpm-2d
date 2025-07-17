//
//  window.h
//  SimpleMPM2D
//
//  Created by Yonghao Yue on 2020/08/10.
//

#ifndef window_h
#define window_h

struct WindowInfo
{
  int width;
  int height;
  double framesize_windowsize_scale_x;
  double framesize_windowsize_scale_y;
};

extern WindowInfo g_Window;

#endif /* window_h */
