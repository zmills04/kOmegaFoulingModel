// clDisplay.h: interface for the clDisplay class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_)
#define AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_


//TODO: Clean this up and implement openGL calls in this cpp/h file


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clProblem.h"
class clDisplay  
{
public:
	time_t StartTimeSys;
	char cStartTimeSys[26];
	unsigned int counter;
	//OpenGL objects
	MSG msg;
	HGLRC         gGlCtx;
	HDC           gHdc;
	BOOL quit = FALSE;

	int down_pixel_x = 0;
	int down_clicked = 0;
	int down_pixel_y = 0;
	int mouse_save = 0;
	double t2, t1, totalElapsedTime = 0.;
	int frameCount = 0, frameRefCount = 90;

	//size of particles and lines
	double psize = POINT_SIZES, lsize = LINE_SIZES;
	//variables used to define window parameters
	cl_double2 Window_Center, Window_Dims, Window_Shift;
	void set_window(int framesPerSec);
	void check_message();
	int enableGLAndGetGLContext(cl_platform_id platform);

	clDisplay(){ counter = 0; };
	virtual ~clDisplay(){};

	static int displayDevices(cl_platform_id platform, cl_device_type deviceType);

	void finish();
	
	void ini(void);      
	
	void init_glew();

	void Render_Domain();

	void ShowTime(void);

	void start(){ShowTime();}

	void step();
};

#endif // !defined(AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_)
