// clDisplay.cpp: implementation of the clDisplay class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clDisplay.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"
#include "oclEnvironment.h"
static HWND   gHwnd;


LRESULT CALLBACK WndProc1(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message)
	{

	case WM_CREATE:
		return 0;

	case WM_CLOSE:
		PostQuitMessage(0);
		return 0;

	case WM_DESTROY:
		return 0;

	case WM_LBUTTONDOWN:
	{
						   if (d.mouse_save)
						   {
							   d.down_pixel_x = GET_X_LPARAM(lParam);
							   d.down_pixel_y = GET_Y_LPARAM(lParam);
							   d.down_clicked = 1;
						   }
		return 0;
	}
	case WM_LBUTTONUP:
	{
						 if (d.mouse_save && d.down_clicked)
						 {
							 d.down_clicked = 0;
							 int spacingx = GET_X_LPARAM(lParam) - d.down_pixel_x;
							 int spacingy = d.down_pixel_y - GET_Y_LPARAM(lParam);

							 double ppx = (double)screenWidth / 2. / d.Window_Dims.x;
							 double ppy = (double)screenHeight / 2. / d.Window_Dims.y;

							 double spacing_x = (double)spacingx / ppx;
							 double spacing_y = (double)spacingy / ppy;

							 double locationx = (double)d.down_pixel_x / ppx +										(d.Window_Center.x - d.Window_Dims.x);

							 double locationy = d.Window_Center.y + d.Window_Dims.y -								(double)d.down_pixel_y / ppy;

							 vtr.saveBox(locationx, locationy, spacing_x, spacing_y);
						 }
						return 0;

	}
	

	case WM_KEYDOWN:
		switch (wParam)
		{
		case VK_UP:
		{
					  glMatrixMode(GL_PROJECTION);
					  glLoadIdentity();
					  d.Window_Center.y -= d.Window_Shift.y;
					  gluOrtho2D(d.Window_Center.x - d.Window_Dims.x, d.Window_Center.x + d.Window_Dims.x, d.Window_Center.y - d.Window_Dims.y, d.Window_Center.y + d.Window_Dims.y);
					  return 0;
		}
		case VK_DOWN:
		{
						glMatrixMode(GL_PROJECTION);
						glLoadIdentity();
						d.Window_Center.y += d.Window_Shift.y;
						gluOrtho2D(d.Window_Center.x - d.Window_Dims.x, d.Window_Center.x + d.Window_Dims.x, d.Window_Center.y - d.Window_Dims.y, d.Window_Center.y + d.Window_Dims.y);
						return 0;
		}
		case VK_RIGHT:
		{
						 glMatrixMode(GL_PROJECTION);
						 glLoadIdentity();
						 d.Window_Center.x -= d.Window_Shift.x;
						 gluOrtho2D(d.Window_Center.x - d.Window_Dims.x, d.Window_Center.x + d.Window_Dims.x, d.Window_Center.y - d.Window_Dims.y, d.Window_Center.y + d.Window_Dims.y);
						 return 0;
		}
		case VK_LEFT:
		{
						glMatrixMode(GL_PROJECTION);
						glLoadIdentity();
						d.Window_Center.x += d.Window_Shift.x;
						gluOrtho2D(d.Window_Center.x - d.Window_Dims.x, d.Window_Center.x + d.Window_Dims.x, d.Window_Center.y - d.Window_Dims.y, d.Window_Center.y + d.Window_Dims.y);
						return 0;
		}
			return 0;
		}
	case WM_CHAR:
		switch (wParam)
		{
		case 'q':
		{
					vfl.save_variables();
					return 0;

		}
		case 'z':
		{
					//vtr.Save_SS();
					return 0;
		}
		case 'y':
		{
					//vtr.NodV.save_txt_from_device("nV");
					return 0;
		}
		case 'w':
		{
					//vtr.P.save_txt_from_device_full("trp"); 
					//vtr.Ploc.save_txt_from_device("Ploc");
					return 0;
		}
		case 'o':
		{
					glMatrixMode(GL_PROJECTION);
					glLoadIdentity();
					d.Window_Dims.x += d.Window_Shift.x;
					d.Window_Dims.y += d.Window_Shift.y;
					d.Window_Shift.x = d.Window_Dims.x * 0.10;
					d.Window_Shift.y = d.Window_Dims.x * 0.10;
					gluOrtho2D(d.Window_Center.x - d.Window_Dims.x, d.Window_Center.x + d.Window_Dims.x, d.Window_Center.y - d.Window_Dims.y, d.Window_Center.y + d.Window_Dims.y);
					return 0;

		}
		case 'i':
		{
					glMatrixMode(GL_PROJECTION);
					glLoadIdentity();
					d.Window_Dims.x -= d.Window_Shift.x;
					d.Window_Dims.y -= d.Window_Shift.y;
					d.Window_Shift.x = d.Window_Dims.x * 0.10;
					d.Window_Shift.y = d.Window_Dims.x * 0.10;
					gluOrtho2D(d.Window_Center.x - d.Window_Dims.x, d.Window_Center.x + d.Window_Dims.x, d.Window_Center.y - d.Window_Dims.y, d.Window_Center.y + d.Window_Dims.y);
					return 0;
		}
		case 'n':
		{
					d.lsize -= 0.2;
					return 0;

		}
		case 't':
		{
					d.lsize += 0.2;
					return 0;

		}
		case 's':
		{
					d.psize -= 0.2;
					return 0;

		}
		case 'm':
		{
					d.mouse_save ^= 1;
					return 0;

		}
		case 'b':
		{
					d.psize += 0.2;
					return 0;

		}
			return 0;
		}
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);

	}
}

void clDisplay::finish()
{
	ShowTime();
	printf("\n");
}

void clDisplay::ini()
{
	counter = 0;
	time(&StartTimeSys);
	strcpy(cStartTimeSys,&(ctime(&StartTimeSys)[11]));
};

void clDisplay::ShowTime()
{
	time_t CurrTimeSys, ElapsTimeSys, FinishTimeSys, RemTimeSys;
	char cFinishTimeSys[26];

	time(&CurrTimeSys);

	ElapsTimeSys = CurrTimeSys - StartTimeSys;
	FinishTimeSys = StartTimeSys + (time_t)((double)ElapsTimeSys * (p.StopTime - p.RestartTime) / (p.TimeN - p.RestartTime));
	int ElapsH = 0, ElapsM = 0, ElapsS = 0;
	if (ElapsTimeSys > 0)
	{
		ElapsH = (int)(ElapsTimeSys / (60 * 60));
		ElapsM = (int)((ElapsTimeSys - ElapsH * 60 * 60) / 60);
		ElapsS = (int)((ElapsTimeSys - ElapsH * 60 * 60 - ElapsM * 60));
	}
	int RemH = -1, RemM = 0, RemS = 0;
	if (FinishTimeSys > 0)
	{
		RemTimeSys = FinishTimeSys - CurrTimeSys;
		RemH = (int)(RemTimeSys / (60 * 60));
		RemM = (int)((RemTimeSys - RemH * 60 * 60) / 60);
		RemS = (int)((RemTimeSys - RemH * 60 * 60 - RemM * 60));
	}

	if (FinishTimeSys > 0)
	{
		strcpy(cFinishTimeSys, &(ctime(&FinishTimeSys)[11]));
	}
	else
	{
		sprintf(cFinishTimeSys, "00:00:00");
	}
	
	//vlb.Calculate_U_mean();
	//vlb.Mass_Cor.read_from_buffer();
	
	//printf("Mass%5.1f T%f U%f Ro%g E%d:%.2d:%.2d R%d:%.2d:%.2d P%5.2f%%\r",
	//	vls.Masses(0), vfd.Reduce_T.reduceSingle() / (double)(p.nX*p.Channel_Height), vtr.Umean_Current, vlb.Mass_Cor(0)*81., ElapsH, ElapsM, ElapsS, RemH, RemM, RemS, 100. * c.Time / c.StopTime);

};
void clDisplay::step()
{
	counter++;
	if (++counter < p.displaySignalFreq)
		return;

	ShowTime();
	counter = 0;
}

void clDisplay::init_glew()
{
	glewInit();
	if (!glewIsSupported("GL_VERSION_2_0 " "GL_ARB_pixel_buffer_object"))
	{
		std::cerr << "Support for necessary OpenGL extensions missing."
			<< std::endl;
		return exit(0);
	}

	glClearColor(1., 1., 1., 1.);

	glViewport(0., 0., screenWidth, screenHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	double ar = (double)p.nY / (double)p.nX;
	d.Window_Center = { { 100., (double)p.Channel_Height/2. } };
	d.Window_Dims = { { (double)p.Channel_Height, (double)p.Channel_Height } };
	d.Window_Shift = { { 2., 2. } };

	gluOrtho2D(d.Window_Center.x - d.Window_Dims.x, d.Window_Center.x + d.Window_Dims.x, d.Window_Center.y - d.Window_Dims.y, d.Window_Center.y + d.Window_Dims.y);
}

int clDisplay::displayDevices(cl_platform_id platform, cl_device_type deviceType)
{
	cl_int status;
	// Get platform name
	char platformVendor[1024];
	status = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(platformVendor),
		platformVendor, NULL);
	CHECK_OPENCL_ERROR(status, "clGetPlatformInfo failed");
	std::cout << "\nSelected Platform Vendor : " << platformVendor << std::endl;
	// Get number of devices available
	cl_uint deviceCount = 0;
	status = clGetDeviceIDs(platform, deviceType, 0, NULL, &deviceCount);
	CHECK_OPENCL_ERROR(status, "clGetDeviceIDs failed");
	cl_device_id* deviceIds = (cl_device_id*)malloc(sizeof(cl_device_id)*
		deviceCount);
	CHECK_ALLOCATION(deviceIds, "Failed to allocate memory(deviceIds)");
	// Get device ids
	status = clGetDeviceIDs(platform, deviceType, deviceCount, deviceIds, NULL);
	CHECK_OPENCL_ERROR(status, "clGetDeviceIDs failed");
	// Print device index and device names
	int afsg = CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
	for (cl_uint i = 0; i < deviceCount; ++i)
	{
		char deviceName[1024];
		status = clGetDeviceInfo(deviceIds[i], CL_DEVICE_NAME, sizeof(deviceName),
			deviceName, NULL);

		cl_int dims;
		status = clGetDeviceInfo(deviceIds[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_int), &dims, NULL);
		size_t workitem_size[3], workgroup_size;
		cl_uint dsize;
		status = clGetDeviceInfo(deviceIds[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(workitem_size), &workitem_size, NULL);
		CHECK_OPENCL_ERROR(status, "clGetDeviceInfo failed");
		status = clGetDeviceInfo(deviceIds[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size, NULL);
		status = clGetDeviceInfo(deviceIds[i], CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(dsize), &dsize, NULL);

		size_t extensionSize;
		clGetDeviceInfo(deviceIds[i], CL_DEVICE_EXTENSIONS, 0, NULL, &extensionSize);
		char *extensions;
		extensions = (char*)malloc(extensionSize);
		status = clGetDeviceInfo(deviceIds[i], CL_DEVICE_EXTENSIONS, extensionSize, extensions, NULL);

		std::cout << "Device " << i << " : " << deviceName
			<< " Device ID is " << deviceIds[i] << std::endl << "Work item size "
			<< workitem_size[0] << ", " << workitem_size[1] << ", " << workitem_size[2] << "\n" << "Work group size "
			<< workgroup_size << "\n" << "Double Vector Size " << dsize << "\n";
		printf("Device Extensions:\n");

		for (int ij = 0; ij < (int)strlen(extensions); ij++)
		{
			if (extensions[ij] == ' ')
				extensions[ij] = '\n';
		}
		printf("%s\n", extensions);
		FREE(extensions);

	}
	free(deviceIds);
	return SDK_SUCCESS;
}


void clDisplay::Render_Domain()
{
//	d.frameCount++;
//
//	vtr.P_vbo.release_gl_objects(p.IOqueue);
//	vls.LSt_vbo.release_gl_objects(p.IOqueue);
//	vls.LSb_vbo.release_gl_objects(p.IOqueue);
//	clFinish(p.IOqueue);
//
//	d.check_message();
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//
//	glEnable(GL_BLEND);
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//
//	//particles
//	vtr.P_vbo.Render_VBO(d.psize);
//
//	//Render Pos of Wall
//	vls.LSt0_vbo.Render_VBO(d.lsize);
//	vls.LSb0_vbo.Render_VBO(d.lsize);
//
//	//Render Current Pos of Interface with Points indicating Nodes
//	vls.LSt_vbo.Render_VBO(d.lsize);
//	vls.LSb_vbo.Render_VBO(d.lsize);
//	vls.LSt_vbo.Render_Points_Const(d.lsize);
//	vls.LSb_vbo.Render_Points_Const(d.lsize);
//
//#ifdef OPENGL_GRIDLINES
//	vls.LinesH.Render_Individual_Lines_Const(d.lsize / 2.);
//	vls.LinesV.Render_Individual_Lines_Const(d.lsize / 2.);
//#endif 
//
//	vtr.P_vbo.acquire_gl_objects(p.IOqueue);
//	vls.LSt_vbo.acquire_gl_objects(p.IOqueue);
//	vls.LSb_vbo.acquire_gl_objects(p.IOqueue);
//	clFlush(p.IOqueue);
//	SwapBuffers(d.gHdc);
//
//	d.t2 = clock() * CLOCKS_PER_SEC;
//	d.totalElapsedTime += (double)(d.t2 - d.t1);
//	d.t1 = clock() * CLOCKS_PER_SEC;
//	if (d.frameCount && d.frameCount > d.frameRefCount)
//	{
//		double fMs = (double)((d.totalElapsedTime / (double)CLOCKS_PER_SEC) /
//			(double)(d.frameCount * c.dTtr) );
//		int framesPerSec = (int)(1.0 / (fMs / CLOCKS_PER_SEC));
//
//		d.set_window(framesPerSec);
//		d.frameCount = 0;
//		d.totalElapsedTime = 0.0;
//	}
//
//	glFinish();
//	clFinish(p.IOqueue);
}

void clDisplay::check_message()
{
	if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
	{
		// handle or dispatch messages
		if (msg.message == WM_QUIT)
		{
			quit = true;
		}
		else
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}
}

void clDisplay::set_window(int framesPerSec)
{
	char title[256];
	sprintf_s(title, 256, "t = %d | %d fps ", (int)p.Time, framesPerSec);
	SetWindowText(gHwnd, title);
}

int clDisplay::enableGLAndGetGLContext(cl_platform_id platform)
{
	//cl_int status;
	//bool ret = false;
	//int  pfmt;
	//PIXELFORMATDESCRIPTOR pfd;
	//ZeroMemory(&pfd, sizeof(pfd));
	//pfd.nSize = sizeof(pfd);
	//pfd.nVersion = 1;
	//pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL |
	//	PFD_DOUBLEBUFFER;
	//pfd.iPixelType = PFD_TYPE_RGBA;
	//pfd.cColorBits = 24;
	//pfd.cDepthBits = 16;
	//pfd.iLayerType = PFD_MAIN_PLANE;


	//WNDCLASS windowclass;

	//windowclass.style = CS_OWNDC;
	//windowclass.lpfnWndProc = WndProc1;
	//windowclass.cbClsExtra = 0;
	//windowclass.cbWndExtra = 0;
	//windowclass.hInstance = NULL;
	//windowclass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	//windowclass.hCursor = LoadCursor(NULL, IDC_ARROW);
	//windowclass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	//windowclass.lpszMenuName = NULL;
	//windowclass.lpszClassName = reinterpret_cast<LPCSTR>("SimpleGL");
	//RegisterClass(&windowclass);

	//gHwnd = CreateWindow(reinterpret_cast<LPCSTR>("SimpleGL"),
	//	reinterpret_cast<LPCSTR>("OpenGL Texture Renderer"),
	//	WS_CAPTION | WS_POPUPWINDOW,
	//	0,
	//	0,
	//	screenWidth,
	//	screenHeight,
	//	NULL,
	//	NULL,
	//	windowclass.hInstance,
	//	NULL);
	//d.gHdc = GetDC(gHwnd);

	//pfmt = ChoosePixelFormat(d.gHdc,
	//	&pfd);

	//int iFormat = ChoosePixelFormat(d.gHdc, &pfd);
	//SetPixelFormat(d.gHdc, iFormat, &pfd);
	//d.gGlCtx = wglCreateContext(d.gHdc);
	//wglMakeCurrent(d.gHdc, d.gGlCtx);

	//cl_context_properties properties[] =
	//{
	//	CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(),
	//	CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(),
	//	CL_CONTEXT_PLATFORM, (cl_context_properties)platform,
	//	0
	//};

	//// Create OpenCL context from device's id
	//p.context = clCreateContext(properties,
	//	1,
	//	&p.device_use,
	//	0,
	//	0,
	//	&status);
	//CHECK_OPENCL_ERROR(status, "clCreateContext failed!!");
	//ShowWindow(gHwnd, SW_SHOW);		// everything went OK, show the window
	//UpdateWindow(gHwnd);
	return SDK_SUCCESS;
}