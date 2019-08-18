#include "oclEnvironment.h"
#include "clDisplay.h"



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
		if (clEnv::instance()->mouse_save)
		{
			clEnv::instance()->down_pixel_x = GET_X_LPARAM(lParam);
			clEnv::instance()->down_pixel_y = GET_Y_LPARAM(lParam);
			clEnv::instance()->down_clicked = 1;
		}
		return 0;
	}
	case WM_LBUTTONUP:
	{
		if (clEnv::instance()->mouse_save && clEnv::instance()->down_clicked)
		{
			clEnv::instance()->down_clicked = 0;
			int spacingx = GET_X_LPARAM(lParam) - clEnv::instance()->down_pixel_x;
			int spacingy = clEnv::instance()->down_pixel_y - GET_Y_LPARAM(lParam);

			double ppx = (double)screenWidth / 2. / clEnv::instance()->Window_Dims.x;
			double ppy = (double)screenHeight / 2. / clEnv::instance()->Window_Dims.y;

			double spacing_x = (double)spacingx / ppx;
			double spacing_y = (double)spacingy / ppy;

			double locationx = (double)(clEnv::instance()->down_pixel_x) / ppx +
				(clEnv::instance()->Window_Center.x - clEnv::instance()->Window_Dims.x);

			double locationy = clEnv::instance()->Window_Center.y + 
				clEnv::instance()->Window_Dims.y - 
				(double)(clEnv::instance()->down_pixel_y) / ppy;

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
			clEnv::instance()->Window_Center.y -= clEnv::instance()->Window_Shift.y;
			gluOrtho2D(clEnv::instance()->Window_Center.x - clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.x + clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.y - clEnv::instance()->Window_Dims.y,
				clEnv::instance()->Window_Center.y + clEnv::instance()->Window_Dims.y);
			return 0;
		}
		case VK_DOWN:
		{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			clEnv::instance()->Window_Center.y += clEnv::instance()->Window_Shift.y;
			gluOrtho2D(clEnv::instance()->Window_Center.x - clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.x + clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.y - clEnv::instance()->Window_Dims.y,
				clEnv::instance()->Window_Center.y + clEnv::instance()->Window_Dims.y);
			return 0;
		}
		case VK_RIGHT:
		{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			clEnv::instance()->Window_Center.x -= clEnv::instance()->Window_Shift.x;
			gluOrtho2D(clEnv::instance()->Window_Center.x - clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.x + clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.y - clEnv::instance()->Window_Dims.y,
				clEnv::instance()->Window_Center.y + clEnv::instance()->Window_Dims.y);
			return 0;
		}
		case VK_LEFT:
		{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			clEnv::instance()->Window_Center.x += clEnv::instance()->Window_Shift.x;
			gluOrtho2D(clEnv::instance()->Window_Center.x - clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.x + clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.y - clEnv::instance()->Window_Dims.y,
				clEnv::instance()->Window_Center.y + clEnv::instance()->Window_Dims.y);
			return 0;
		}
		return 0;
		}
	case WM_CHAR:
		switch (wParam)
		{
		case 'q':
		{
			vfl.saveVariables();
			return 0;

		}
		case 'z':
		{
			vtr.wallShear.saveDebug();
			return 0;
		}
		//case 'y':
		//{
		//	vtr.nodV.saveFromDevice(true, trStructBase::saveTxtFl);
		//	return 0;
		//}
		case 'w':
		{
			vtr.P.saveFromDevice(true, trStructBase::saveTxtFl);
			vtr.parSort.Ploc.save_txt_from_device();
			return 0;
		}
		case 'o':
		{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			clEnv::instance()->Window_Dims.x += clEnv::instance()->Window_Shift.x;
			clEnv::instance()->Window_Dims.y += clEnv::instance()->Window_Shift.y;
			clEnv::instance()->Window_Shift.x = clEnv::instance()->Window_Dims.x * 0.10;
			clEnv::instance()->Window_Shift.y = clEnv::instance()->Window_Dims.x * 0.10;
			gluOrtho2D(clEnv::instance()->Window_Center.x - clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.x + clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.y - clEnv::instance()->Window_Dims.y,
				clEnv::instance()->Window_Center.y + clEnv::instance()->Window_Dims.y);
			return 0;

		}
		case 'i':
		{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			clEnv::instance()->Window_Dims.x -= clEnv::instance()->Window_Shift.x;
			clEnv::instance()->Window_Dims.y -= clEnv::instance()->Window_Shift.y;
			clEnv::instance()->Window_Shift.x = clEnv::instance()->Window_Dims.x * 0.10;
			clEnv::instance()->Window_Shift.y = clEnv::instance()->Window_Dims.x * 0.10;
			gluOrtho2D(clEnv::instance()->Window_Center.x - clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.x + clEnv::instance()->Window_Dims.x,
				clEnv::instance()->Window_Center.y - clEnv::instance()->Window_Dims.y,
				clEnv::instance()->Window_Center.y + clEnv::instance()->Window_Dims.y);
			return 0;
		}
		case 'n':
		{
			clEnv::instance()->lsize -= 0.2;
			return 0;

		}
		case 't':
		{
			clEnv::instance()->lsize += 0.2;
			return 0;

		}
		case 's':
		{
			clEnv::instance()->psize -= 0.2;
			return 0;

		}
		case 'm':
		{
			clEnv::instance()->mouse_save ^= 1;
			return 0;

		}
		case 'b':
		{
			clEnv::instance()->psize += 0.2;
			return 0;

		}
		return 0;
		}
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);

	}
}



clEnv *clEnv::s_instance = nullptr;
HWND clEnv::gHwnd = nullptr;


void clEnv::iniCLGLContext()
{
	cl_int status;
	BOOL ret = FALSE;
	int  pfmt;
	PIXELFORMATDESCRIPTOR pfd;
	ZeroMemory(&pfd, sizeof(pfd));
	pfd.nSize = sizeof(pfd);
	pfd.nVersion = 1;
	pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL |
		PFD_DOUBLEBUFFER;
	pfd.iPixelType = PFD_TYPE_RGBA;
	pfd.cColorBits = 24;
	pfd.cDepthBits = 16;
	pfd.iLayerType = PFD_MAIN_PLANE;


	WNDCLASS windowclass;

	windowclass.style = CS_OWNDC;
	windowclass.lpfnWndProc = WndProc1;
	windowclass.cbClsExtra = 0;
	windowclass.cbWndExtra = 0;
	windowclass.hInstance = NULL;
	windowclass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	windowclass.hCursor = LoadCursor(NULL, IDC_ARROW);
	windowclass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	windowclass.lpszMenuName = NULL;
	windowclass.lpszClassName = reinterpret_cast<LPCSTR>("SimpleGL");
	RegisterClass(&windowclass);

	gHwnd = CreateWindow(reinterpret_cast<LPCSTR>("SimpleGL"),
		reinterpret_cast<LPCSTR>("OpenGL Texture Renderer"),
		WS_CAPTION | WS_POPUPWINDOW,
		0,
		0,
		screenWidth,
		screenHeight,
		NULL,
		NULL,
		windowclass.hInstance,
		NULL);
	gHdc = GetDC(gHwnd);

	pfmt = ChoosePixelFormat(gHdc,
		&pfd);

	int iFormat = ChoosePixelFormat(gHdc, &pfd);
	SetPixelFormat(gHdc, iFormat, &pfd);
	gGlCtx = wglCreateContext(gHdc);
	wglMakeCurrent(gHdc, gGlCtx);

	cl_context_properties properties[] =
	{
		CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(),
		CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(),
		CL_CONTEXT_PLATFORM, (cl_context_properties)platform,
		0
	};

	// Create OpenCL context from device's id
	context = clCreateContext(properties, 1, &device, nullptr, nullptr, &status);
	ERROR_CHECKING(status, "create context failed", ERROR_OCL_INITIALIZATION);
	ShowWindow(gHwnd, SW_SHOW);		// everything went OK, show the window
	UpdateWindow(gHwnd);
}



void clEnv::initGlew()
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
	Window_Center = { { 100., (double)p.Channel_Height / 2. } };
	Window_Dims = { { (double)p.Channel_Height, (double)p.Channel_Height } };
	Window_Shift = { { 2., 2. } };

	gluOrtho2D(	Window_Center.x - Window_Dims.x, 
				Window_Center.x + Window_Dims.x,
				Window_Center.y - Window_Dims.y,
				Window_Center.y + Window_Dims.y);
}

void clEnv::checkMessage()
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

void clEnv::setWindow(int framesPerSec)
{
	char title[256];
	sprintf_s(title, 256, "t = %d | %d fps ", (int)p.Time, framesPerSec);
	SetWindowText(gHwnd, title);
}


void clEnv::renderDomain()
{
	frameCount++;
	
	vtr.glParticles.P_vbo.release_gl_objects(IOqueue);
	vtr.glParticles.LSt_vbo.release_gl_objects(IOqueue);
	vtr.glParticles.LSb_vbo.release_gl_objects(IOqueue);
	clFinish(IOqueue);
	
	checkMessage();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	//particles
	vtr.glParticles.P_vbo.Render_VBO(psize);
	
	//Render Pos of Wall
	vtr.glParticles.LSt0_vbo.Render_VBO(lsize);
	vtr.glParticles.LSb0_vbo.Render_VBO(lsize);
	
	//Render Current Pos of Interface with Points indicating Nodes
	vtr.glParticles.LSt_vbo.Render_VBO(lsize);
	vtr.glParticles.LSb_vbo.Render_VBO(lsize);
	vtr.glParticles.LSt_vbo.Render_Points_Const(lsize);
	vtr.glParticles.LSb_vbo.Render_Points_Const(lsize);
	if(vtr.glParticles.glGridFlag)
	{
		vtr.glParticles.LinesH.Render_Individual_Lines_Const(lsize / 2.);
		vtr.glParticles.LinesV.Render_Individual_Lines_Const(lsize / 2.);
	} 
	
	vtr.glParticles.P_vbo.acquire_gl_objects(IOqueue);
	vtr.glParticles.LSt_vbo.acquire_gl_objects(IOqueue);
	vtr.glParticles.LSb_vbo.acquire_gl_objects(IOqueue);
	clFlush(IOqueue);
	SwapBuffers(gHdc);
	
	t2 = clock() * CLOCKS_PER_SEC;
	totalElapsedTime += (double)(t2 - t1);
	t1 = clock() * CLOCKS_PER_SEC;
	if (frameCount && frameCount > frameRefCount)
	{
		double fMs = (double)((totalElapsedTime / (double)CLOCKS_PER_SEC) /
			(double)(frameCount * p.dTtr) );
		int framesPerSec = (int)(1.0 / (fMs / CLOCKS_PER_SEC));
	
		setWindow(framesPerSec);
		frameCount = 0;
		totalElapsedTime = 0.0;
	}
	
	glFinish();
	clFinish(IOqueue);
}

void clEnv::displayDevices()
{
	cl_int status;
	// Get platform name
	char platformVendor[1024];
	status = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(platformVendor),
		platformVendor, NULL);

	ERROR_CHECKING(status, "clGetPlatformInfo failed", ERROR_OCL_INITIALIZATION);
	std::cout << "\nSelected Platform Vendor : " << platformVendor << std::endl;
	// Get number of devices available
	cl_uint deviceCount = 0;
	status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &deviceCount);
	ERROR_CHECKING(status, "clGetDeviceIDs failed", ERROR_OCL_INITIALIZATION);
	cl_device_id* deviceIds = (cl_device_id*)malloc(sizeof(cl_device_id) *
		deviceCount);
	CHECK_ALLOCATION(deviceIds, "Failed to allocate memory(deviceIds)", ERROR_OCL_INITIALIZATION);
	// Get device ids
	status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, deviceCount, deviceIds, NULL);
	ERROR_CHECKING(status, "clGetDeviceIDs failed", ERROR_OCL_INITIALIZATION);
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
		ERROR_CHECKING(status, "clGetDeviceInfo failed", ERROR_OCL_INITIALIZATION);
		status = clGetDeviceInfo(deviceIds[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size, NULL);
		status = clGetDeviceInfo(deviceIds[i], CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(dsize), &dsize, NULL);

		size_t extensionSize;
		clGetDeviceInfo(deviceIds[i], CL_DEVICE_EXTENSIONS, 0, NULL, &extensionSize);
		char* extensions;
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
}

void clEnv::buildProgram(cl_program &program, std::string &KerString, 
	programType ptype_, std::string name_)
{
	int status;
	size_t lengths[1] = { KerString.size() };
	const char* sources[1] = { KerString.data() };

	program = clCreateProgramWithSource(context, 1, sources, lengths, &status);
	ERROR_CHECKING(status, "Error creating program from source", ERROR_OCL_INITIALIZATION);

	char oclVersion[1024];
	size_t oclVersionSize;
	status = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(oclVersion), oclVersion, &oclVersionSize);
	ERROR_CHECKING(status, "Error getting platform information", ERROR_OCL_INITIALIZATION);

	std::string platver(oclVersion, oclVersionSize);
	if (platver.find("OpenCL 2.") != std::string::npos)
	{
		status = clBuildProgram(program, 1, &device, "-cl-mad-enable -cl-std=CL1.2", nullptr, nullptr);
	}
	else
	{
		status = clBuildProgram(program, 1, &device, "-cl-std=CL2.0 -cl-mad-enable", nullptr, nullptr);
	}

	printBuildInfo(program, status, name_);
	ERROR_CHECKING(status, "Build program failed", ERROR_OCL_INITIALIZATION);
	programFl[ptype_] = true;
}

// Initializes context and sets it as the default
void clEnv::iniContext()
{
	const cl_context_properties properties[] =
	{
		CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties> (platform), 0, 0
	};
	int status;
	context = clCreateContext(properties, 1, &device, nullptr, nullptr, &status);
	ERROR_CHECKING(status, "create context failed", ERROR_OCL_INITIALIZATION);
	contextFlag = true;
//	staticBaseVar::setArrayContext(&context);
}



// Initializes device and sets it as default. If DEVICE_ID defined in constdef.h is greater than number devices,
// this function will default to device 0
void clEnv::iniDevice(cl_uint deviceID)
{
	cl_uint deviceIdCount = 0;
	clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &deviceIdCount);

	ERROR_CHECKING(deviceIdCount == 0, "No OpenCL devices exist on system", ERROR_OCL_INITIALIZATION);
	std::vector<cl_device_id> deviceIds(deviceIdCount);
	clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, deviceIdCount, deviceIds.data(), nullptr);

#ifdef PRINT_OPENCL_INFO
	clDisplay::displayDevices(platform, CL_DEVICE_TYPE_GPU);
#endif


	if (deviceID >= deviceIdCount)
	{
		deviceID = 0;
		std::cout << "Warning: Defined OpenCL device does not exist, defaulting to device 0" << std::endl;
	}

	device = deviceIds[deviceID];
	deviceFlag = true;
}

// Initializes platform and sets it as the default
void clEnv::iniPlatform()
{
	cl_uint platformIdCount = 0;
	clGetPlatformIDs(0, nullptr, &platformIdCount);
	std::vector<cl_platform_id> platformIds(platformIdCount);
	clGetPlatformIDs(platformIdCount, platformIds.data(), nullptr);
	platform = platformIds[0];
}

void clEnv::iniQueues()
{
	queProps[0] = CL_QUEUE_PROPERTIES;
	queProps[1] = 0;
	queProps[2] = 0;
	queProps[3] = 0;

	cl_int status;

	TRqueue = clCreateCommandQueueWithProperties(context, device, queProps, &status);
	ERROR_CHECKING(status, "Error Creating TRqueue", ERROR_OCL_INITIALIZATION);

	LBqueue = clCreateCommandQueueWithProperties(context, device, queProps, &status);
	ERROR_CHECKING(status, "Error Creating LBqueue", ERROR_OCL_INITIALIZATION);

	IOqueue = clCreateCommandQueueWithProperties(context, device, queProps, &status);
	ERROR_CHECKING(status, "Error Creating IOqueue", ERROR_OCL_INITIALIZATION);

	FDqueue = clCreateCommandQueueWithProperties(context, device, queProps, &status);
	ERROR_CHECKING(status, "Error Creating FDqueue", ERROR_OCL_INITIALIZATION);

	queueFlag = true;
}

void clEnv::printBuildInfo(cl_program &program, int status, std::string name_)
{
	size_t log_size;
	int cl_err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
	char *log = (char *)malloc(log_size);
	cl_err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
	
	std::ofstream ofs("x64" SLASH "clSource" SLASH + name_ + "_build_output.txt", std::ofstream::out);
	if (log_size > 1)
	{
		ofs << log;
		if (status == 0)
			printf("Program compiled with warnings see build log for info\n");
		else
			printf("Error compiling program see build log for info\n");
	}
	else
	{
		ofs << "Source Built Successfully";
		std::cout << "Successful Build of " << name_ << " Source" << std::endl;
	}
	ofs.close();
	FREE(log);
}

void clEnv::Release_Event(cl_event evt, std::string evt_name)
{
	cl_int ret_val;
	for (int i = 0; i < 10000; i++)
	{
		cl_int ret = clGetEventInfo(evt, 
			CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &ret_val, NULL);

		if (ret_val == CL_COMPLETE)
		{
			clReleaseEvent(evt);
			return;
		}

		delay_func(0.01);
	}

	ERROR_CHECKING(ret_val, "Event " + evt_name + " was unable to be released", ret_val);
}