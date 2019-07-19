#include "oclEnvironment.h"
#include "clDisplay.h"



clEnv *clEnv::s_instance = nullptr;




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