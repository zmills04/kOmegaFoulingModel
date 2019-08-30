// oclEnvironment.h: Handles Initialization of OpenCL, and contains all
// opencl related variables (i.e. queues, context, waitlists) 
// Program Sources are generated elsewhere, but are compiled here
//
// (c) Zach Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_OCLENVIRONMENT_H__INCLUDED_)
#define AFX_OCLENVIRONMENT_H__INCLUDED_

#pragma once


#include "StdAfx.h"
#include "logger.h"


#define IOQUEUE			*clEnv::instance()->getIOqueue()
#define IOQUEUE_REF		clEnv::instance()->getIOqueue()
#define LBQUEUE			*clEnv::instance()->getLBqueue()
#define LBQUEUE_REF		clEnv::instance()->getLBqueue()
#define FDQUEUE			*clEnv::instance()->getFDqueue()
#define FDQUEUE_REF		clEnv::instance()->getFDqueue()
#define TRQUEUE			*clEnv::instance()->getTRqueue()
#define TRQUEUE_REF		clEnv::instance()->getTRqueue()
#define CLCONTEXT		*clEnv::instance()->getContext()
#define CLCONTEXT_REF	clEnv::instance()->getContext()


class clEnv
{
public:
	enum programType { ReduceType = 0, BaseType = 1, BiCGSType = 2 };
	
	static clEnv *instance()
	{
		if (!s_instance)
			s_instance = new clEnv;
		return s_instance;
	}

	void setTurbFlag(bool tflag_) {	turbFlag = tflag_; }
	void setTempFlag(bool tflag_) {	tempFlag = tflag_; }
	void setTrFlag(bool tflag_) { trFlag = tflag_; }
	void setRestartFlag(bool rflag_) { restartFlag = rflag_; }
	bool getRestartFlag() { return restartFlag; }


	// Initializes context and sets it as the default
	void iniContext();

	// Initializes platform and sets it as the default
	void iniPlatform();

	// Initializes device and sets it as default. If DEVICE_ID defined in constdef.h is greater than number devices,
	// this function will default to device 0
	void iniDevice(cl_uint deviceID);

	void iniQueues();

	void buildProgram(cl_program &program, std::string &KerString, programType ptype_, std::string name_);

	void printBuildInfo(cl_program &program, int status, std::string name_);

	// Waits till events have completed before releasing cl_event object and returning
	void Release_Event(cl_event evt, std::string evt_name);

	void finishQueues()
	{
		flushQueues();
		clFinish(IOqueue);
		clFinish(LBqueue);
		clFinish(FDqueue);
		clFinish(TRqueue);
	}

	void flushQueues()
	{
		clFlush(IOqueue);
		clFlush(LBqueue);
		clFlush(FDqueue);
		clFlush(TRqueue);
	}

	void iniOpenCL(int deviceid_, bool useOpenGL)
	{
		iniPlatform();
		LOGMESSAGE("openCL platform initialized");
		iniDevice((cl_uint)deviceid_);
		LOGMESSAGE("Device number " + std::to_string(deviceid_) + " initialized");
		if (!useOpenGL)
		{
			iniContext();
			LOGMESSAGE("openCL context initialized");
		}
		else // GL implementation not working at the moment
		{
			iniCLGLContext();
			LOGMESSAGE("openCL-GL context initialized");
		}
		iniQueues();
		LOGMESSAGE("openCL queues initialized");

	}

	void iniCLGLContext();

	void displayDevices();



//////////////OpenGL Variables and Functions///////////////////////
	MSG msg;
	HGLRC         gGlCtx;
	static HWND gHwnd;
	HDC           gHdc;
	bool quit = false;

	int down_pixel_x = 0;
	int down_clicked = 0;
	int down_pixel_y = 0;
	int mouse_save = 0;
	double t2, t1, totalElapsedTime = 0.;
	int frameCount = 0, frameRefCount = 90;
	double psize = POINT_SIZES, lsize = LINE_SIZES;
	cl_double2 Window_Center, Window_Dims, Window_Shift;
	void setWindow(int framesPerSec);
	void checkMessage();
	void initGlew();
	void renderDomain();

///////////////////////////////////////////////////////////////////



	cl_command_queue_properties* getQueueProperties() { return queProps; }
	cl_context* getContext() { return &context; }
	cl_device_id* getDevice() { return &device; }
	cl_platform_id* getPlatform() { return &platform; }
	cl_command_queue* getLBqueue() { return &LBqueue; }
	cl_command_queue* getTRqueue() { return &TRqueue; }
	cl_command_queue* getIOqueue() { return &IOqueue; }
	cl_command_queue* getFDqueue() { return &FDqueue; }
	cl_event* getLBevt() { return &LBevent; }
	cl_event* getTRevt() { return &TRevent; }
	cl_event* getIOevt() { return &IOevent; }
	cl_event* getFDevt() { return &FDevent; }

	cl_command_queue TRqueue, IOqueue, FDqueue, LBqueue;
	cl_event TRevent, LBevent, IOevent, FDevent;
	cl_command_queue_properties queProps[4];
	bool turbFlag, trFlag, tempFlag, contextFlag, deviceFlag, queueFlag;
	bool programFl[3];
	bool restartFlag;

	cl_context context;
	cl_device_id device;
	cl_platform_id platform;
	char platformInfo;
	cl_program **programs;

private:
	clEnv()
	{
		turbFlag = false;
		tempFlag = false;
		trFlag = false;
		contextFlag = false;
		deviceFlag = false;
		queueFlag = false;
		programFl[0] = false;
		programFl[1] = false;
		programFl[2] = false;
	}



	static clEnv *s_instance;
	~clEnv()
	{
		if (queueFlag)
		{
			clReleaseCommandQueue(FDqueue);
			clReleaseCommandQueue(LBqueue);
			clReleaseCommandQueue(TRqueue);
			clReleaseCommandQueue(IOqueue);
		}
		if (contextFlag)
			clReleaseContext(context);
		if (deviceFlag)
			clReleaseDevice(device);

		if (programFl[0])
		{
			clReleaseProgram(*programs[0]);
		}
		if (programFl[1])
		{
			clReleaseProgram(*programs[1]);
		}
		if (programFl[2])
		{
			clReleaseProgram(*programs[2]);
		}
	};

};




//void Release_Event(cl_event evtame)
//{
//	//cl_int ret_val;
//	//for (int i = 0; i < 10000; i++)
//	//{
//	//	cl_int ret = clGetEventInfo(evtame, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &ret_val, NULL);
//
//	//	if (ret_val == CL_COMPLETE)
//	//	{
//	//		clReleaseEvent(evt);
//	//		return;
//	//	}
//
//	//	delay(0.01);
//	//}
//	//FILE *stream;
//	//stream = fopen("error.txt", "w+");
//	////	cl_int refcount;
//	////	cl_int ret = clGetEventInfo(evt, CL_EVENT_REFERENCE_COUNT, sizeof(cl_int), &refcount, NULL);
//	////	ret = clGetEventInfo(evt, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &ret_val, NULL);
//	//if (ret_val == CL_QUEUED)
//	//{
//	//	fprintf(stream, "Error releasing event %s, with return CL_QUEUED\n", evt_name);
//	//}
//	//else if (ret_val == CL_SUBMITTED)
//	//{
//	//	fprintf(stream, "Error releasing event %s, with return CL_SUBMITTED\n", evt_name);
//	//}
//	//if (ret_val == CL_RUNNING)
//	//{
//	//	fprintf(stream, "Error releasing event %s, with return CL_RUNNING\n", evt_name);
//	//}
//	//else
//	//{
//	//	fprintf(stream, "Error unknown releasing event %s\n", evt_name);
//	//}
//
//	//fclose(stream);
//
//	//vfl.save_variables();
//
//	//exit(0);
//}
//
//
//
//
//void AddKernels()
//{
//	//		double Mu = MU_NUMBER;// (1. / 3.)*(vlb.tau - 0.5);
//	//		char char_temp[100];
//	//		std::string defines;
//	//		sprintf(char_temp, "#define FULLSIZEX %d\n", p.nX);
//	//		defines = char_temp;
//	//		sprintf(char_temp, "#define FULLSIZEY %d\n", p.Channel_Height);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define CHANNEL_LENGTH %d\n", p.nX);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define CHANNEL_HEIGHT %d\n", p.nY);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define CHANNEL_LENGTH_FULL %d\n", p.XsizeFull);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define DIST_SIZE %d\n", p.XsizeFull*p.nY);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define FULLSIZEXY %d\n", p.Channel_Height*p.nX);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define FULLSIZEY_UT %d\n", p.Channel_Height + 1);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define DOMAIN_SIZE_Y %d\n", p.nY);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define NUMBLOCKS %d\n", numBlocks);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define INI_MASS %g\n", (double)(p.nX*p.Channel_Height));
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define tau0 (%12.11g)\n", LB_TAU_NUMBER);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define omega_f (%12.11g)\n", 2.*LB_TAU_NUMBER);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define omega_f_float (%12.11g)\n", 2.*LB_TAU_NUMBER);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define omega_1f (%12.11g)\n", OMEGA_1F);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define omega_diff (%12.11g)\n", 2.*LB_TAU_NUMBER - OMEGA_1F);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define BETA_INV (%12.11g)\n", 1. / LB_TAU_NUMBER);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define BETA_VAL (%12.11g)\n", LB_TAU_NUMBER);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define BETA_VAL_FLOAT (%12.11g)\n", LB_TAU_NUMBER);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define UX_INLET_VAL (%f)\n", UX_INLET);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define TFD_X_IN_VAL (%f)\n", TFD_X_IN);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define FTERM_VAL (%12.11g)\n", FORCE_TERM_VAL);
//	//		defines.append(char_temp);
//	//#ifdef USE_FLOATS
//	//		sprintf(char_temp, "#define USE_FLOATS\n");
//	//		defines.append(char_temp);
//	//#endif
//	//#ifdef SOA_STORAGE
//	//		sprintf(char_temp, "#define SOA_STORAGE\n");
//	//		defines.append(char_temp);
//	//#endif
//	//
//	//#ifdef INLET_OUTLET_BC
//	//		sprintf(char_temp, "#define INLET_OUTLET_BC\n");
//	//		defines.append(char_temp);
//	//#endif
//	//
//	//
//	//#ifdef DISABLE_TURBULENT_VISC
//	//		sprintf(char_temp, "#define DISABLE_TURBULENT_VISC\n");
//	//		defines.append(char_temp);
//	//#endif
//	//		double ypl = 11.;
//	//		for (int i = 0; i < 11; i++)
//	//		{
//	//			ypl = log(max(ROUGHNESS_FACTOR*ypl, 1)) / 0.41;
//	//		}
//	//		// specific to flow in straight channel
//	//		// based on openfoam's implementation
//	//		// see http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2016/FangqingLiu/openfoamFinal.pdf and openfoam documentation
//	//
//	//		double Returb = sqrt(RE_NUMBER*3.);
//	//		double ypluswall = Returb / PIPE_RADIUS / 2.;
//	//		double Utau = Returb*MU_NUMBER / PIPE_RADIUS;
//	//		double kwall;
//	//		if (ypluswall < ypl)
//	//		{
//	//			double Cf = (1. / ((ypluswall + 11.)*(ypluswall + 11.)) + 2.*ypluswall / 11. / 11. / 11. - 1. / 11. / 11.);
//	//			kwall = Utau*Utau*(2400. / 1.9 / 1.9 * Cf);
//	//		}
//	//		else
//	//		{
//	//			exit(-158234);
//	//		}
//	//		double nutwall = MU_NUMBER * (0.41 * ypluswall / log(ROUGHNESS_FACTOR*ypluswall) - 1.);
//	//		double omega_vis = 6.*MU_NUMBER / (0.075*ypluswall*ypluswall), omega_log = Utau / (0.3*ypluswall*0.41);
//	//		double omegawall = sqrt(omega_vis*omega_vis + omega_log*omega_log);
//	//		vlb.kWallVal = kwall;
//	//		vlb.oWallVal = omegawall;
//	//
//	//		sprintf(char_temp, "#define RE_TURBULENT (%12.11g)\n", Returb);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define UTAU_VAL (%12.11g)\n", Utau);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define YPLUS_WALL_NODE (%12.11g)\n", ypluswall);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define K_WALL_VALUE (%12.11g)\n", kwall);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define NUT_WALL_VALUE (%12.11g)\n", nutwall);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define OMEGA_WALL_VALUE (%12.11g)\n", omegawall);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define PR_TURB_NUMBER (%12.11g)\n", PR_TURB_NUMBER);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define PR_NUMBER (%12.11g)\n", PR_NUMBER);
//	//		defines.append(char_temp);
//	//
//	//
//	//		////////////  Workgroup Sizes ////////////////////////////////////////
//	//		sprintf(char_temp, "#define BLOCK_SIZE %d\n", BLOCKSIZE);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_RED %d\n", WORKGROUPSIZE_RED);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZEX_LB %d\n", WORKGROUPSIZEX_LB);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZEY_LB %d\n", WORKGROUPSIZEY_LB);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_INLET %d\n", WORKGROUPSIZE_INLET);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_IBB %d\n", WORKGROUPSIZE_IBB);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZEX_FD %d\n", WORKGROUPSIZEX_FD);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZEY_FD %d\n", WORKGROUPSIZEY_FD);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_NU %d\n", WORKGROUPSIZE_NU);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_RERELEASE %d\n", (int)WORKGROUPSIZE_RERELEASE);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_TR %d\n", WORKGROUPSIZE_TR);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_TR_WALL %d\n", WORKGROUPSIZE_TR_WALL);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_TR_GL %d\n", WORKGROUPSIZE_TR_GL);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_TR_SHEAR %d\n", WORKGROUPSIZE_TR_SHEAR);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_SORT %d\n", WORKGROUPSIZE_SORT);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_PLOC %d\n", WORKGROUPSIZE_PLOC);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_FL_SHIFT %d\n", WORKGROUPSIZE_FL_SHIFT);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_FL_RAMP %d\n", NUM_INLET_OUTLET_NODES);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_UPDATEM %d\n", WORKGROUPSIZE_UPDATEM);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_UPDATEFD %d\n", WORKGROUPSIZE_UPDATEFD);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_UPDATEWALL %d\n", WORKGROUPSIZE_UPDATEWALL);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_UPDATE_GL %d\n", WORKGROUPSIZE_UPDATE_GL);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_SHIFT_PAR %d\n", WORKGROUPSIZE_SHIFT_PAR);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_TR_DISTS %d\n", WORKGROUPSIZE_TR_DISTS);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define WORKGROUPSIZE_TR_WALL_REFLECT %d\n", WORKGROUPSIZE_TR_WALL_REFLECT);
//	//		defines.append(char_temp);
//	//		////////////////////////////////////////////////////////////////////
//	//#ifdef OPENCL_VERSION_1_2
//	//		sprintf(char_temp, "#define OPENCL_VERSION_1_2\n");
//	//		defines.append(char_temp);
//	//#endif
//	//#ifdef TR_SOLVER
//	//		sprintf(char_temp, "#define START_THERMO_VEL %g\n", (double)START_THERMO_VEL);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define MAX_BL_PER_NODE %d\n", MAX_BL_PER_NODE);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define TRC_NUM_TRACERS %d\n", TRC_NUM_TRACERS);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define NUM_PAR_SIZES %d\n", vtr.Nd);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define FULLSIZEX_TR %d\n", TrDomainSize.x);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define FULLSIZEY_TR %d\n", TrDomainSize.y);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define MU_NUMBER %g\n", Mu);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define DTTR %g\n", c.dTtr);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define DTTR_WALL %g\n", c.dTtr_wall);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define FOUL_SIZE_SWITCH_SOOT2 %g\n", FOUL_SIZE_SWITCH_SOOT2);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define KAIR %g\n", (THERMAL_CONDUCTIVITY_AIR * DELTA_F / DELTA_T));
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define KSOOT %g\n", (THERMAL_CONDUCTIVITY_FOUL * DELTA_F / DELTA_T));
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define ALPHA_DEN_AIR %g\n", vfd.Alpha_den_air);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define ALPHA_DEN_SOOT %g\n", vfd.Alpha_den_soot);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define ALPHA_TEST_AIR %g\n", (THERMAL_CONDUCTIVITY_AIR * DELTA_F / DELTA_T) / vfd.Alpha_den_air);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define ALPHA_TEST_SOOT %g\n", (THERMAL_CONDUCTIVITY_FOUL * DELTA_F / DELTA_T) / vfd.Alpha_den_soot);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define PAR_TIMER_START %d\n", MIN(TIME_BEFORE_RERELEASE, 20000));
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define DEP_TIMER_START %d\n", MIN(TIME_BEFORE_STICK, 20000));
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define X_START_VAL_INT %d\n", X_RELEASE_POS);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define X_MAX_VAL_INT %d\n", X_STOP_POS);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define X_START_VAL %g\n", (double)X_RELEASE_POS);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define X_MIN_VAL %g\n", (double)X_RELEASE_POS - 0.05);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define X_MAX_VAL %6.1f\n", (double)X_STOP_POS);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define Y_MIN_VAL %g\n", 0.5);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define Y_MAX_VAL %g\n", (double)p.nY - 0.5);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define XVALS_HEIGHT %d\n", CHANNEL_HEIGHT + 2);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define NUM_PAR_GL_DIV %d\n", NUM_PAR_GL_DIV);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define ALPHA_FLUID %g\n", vfd.Alpha_fluid);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define ALPHA_FOUL %g\n", vfd.Alpha_foul);
//	//		defines.append(char_temp);
//	//		if (MAX_NUM_BL_ROLL != -1)
//	//		{
//	//			sprintf(char_temp, "#define MAX_NUM_BL_ROLL %d\n", MAX_NUM_BL_ROLL);
//	//			defines.append(char_temp);
//	//		}
//	//
//	//		int count = 1;
//	//		int tempvar = 2;
//	//		while (tempvar < WORKGROUPSIZE_SORT)
//	//		{
//	//			count++;
//	//			tempvar *= 2;
//	//		}
//	//		sprintf(char_temp, "#define SORT_NUM_MERGES %d\n", count);
//	//		defines.append(char_temp);
//	//
//	//		sprintf(char_temp, "constant double Par_multiplier[%d] = { %20.16g", vtr.Nd, (vtr.V_par[0] / LS_SPACING));
//	//		defines.append(char_temp);
//	//		for (int i = 1; i < vtr.Nd; i++)
//	//		{
//	//			sprintf(char_temp, ", %20.16g", (vtr.V_par[i] / LS_SPACING));
//	//			defines.append(char_temp);
//	//		}
//	//		sprintf(char_temp, "};\n");
//	//		defines.append(char_temp);
//	//
//	//
//	//		sprintf(char_temp, "#define PERCENT_USED_IN_SMOOTHING %g\n", PERCENT_USED_IN_SMOOTHING);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define PERCENT_NOT_USED_IN_SMOOTHING %g\n", (1. - PERCENT_USED_IN_SMOOTHING));
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define NEIGHS_PER_SIDE_SMOOTHING %d\n", NEIGHS_PER_SIDE_SMOOTHING);
//	//		defines.append(char_temp);
//	//
//	//
//	//#endif
//	//#ifdef TFD_SOLVER
//	//		sprintf(char_temp, "#define TIN %g\n", vfd.ROE_INX);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define TMIN %g\n", vfd.T_Actual_Min);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define TDIFF %g\n", vfd.T_Actual_Diff);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define FULLSIZEXM1 %d\n", p.nX - 1);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define NU_MULTIPLIER %g\n", 4.*vlb.Pipe_radius);
//	//		defines.append(char_temp);
//	//		sprintf(char_temp, "#define LOCAL_SOLVE_SIZE_FD %d\n", WORKGROUPSIZEX_FD*WORKGROUPSIZEY_FD);
//	//		defines.append(char_temp);
//	//#endif
//
//#ifdef	DEBUG_TURBARR
//	appendKerString("DEBUG_ARRAYS");
//#endif
//
//#ifdef USE_OPENGL
//	appendKerString("USE_OPENGL");
//#endif
//
//	std::ifstream in_heading("source\\Heading.cl");
//	std::string result_heading((std::istreambuf_iterator<char>(in_heading)), std::istreambuf_iterator<char>());
//	KerString.append(result_heading);
//	KerString.append("\n");
//
//	std::ifstream in_struct("source\\Structures.cl");
//	std::string result_struct((std::istreambuf_iterator<char>(in_struct)), std::istreambuf_iterator<char>());
//	KerString.append(result_struct);
//	KerString.append("\n");
//
//	std::ifstream in_out("source\\Output_kernels.cl");
//	std::string result_out((std::istreambuf_iterator<char>(in_out)), std::istreambuf_iterator<char>());
//	KerString.append(result_out);
//
//	std::ifstream in_lb("source\\LB_kernels.cl");
//	std::string result_lb((std::istreambuf_iterator<char>(in_lb)), std::istreambuf_iterator<char>());
//	KerString.append(result_lb);
//	KerString.append("\n");
//
//	if (tempFlag)
//	{
//		std::ifstream in_t("source\\T_kernels.cl");
//		std::string result_t((std::istreambuf_iterator<char>(in_t)), std::istreambuf_iterator<char>());
//		KerString.append(result_t);
//		KerString.append("\n");
//	}
//
//	if (turbFlag)
//	{
//		std::ifstream in_ko("source\\kOmega_kernels_testing.cl");
//		std::string result_ko((std::istreambuf_iterator<char>(in_ko)), std::istreambuf_iterator<char>());
//		KerString.append(result_ko);
//		KerString.append("\n");
//	}
//
//	if (trFlag)
//	{
//		std::ifstream in_tr("source\\TR_kernels.cl");
//		std::string result_tr((std::istreambuf_iterator<char>(in_tr)), std::istreambuf_iterator<char>());
//		KerString.append(result_tr);
//
//		std::ifstream in_sk("source\\Sort_kernels.cl");
//		std::string result_sk((std::istreambuf_iterator<char>(in_sk)), std::istreambuf_iterator<char>());
//		KerString.append(result_sk);
//
//		std::ifstream in_fl("source\\FL_kernels.cl");
//		std::string result_fl((std::istreambuf_iterator<char>(in_fl)), std::istreambuf_iterator<char>());
//		KerString.append(result_fl);
//
//	}
//
//	std::ofstream outfile("temp_file.cl");
//	outfile << KerString;
//	outfile.close();
//}








#endif