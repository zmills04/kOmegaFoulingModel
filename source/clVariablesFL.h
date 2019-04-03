// clVariablesFL.h: Fouling layer class.
//
// (c) Zach Mills, 2015 
//////////////////////////////////////////////////////////////////////

//TODO: need to save necessary variables (both as txt and bin files)
//	need to clean up unncessary functions and variables and comment
//	need to set up load functions for restarting simulation

#if !defined(AFX_CLVARIABLESFL_H__INCLUDED_)
#define AFX_CLVARIABLESFL_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clProblem.h"



class clVariablesFL
{
public:
	// Func Pointer for calling loadParams
	std::function<void(void)> loadParamsPtr;
	

	clVariablesFL()
	{
		loadParamsPtr = std::bind(&clVariablesFL::loadParams, this);
	};
	virtual ~clVariablesFL()
	{
		//Release_Objects(); 
	};

	Kernel update_FL_kernel[4];
	Kernel update_LS_kernel[5];
	Kernel update_LB_kernel[2];
	Kernel update_FD_kernel;
	Kernel update_TR_kernel[5];
	Kernel update_GL_kernel;

	Array1DFI FI;
	Array1DRI RI;
	Array2Dd BLdep_tot, BLdep_tot_temp;
	Array1Dd IO_ind_dist;
	Array1Dd Sum_M_temp;
	Array1Dd Debug_out;

	bool flSolverFlag, saveOnStartFlag;
	bool restartRunFlag;
	cl_uint4 IO_end;
	cl_uint Num_active_nodes;
	cl_uint Num_IO_indicies;
	int FL_timer, flTimePerUpdate;
	int Smooth_timer, flTimePerSmooth, neighsPerSideSmoothing;
	int SaveStepNum = 0;
	double smoothingPct;

	void setSourceDefines();
	void testRestartRun();
	void saveParams();
	void loadParams();
	void ini();
	void ini_foulI();
	void ini_IO_vars();//done
	cl_double2 get_center(cl_double2 P0, cl_double2 P1);
	void update_FL();
	void update_FD(cl_event *wait);
	void update_LS(cl_event *evt);
	void update_LB();
	void update_TR(cl_event *wait_fill, int Num_Wnodes_temp);
	void update_GL();
	void update();
	void setKernelArgs();
	void Release_Objects();
	void allocateBuffers();
	void allocateArrays();
	void ini_group_sizes();
	void save2file();
	void save_variables();
	void createKernels();
	//saves variables necessary for debugging
	void savedebug();
	void freeHostMem();
	BOOL test_bounds();
	BOOL Restart_Run();
	double Gaussian_Kernel(double input);
	void saveRestartFiles();
	void UpdateRestart();
	void RenameDebug_Files(int dirnumber);
	void CallRename(char* file, const char *fol);
};


#endif