// clProblem.h: Class which tracks simulation time, and output
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLPROBLEM_H__INCLUDED_)
#define AFX_CLPROBLEM_H__INCLUDED_

#pragma once

#include "StdAfx.h"
#include "Kernels.h"
#include "clDisplay.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"



#define GetCurrentDir _getcwd
// yaml parameter file parsing/emitting
#define systemParamNum		1
#define fluidParamNum		2
#define thermalParamNum		3
#define kOmegaParamNum		4
#define lsParamNum			5
#define trParamNum			6
#define flParamNum			7

class clProblem 
{
public:
	// Func Pointer for calling loadParams
	std::function<void(void)> loadParamsPtr;
	
	clProblem() : Timeout("time")
	{
		loadParamsPtr = std::bind(&clProblem::loadParams, this);
	};
	~clProblem()
	{
	};

//////////////////////////////////////////////////////////////////////////////////////
//////////             System and Run parameters, functions for              /////////
//////////               Reading from parameter yaml file, and               /////////
//////////            Function to emit to yaml file for restarts             /////////
//////////////////////////////////////////////////////////////////////////////////////
	const static double Pi, R, Kb;
	const static double eps;
	static double DELTA_L, DELTA_T, DELTA_M, DELTA_F, DELTA_P; //delta_p is dX^2*dM/dT^2 (units of energy) Not sure why P is used to identify this.
	 
	// X location of Wavy Section Start, period length of waves, and number of period
	int xWavyStart, wavyPeriodLen, numWavyPeriods;
	
	static double Pipe_radius;
	static int Channel_Height, nX, nY, XsizeFull;
	static unsigned int FullSize, distSize;
	static cl_int2 nn;
	static bool useOpenGL;
	static int DeviceID;
	static unsigned int Time, TimeN;
	static double dTlb, dTtfd, dTtr, dTtr_wall;
	static unsigned int StopTime, StartTime, RestartTime;
	static int trSteps, tfdSteps, trSteps_wall;
	static unsigned int Time_Btw_Clump;
	static unsigned int Next_Clump_Time;
	static unsigned int saveBinStep, saveDumpStep; // per number of dump steps (usually both set to 1)
	static unsigned int saveDumpStepStart; // dump step number to start saving (usually 0)
	
	// nextDumpSave, nextDumpBin indicates next dumpstep to save dump/bin files 
	// currentDumpStep range [0, DumpStepNum].
	static int currentDumpStep, nextDumpSave, nextDumpBin; 
	
	static unsigned int saveAvgStartTime, saveTimeStepAvg, 
		saveStepNumAvg, nextSaveAvgTime;
	static unsigned int saveNuStartTime, saveTimeStepNu,
		saveStepNumNu, nextSaveNuTime;
	static unsigned int saveSSStartTime, saveTimeStepSS,
		saveStepNumSS, nextSaveSSTime;
	static unsigned int saveIOStartTime, saveTimeStepIO,
		saveStepNumIO, nextSaveIOTime;
	static unsigned int dumpStartTime, dumpTimeStep,
		dumpStepNum, nextDumpStepTime; // just output, no saving in results folder
	
	static bool avgNumStepDef, nuNumStepDef, ssNumStepDef,
		ioNumStepDef, dumpNumStepDef;
	
	static unsigned int displaySignalFreq;
	static cl_short IBB_Flag[8];

	static bool restartRunFlag, flOutputDump;


	YAML::Parser parser;
	YAML::Node	yamlIn;
	YAML::Emitter *yamlOut;
	
	// Sets defines for opencl kernels
	void setSourceDefines();

	// Helper function to write parameter values to yaml emitter
	template <typename T>
	void setParameter(std::string pname_, T val_);

	/// writes system parameters to yaml file for restart
	void saveSystemParams();

	// Writes parameters from all methods to yaml file for restart
	void saveParameters();

	/// Calls all loadParameter functions throughout code
	void loadRunParameters(std::string runparams_);

	/// Loads system parameters (time and domain size params)
	void loadParams();

	/// Iterates through yaml file to find document associated with loadFunction
	/// and calls load function. Creates and passes empty doc it doesnt find document
	void parseAndLoadParams(std::string &finrunparams_, int paramnum_,
		std::function<void(void)> &loadparamptr_);

	// creates file with just the paramNameNumber
	// to pass to allow a loadParams function to load defaults
	void parseEmptyFile(int addnum_);

	// Helper function to parse parameter file for various save time info sets
	void getSaveStepInfo(std::string varname_, unsigned int &starttime_,
		unsigned int &numsteps_, unsigned int &timestep_, bool &numdefined_,
		const unsigned int definedstart_, const unsigned int definednumsteps_);

	// Helper function to search for parameter in yaml document
	// returns value from yaml file if found otherwise, default returned
	template <typename T>
	T getParameter(std::string pname_, T defval)
	{
		T retval_;
		if (const YAML::Node *pName = yamlIn.FindValue(pname_))
		{
			*pName >> retval_;
		}
		else
		{
			retval_ = defval;
		}
		return retval_;
	}

	template <typename T>
	bool getParameter(std::string pname_, T& pval_, T defval)
	{
		bool retval_;
		if (const YAML::Node *pName = yamlIn.FindValue(pname_))
		{
			*pName >> pval_;
			retval_ = true;
		}
		else
		{
			pval_ = defval;
			retval_ = false;
		}
		return retval_;
	}

	// returns 0 if neither parameters provided in file, 1 if pname1_ is found,
	// 2 if pname2_ is found, and throws error if both are found.
	template <typename T>
	int getParameter(std::string pname1_, std::string pname2_, T &val1_, T &val2_)
	{
		int retval_ = 0;
		if (const YAML::Node *pName = yamlIn.FindValue(pname1_))
		{
			*pName >> val1_;
			retval_ = 1;
		}
		if (const YAML::Node *pName = yamlIn.FindValue(pname2_))
		{
			*pName >> val2_;

			ERROR_CHECKING(retval_ == 1, "Cannot define both " + pname1_ + " and "\
				+ pname2_ + "in parameters file", ERROR_LOADING_PARAMS);

			retval_ = 2;
		}

		return retval_;
	}


	// searches for parameter and returns value from yaml file
	// throws error and exits if parameter not found. (for non-default parameters)
	template <typename T>
	T getParameter(std::string pname_);

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

	std::string outDir;
	Array1Dd Timeout;
	int sort_called = 0;
	int update_walls = 0;


	// Simulation control
	void ini();
	void step();
	void finish();

	// Update output functions
	//void updateAvgs();
	//void updateNu();
	//void updateSS();
	//void updateIODist();
	
	// Save output functions
	//void SaveAvgs();
	//void SaveNu();
	//void SaveSS();
	//void SaveIO();
	void DumpStep();
	void RenameOutputFiles();
	void freeHostArrays();




	/// Convenience Functions for file operations
	//void CopyFile(std::string NameSrc, std::string NameDest);
	//void MakeDir(std::string &NewDir);
	//void cleanfiles();


	// Empty functions to allow for templated function pointer wrapper
	// without having to use template specialization or SFINAE
	void createKernels() {}
	void setKernelArgs() {}
};

#endif // !defined(AFX_CLPROBLEM_H__INCLUDED_)
