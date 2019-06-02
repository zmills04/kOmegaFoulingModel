// Class used to generate openCL source code. Mostly builds source code
// from kernel files (unlike Sparse and Reduce generators). iniList
// contains pointers to function creating kernels and setting arguments,
// so that vlb, vfd, etc only need a single initialization function.
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SOURCEGENERATOR_H__INCLUDED_)
#define AFX_SOURCEGENERATOR_H__INCLUDED_

#pragma once

#include "StdAfx.h"


#ifndef CHECK_KERNEL_STATUS
#define STATUS_CHECK_RUN(codeval, msg) do {} while(0)
#else
#define STATUS_CHECK_RUN(codeval, msg)	if(codeval) { Gen_Error_Msg(codeval, msg); }
#endif

#define STATUS_CHECK_SETUP(codeval, msg)	if(codeval) { Gen_Error_Msg(codeval, msg); }


#define GetSourceProgram		*sourceGenerator::SourceInstance()->getProgram()
#define GetSourceProgramAdd		sourceGenerator::SourceInstance()->getProgram()
#define SOURCEINSTANCE			sourceGenerator::SourceInstance()

template <typename T>
struct argSetter
{
	Kernel *ker;
	int ind;
	size_t varsize;
	T var;
	argSetter(Kernel *ker_, const int ind_, T var_)
	{
		ker = ker_;
		ind = ind_;
		varsize = sizeof(T);
		var = var;
	}
	argSetter(const size_t locmemsize, Kernel *ker_, const int ind_)
	{
		ker = ker_;
		varsize = locmemsize;
		ind = ind_;
		var = nullptr;
	}
};



class sourceGenerator
{
private:
	static sourceGenerator *s_instance;
	sourceGenerator(const std::string type_, const clEnv::programType ptype_)
	{
		type = type_;
		pType = ptype_;
		setupGenerator();
		addFile2Kernel("Heading.cl");
		addFile2Kernel("Structures.cl");
		addFile2Kernel("outputKernels.cl");
	}

	sourceGenerator()
	{
		setupGenerator();
	}

	~sourceGenerator()
	{
		if (programFl)
		{
			clReleaseProgram(program);
		}
	}

public:
	friend class ReduceGenerator;
	friend class BiCGStabGenerator;
	static sourceGenerator *SourceInstance()
	{
		if (!s_instance)
			s_instance = new sourceGenerator("sourceGenerator", clEnv::BaseType);
		return s_instance;
	}

	void setupGenerator()
	{
		programSrc = "";
		context = clEnv::instance()->getContext();
		device = clEnv::instance()->getDevice();
		ioQue = clEnv::instance()->getIOqueue();
		programFl = false;
	}



	cl_program* getProgram() { return &program; }
	std::string& getDefineStr() { return defineStr; }

	std::string type;
	clEnv::programType pType;





	// This will mostly store info for derived classes, since this class
	// is generating kernels mostly from pre-written kernels
	typedef std::map<unsigned int, std::string> KernelMap;

	// Return value for hashing function
	enum hashResult { newKernel, oldKernel };

	// Appends a string containing a kernel function
	// to the source string to be built
	void addString2Kernel(const std::string &kerStr);

	void addString2Defines(const std::string &defstr_);

	// Appends a kernel from a file to the source string
	void addFile2Kernel(const std::string &fname_);

	// Helper functions which generate a string containing 
	// a #define statement with the name name_ and value of val
	void addDefine(std::string &kerstr, const std::string varname, double value);

	void addDefine(std::string &kerstr, const std::string varname, int value);

	void addDefine(std::string &kerstr, const std::string varname, std::string value);

	void addDefine(std::string &kerstr, const std::string varname);

	// Appends kernel initialization functions to the lists, which will
	// call functions after program is compiled
	void addIniFunction(std::function<void(void)> &createptr_, std::function<void(void)> &setargptr_)
	{
		setArgumentList.push_back(setargptr_);
		createKernelList.push_back(createptr_);
	}

	// Compiles openCL kernels
	void buildSource()
	{
		// need to add defines to beginning of base programs
		if (pType == clEnv::BaseType)
		{
			defineStr.append(programSrc);
			programSrc = defineStr;
		}
		std::ofstream outfile("x64" SLASH "clSource" SLASH + type + ".cl");
		outfile << programSrc;
		outfile.close();
		clEnv::instance()->buildProgram(program, programSrc, pType, type);
		programFl = true;
		callIniKernels();
	}

	// Helper function to replace strings when generating kernels
	void findAndReplace(std::string & data, std::string toSearch,
		std::string replaceStr);

	// Generates error message in file (if LOG_ERROR_IN_FILE is defined),
	// prints to cli and exits with code Errcode.
	void Gen_Error_Msg(int Errcode, std::string mesg);

	// Iterates through iniList and calls initialization functions
	virtual void callIniKernels();

	// Hashing function used to check if a kernel has already
	// been created mostly used by derived classes
	static unsigned int rsHash(const std::string &key);

	// Check is a kernel has already been built, and if so returns the
	// name of the function instead of appending the string to the 
	// source string. 
	hashResult checkHash(std::string &kerstr, std::string &kername);

	// kernel map contains the hashes of kernels already created 
	// using the generator (created and appended to the source
	// string, not compiled)
	KernelMap kernel_map;

	// Contains pointers to functions which will create kernels and
	// set arguments
	std::list<std::function<void(void)>> createKernelList;
	std::list<std::function<void(void)>> setArgumentList;



	// These openCL variables are used to create
	// kernels
	static cl_context *context;
	static cl_device_id *device;
	static cl_command_queue *ioQue;


	// String containing kernels to be built
	std::string programSrc;
	std::string defineStr;

	// clEnv will keep a pointer to all programs as well
	cl_program program;
	bool programFl;
};














#endif