// StdAfx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
// (c) Alexander Alexeev, 2006 
//

#ifndef _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif
//#define _CRT_SECURE_NO_DEPRECATE
#pragma warning(disable : 4996)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 6054)

#if !defined(AFX_STDAFX_H__24F045E1_A23A_42BE_BCBB_F3F2F6EE9AF1__INCLUDED_)
#define AFX_STDAFX_H__24F045E1_A23A_42BE_BCBB_F3F2F6EE9AF1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(SLASH)
#if defined(_WIN32)
#define SLASH "\\"
#else 
#define SLASH "/"
#endif //_WIN32
#endif //SLASH

#if !defined(_WIN32)
#define _OPEN_MP
#endif //_WIN32

#if defined(_WIN32)
#define _OPERATION_SYSTEM  "WIN32"
#else 
#define _OPERATION_SYSTEM  "UNIX"
#endif //_WIN32

// TODO: reference additional headers your program requires here
// TODO: remove unnecessary functions, defines and includes
// TODO: add comments
// TODO: move functions to separate file



#include "ConstDef.h"


#ifdef OPENCL_VERSION_1_2
#define CL_HPP_TARGET_OPENCL_VERSION	120
#define CL_HPP_CL_1_2_DEFAULT_BUILD
#endif



//#ifdef USE_OPENGL
#include <Windows.h>
#include <GL/glew.h>
#include <CL/cl_gl.h>
#include <windowsx.h>
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")
//#endif
#include <limits.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>


#include <cassert>
#include <clsparse_export.h>
#include <vector>
#include <stdio.h>
#include <map>
#include <list>
#include <fstream>
#include <sstream>
#if defined(_WIN32)
#include <direct.h>
#endif //_WIN32
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <string>
#include <CL/opencl.h>
#include <CL\cl2.hpp>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <thread>
#include <clSPARSE.h>
#include <clSPARSE-error.h>
#pragma warning( push )
#pragma warning( disable : 4146)
#include <yaml-cpp/yaml.h>
#pragma warning( pop )

typedef unsigned int localSize_t;
typedef cl_uint globalSize_t;

typedef struct cl_bool2
{
	cl_bool2(bool b1_, bool b2_)
	{
		x = b1_;
		y = b2_;
	}
	bool x;
	bool y;
}cl_bool2;

#pragma OPENCL_EXTENSION cl_amd_fp64: enable
#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_subgroups : enable

#define LOGMESSAGE(mess_)	Logger::Instance()->writeToLogFile(mess_)
#define LOGPARTIAL(mess_)	Logger::Instance()->writePartialToLogFile(mess_)

#define PRINT_ERROR_FUNC(msg, error_num)\
	std::string ERRMESG_ = "Error in ";\
	ERRMESG_.append(__FILE__);\
	ERRMESG_.append(" at line ");\
	ERRMESG_.append(std::to_string(__LINE__));\
	ERRMESG_.append("\nMessage: ");\
	ERRMESG_.append(msg);\
	LOGMESSAGE(ERRMESG_);\
	std::cout << "Error in " <<  __FILE__  << " at line " << __LINE__ << "\n";\
	std::cout << "See log file for more information\n";\
	exit(error_num);

#define ERROR_CHECKING(test, msg, error_num)\
	if ((test))\
	{\
		PRINT_ERROR_FUNC(msg, error_num);\
	}

#define ERROR_CHECKING_OCL(test, msg, error_num) \
	if ((test)) \
	{\
		PRINT_ERROR_FUNC("OpenCL Error Number: " + std::to_string(static_cast<int>(test)) + "\n" + msg, error_num); \
	}


#define FREE(ptr) \
{ \
if (ptr != NULL) \
{ \
	free(ptr); \
	ptr = NULL; \
} \
}

#define FREE_OCL_BUFFER(buf) if(buf != nullptr) {clReleaseMemObject(buf);}

#define FREE_OCL_BUFFER_ARRAY(bufarr_, size_) \
if (bufarr_ != nullptr) \
{ \
	for (unsigned int i = 0; i < size_; i++) \
	{ \
		clReleaseMemObject(bufarr_[i]); \
	} \
	delete[] bufarr_; \
}

enum FileMode { FileIn, FileOut, FileAppend, BinaryIn, BinaryOut };

enum statKernelType {
	uniformKer, triangularKer, parabolicKer,
	quarticKer, triWeightKer, triCubeKer, gaussianKer, cosineKer,
	logisticKer, sigmoidKer
};

//#ifndef PRINT_ERROR_MESSAGES
//#define PRINT_ERROR					do {} while(0)
//#define PRINT_ERROR_W_MSG(msg_)		do {} while(0)
//#else
//#define PRINT_ERROR					FILE *stream;\
//									stream = fopen("error.txt", "w+");\
//									fprintf(stream, "Error in %s at line %d\n", __FILE__, __LINE__);\
//									fclose(stream);\
//									printf("Error in %s at line %d\n", __FILE__, __LINE__);\
//									exit(__LINE__);
//
//#define PRINT_ERROR_W_MSG(msg_)		FILE *stream;\
//									stream = fopen("error.txt", "w+");\
//									fprintf(stream, "Error in %s at line %d\n%s\n", __FILE__, __LINE__, msg_);\
//									fclose(stream);\
//									printf("Error in %s at line %d\n%s\n", __FILE__, __LINE__, msg_);\
//									exit(__LINE__);
//#endif

//typedef struct Nodetemp
//{
//	cl_double2 dX;		//Initial spacings between lattice points and walls
//	cl_double2 dX_cur;	//current spacings between lattice points and walls
//	cl_int4 neigh;		//Indicies of 4 lattice points enclosing square (U and T array are of size vlb.nX*(vlb.Channel_Height+1))
//	//where additional lattice point is used to represent lattice points in wall
//	cl_short BLind[MAX_BL_PER_NODE];	//Inidices of three closest boundary links (-1 used to represent when less than three BL's near)
//	cl_short type;		//type of of node (determined by types of nodes neighbors are)
//	cl_char Wall_Flag;  //0 is no walls nearby, 1 if bottom wall nearby and 2 if top wall nearby
//	cl_char Pad[(14 - MAX_BL_PER_NODE * 2)];
//} nodet;

//typedef struct Ramp_info
//{
//	double Yval_start[4];
//	cl_uint IO_end[4]; //Index of LS nodes at beginning and end of TR area
//	cl_int Ind_dir[4]; //direction traversed from IO_end when creating ramp
//} rampI;

//typedef struct Ramp_info
//{
//	double Ybegin;
//	double Coeff;
//	cl_uint IOind; //Index of LS nodes at beginning and end of TR areatypename
//	cl_uint Cind; //direction traversed from IO_end when creating ramp
//	cl_char pad[8];
//} rampI;



#define ARRAY_ERROR_CHECKING
//#define STATUS_CHECK_SETUP(codeval, msg)	if(codeval) { Gen_Error_Msg(codeval, msg); }


#define Subtract2(A, B) { {A.x - B.x, A.y - B.y} } 
#define Add2(A, B)		{ {A.x + B.x, A.y + B.y}}
#define GETLEN(A)		(sqrt(A.x*A.x + A.y*A.y))
#define Multiply2(A, B)	{{A.x*B, A.y*B}}
#define Divide2(A, B)	{{A.x/B, A.y/B}}
#define VNORM(A)		Divide2(A, GETLEN(A))
#define DOT_PROD(A, B)	(A.x*B.x + A.y*B.y)
#define MOD2(A, B)		{{ MOD(A.x, B.x), MOD(A.y, B.y)}}
#define DOT_PROD_3D(A, B)	(A.x*B.x + A.y*B.y + A.z*B.z)

#define MOD(A, B)		(((A) % (B) + (B)) % (B))
#define MODFAST(A, B)	(((A) >= (B)) ? ((A) - (B)) : (((A) < 0) ? ((A) + (B)) : (A)))
#define MAX(A, B)		(((A) > (B)) ? (A) : (B))
#define MIN(A, B)		(((A) < (B)) ? (A) : (B))
#define COMPINT2(A, B)		((A.x == B.x) && (A.y == B.y))

#ifdef _DEBUG
#define LB_DEBUG
//#define LS_DEBUG
//#define TR_DEBUG
#endif //_DEBUG

#define VERSION		6.00


#if !defined(__EXTERN_VARIABLES__)
#define __EXTERN_VARIABLES__

//class sourceGenerator;
//class BiCGStabGenerator;
//class ReduceGenerator;
class clVariablesLS;
class clVariablesLB;
class clVariablesTR;
class clVariablesFD;
class clVariablesFL;
class clProblem;
class clDisplay;

extern clVariablesLS vls;
extern clVariablesLB vlb;
extern clVariablesTR vtr;
extern clVariablesFD vfd;
extern clVariablesFL vfl;
extern clProblem p;
extern clDisplay d;

#endif // !defined(__EXTERN_VARIABLES__)

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.



#endif // !defined(AFX_STDAFX_H__24F045E1_A23A_42BE_BCBB_F3F2F6EE9AF1__INCLUDED_)



