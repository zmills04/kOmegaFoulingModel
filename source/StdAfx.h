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
#include <clSPARSE.h>
#include <clSPARSE-error.h>
#pragma warning( push )
#pragma warning( disable : 4146)
#include <yaml-cpp/yaml.h>
#pragma warning( pop )

typedef unsigned int localSize_t;
typedef cl_uint globalSize_t;
typedef cl_uint2 cl_bool2;

#pragma OPENCL_EXTENSION cl_amd_fp64: enable
#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_subgroups : enable
//#pragma OPENCL EXTENSION CL_KHR_gl_sharing : enable


#define SDK_SUCCESS 0
#define SDK_FAILURE 1
#define SDK_EXPECTED_FAILURE 2


#define CHECK_ALLOCATION(actual, msg) \
if (actual == NULL) \
{ \
	std::cout << "Location : " << __FILE__ << ":" << __LINE__ << std::endl; \
	return SDK_FAILURE; \
}

#define CHECK_ERROR(actual, reference, msg) \
if (actual != reference) \
{ \
	std::cout << "Location : " << __FILE__ << ":" << __LINE__ << std::endl; \
	return SDK_FAILURE; \
}

#define CHECK_ERROR_WITH_EXIT(actual, reference, msg) \
if (actual != reference) \
{ \
	std::cout << "Location : " << __FILE__ << ":" << __LINE__ << std::endl; \
	exit(actual);\
}

#define PRINT_ERROR_FUNC(msg, error_num) \
	std::fstream stream;\
	stream.open("error.txt", std::ios_base::out);\
	stream << "Error in " <<  __FILE__  << " at line " << __LINE__ << "\n";\
	stream << "Message: " << msg << "\n";\
	stream.close();\
	std::cout << "Error in " <<  __FILE__  << " at line " << __LINE__ << "\n";\
	std::cout << "See error.txt for more information\n";\
	exit(error_num);


#define ERROR_CHECKING(test, msg, error_num) if ((test)) { PRINT_ERROR_FUNC(msg, error_num); }

#define ERROR_CHECKING_OCL(test, msg) \
if ((test)) \
{ \
	FILE *stream;\
	stream = fopen("error.txt", "w+");\
	fprintf(stream, "Error in %s at line %d\n", __FILE__, __LINE__);\
	fprintf(stream, "OpenCL Error Number: %d\n", test);\
	fprintf(stream, "Message: " msg "\n");\
	fclose(stream);\
	printf("Error in %s at line %d\n. See error.txt for more information\n", __FILE__, __LINE__);\
	exit(test);\
}

#define FREE(ptr) \
{ \
if (ptr != NULL) \
{ \
	free(ptr); \
	ptr = NULL; \
} \
}

#define RELEASE_EVENT(ptr) \
{ \
if (ptr != NULL) \
{ \
	clReleaseEvent(ptr); \
} \
}

#ifdef _WIN32
#define ALIGNED_FREE(ptr) \
{ \
if (ptr != NULL) \
{ \
	_aligned_free(ptr); \
	ptr = NULL; \
} \
}
#endif
#define CHECK_OPENCL_ERROR(actual, msg) \
if (checkVal(actual, CL_SUCCESS, msg)) \
{ \
	std::cout << "Location : " << __FILE__ << ":" << __LINE__ << std::endl; \
	return SDK_FAILURE; \
}

#define OPENCL_EXPECTED_ERROR(msg) \
{ \
	return SDK_EXPECTED_FAILURE; \
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




template<typename T>
static int checkVal(
	T input,
	T reference,
	std::string message, bool isAPIerror = true)
{
	if (input == reference)
	{
		return SDK_SUCCESS;
	}
	else
	{
		if (isAPIerror)
		{
			std::cout << "Error: " << message << " Error code : ";
		}
		return SDK_FAILURE;
	}
}


enum FileMode { FileIn, FileOut, FileAppend, BinaryIn, BinaryOut };



#ifndef PRINT_ERROR_MESSAGES
#define PRINT_ERROR					do {} while(0)
#define PRINT_ERROR_W_MSG(msg_)		do {} while(0)
#else
#define PRINT_ERROR					FILE *stream;\
									stream = fopen("error.txt", "w+");\
									fprintf(stream, "Error in %s at line %d\n", __FILE__, __LINE__);\
									fclose(stream);\
									printf("Error in %s at line %d\n", __FILE__, __LINE__);\
									exit(__LINE__);

#define PRINT_ERROR_W_MSG(msg_)		FILE *stream;\
									stream = fopen("error.txt", "w+");\
									fprintf(stream, "Error in %s at line %d\n%s\n", __FILE__, __LINE__, msg_);\
									fclose(stream);\
									printf("Error in %s at line %d\n%s\n", __FILE__, __LINE__, msg_);\
									exit(__LINE__);
#endif

typedef struct Nodetemp
{
	cl_double2 dX;		//Initial spacings between lattice points and walls
	cl_double2 dX_cur;	//current spacings between lattice points and walls
	cl_int4 neigh;		//Indicies of 4 lattice points enclosing square (U and T array are of size vlb.nX*(vlb.Channel_Height+1))
	//where additional lattice point is used to represent lattice points in wall
	cl_short BLind[MAX_BL_PER_NODE];	//Inidices of three closest boundary links (-1 used to represent when less than three BL's near)
	cl_short type;		//type of of node (determined by types of nodes neighbors are)
	cl_char Wall_Flag;  //0 is no walls nearby, 1 if bottom wall nearby and 2 if top wall nearby
	cl_char Pad[(14 - MAX_BL_PER_NODE * 2)];
} nodet;

//typedef struct Ramp_info
//{
//	double Yval_start[4];
//	cl_uint IO_end[4]; //Index of LS nodes at beginning and end of TR area
//	cl_int Ind_dir[4]; //direction traversed from IO_end when creating ramp
//} rampI;

typedef struct Ramp_info
{
	double Ybegin;
	double Coeff;
	cl_uint IOind; //Index of LS nodes at beginning and end of TR areatypename
	cl_uint Cind; //direction traversed from IO_end when creating ramp
	cl_char pad[8];
} rampI;



//inline std::istream& operator>>(std::istream& is, const par& el)
//{
//	return is;
//}
//
//inline std::istream& operator>>(std::istream& is, const bLinks& el)
//{
//	return is;
//}


inline std::istream& operator>>(std::istream& is, const rampI& el)
{
	return is;
}

//inline std::istream& operator>>(std::istream& is, const foulI& el)
//{
//	return is;
//}
//

//inline std::istream& operator>>(std::istream& is, const Pparam& el)
//{
//	return is;
//}

inline std::istream& operator>>(std::istream& is, const nodet& el)
{
	return is;
}

//
//inline std::istream& operator>>(std::istream& is, const nodeI& el)
//{
//	return is;
//}
//
//inline std::istream& operator>>(std::istream& is, const nodeV& el)
//{
//	return is;
//}
//
//inline std::istream& operator>>(std::istream& is, const nodeC& el)
//{
//	return is;
//}

//inline std::istream& operator>>(std::istream& is, const Nact& el)
//{
//	return is;
//}



//inline std::ostream& operator<<(std::ostream& os, const BLbound& el)
//{
//	os << "BL_bounds: Bottom (min, max) = (" << el.MIN_BL_BOT << ", " << el.MAX_BL_BOT << ")\n";
//	os << "Top (min, max) = (" << el.MIN_BL_TOP << ", " << el.MAX_BL_TOP << ")\n\n";
//	return os;
//}
//
//inline std::ostream& operator<<(std::ostream& os, const Trparam& el)
//{
//	os << "Tr Params:\n";
//	os << "\tTop Location = " << el.Top_location << "\n";
//	os << "\tBottom Location = " << el.Bottom_location << "\n";
//	os << "\tUmax Value = " << el.umax_val << "\n";
//	os << "\tWall spacing at inlet = " << el.bval << "\n";
//	os << "\tBottom Wall Location at Inlet = " << el.offset_y << "\n";
//	os << "\tX position of Release = " << el.X_release << "\n";
//	os << "\tBL at bottom where released = " << el.BL_rel_bot << "\n";
//	os << "\tBL at top where released = " << el.BL_rel_top << "\n";
//	os << "\tUvals_start = " << el.Uvals_start << "\n";
//	return os;
//}
//
//inline std::ostream& operator<<(std::ostream& os, const par& el)
//{
//#ifdef DUMP_PAR_FULL
//#ifdef ARRAY_DEBUG
//	os << el.type << "\t" << el.Num_rep << "\t" << el.Dep_Flag << "\t"
//		<< el.Dep_timer << "\t" << el.timer << "\t" << el.loc << "\t"
//		<< std::setw(10) << el.pos.x << "\t" << std::setw(10) << el.pos.y;
//#else
//	os << el.type << "\t" << el.Num_rep << "\t" << el.Dep_Flag << "\t"
//		<< el.Dep_timer << "\t" << el.timer << "\t" << el.loc << "\t"
//		<< el.pos.x << "\t" << el.pos.y;
//#endif
//
//#else
//
//#ifdef ARRAY_DEBUG
//	os << std::setw(10) << el.pos.x << "\t" << std::setw(10) << el.pos.y;
//#else
//	os << el.pos.x << "\t" << el.pos.y;
//#endif
//#endif
//	return os;
//}
//
//
//// Saves full BL structure
//
//inline std::ostream& operator<<(std::ostream& os, const bLinks& el)
//{
//#ifdef DUMP_BLINKS_FULL
//#ifdef ARRAY_DEBUG
//	os << std::setw(10) << el.vP0.x << "\t" << std::setw(10) << el.vP0.y <<
//		"\t" << std::setw(10) << el.vP1.x << "\t" << std::setw(10) << el.vP1.y <<
//		"\t" << std::setw(10) << el.vTvec.x << "\t" << std::setw(10) << el.vTvec.y <<
//		"\t" << std::setw(10) << el.vNvec.x << "\t" << std::setw(10) << el.vNvec.y <<
//		"\t" << std::setw(10) << el.Tau << "\t" << std::setw(10) << el.blLen <<
//		"\t" << el.Node_loc << "\t" << el.Color_ind << "\t" << el.P0ind << "\t" <<
//		el.P1ind << "\t" << el.dir;
//#else
//	os << el.vP0.x << "\t" << el.vP0.y << "\t" << el.vP1.x << "\t" << el.vP1.y <<
//		"\t" << el.vTvec.x << "\t" << el.vTvec.y << "\t" << el.vNvec.x << "\t" <<
//		el.vNvec.y << "\t" << el.Tau << "\t" << el.blLen << "\t" << el.Node_loc <<
//		"\t" << el.Color_ind << "\t" << el.P0ind << "\t" << el.P1ind << "\t" <<
//		el.dir;
//#endif
//
//#else
//
//#ifdef ARRAY_DEBUG
//	os << std::setw(10) << el.Tau;
//#else
//	os << el.Tau;
//#endif
//
//#endif
//	return os;
//}


inline std::ostream& operator<<(std::ostream& os, const rampI& el)
{
#ifdef ARRAY_DEBUG
	os << std::setw(10) << el.Ybegin << "\t" << std::setw(10) << el.Coeff <<
		"\t" << el.IOind << "\t" << el.Cind;
#else
	os << el.Ybegin << "\t" << el.Coeff << "\t" << el.IOind <<
		"\t" << el.Cind;
#endif
	return os;
}


//
//inline std::ostream& operator<<(std::ostream& os, const foulI& el)
//{
//#ifdef ARRAY_DEBUG
//	os << std::setw(10) << el.WeightsL.x << "\t" << std::setw(10) <<
//		el.WeightsL.y << "\t" << std::setw(10) << el.WeightsL.z << "\t" << std::setw(10) <<
//		el.WeightsL.w << "\t" << std::setw(10) << el.WeightsR.x << "\t" << std::setw(10) <<
//		el.WeightsR.y << "\t" << std::setw(10) << el.WeightsR.z << "\t" << std::setw(10) <<
//		el.WeightsR.w << "\t" << std::setw(10) << el.vN.x << "\t" << std::setw(10) <<
//		el.vN.y << "\t" << std::setw(10) << el.disp << "\t";
//	os << el.BL_ind << "\t" << el.C_ind;
//#else
//	os << el.WeightsL.x << "\t" << el.WeightsL.y << "\t" << el.WeightsL.z << "\t" << el.WeightsL.w << "\t" << el.WeightsR.x << "\t" << el.WeightsR.y
//		<< "\t" << el.WeightsR.z << "\t" << el.WeightsR.w << "\t" << el.vN.x << "\t" << el.vN.y << "\t" << el.disp << "\t";
//	os << el.BL_ind << "\t" << el.C_ind;
//#endif
//	return os;
//}
//
//
//inline std::ostream& operator<<(std::ostream& os, const Pparam& el)
//{
//#ifdef ARRAY_DEBUG
//	os << std::setw(10) << el.Dp << "\t" << std::setw(10) << el.Mp << "\t" <<
//		std::setw(10) << el.Q_A_prime.x << "\t" << std::setw(10) << el.Q_A.x <<
//		"\t" << std::setw(10) << el.tau_crit.x << "\t" << std::setw(10) <<
//		el.tau_crit.y << "\t" << std::setw(10) << el.Kth << "\t" << std::setw(10) <<
//		el.D_dist << "\t" << std::setw(10) << el.L_coeff << "\t" << std::setw(10) <<
//		el.D_coeff;
//#else
//	os << el.Dp << "\t" << el.Mp << "\t" << el.Q_A_prime.x << "\t" << el.Q_A.x <<
//		"\t" << el.tau_crit.x << "\t" << el.tau_crit.y << "\t" << el.Kth << "\t" <<
//		el.D_dist << "\t" << el.L_coeff << "\t" << el.D_coeff;
//#endif
//	return os;
//}



inline std::ostream& operator<<(std::ostream& os, const nodet& el)
{
#ifdef ARRAY_DEBUG
	os << std::setw(10) << el.dX.x << "\t" << std::setw(10) << el.dX.y << "\t"
		<< std::setw(10) << el.dX_cur.x << "\t" << std::setw(10) << el.dX_cur.y <<
		"\t" << el.neigh.x << "\t" << el.neigh.y << "\t" << el.neigh.z << "\t" <<
		el.neigh.w << "\t" << el.type << "\t" << el.Wall_Flag << "\t";
#else
	os << el.dX.x << "\t" << el.dX.y << "\t" << el.dX_cur.x << "\t" << el.dX_cur.y <<
		"\t" << el.neigh.x << "\t" << el.neigh.y << "\t" << el.neigh.z << "\t" <<
		el.neigh.w << "\t" << el.type << "\t" << el.Wall_Flag << "\t";
#endif
	for (int k = 0; k < MAX_BL_PER_NODE - 1; k++)
	{
		os << el.BLind[k] << "\t";
	}
	os << el.BLind[MAX_BL_PER_NODE - 1];
	return os;
}

//
//inline std::ostream& operator<<(std::ostream& os, const nodeI& el)
//{
//	for (int k = 0; k < MAX_BL_PER_NODE; k++)
//	{
//		os << el.BLind[k] << "\t";
//	}
//	os << el.Wall_Flag;
//	return os;
//}
//
//
//
//
//inline std::ostream& operator<<(std::ostream& os, const nodeV& el)
//{
//#ifdef ARRAY_DEBUG
//	os << std::setw(10) << el.Temps.x << "\t" << std::setw(10) << el.Temps.y << "\t"
//		<< std::setw(10) << el.Temps.z << "\t" << std::setw(10) << el.Temps.w <<
//		"\t";
//	os << std::setw(10) << el.U00.x << "\t" << std::setw(10) << el.U00.y << "\t"
//		<< std::setw(10) << el.U10.x << "\t" << std::setw(10) << el.U10.y <<
//		"\t";
//	os << std::setw(10) << el.U01.x << "\t" << std::setw(10) << el.U01.y << "\t"
//		<< std::setw(10) << el.U11.x << "\t" << std::setw(10) << el.U11.y;
//#else
//	os << el.Temps.x << "\t" << el.Temps.y << "\t" << el.Temps.z << "\t" <<
//		el.Temps.w << "\t" << el.U00.x << "\t" << el.U00.y << "\t" <<
//		el.U10.x << "\t" << el.U10.y << "\t" << el.U01.x << "\t" <<
//		el.U01.y << "\t" << el.U11.x << "\t" << el.U11.y;
//#endif
//	return os;
//}
//
//
//
//
//inline std::ostream& operator<<(std::ostream& os, const nodeC& el)
//{
//#ifdef ARRAY_DEBUG
//	os << std::setw(10) << el.CoeffT00.x << "\t" << std::setw(10) << el.CoeffT00.y << "\t"
//		<< std::setw(10) << el.CoeffT00.z << "\t" << std::setw(10) << el.CoeffT00.w <<
//		"\t";
//	os << std::setw(10) << el.CoeffT01.x << "\t" << std::setw(10) << el.CoeffT01.y << "\t"
//		<< std::setw(10) << el.CoeffT01.z << "\t" << std::setw(10) << el.CoeffT01.w <<
//		"\t";
//	os << std::setw(10) << el.CoeffT10.x << "\t" << std::setw(10) << el.CoeffT10.y << "\t"
//		<< std::setw(10) << el.CoeffT10.z << "\t" << std::setw(10) << el.CoeffT10.w <<
//		"\t";
//	os << std::setw(10) << el.CoeffT11.x << "\t" << std::setw(10) << el.CoeffT11.y << "\t"
//		<< std::setw(10) << el.CoeffT11.z << "\t" << std::setw(10) << el.CoeffT11.w <<
//		"\t";
//	os << std::setw(10) << el.CoeffU00.x << "\t" << std::setw(10) << el.CoeffU00.y << "\t"
//		<< std::setw(10) << el.CoeffU00.z << "\t" << std::setw(10) << el.CoeffU00.w <<
//		"\t";
//	os << std::setw(10) << el.CoeffU01.x << "\t" << std::setw(10) << el.CoeffU01.y << "\t"
//		<< std::setw(10) << el.CoeffU01.z << "\t" << std::setw(10) << el.CoeffU01.w <<
//		"\t";
//	os << std::setw(10) << el.CoeffU10.x << "\t" << std::setw(10) << el.CoeffU10.y << "\t"
//		<< std::setw(10) << el.CoeffU10.z << "\t" << std::setw(10) << el.CoeffU10.w <<
//		"\t";
//	os << std::setw(10) << el.CoeffU11.x << "\t" << std::setw(10) << el.CoeffU11.y << "\t"
//		<< std::setw(10) << el.CoeffU11.z << "\t" << std::setw(10) << el.CoeffU11.w <<
//		"\t";
//#else
//	os << el.CoeffT00.x << "\t" << el.CoeffT00.y << "\t" << el.CoeffT00.z <<
//		"\t" << el.CoeffT00.w << "\t" << el.CoeffT01.x << "\t" << el.CoeffT01.y <<
//		"\t" << el.CoeffT01.z << "\t" << el.CoeffT01.w << "\t" << el.CoeffT10.x <<
//		"\t" << el.CoeffT10.y << "\t" << el.CoeffT10.z << "\t" << el.CoeffT10.w <<
//		"\t" << el.CoeffT11.x << "\t" << el.CoeffT11.y << "\t" << el.CoeffT11.z <<
//		"\t" << el.CoeffT11.w << "\t" << el.CoeffU00.x << "\t" << el.CoeffU00.y <<
//		"\t" << el.CoeffU00.z << "\t" << el.CoeffU00.w << "\t" << el.CoeffU01.x <<
//		"\t" << el.CoeffU01.y << "\t" << el.CoeffU01.z << "\t" << el.CoeffU01.w <<
//		"\t" << el.CoeffU10.x << "\t" << el.CoeffU10.y << "\t" << el.CoeffU10.z <<
//		"\t" << el.CoeffU10.w << "\t" << el.CoeffU11.x << "\t" << el.CoeffU11.y <<
//		"\t" << el.CoeffU11.z << "\t" << el.CoeffU11.w << "\t";
//#endif
//	os << el.neigh.x << "\t" << el.neigh.y << "\t" << el.neigh.z << "\t" <<
//		el.neigh.w;
//	return os;
//}
//
//
//
//inline std::ostream& operator<<(std::ostream& os, const Nact& el)
//{
//	os << el.ii00.x << "\t" << el.ii00.y << "\t" << el.ii10.x << "\t" <<
//		el.ii10.y << "\t" << el.ii01.x << "\t" << el.ii01.y << "\t" <<
//		el.ii11.x << "\t" << el.ii11.y;
//	return os;
//}







inline std::ostream& operator<<(std::ostream& os, const cl_int2& el)
{
	os << el.x << "\t" << el.y;
	return os;
}
//std::ostream& operator<<(std::ostream& os, const cl_int3& el)
//{
//	os << el.x << "\t" << el.y << "\t" << el.z;
//	return os;
//}
inline std::ostream& operator<<(std::ostream& os, const cl_int4& el)
{
	os << el.x << "\t" << el.y << "\t" << el.z << "\t" << el.w;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_int8& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_int16& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const cl_uint2& el)
{
	os << el.x << "\t" << el.y;
	return os;
}
//std::ostream& operator<<(std::ostream& os, const cl_uint3& el)
//{
//	os << el.x << "\t" << el.y << "\t" << el.z;
//	return os;
//}
inline std::ostream& operator<<(std::ostream& os, const cl_uint4& el)
{
	os << el.x << "\t" << el.y << "\t" << el.z << "\t" << el.w;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_uint8& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_uint16& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const cl_ushort2& el)
{
	os << el.x << "\t" << el.y;
	return os;
}
//std::ostream& operator<<(std::ostream& os, const cl_uint3& el)
//{
//	os << el.x << "\t" << el.y << "\t" << el.z;
//	return os;
//}
inline std::ostream& operator<<(std::ostream& os, const cl_ushort4& el)
{
	os << el.x << "\t" << el.y << "\t" << el.z << "\t" << el.w;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_ushort8& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_ushort16& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const cl_double2& el)
{
	os << el.x << "\t" << el.y;
	return os;
}
//std::ostream& operator<<(std::ostream& os, const cl_double3& el)
//{
//	os << el.x << "\t" << el.y << "\t" << el.z;
//	return os;
//}
inline std::ostream& operator<<(std::ostream& os, const cl_double4& el)
{
	os << el.x << "\t" << el.y << "\t" << el.z << "\t" << el.w;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_double8& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_double16& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const cl_float2& el)
{
	os << el.x << "\t" << el.y;
	return os;
}
//std::ostream& operator<<(std::ostream& os, const cl_float3& el)
//{
//	os << el.x << "\t" << el.y << "\t" << el.z;
//	return os;
//}
inline std::ostream& operator<<(std::ostream& os, const cl_float4& el)
{
	os << el.x << "\t" << el.y << "\t" << el.z << "\t" << el.w;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_float8& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}
inline std::ostream& operator<<(std::ostream& os, const cl_float16& el)
{
	os << el.lo << "\t" << el.hi;
	return os;
}




inline std::istream& operator>>(std::istream& is, const cl_int2& el)
{
	int *xval_ = (int*)&el.x;
	int *yval_ = (int*)&el.y;
	is >> *xval_;
	is >> *yval_;
	return is;
}
//std::istream& operator>>(std::istream& is, const cl_int3& el)
//{
//	int *xval_ = (int*)&el.x;
//	int *yval_ = (int*)&el.y;
//	int *zval_ = (int*)&el.z;
//	is >> *xval_;
//	is >> *yval_;
//	is >> *zval_;
//	return is;
//}
inline std::istream& operator>>(std::istream& is, const cl_int4& el)
{
	int *xval_ = (int*)&el.x;
	int *yval_ = (int*)&el.y;
	int *zval_ = (int*)&el.z;
	int *wval_ = (int*)&el.w;
	is >> *xval_;
	is >> *yval_;
	is >> *zval_;
	is >> *wval_;
	return is;
}
inline std::istream& operator>>(std::istream& is, const cl_int8& el)
{
	int *s1val_ = (int*)&el.s0;
	int *s2val_ = (int*)&el.s1;
	int *s3val_ = (int*)&el.s2;
	int *s4val_ = (int*)&el.s3;
	int *s5val_ = (int*)&el.s4;
	int *s6val_ = (int*)&el.s5;
	int *s7val_ = (int*)&el.s6;
	int *s8val_ = (int*)&el.s7;
	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	return is;
}

inline std::istream& operator>>(std::istream& is, const cl_int16& el)
{
	int *s1val_ = (int*)&el.s0;
	int *s2val_ = (int*)&el.s1;
	int *s3val_ = (int*)&el.s2;
	int *s4val_ = (int*)&el.s3;
	int *s5val_ = (int*)&el.s4;
	int *s6val_ = (int*)&el.s5;
	int *s7val_ = (int*)&el.s6;
	int *s8val_ = (int*)&el.s7;
	int *s9val_ = (int*)&el.s8;
	int *s10val_ = (int*)&el.s9;
	int *s11val_ = (int*)&el.sa;
	int *s12val_ = (int*)&el.sb;
	int *s13val_ = (int*)&el.sc;
	int *s14val_ = (int*)&el.sd;
	int *s15val_ = (int*)&el.se;
	int *s16val_ = (int*)&el.sf;

	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	is >> *s9val_;
	is >> *s10val_;
	is >> *s11val_;
	is >> *s12val_;
	is >> *s13val_;
	is >> *s14val_;
	is >> *s15val_;
	is >> *s16val_;
	return is;
}






inline std::istream& operator>>(std::istream& is, const cl_uint2& el)
{
	unsigned int *xval_ = (unsigned int*)&el.x;
	unsigned int *yval_ = (unsigned int*)&el.y;
	is >> *xval_;
	is >> *yval_;
	return is;
}
//std::istream& operator>>(std::istream& is, const cl_uint3& el)
//{
//	unsigned int *xval_ = (unsigned int*)&el.x;
//	unsigned int *yval_ = (unsigned int*)&el.y;
//	unsigned int *zval_ = (unsigned int*)&el.z;
//	is >> *xval_;
//	is >> *yval_;
//	is >> *zval_;
//	return is;
//}
inline std::istream& operator>>(std::istream& is, const cl_uint4& el)
{
	unsigned int *xval_ = (unsigned int*)&el.x;
	unsigned int *yval_ = (unsigned int*)&el.y;
	unsigned int *zval_ = (unsigned int*)&el.z;
	unsigned int *wval_ = (unsigned int*)&el.w;
	is >> *xval_;
	is >> *yval_;
	is >> *zval_;
	is >> *wval_;
	return is;
}
inline std::istream& operator>>(std::istream& is, const cl_uint8& el)
{
	unsigned int *s1val_ = (unsigned int*)&el.s0;
	unsigned int *s2val_ = (unsigned int*)&el.s1;
	unsigned int *s3val_ = (unsigned int*)&el.s2;
	unsigned int *s4val_ = (unsigned int*)&el.s3;
	unsigned int *s5val_ = (unsigned int*)&el.s4;
	unsigned int *s6val_ = (unsigned int*)&el.s5;
	unsigned int *s7val_ = (unsigned int*)&el.s6;
	unsigned int *s8val_ = (unsigned int*)&el.s7;
	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	return is;
}

inline std::istream& operator>>(std::istream& is, const cl_uint16& el)
{
	unsigned int *s1val_ = (unsigned int*)&el.s0;
	unsigned int *s2val_ = (unsigned int*)&el.s1;
	unsigned int *s3val_ = (unsigned int*)&el.s2;
	unsigned int *s4val_ = (unsigned int*)&el.s3;
	unsigned int *s5val_ = (unsigned int*)&el.s4;
	unsigned int *s6val_ = (unsigned int*)&el.s5;
	unsigned int *s7val_ = (unsigned int*)&el.s6;
	unsigned int *s8val_ = (unsigned int*)&el.s7;
	unsigned int *s9val_ = (unsigned int*)&el.s8;
	unsigned int *s10val_ = (unsigned int*)&el.s9;
	unsigned int *s11val_ = (unsigned int*)&el.sa;
	unsigned int *s12val_ = (unsigned int*)&el.sb;
	unsigned int *s13val_ = (unsigned int*)&el.sc;
	unsigned int *s14val_ = (unsigned int*)&el.sd;
	unsigned int *s15val_ = (unsigned int*)&el.se;
	unsigned int *s16val_ = (unsigned int*)&el.sf;

	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	is >> *s9val_;
	is >> *s10val_;
	is >> *s11val_;
	is >> *s12val_;
	is >> *s13val_;
	is >> *s14val_;
	is >> *s15val_;
	is >> *s16val_;
	return is;
}

inline std::istream& operator>>(std::istream& is, const cl_ushort2& el)
{
	unsigned int* xval_ = (unsigned int*)& el.x;
	unsigned int* yval_ = (unsigned int*)& el.y;
	is >> *xval_;
	is >> *yval_;
	return is;
}
//std::istream& operator>>(std::istream& is, const cl_uint3& el)
//{
//	unsigned int *xval_ = (unsigned int*)&el.x;
//	unsigned int *yval_ = (unsigned int*)&el.y;
//	unsigned int *zval_ = (unsigned int*)&el.z;
//	is >> *xval_;
//	is >> *yval_;
//	is >> *zval_;
//	return is;
//}
inline std::istream& operator>>(std::istream& is, const cl_ushort4& el)
{
	unsigned int* xval_ = (unsigned int*)& el.x;
	unsigned int* yval_ = (unsigned int*)& el.y;
	unsigned int* zval_ = (unsigned int*)& el.z;
	unsigned int* wval_ = (unsigned int*)& el.w;
	is >> *xval_;
	is >> *yval_;
	is >> *zval_;
	is >> *wval_;
	return is;
}
inline std::istream& operator>>(std::istream& is, const cl_ushort8& el)
{
	unsigned int* s1val_ = (unsigned int*)& el.s0;
	unsigned int* s2val_ = (unsigned int*)& el.s1;
	unsigned int* s3val_ = (unsigned int*)& el.s2;
	unsigned int* s4val_ = (unsigned int*)& el.s3;
	unsigned int* s5val_ = (unsigned int*)& el.s4;
	unsigned int* s6val_ = (unsigned int*)& el.s5;
	unsigned int* s7val_ = (unsigned int*)& el.s6;
	unsigned int* s8val_ = (unsigned int*)& el.s7;
	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	return is;
}

inline std::istream& operator>>(std::istream& is, const cl_ushort16& el)
{
	unsigned int* s1val_ = (unsigned int*)& el.s0;
	unsigned int* s2val_ = (unsigned int*)& el.s1;
	unsigned int* s3val_ = (unsigned int*)& el.s2;
	unsigned int* s4val_ = (unsigned int*)& el.s3;
	unsigned int* s5val_ = (unsigned int*)& el.s4;
	unsigned int* s6val_ = (unsigned int*)& el.s5;
	unsigned int* s7val_ = (unsigned int*)& el.s6;
	unsigned int* s8val_ = (unsigned int*)& el.s7;
	unsigned int* s9val_ = (unsigned int*)& el.s8;
	unsigned int* s10val_ = (unsigned int*)& el.s9;
	unsigned int* s11val_ = (unsigned int*)& el.sa;
	unsigned int* s12val_ = (unsigned int*)& el.sb;
	unsigned int* s13val_ = (unsigned int*)& el.sc;
	unsigned int* s14val_ = (unsigned int*)& el.sd;
	unsigned int* s15val_ = (unsigned int*)& el.se;
	unsigned int* s16val_ = (unsigned int*)& el.sf;

	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	is >> *s9val_;
	is >> *s10val_;
	is >> *s11val_;
	is >> *s12val_;
	is >> *s13val_;
	is >> *s14val_;
	is >> *s15val_;
	is >> *s16val_;
	return is;
}




inline std::istream& operator>>(std::istream& is, const cl_double2& el)
{
	double *xval_ = (double*)&el.x;
	double *yval_ = (double*)&el.y;
	is >> *xval_;
	is >> *yval_;
	return is;
}
//std::istream& operator>>(std::istream& is, const cl_double3& el)
//{
//	double *xval_ = (double*)&el.x;
//	double *yval_ = (double*)&el.y;
//	double *zval_ = (double*)&el.z;
//	is >> *xval_;
//	is >> *yval_;
//	is >> *zval_;
//	return is;
//}
inline std::istream& operator>>(std::istream& is, const cl_double4& el)
{
	double *xval_ = (double*)&el.x;
	double *yval_ = (double*)&el.y;
	double *zval_ = (double*)&el.z;
	double *wval_ = (double*)&el.w;
	is >> *xval_;
	is >> *yval_;
	is >> *zval_;
	is >> *wval_;
	return is;
}
inline std::istream& operator>>(std::istream& is, const cl_double8& el)
{
	double *s1val_ = (double*)&el.s0;
	double *s2val_ = (double*)&el.s1;
	double *s3val_ = (double*)&el.s2;
	double *s4val_ = (double*)&el.s3;
	double *s5val_ = (double*)&el.s4;
	double *s6val_ = (double*)&el.s5;
	double *s7val_ = (double*)&el.s6;
	double *s8val_ = (double*)&el.s7;
	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	return is;
}

inline std::istream& operator>>(std::istream& is, const cl_double16& el)
{
	double *s1val_ = (double*)&el.s0;
	double *s2val_ = (double*)&el.s1;
	double *s3val_ = (double*)&el.s2;
	double *s4val_ = (double*)&el.s3;
	double *s5val_ = (double*)&el.s4;
	double *s6val_ = (double*)&el.s5;
	double *s7val_ = (double*)&el.s6;
	double *s8val_ = (double*)&el.s7;
	double *s9val_ = (double*)&el.s8;
	double *s10val_ = (double*)&el.s9;
	double *s11val_ = (double*)&el.sa;
	double *s12val_ = (double*)&el.sb;
	double *s13val_ = (double*)&el.sc;
	double *s14val_ = (double*)&el.sd;
	double *s15val_ = (double*)&el.se;
	double *s16val_ = (double*)&el.sf;

	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	is >> *s9val_;
	is >> *s10val_;
	is >> *s11val_;
	is >> *s12val_;
	is >> *s13val_;
	is >> *s14val_;
	is >> *s15val_;
	is >> *s16val_;
	return is;
}





inline std::istream& operator>>(std::istream& is, const cl_float2& el)
{
	float *xval_ = (float*)&el.x;
	float *yval_ = (float*)&el.y;
	is >> *xval_;
	is >> *yval_;
	return is;
}
//std::istream& operator>>(std::istream& is, const cl_float3& el)
//{
//	float *xval_ = (float*)&el.x;
//	float *yval_ = (float*)&el.y;
//	float *zval_ = (float*)&el.z;
//	is >> *xval_;
//	is >> *yval_;
//	is >> *zval_;
//	return is;
//}
inline std::istream& operator>>(std::istream& is, const cl_float4& el)
{
	float *xval_ = (float*)&el.x;
	float *yval_ = (float*)&el.y;
	float *zval_ = (float*)&el.z;
	float *wval_ = (float*)&el.w;
	is >> *xval_;
	is >> *yval_;
	is >> *zval_;
	is >> *wval_;
	return is;
}
inline std::istream& operator>>(std::istream& is, const cl_float8& el)
{
	float *s1val_ = (float*)&el.s0;
	float *s2val_ = (float*)&el.s1;
	float *s3val_ = (float*)&el.s2;
	float *s4val_ = (float*)&el.s3;
	float *s5val_ = (float*)&el.s4;
	float *s6val_ = (float*)&el.s5;
	float *s7val_ = (float*)&el.s6;
	float *s8val_ = (float*)&el.s7;
	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	return is;
}

inline std::istream& operator>>(std::istream& is, const cl_float16& el)
{
	float *s1val_ = (float*)&el.s0;
	float *s2val_ = (float*)&el.s1;
	float *s3val_ = (float*)&el.s2;
	float *s4val_ = (float*)&el.s3;
	float *s5val_ = (float*)&el.s4;
	float *s6val_ = (float*)&el.s5;
	float *s7val_ = (float*)&el.s6;
	float *s8val_ = (float*)&el.s7;
	float *s9val_ = (float*)&el.s8;
	float *s10val_ = (float*)&el.s9;
	float *s11val_ = (float*)&el.sa;
	float *s12val_ = (float*)&el.sb;
	float *s13val_ = (float*)&el.sc;
	float *s14val_ = (float*)&el.sd;
	float *s15val_ = (float*)&el.se;
	float *s16val_ = (float*)&el.sf;

	is >> *s1val_;
	is >> *s2val_;
	is >> *s3val_;
	is >> *s4val_;
	is >> *s5val_;
	is >> *s6val_;
	is >> *s7val_;
	is >> *s8val_;
	is >> *s9val_;
	is >> *s10val_;
	is >> *s11val_;
	is >> *s12val_;
	is >> *s13val_;
	is >> *s14val_;
	is >> *s15val_;
	is >> *s16val_;
	return is;
}



#define ARRAY_ERROR_CHECKING
#define STATUS_CHECK_SETUP(codeval, msg)	if(codeval) { Gen_Error_Msg(codeval, msg); }


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

class sourceGenerator;
class BiCGStabGenerator;
class ReduceGenerator;
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


#include "HelperFuncs.h"
#include "Kernels.h"
#include "Array.h"
#include "oclEnvironment.h"

