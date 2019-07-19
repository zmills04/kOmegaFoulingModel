

#if !defined(AFX_HELPERFUNCS_H__INCLUDED_)
#define AFX_HELPERFUNCS_H__INCLUDED_

#pragma once


#include "StdAfx.h"

enum FuncPtrType { ptr2CreateKernels, ptr2SetKernelArgs, ptr2LoadParams };

template <class T>
void functionPointerWrapper(void* pt2Object, FuncPtrType fptype_);

void delay_func(double sec);

bool fileopen(std::fstream &stream, std::string &FileName, FileMode mode_, const int att);
bool fileopen(std::fstream &stream, std::string &FileName, FileMode Mode);

void MakeDir(std::string& NewDir);

void CopyFile(std::string NameSrc, std::string NameDest);

void RenameFile(std::string SourceFileName, std::string DestFileName);

double boxMuller(double meanval, double stdval);





void operator >> (const YAML::Node& node, cl_uint4 vec_);
YAML::Emitter& operator << (YAML::Emitter& out, const cl_uint4 vec_);














#endif