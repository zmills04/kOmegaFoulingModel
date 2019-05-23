// particleDisplay.h: Class containing variables and functions used
// to make openGL calls which display particles on screen.
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PARTICLEDISPLAY_H__INCLUDED_)
#define AFX_PARTICLEDISPLAY_H__INCLUDED_

#pragma once

#include "StdAfx.h"



class particleDisplay
{
public:

	particleDisplay() : P_vbo("pVBO"), parColorList("parColorList"),
		Tcrit_color("TcritColor")
	{}

	~particleDisplay()
	{}

	Kernel TR_GL_kernel;
	
	int numParGL;
	int OpenGL_timer;
	int numParPerGLObj; // number of tracers each gl object represents
	bool glGridFlag;
	double pointSizes; // initial size of particles
	double lineSizes; // initial thickness of lines
	int glScreenWidth; // width of opengl window
	int glScreenHeight; // height of opengl window


	Array1DGL P_vbo;
	Array1Df parColorList;
	Array1Dd Tcrit_color;

	void allocateArrays();
	void allocateBuffers();
	void createKernels();
	void ini();
	void loadParams();
	//void save2file();
	//void saveDebug();
	void saveParams();
	//void saveRestartFiles();
	void setKernelArgs();
	void setSourceDefines();
	//void testRestartRun();


	void iniParticleColors();


};














#endif