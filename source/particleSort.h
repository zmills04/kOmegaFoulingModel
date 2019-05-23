// particleProperties.h: Class storing kernels, parameters and arrays
// used to sort particles.
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////


// TODO: make sure that UVALS_START location doesnt change due to
//			growing fouling layer. If so, need to implement methods
//			to handle if this happens

#if !defined(AFX_PARTICLESORT_H__INCLUDED_)
#define AFX_PARTICLESORT_H__INCLUDED_

#include "StdAfx.h"


#pragma once


class particleSort
{
public:
	particleSort() : Ploc("trPloc"), Umax_val("trUmaxVal"), Ptemp("ParTemp"),
		trP("TrParam")
	{}

	~particleSort()
	{
	}


	// Sorts particles
	// [0] local merge
	// [1] and [2] global merge (alternating between buffers)
	// [3] updates location information
	Kernel sortKernel[3];

	// Re-releases particles at inlet after sort
	Kernel trReReleaseKernel;

	// Calculates Umax for Re-release of particles
	// (for generating the random distribition.
	// I believe that the method is discussed in thesis)
	Kernel getUmaxKernel;
	
	// Called after particles are clumped together.
	// Checks to ensure that clumping process did not shift
	// particle outside of domain, and moves it if it has.
	Kernel clumpParticleKernels;

	int localRange;
	bool oddMergesFlag;
	int sortTimer;
	int numStepsBtwSort;
	cl_int4 Ploc_inds;
	cl_uint numMerges;
	double inletConcDtDivDx;   //inlet concentration in #particles*dt/dx units (multiply by Umean, which is not know a priori to get actual concentration at inlet)
	int avgParPerRelease;
	// lsc.C, and BL info that can be set as defines,
	// and might be useful for other things such as 
	// calculating kernel sizes, etc.
	int BL_rel_bot, BL_rel_top;
	int BL_stop_bot, BL_stop_top;
	int LSC_rel_bot, LSC_rel_top;
	int LSC_stop_bot, LSC_stop_top;


	//tracks the location of particles after sort
	Array1Dv2i Ploc;

	//parameters used to re-release particles
	TrParam trP;

	// Single element array used for storing max inlet velocity
	Array1Dd Umax_val;

	// Temp array for particles used during clumping (host mem) and 
	// during sort (device mem)
	Par Ptemp;

	// Array used for tracking new location of sorted particles
	cl_mem sortLocs1, sortLocs2;

	////////////////////////////
	//    Basic Functions     //
	////////////////////////////
	void allocateArrays();
	void allocateBuffers();
	void createKernels();
	void ini();
	void loadParams();
	void save2file();
	void saveDebug();
	void saveParams();
	void saveRestartFiles();
	void setKernelArgs();
	void setSourceDefines();
	bool testRestartRun();


	void initialSort();
	void iniTrp();
	void operator()(cl_event *TR_prev_evt = nullptr, int numevt = 0);
	void reReleasePar();
	void sortParticlesForClumping();
	void updateTrp();
	void clumpParticles();

};



#endif