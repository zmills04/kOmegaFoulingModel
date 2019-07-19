// particleSort.h: Class storing kernels, parameters and arrays
// used to sort particles, and re-release deposited particles.
// This will be a subclass contained in clVariablesTR.
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////


// TODO: make sure that UVALS_START location doesnt change due to
//			growing fouling layer. If so, need to implement methods
//			to handle if this happens

#if !defined(AFX_PARTICLESORT_H__INCLUDED_)
#define AFX_PARTICLESORT_H__INCLUDED_

#include "StdAfx.h"
#include "particleStructs.h"

#pragma once


class particleSort
{
public:

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////    CONSTRUCTOR/DESTRUCTOR/ENUMS/FUNC_PTRS     ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

	particleSort() : Ploc("trPloc"), Umax_val("trUmaxVal"), Ptemp("ParTemp"),
		trP("TrParam"), clMemIniFlag(false)
	{}

	~particleSort()
	{
		// need to release buffers if they have been initialized
		if (clMemIniFlag)
		{
			clReleaseMemObject(sortLocs1);
			clReleaseMemObject(sortLocs2);
		}
		clMemIniFlag = false;
	}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////                                               ///////////////
//////////////           KERNELS/REDUCTIONS/SOLVERS          ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	
    // Sorting Kernels:
	Kernel localMergeKernel; // Performs initial local merge
	DualKernel globalMergeKernel; // performs global merge 
								  // (alternates between buffers)
	Kernel updateLocationKernel;

	// Re-releases particles at inlet after sort
	Kernel trReReleaseKernel;

	// Calculates Umax for Re-release of particles
	// (for generating the random distribition.
	// I believe that the method is discussed in thesis)
	// For OpenCL version 1.2 the implementation
	// is basically the same as the final reduce step for generic reduce
	// kernels. For OpenCL 2.0, implementation uses built-in max reduce
	// function. Current implementation only allows for 1 work group,
	// so size of channel height must be <= max WG size for opencl 2.0 or
	// twice the max WG size for opencl 1.2 (throws error if not).
	// (currently, max WG size hardcoded as 256).
	Kernel getUmaxKernel;

	// Called after particles are clumped together.
	// Checks to ensure that clumping process did not shift
	// particle outside of domain, and moves it if it has.
	Kernel clumpParticleKernels;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////                   ARRAYS                      ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////	
//////////////                   Data Arrays                 ///////////////
////////////////////////////////////////////////////////////////////////////
		
	//tracks the location of particles after sort
	Array1Dv2i Ploc;

////////////////////////////////////////////////////////////////////////////	
//////////////                   Method Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////

	//parameters used to re-release particles
	TrParam trP;

	// Single element array used for storing max inlet velocity
	// Array is used to simplify reading from device, etc.
	Array1Dd Umax_val;

	// Temp array for particles used during clumping (host mem) and 
	// during sort (device mem)
	Par Ptemp;


////////////////////////////////////////////////////////////////////////////	
//////////////                   Output Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////	




////////////////////////////////////////////////////////////////////////////	
//////////////                   Display Arrays              ///////////////
////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////                 OPENCL BUFFERS                ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

	// Array used for tracking new location of sorted particles
	cl_mem sortLocs1, sortLocs2;


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////                   VARIABLES                   ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////	
//////////////             Run Parameter Variables           ///////////////
////////////////////////////////////////////////////////////////////////////	
	
	// maximum number of particles to realease in a single step. This
	// allows for staggered release at the beginning of simulation,
	// and after a clumping.
	int avgParPerRelease;

	// note: these should be in terms of lb times steps, and not
	// vtr steps, which will not necessarily occur every lb step.
	int sortTimer;		// Time steps until next sort is called		 
	int numStepsBtwSort;// Number of time steps between sort calls

	bool clumpFlag; // flag for whether or not to clump particles
					// during simulation
	
////////////////////////////////////////////////////////////////////////////	
//////////////                Method Variables               ///////////////
////////////////////////////////////////////////////////////////////////////	

	int localRange;			// used by merge sort kernels
	bool oddMergesFlag;		// flag indicating if an odd number of global
							// merges are being used, which results in
							// final sort step to end with data in temp Ptemp
	cl_int4 Ploc_inds;		// Used to quickly set first two elements of Ploc
	cl_uint numMerges;		// number of global merges
	double inletConcDtDivDx;// inlet concentration in (# of particles)*dt/dx units
							// value is multiplied by Umean to get total # of 
							// to be released after a given sort step (where # of
							// particles is the represented particles, not number
							// of tracers)

	// lsc.C, and BL info that can be set as defines,
	// and might be useful for other things such as 
	// calculating kernel sizes, etc.
	int BL_rel_bot, BL_rel_top;
	int BL_stop_bot, BL_stop_top;
	int LSC_rel_bot, LSC_rel_top;
	int LSC_stop_bot, LSC_stop_top;

	bool clMemIniFlag; // flag indicating that cl_mem was initialized


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////                  BASE FUNCTIONS               ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	

	// Allocates host arrays containing data
	void allocateArrays();

	// Allocates device buffers, and copies host contents to them
	void allocateBuffers();

	// Creates kernels from compiled openCL source code. Pointers to 
	// these functions are passed to sourceGenerator to allow for 
	// them to be called after compilation
	void createKernels();

	// Frees arrays on host which are no longer needed to save memory
	void freeHostArrays();

	// Inititialization function
	void ini();

	// Loads parameters passed in yaml parameter file, (also reads in 
	// restart variables when a run is restarted)
	void loadParams();

	// Copies saved files from main folder into results folder to ensure
	// that next files do not save 
	void renameSaveFiles();

	// Writes output data to file(specific arrays, not all of them)
	void save2file();

	// Writes additional data to file for debugging purposes
	void saveDebug();

	// Saves parameters to yaml file used for restarting runs
	void saveParams();

	// Saves bin files necessary to restart run
	void saveRestartFiles();

	// Saves time output data (i.e. avg velocity, shear, Nu, etc)
	void saveTimeData();

	// Sets arguments to kernels.Pointer to this function is passed to 
	// sourceGenerator along with createKernels so that after kernels 
	// are created, all arguments are set.
	// Note: parameters cannot be set in any other initialization functions,
	// because kernel creation is last initialization step.
	void setKernelArgs();

	// Prepends macro definitions specific to the class method to opencl source
	void setSourceDefines();

	// Tests to see if run is new, or if it is a restart.If information is 
	// missing, the run cannot restart.
	// Note: The run will be able to load temp and velocity data while still
	// being a new run for remaining methods
	bool testRestartRun();

	// Calls kernels to save time data (umean, avg density, etc) to array
	// on device, which will eventually be saved if it reaches its max 
	// size, or a save step is reached
	void updateTimeData();





	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////	
	//////////////                                               ///////////////
	//////////////            CLASS SPECIFIC FUNCTIONS           ///////////////
	//////////////                                               ///////////////
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////	
	//////////////            Initialization Functions           ///////////////
	////////////////////////////////////////////////////////////////////////////
	
	// Performs initial sort of particles, and also called after clumping to
	// sort particles.
	void initialSort();

	// Initializes Trp structure, which contains variables used by sorting
	// kernels.
	void iniTrp();


	////////////////////////////////////////////////////////////////////////////	
	//////////////              Updating Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	// Updates data in Trp structure for use in sorting kernels
	void updateTrp();

	////////////////////////////////////////////////////////////////////////////	
	//////////////               Solving Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	// combines tracers in close proximity to one another in order to reduce
	// number of tracers currently in simulation domain, while keeping the number
	// of particles constant. (This is necessary because large numbers of 
	// particles will collect in the recirulation regions and take a long time
	// to deposit.
	void clumpParticles();

	// Sorts particles
	void operator()(cl_event* TR_prev_evt = nullptr, int numevt = 0);
	
	// Releases particles which have been deposited or removed by clumping 
	void reReleasePar();

	// Sorts particles before clumping is performed
	void sortParticlesForClumping();


	////////////////////////////////////////////////////////////////////////////	
	//////////////                Output Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////	
	//////////////                Display Functions              ///////////////
	////////////////////////////////////////////////////////////////////////////

};



#endif