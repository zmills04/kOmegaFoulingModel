// shearStress.h: Class containing variables and methods used to 
// to calculate the shear stress in the channel, and remove 
// deposited particles based on the magnitude of the local shear.
// This will be a subclass in clVariablesTR.
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SHEARSTRESS_H__INCLUDED_)
#define AFX_SHEARSTRESS_H__INCLUDED_

#include "StdAfx.h"
#include "Kernels.h"
#include "Array.h"
#include "particleStructs.h"


#pragma once



// TODO: ordering of vtr.BL is [vls.BL(0) -> vls.BL(vls.nBL/2), vls.BL(nBL-1) -> vls.BL(vls.nBL/2+1)]
//		make sure that methods used here and in trKernelsShear reflect that.

class shearStress
{
public:
	

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////    CONSTRUCTOR/DESTRUCTOR/ENUMS/FUNC_PTRS     ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	shearStress() : sInds("ssSind"), blInds("ssblInds"),
		ssWeights("ssWeights"),	shearCoeffs("ssShearCoeffs"), 
		blIndsLoc("ssblIndsLoc"), Tau("ssTau"), ssOutput("ShearStress")
	{}

	~shearStress()
	{}

	enum saveShearLocs { bothWalls, topWall, botWall };

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////                                               ///////////////
//////////////           KERNELS/REDUCTIONS/SOLVERS          ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

	// Shear calculation kernels
	DualKernel nodeShearKernel;// Calculates shear at wall boundary nodes
	Kernel wallShearKernel;	// interpolates to wall boundary locations

	// Shear removal of deposited particles
	// variable at option index is cl_int2 {{ offset (index of first BL 
	// deposited particle), max_el (total number of BL deposited particles) }}
	Kernel trShearRemovalKernel;

	// Updates Shear Stress coefficients (called in vfl)
	// two step process
	Kernel updateSSKernel[2]; // option index for [1] is number of wall
							  // boundary nodes

	// Prepares shear values to be written to file
	Kernel saveShearKernel;



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

	//Shear stresses calculated at boundary nodes
	Array1Dd Tau;

////////////////////////////////////////////////////////////////////////////	
//////////////                   Method Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////

	// Index of first node used in averaging to get shear at wall
	Array1Dv4i sInds;

	// Index of BL array each index in SS arrays corresponds to	
	Array1Di blInds;

	// Weights for averaging shear when interpolating from SS at nodes to SS at wall
	Array1Dv4d ssWeights;

	// Used to update SS coefficients
	Array1Dv2i blIndsLoc;

	// locations of nodes where shear is calculated 
	// (absolute location i.e. i+FullSizeX*j)
	// Array1Di shearInds; // use vls.ssArr;

	// Coefficients used in calc. of shear at nodes
	Array1Dv8d shearCoeffs;


////////////////////////////////////////////////////////////////////////////	
//////////////                   Output Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////	

	//Array to store Shear Stresses for output
	TimeData<double> ssOutput;


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
	int maxOutLinesSS;
	bool calcSSFlag;		// flag specifying if shear should be saved
	saveShearLocs ssLocSave; // specifies regions to save (top wall, bot wall or both)

	// Total number of nodes at which shear stress is calculated (fluctuates
	// as surface fouls). Based on size of vls.ssArr
	int shearNodeSize;		// actual size of arrays
	int shearNodeDynSize;	// full size of arrays (sizes based on vls.ssArr)

	// Wall BLs shear calculations start and stop and total number of BLs to
	// calculate (for top and bottom)
	int shearBLStartBot, shearBLStartTop;
	int shearBLStopBot, shearBLStopTop;
	int shearBLSizeTop, shearBLSizeBot, shearBLSizeTotal;

////////////////////////////////////////////////////////////////////////////	
//////////////                Method Variables               ///////////////
////////////////////////////////////////////////////////////////////////////	
	
	// number of shear values saved in ssOutput each save step
	int ssOutputSize;	


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
	
	// Initializes shear coefficients used to calculate shear at boundary nodes
	// and, following this, at boundary links
	void iniShearCoeffs();


	////////////////////////////////////////////////////////////////////////////	
	//////////////              Updating Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////
	
	// Updates arrays used in the calculation of shear stress at boundary nodes
	// and, following this, at boundary links
	void updateShearArrays(bool reSizeFlag);

	// This function is called after sort, and it will update the kernel
	// arguments used in trShearRemovalKernel until the next sort.
	void updateParRemArgs();


	////////////////////////////////////////////////////////////////////////////	
	//////////////               Solving Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////
	
	// Calls kernels to calculate shear at boundary links
	void calculateShear(cl_command_queue *que_ = nullptr,
		cl_event* waitevt = nullptr, int numwait = 0, cl_event *evt_ = nullptr);

	// Tests particles to see if shear is sufficient for removal
	void shearRemoval(cl_command_queue *que_ = nullptr,
		cl_event* waitevt = nullptr, int numwait = 0, cl_event *evt_ = nullptr);


	////////////////////////////////////////////////////////////////////////////	
	//////////////                Output Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////	
	//////////////                Display Functions              ///////////////
	////////////////////////////////////////////////////////////////////////////



};



#endif