// Class containing data/methods for kOmega turbulence model.
// clVariablesLB will contain an instance of this class rather than
// there being a separate global instance in lbls.cpp.
//
// (c) Zachary Grant Mills 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_KOMEGA_H__INCLUDED_)
#define AFX_KOMEGA_H__INCLUDED_

#pragma once

#include "SparseMatrix.h"
#include "Reducer.h"
#include "BiCGStabSolver.h"



class kOmega
{
public:

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////    CONSTRUCTOR/DESTRUCTOR/ENUMS/FUNC_PTRS     ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

	kOmega() : WallD("wallD"),
		Diff_Omega("diffOmega"), Diff_K("diffK"), Fval_array("Fvals"),
		Nut_array("lbnut"), dKdO_array("dKdO"), Sxy_array("lbsxy"),
		Kappa_array("lbkappa"), Omega_array("lbomega"),
		KappaInds(M_SOLID_NODE, M_FLUID_NODE, SOLID_BOUNDARY_NODE),
		OmegaInds(M_SOLID_NODE, M_FLUID_NODE, SOLID_BOUNDARY_NODE)
	{}

	~kOmega() {}




	// Used to specify the type of debugging being used for kOmega functions
	enum { lbDbgSave, koDbgSave1, koDbgSave2, koDbgSave, DbgSave };



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////            KERNELS/SOLVERS/REDUCERS           ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


	Kernel kOmegaUpdateDiffCoeffs;	// Step 1 updating coefficients for kOmega Eqns
	Kernel kOmegaUpdateCoeffs;		// Step 2 updating coefficients for kOmega Eqns
	Kernel updateWallDKernel;		// updates wall distance array

	// Reduction kernels
	Reducer<double, double, ReduceGenerator::Min> minKappa, minOmega;
	Reducer<double, double> sumKappa, sumOmega;

	// Solvers for Kappa and Omega
	BiCGStabSolver Kappa;
	BiCGStabSolver Omega;


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

	// Arrays that contain actual kappa and omega values, which
	// directly addressed during initialization/reading bin files
	// pointers to these arrays are passed to Solvers, which
	// allows for arrays to be accessed from those classes
	Array2Dd Kappa_array, Omega_array; 
	
	Array2Dd Nut_array; // turbulent viscosity 			



////////////////////////////////////////////////////////////////////////////	
//////////////                 Method Arrays                 ///////////////
////////////////////////////////////////////////////////////////////////////

	// CSR indicies for kappa and omega, which are passed to 
	// solvers
	CSR_Inds_Periodic KappaInds, OmegaInds;
	
	// Arrays for kOmega solver coefficients,
	// which are used to calculate the soln matrix
	Array2Dd WallD;	// wall distance of each node
	Array2Dd Diff_Omega, Diff_K; // Diffusion terms
	Array2Dd Fval_array; // F term 
	Array3Dd dKdO_array; // derivatives of kappa and omega
	Array3Dd Sxy_array; // strain rate Tensor terms xx,xy and yy

	// Used for debugging kOmega
#ifdef DEBUG_TURBARR
	Array3Dd lbDbgArr, koDbgArr1, koDbgArr2;
#endif

////////////////////////////////////////////////////////////////////////////	
//////////////                 Output Arrays                 ///////////////
////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////	
//////////////                Display Arrays                 ///////////////
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
	// Switches for initialization and to turn on turb. model solver
	bool perturbVelField;
	bool kOmegaSolverFlag, iniTurbVelocity;

	// Flags used to indicate if simulation is being restarted from
	// a checkpoint
	bool kOmegaLoadedFlag;

	// Parameters for perturbing initial velocity distribution
	double perturbDUPlus, perturbEpsilon;
	double perturbBeta, perturbAlpha, perturbSigma;

	// Parameters used by kOmega Solver
	double roughnessFactor, UtauVal, ReTurbVal;

	// Initial Values to fill k and Omega Arrays with
	double kIniVal, omegaIniVal;
 


	// Solver parameters for Kappa and Omega eqns
	double kappaMaxRelTol, kappaMaxAbsTol, omegaMaxRelTol, omegaMaxAbsTol;
	int kappaMaxIters, omegaMaxIters;


	// Search radius for BLs when calculating wall distances
	int wallDSearchRar;

	// Characteristic Values used to initialize Kappa/Omega
	double TurbIntensity, TurbLScale_inv;

	// Wall function values
	double kappaWallVal, omegaWallVal, nutWallVal, yPlusWallVal;


////////////////////////////////////////////////////////////////////////////	
//////////////                Method Variables               ///////////////
////////////////////////////////////////////////////////////////////////////	
	




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

	// initialization function for time data array
	void iniTimeData();

	// Loads parameters passed in yaml parameter file, (also reads in 
	// restart variables when a run is restarted)
	void loadParams();

	// Copies saved files from main folder into results folder to ensure
	// that next files do not save 
	void renameSaveFiles();

	// Writes output data to file(specific arrays, not all of them)
	void save2file();

	// Writes additional data to file for debugging purposes
	void saveDebug(int saveFl);

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

	// called by FL class to update data
	void update();

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

	// Iterates through BL indicies provided by first two arguments and 
	// returns the minimum distance between nLoc and the wall.
	double findMinDist(cl_int2 botSearchInds, cl_int2 topSearchInds, 
		cl_double2& nLoc);

////////////////////////////////////////////////////////////////////////////	
//////////////            Initialization Functions           ///////////////
////////////////////////////////////////////////////////////////////////////

	// initializes k and omega arrays, solvers and reduce classes
	void iniKOmegaArrays();

	// Initializes wall distances array
	void iniWallD();

	// Calculates initial values for k and omega if not provided in yaml
	// params file, calculates nut from these values and fills arrays
	// with these values
	void setInitialValues();
////////////////////////////////////////////////////////////////////////////	
//////////////              Updating Functions               ///////////////
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////	
//////////////               Solving Functions               ///////////////
////////////////////////////////////////////////////////////////////////////

	void Solve();

////////////////////////////////////////////////////////////////////////////	
//////////////                Output Functions               ///////////////
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////	
//////////////                Display Functions              ///////////////
////////////////////////////////////////////////////////////////////////////
};




#endif //AFX_KOMEGA_H__INCLUDED_
