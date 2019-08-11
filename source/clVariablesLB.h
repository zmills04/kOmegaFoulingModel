// clVariables.h: interface for the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLVARIABLESLB_H__INCLUDED_)
#define AFX_CLVARIABLESLB_H__INCLUDED_

#pragma once

#include "SparseMatrix.h"
#include "Reducer.h"
#include "BiCGStabSolver.h"
#include "kOmega.h"



class clVariablesLB
{
public:
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////    CONSTRUCTOR/DESTRUCTOR/ENUMS/FUNC_PTRS     ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	
	// Func Pointer for calling loadParams(bool dynamicFlag) 
	std::function<void(void)> loadParamsPtr;

	// kOmega data/functions/kernels
	kOmega kOmegaClass;
	
	clVariablesLB() : FA("FA"), FB("FB"), Ro_array("lbro"),
		Ux_array("lbux"), Uy_array("lbuy"), Inlet_Vel("inletVels"),
		lbOut("lbOut")//, NodeType("nodeType")
	{
		loadParamsPtr = std::bind(&clVariablesLB::loadParams, this);
	};

	virtual ~clVariablesLB()
	{ 
	};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////                                               ///////////////
//////////////           KERNELS/REDUCTIONS/SOLVERS          ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

	DualKernel collisionKernel;	// Kernel performs collision and propogation

#ifndef IN_KERNEL_IBB
	DualKernel ibbKernel;		// Kernel applies IBB BC's
#endif

	DualKernel lbOutflow;		// Applies Outflow BC, if inflow/outflow is used
	Kernel calcRoJKernel;		// calculates rho, Ux, Uy from Distribution
								// used when restarting run

	Reducer<double, double> sumUx, sumUy, sumRo;
	
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
	
	Array3Dd FA, FB; // LB distributions
	
	Array2Dd Ro_array;	// Density array
	Array2Dd Ux_array, Uy_array; // Velocity array

	
////////////////////////////////////////////////////////////////////////////	
//////////////                   Method Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////
	
	//Array2Di NodeType; // indicator for node type (solid, fluid, fouling layer)
					   // and whether or not a distribution will bounce back from
					   // a given direction

	// IBB stuff for lbm not currently implemented
	Array1Dv2i IBB_loc;		// locations used in IBB step
	Array1Dv2d IBB_coeff;	// Coefficients used in IBB step
	
	Array1Dd Inlet_Vel;	// Used by inlet velocity BC (not implemented)
	

////////////////////////////////////////////////////////////////////////////	
//////////////                   Output Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////	
	TimeData<double> lbOut; // stores time data (i.e. umean, total density, etc)



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

	// Flow properties (Ux_inlet is for inlet/outlet BC, Fval for periodic)
	double Re, Ux_inlet, Fval;
	bool reDefined; // false if parameter file contained Fval, else true

	// Fluid properties
	double tau, MuVal, RhoVal;
	bool tauDefined; // false if parameter file contained MuVal, else true

	// Read Fluid Properties (for calculating deltaT, deltaM, etc)
	double nuAir, rhoAir;

	// Switches for initialization and to turn on turb. model solver
	bool iniPoiseuille, runLBFDFirst;
	bool saveMacroStart;

	// Flags used to indicate if simulation is being restarted from
	// a checkpoint
	bool fLoadedFlag, restartRunFlag;

	// Flag indicating whether or not to use turbulent velocity BC's
	// at boundary nodes
	bool turbVelBCFlag;

	// Flag indicating if averages for flow variables should be saved
	bool saveAvgInfo;

	// y-Size of lbOut (max number of save points before needing to output)
	int lbOutSize;

	// Times for pre-particle release LBM/kOmega/Temp solvers
	unsigned int tlbmIniDumpStep, tlbmNumIniSteps_LB, tlbmNumIniSteps_TLBM;

	double geomPenaltyFactor; // mean velocity calculated from set Re
							  // is scaled by this amount, when calculating values
							  // used to initialize flow variables. This is because the
							  // additional pressure drop induced by geometry will cause 
							  // velocities to be lower than those in a straight
							  // channel with same pressure drop (and eqns used to calculate
							  // initial values are based parallel plate geometry).

	// Variables used in the iterative flowrate solver (not implemented)
	// Note: not all of these are read from/written to parameter file
	// at the moment (need to implement this if using flow solver)
	unsigned int numIntervalsPerAvg, timeBetweenIntervals;
	unsigned int pauseBtwAdj, Next_Update_Time;
	double flowrateMaxPercentDiff, Fval_prev, Fval_orig, Fval_Diff;
	double Um_keep, Um_prev;

////////////////////////////////////////////////////////////////////////////	
//////////////                Method Variables               ///////////////
////////////////////////////////////////////////////////////////////////////	
	
	int alter;			// used to alternate between two collision kernels
	int IBB_array_len;	// length of IBB array

	
	// Number of distributions
	int nL = 9;

	// Reverse direction for bounce back
	const static int rev[9];

	// Lattice Velocities
	const static int Cx[9];
	const static int Cy[9];
	const static cl_int2 Cxy[9];
	const static cl_double2 Cxy_double[9];

	// Distribution weights
	const static double Weigh[9];
	
	// Cs^2 (speed of sound squared = 1/3)
	const static double Cs2;

	// Boundary flags
	const static int boundsArr[9];
#ifdef IN_KERNEL_IBB
	const static int boundsArrT1[9];
	const static int boundsArrT2[9];
#endif

	
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

	// Loads parameters passed in yaml parameter files, (also reads in 
	// restart variables when a run is restarted). 
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

////////////////////////////////////////////////////////////////////////////	
//////////////            Initialization Functions           ///////////////
////////////////////////////////////////////////////////////////////////////

	// Calculates Feq values for a given velocity and density
	void calcFeqVals(double* Feq_vals, double rho_val, cl_double2 u_val);

	// Solves LB/kOmega/Temp equations to get initial distributions
	// before releasing particles
	void getInitialFlow();
	void getIntialFlowAndTemp();

	//initializes F distributions to Weight values
	void iniDists();

	//initializes F distributions to parabolic flow
	void iniDists(double Umaxval);

	//creates IBB arrays from arrays in clVariablesLS
	//void iniIBB();

	//Initializes velocities at inlet nodes
	void iniInletVels();

	// Initializes and nodeType Array
	//void iniNodeType();

	// calculates rho and U values for arrays when distributions
	// have been loaded from a bin file.
	void iniRhoUFromDists();

////////////////////////////////////////////////////////////////////////////	
//////////////              Updating Functions               ///////////////
////////////////////////////////////////////////////////////////////////////

	//Updates IBB arrays
	//void updateIBBArrays(bool reSizeFlag);

////////////////////////////////////////////////////////////////////////////	
//////////////              Solving Functions                ///////////////
////////////////////////////////////////////////////////////////////////////

	// Calculates current mean velocity of fluid
	double calcUmean();

	//Enqueues LB kernels
	void Solve();

////////////////////////////////////////////////////////////////////////////	
//////////////                Output Functions               ///////////////
////////////////////////////////////////////////////////////////////////////
	
	// Saves each FA/FB distribution as an individual text file for
	// debugging purposes
	void saveDistributions();


////////////////////////////////////////////////////////////////////////////	
//////////////                Display Functions              ///////////////
////////////////////////////////////////////////////////////////////////////	


};

#endif // !defined(AFX_CLVARIABLESLB_H__INCLUDED_)