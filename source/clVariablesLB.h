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



class clVariablesLB
{
public:
	// Func Pointer for calling loadParams
	std::function<void(void)> loadParamsPtr;

	clVariablesLB() : FA("FA"), FB("FB"), NodeType("nodeType"),
		Ro_array("lbro"), Ux_array("lbux"), Uy_array("lbuy"),
		Inlet_Vel("inletVels"), dXCoeffs("dXCoeffs"), WallD("wallD"),
		Diff_Omega("diffOmega"), Diff_K("diffK"), Fval_array("Fvals"),
		Nut_array("lbnut"), dKdO_array("dKdO"), Sxy_array("lbsxy"),
		Kappa_array("lbkappa"), Omega_array("lbomega")
	{
		loadParamsPtr = std::bind(&clVariablesLB::loadParams, this);


	};
	virtual ~clVariablesLB()
	{ 
	};

	enum dxCoeffInd { dxind_e, dxind_w, dxind_c, dyind_n, dyind_s, dyind_c, dx2ind_e, dx2ind_w, dx2ind_c, dy2ind_n, dy2ind_s, dy2ind_c, };
	enum { lbDbgSave, koDbgSave1, koDbgSave2, koDbgSave, DbgSave };
////////////////////////////////////////////////////////////////////////////	
//////////////           Kernels/Reductions/Solvers          ///////////////
////////////////////////////////////////////////////////////////////////////	
	DualKernel Collision_kernel;	// Kernel performs collision and propogation
	DualKernel IBB_kernel_Fluid;	// Kernel applies IBB BC's
	DualKernel LB_Outflow;			// Applies Outflow BC, if inflow/outflow is used
	Kernel kOmegaUpdateDiffCoeffs;	// Step 1 updating coefficients for kOmega Eqns
	Kernel kOmegaUpdateCoeffs;		// Step 2 updating coefficients for kOmega Eqns
	Kernel calcRoJKernel;

	Reducer<double, ReduceGenerator::Min> minKappa, minOmega;
	Reducer<double> sumUx, sumUy, sumRo, sumKappa, sumOmega;

	BiCGStabSolver Kappa;
	BiCGStabSolver Omega;
	

////////////////////////////////////////////////////////////////////////////	
//////////////                   Arrays                      ///////////////
////////////////////////////////////////////////////////////////////////////	
	Array3Dd FA, FB; // LB distributions
	Array2Di NodeType; // indicator for node type (solid, fluid, fouling layer)
					   // and whether or not a distribution will bounce back from
					   // a given direction

	// IBB stuff for lbm not currently implemented
	Array1Dv2i IBB_loc;		// locations used in IBB step
	Array1Dv2d IBB_coeff;	// Coefficients used in IBB step
	
	// Macroscopic Variables for fluid
	Array2Dd Ro_array;	// Density array
	Array2Dd Ux_array, Uy_array; // Velocity array
	

	Array1Dd Inlet_Vel;	// Parabolic velocity values at each node at inlet
	
	CSR_Inds_Periodic KappaInds, OmegaInds;
	Array2Dd Kappa_array, Omega_array; // arrays used to read in kOmega values
										// from bin files and or intialize arrays

	// Arrays for kOmega solver coefficients,
	// which are used to calculate the soln matrix
	Array3Dd dXCoeffs;
	Array2Dd WallD, Diff_Omega, Diff_K, Fval_array;
	Array2Dd Nut_array;
	Array3Dd dKdO_array;
	Array3Dd Sxy_array;

	// Used for debugging kOmega
#ifdef DEBUG_TURBARR
	Array3Dd lbDbgArr, koDbgArr1, koDbgArr2;
#endif



////////////////////////////////////////////////////////////////////////////	
//////////////                   Variables                   ///////////////
////////////////////////////////////////////////////////////////////////////	

	// Characteristic Values used to initialize Kappa/Omega
	double TurbIntensity, TurbLScale_inv;
	
	// Wall function values
	double kappaWallVal, omegaWallVal, nutWallVal, yPlusWallVal;
	
	int alter;			// used to alternate between two collision kernels
	int Save_loc;		// Tracks location to save RoJ_out data at each time step
	int IBB_array_len;	// length of IBB array

////////////////////////////////////////////////////////////////////////////	
//////////////               Run Parameters                  ///////////////
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
	bool perturbVelField, iniPoiseuille, runLBFDFirst;
	bool kOmegaSolverFlag, saveMacroStart, iniTurbVelocity;

	// Flags used to indicate if simulation is being restarted from
	// a checkpoint
	bool fLoadedFlag, kOmegaLoadedFlag, restartRunFlag;

	// Parameters for perturbing initial velocity distribution
	double perturbDUPlus, perturbEpsilon, perturbBeta, perturbAlpha, perturbSigma;

	// Parameters used by kOmega Solver
	double roughnessFactor, UtauVal, ReTurbVal;

	// Times for pre-particle release LBM/kOmega/Temp solvers
	unsigned int tlbmIniDumpStep, tlbmNumIniSteps_LB, tlbmNumIniSteps_TLBM;

	// Variables used in the iterative flowrate solver
	unsigned int numIntervalsPerAvg, timeBetweenIntervals;
	unsigned int pauseBtwAdj, Next_Update_Time;
	double flowrateMaxPercentDiff, Fval_prev, Fval_orig, Fval_Diff;
	double Um_keep, Um_prev;

	// Solver parameters for Kappa and Omega eqns
	double kappaMaxRelTol, kappaMaxAbsTol, omegaMaxRelTol, omegaMaxAbsTol;
	int kappaMaxIters, omegaMaxIters;
	
////////////////////////////////////////////////////////////////////////////	
//////////////               Member Functions                ///////////////
////////////////////////////////////////////////////////////////////////////	


	//allocates FA,FB and macroscopic variable arrays/buffers
	void allocateArrays(); 

	// calculates coefficients for differentiation of K and Omega  
	void calcDxTerms();

	// Calculates Feq values for a given velocity and density
	void calcFeqVals(double *Feq_vals, double rho_val, cl_double2 u_val);

	// Calculates turbulent viscosity array from values in kappa and omega
	void calcNutArray();

	//Creates kernels and sets global and local sizes
	void createKernels();

	//Frees memory allocated on host that is unnecessary after its written to device
	void freeHostMem();

	// Solves LB/kOmega/Temp equations to get initial distributions
	// before releasing particles
	void getInitialFlow();
	void getIntialFlowAndTemp();

	//initializes LB system
	void ini(); 
	
	//initializes F distributions to Weight values
	void iniDists();   

	//initializes F distributions to parabolic flow
	void iniDists(double Umaxval);  
	
	void iniLBM();

	//creates IBB arrays from arrays in clVariablesLS
	void iniIBB(); 
	
	//Initializes velocities at inlet nodes
	void iniInletVels();

	// Initializes and nodeType Array
	void iniNodeType();

	// Initializes wall distances array
	void iniWallD();

	// initializes kOmega solver variables
	void iniKOmega();

	//Reads parameters from fluid document in YAML input
	void loadParams();

	//Reads output array from device to host and resets offset value
	void resetOutputs() { Save_loc = 0; }
	
	// Saves macroscopic variables
	// Saves F arrays to bins if CREATE_BIN_FILES is defined
	void save2file();
	
	// Saves Checkpoint files (bin files)
	void saveRestartFiles();

	//saves variables necessary for debugging
	void saveDebug(int saveFl = 4);
	
	// Saves each FA/FB distribution as an individual text file for
	// debugging purposes
	void saveDistributions();

	// Saves parameters in restart file
	void saveParams();

	//Sets arguments for opencl kernels
	void setKernelArgs();

	// Adds defines for kOmega Kernels
	void setSourceDefinesKOmega();
	
	// Adds defines for LBM kernels
	void setSourceDefinesLBM();

	//Enqueues LB kernels
	void SolveLB();
	
	//Enqueues kOmega kernels
	void SolveKOmega();

	// Loads checkpoint files, reads restartFlag parameter
	// and sets appropriate flags
	void testRestartRun();
	
	//Updates IBB arrays
	void updateIBBArrays();

	//Calls Update_output_kernels to get data for tout.txt
	void updateOutputs();

	/*static void functionPointerWrapper(void* pt2Object, FuncPtrType fptype_);*/



////////////////////////////////////////////////////////////////////////////	
//////////////                   Constants                   ///////////////
////////////////////////////////////////////////////////////////////////////	
	int nL = 9;

	const static int rev[9];
	const static int per[9];

	const static int Cx[9];
	const static int Cy[9];

	const static cl_int2 Cxy[9];

	const static cl_double2 Cxy_double[9];

	const static double Weigh[9];
	const static double Cs2;
	const static int LF[LB_NUMBER_CONNECTIONS];

	const static int boundsArr[9];
	const static int boundsArrT1[9];
	const static int boundsArrT2[9];

};

#endif // !defined(AFX_CLVARIABLESLB_H__INCLUDED_)