// clVariables.h: interface for the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLVARIABLESFD_H__INCLUDED_)
#define AFX_CLVARIABLESFD_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clProblem.h"
#include "BiCGStabSolver.h"
#include "SparseMatrix.h"
#include "Reducer.h"

class clVariablesFD
{
public:

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////    CONSTRUCTOR/DESTRUCTOR/ENUMS/FUNC_PTRS     ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

	// Func Pointer for calling loadParams
	std::function<void(void)> loadParamsPtr;

	clVariablesFD() : tempArray("temp"), Alphat("alphat"), 
		alphaArray("alpha"), Nu("Nu")
	{
		loadParamsPtr = std::bind(&clVariablesFD::loadParams, this);
	};
	virtual ~clVariablesFD(){ };

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////                                               ///////////////
//////////////           KERNELS/REDUCTIONS/SOLVERS          ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	   

	Kernel calcNu;			// Calculates Nusselt along wall
	Kernel updateNuVars;		// Updates variables used in calculation of Nu 
	Kernel TempUpdateCoeffs;	// Updates Temperature coefficients from SS
								// to transient
	Kernel SetSSCoeffs;			// Sets SS coefficients for thermal solver
	
	Reducer<double> sumTemp;	

	// Sparse Methods/Classes
	CSR_Inds TempInds;
	BiCGStabSolver Temp;
	

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

	Array2Dd Alphat;		// turbulent_alpha (updated with turbulent
							// diffusion coefficients)

	Array2Dd tempArray;		// Temp array, the pointer of which is passed
							// to the thermal BiCGStab solver

	Array3Dd alphaArray;	// Alpha values for each direction at each node
							// Updated each time step (handles 

	//dX arrays, which were previously implemented in this function moved
	//to clVariablesLS, since they will also be used for turb parameters

////////////////////////////////////////////////////////////////////////////	
//////////////                   Method Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////	
//////////////                   Output Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////	
	Array2Dd Nu;

	//buffers containing Nu(x) along both walls on host and device
#ifdef USE_ORIG_NU
	Array1Dv2d LBdist_bot, LBdist_top;
	Array1Dv2i BLinds, LBnodey;
	Array1Dd LBdist_temp;
#endif

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
	bool restartRunFlag;	//Flag indicating if all conditions are met
							//for thermal solver for restarting run
	bool thermalSolverFlag;	//Flag indicating if temp being solved
	bool tempLoadedFlag;	//Flag indicating if T array loaded from file
	bool calcNuFlag;		//Flag indicating if Nu will be calculated
	bool saveMacroStart;	//Flag indicating if T should be saved after
							//initialization
	bool calcSSTempFlag;	//Flag indicating if steady state temp should be
							//solved during intialization.
	double PrTurbNum, PrNum;//Turbulent Pr and Pr
	double Alpha_fluid;		//Thermal diffusivity of fluid in LB units
	double Alpha_foul;		//Thermal diff of fouling layer in LB units
	double T_Actual_Max;	//Actual temperature at inlet
	double T_Actual_Min;	//Actual wall temp
	double T_Actual_Diff;	//Tinlet - Twall
	double ROE0;			//Dimensionless temp to initialize domain to
	double ROE_INX;			//Dimensionless temp at inlet
	int tempMaxIters;		//Max iterations for BiCGStab solver for temp
	double tempMaxRelTol, tempMaxAbsTol;	//Tolerance for BiCGStab solver
	double kSoot, kAir;		//Thermal conductivity of soot and air
							//in SI units 
	double rhoSoot, cpSoot; //Density and Cp of soot in SI units

////////////////////////////////////////////////////////////////////////////	
//////////////                Method Variables               ///////////////
////////////////////////////////////////////////////////////////////////////	
	
	//Pre-calculated values used to simplify calculation of Alpha array
	double Alpha_den_air, Alpha_den_soot, Alpha_num_soot, Alpha_num_air;

	int Save_loc_Nu;	//offset for current Nu values to save in Nu_buffer
	int numNuNodes;	//Number of nodes where Nu is calculated
	
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

	//Initializes alpha array by filling it with air thermal diffusivity
	void iniAlpha();

	void iniNuKernel();

	void createNuArrays();

	void reIniNuKernel();


	////////////////////////////////////////////////////////////////////////////	
	//////////////              Updating Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	//Fills Alpha Array with Values based on interface thickness
	void updateAlphaArray();


	////////////////////////////////////////////////////////////////////////////	
	//////////////               Solving Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	// Calculates average temp in Domain
	double calcAverageT();

	// Solves for SS temperature to initialize run
	void solveSSThermal();

	// Enqueues kernels to solve T for current time step
	void Solve();


	////////////////////////////////////////////////////////////////////////////	
	//////////////                Output Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////
	
	//Enqueues kernel to update Nu along boundary
	void UpdateNu();


	////////////////////////////////////////////////////////////////////////////	
	//////////////                Display Functions              ///////////////
	////////////////////////////////////////////////////////////////////////////


	






	





	







	////Creates coefficients used for updating A and B matricies
	//void Create_Coeffs();

	//void Find_Neighbors();

	////Initializes dX and alpha arrays
	//void ini_bounds();

	//void Update_Transient_Coeffs();


	////Used in calculation of dX and Alph arrays, called by bcFindNodes_C0 and bcFindNodes
	//void bcFindDirection(int dir, int bl, cl_double2 vC0, cl_double2 vC1,
	//	cl_double2 vCn, cl_int2 vCmin, cl_int2 vCmax, int fflag);

	////used in finding boundary nodes
	//bool bcFindIntersectionNormal(cl_double2 &vC, double &dist, cl_double2 vL0,
	//	cl_double2 vP0, cl_double2 vP1, cl_double2 &vN);

	////Used to fill dX and alpha arrays
	//void bcFindNodes(int bl);  //finds boundary nodes based on vls.C

	//void bcFindNodes_C0(int bl); //boundary nodes based on vls.C0

	////Sets bounday nodes for wall (dX)
	//void bcSetBoundaryNode(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
	//	int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vC2, int xyz);

	////sets boundary nodes for interface (dX_cur)
	//void bcSetInterfaceNode(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
	//	int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vC2, int xyz);

	
	
	
	
	
	
	
	
	////Node spacing arrays (1 between nodes, <= 1 near wall/interface) 
	//Array3Dd dX_full;	//dX arrays in full domain
	//Array3Dd dX_cur_full;	//dX_cur in full domain
	//Array2Dd dX;	//dX(i,j) is initial distance from node i to E (j=0), W (j=1),
	//				// N (j=2), and S (j=3) neighbors/wall
	//Array2Dd dX_cur;//Current distance from wall to neighbor/wall/interface

	
	//Array1Dv4i Neighs;






	

};


#endif // !defined(AFX_CLVARIABLESFD_H__INCLUDED_)