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
	// Func Pointer for calling loadParams
	std::function<void(void)> loadParamsPtr;

	clVariablesFD() : Temp_array("temp"), Alphat("alphat")
	{
		loadParamsPtr = std::bind(&clVariablesFD::loadParams, this);
	};
	virtual ~clVariablesFD(){ };

	// Kernels
	Kernel Nu_kernel, TempUpdateCoeffs;
	Kernel SetSSCoeffs;
	Reducer<double> sumTemp;

	// Sparse Methods/Classes
	CSR_Inds TempInds;
	BiCGStabSolver Temp;
	
	bool thermalSolverFlag, tempLoadedFlag, calcNuFlag, 
		saveMacroStart, calcSSTempFlag;
	double PrTurbNum, PrNum;
	double Alpha_fluid;		//Thermal diffusivity of fluid
	double Alpha_foul;		//Thermal diffusivity of fouling layer
	double T_Actual_Max;	//Actual temperature at inlet
	double T_Actual_Min;	//Actual wall temp
	double T_Actual_Diff;
	double ROE0;			//Initial temperature Temp is set to
	double ROE_INX;
	bool restartRunFlag;	//Set to 1 if run is being restarted	
	int tempMaxIters; 
	double tempMaxRelTol, tempMaxAbsTol;
	double kSoot, kAir, rhoSoot, cpSoot;

	Array2Dd Alphat;				//alpha+turbulent_alpha (updated with turbulent diffusion coefficients)
	Array2Dd Temp_array;
	
	//allocates device/host buffers
	void allocateArrays();

	//Calculates average temp in Domain
	double calcAverageT();

	//Initializes kernels for thermal solver
	void createKernels();

	//Frees memory on host thats not necessary once copied to device buffers
	void freeHostMem();

	//initializes Method
	void ini();

	//Initializes alpha array by filling it with air thermal diffusivity
	void iniAlpha();
	
	void loadParams();

	//Releases buffers
	void releaseObjects();

	//Attempts to load bin files and returns true if all necessary arrays are loaded
	BOOL restartRun();

	//Saves T distribution
	void save2file();

	void saveParams();

	void saveRestartFiles();

	//saves variables necessary for debugging
	void saveDebug();

	//Sets arguments of OpenCL kernels
	void setKernelArgs();

	void setSourceDefines();

	// Solves for SS temperature to initialize run
	void solveSSThermal();
	
	//Enqueues kernels to solve T for current time step
	void solveTemp();

	// Loads checkpoint files, reads restartFlag parameter
	// and sets appropriate flags
	void testRestartRun();

	//Fills Alpha Array with Values based on interface thickness
	void updateAlphaArray();

	//Writes data to buffers allocated on the host.
	int writeToBuffers();




















	/// Nu Calculations 

	Array2Dd Nu;
	//buffers containing Nu(x) along both walls on host and device
#ifdef USE_ORIG_NU
	Array1Dv2d LBdist_bot, LBdist_top;
	Array1Dv2i BLinds, LBnodey;
	Array1Dd LBdist_temp;
#endif
	double Alpha_den_air, Alpha_den_soot, Alpha_num_soot, Alpha_num_air;
	int Save_loc_Nu; //offset for current Nu values to save in Nu_buffer

	cl_int num_nu_nodes;	//Number of nodes where Nu is calculated


	//Enqueues kernel to update Nu along boundary
	void UpdateNu();

	void Create_Nu_arrays();

	void Reinitialize_Nu_kernel();

	void initialize_Nu_kernel();

	//Reads Nu from device and resets Save_loc_Nu to 0
	void reset_Nu_out();










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
	//BOOL bcFindIntersectionNormal(cl_double2 &vC, double &dist, cl_double2 vL0,
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