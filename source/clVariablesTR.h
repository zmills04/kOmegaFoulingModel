// clVariables.h: interface for the clVariables class.
// Due to large number of methods and variables, the methods and
// variables associated with particle sorting and those associated
// with the shear stress calculation have been implemented as their
// own class in separate source files and instances of those classes
// have been included in this class.
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLVARIABLESTR_H__INCLUDED_)
#define AFX_CLVARIABLESTR_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clProblem.h"
#include "particleProperties.h"
#include "particleSort.h"
#include "shearStress.h"
#include "particleDisplay.h"
#include "particleStructs.h"

typedef struct BLbound
{
	cl_int MIN_BL_BOT;
	cl_int	MAX_BL_BOT;
	cl_int MIN_BL_TOP;
	cl_int MAX_BL_TOP;

} BLbound;


class clVariablesTR
{
public:
	// Func Pointer for calling loadParams
	std::function<void(void)> loadParamsPtr;
	
	clVariablesTR()
	{
		loadParamsPtr = std::bind(&clVariablesTR::loadParams, this);
	};
	virtual ~clVariablesTR()
	{
		releaseObjects();
	};



	//////////////////////////////////
	/////    Particle Solver     /////
	//////////////////////////////////
	
	// Updates interpolation nodes with [0] temperature
	// and [1] velocity information (non wall-adjacent nodes)
	Kernel TR_Node_kernel[2];

	// Updates interpolation nodes with [0] temperature
	// and [1] velocity information (wall-adjacent nodes)
	Kernel TR_Wall_Node_kernel[2];

	// Updates particle locations in non-wall adjacent nodes 
	Kernel TR_Update_Par_kernel; 

	// Updates particle locations in wall adjacent nodes
	// performs calculations to deposit/reflect particles
	Kernel TR_Wall_Par_kernel[5];
	
	// Calculates particle distributions at inlet and outlet
	// and prepares it to be written to file
	Kernel Save_Dists_kernel;

	// Sections of vtr split into own classes
	particleSort parSort;
	particleProperties parP;
	shearStress wallShear;
	particleDisplay glParticles;



	// TODO: sort and label (comment) all these variables
	bool restartRunFlag;
	bool trSolverFlag, calcIOFlag, saveMacroStart;
	cl_int2 TrDomainSize;
	int Save_IO_Loc;
	cl_uint2 IO_inds_info;
	double Umean_Current;
	int nActiveNodes;
	int nN;

	int Num_nodes_at_release;
	double X_release;

	double rmax;

	int Num_wall_nodes, Num_wall_nodes_max;




	//////////Simulation Parameters/////////////////
	int indRadiusSearch, maxBLPerNode;
	int parReleaseTime, stopDistX, timeBeforeReRelease;
	int timeBeforeStick, maxNumBLRoll, xReleasePos;
	int xStopPos, reduceDepStop1, reduceDepStop2, startThermoVel;
	double amtReduceDep, xMinVal, massFluxInlet;
	double cutoffRadius, foulSizeSwitchSoot;



	//Particles
	//Array1DP P;	
	Par P;


	//temporary node array used to create NodI and NodC
	// Memory should be deallocated after use
	//Array2DNT *Nod;
	//NodeT* Nod;

	//contains list of BL's near each node and a flag if
	//node contains wall (0 for no wall, 1 for bottom and 2 for top)
	//Array2DNI NodI;
	NodeI NodI;

	//Current index to save BL number when updating NodI array
	Array2Di BLind_ind;

	//Contains the coefficients used to calculate NodV
	//Array2DNC NodC;
	NodeC NodC;

	//Contains velocities and temps used for interpolation at each 
	//Array2DNV NodV; 
	NodeV NodV;

	//rng streams for device (one for each particle)
	Array1Dv2u RandList;

	//Boundary Link info
	//Array1DBL BL;	
	BLinks BL;

	// Number of particles deposited at a node since last update step
	// Will be added to vfl array storing total deposit (after some smoothing) 
	Array2Du BL_dep;

	// particle velocities (info must be kept between updating of location and
	// deposit kernels
	Array1Dv2d PV;		

	//indicies of neighbors surrounding each node
	Array2Di Node_neigh;	
	
	//indicies of wall nodes
	Array1Di Winds;

	///Number of nodes along the wall
	Array1Di Num_W_nodes;	

	//indicates if particle has been updated or not (to avoid two
	//updates when it crosses into a new node)
	Array1Ds Update_flag;	

	//Indicies of BL's that are at start and end of tracer area
	BLbound Bounds;				

	//List of indicies of NodI, NodC etc that are active for tracers
	Array1Di TR_indicies;		
	
	//List of active nodes
	Array1Dv2i Active_indicies;	

	///indicates if a node in domain is active or not
	Array2Di Active_flag;	
	
	// For saving inlet/outlet distributions
	Array1Di IO_inds;
	Array1Du IO_dists_save;

	// Stores information about particles reflecting
	// when updating particles near walls
	Array1Di Reflect_info;
	Array1Di Reflect_inds;


	void allocateArrays();
	void allocateBuffers();
	void setSourceDefines();
	void saveParams();
	void loadParams();
	void testRestartRun();
	void createKernels();
	void setKernelArgs();
	void freeHostMem();
	void releaseObjects();
	void save2file();


	void ini();
	void iniNode();
	void iniParticles();
	int getInd(cl_int2 ij);
	double weightKernelFunc(double Sval);
	double rand1();
	void iniBLinks();
	int getType();
	void iniRand();
	bool testCross(cl_double2 vL0, cl_double2 vLd, cl_double2 vP0, cl_double2 vP1);
	void testNode(int i, int j, cl_double2 vL0, cl_double2 vL1, int blind, int bot_flag);
	void splitNodeArrays();
	void findNodeNeighs();
	void checkNeighNode(int i, int j, int *cur_ind, int nod_loc);
	
	void saveDebug();
	

	void updateTransientCoeffs();
	double getOffset(double xval);
	void saveBox(double x1, double y1, double dx, double dy);
	legacyPar AddPars(legacyPar p1, legacyPar p2);
	void updateIODists();
	void resetIODists();
	void saveRestartFiles();
	void saveRestartFilesIni();    ///files that dont change after initialization (only need bin saved once)
	void callUpdateWallParticles(cl_event *fill_evt = nullptr);
};

///////////////////////////////////////////*/
#endif // !defined(AFX_CLVARIABLESTR_H__INCLUDED_)