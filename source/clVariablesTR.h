// clVariables.h: interface for the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLVARIABLESTR_H__INCLUDED_)
#define AFX_CLVARIABLESTR_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clProblem.h"

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
		//Release_Objects();
	};

	Kernel TR_ReRelease_kernel;
	Kernel TR_Shear_Removal_kernel;
	Kernel TR_Update_Par_kernel;
	Kernel TR_Wall_Par_kernel[5];
	Kernel TR_GL_kernel;
	Kernel TR_Node_kernel[2], TR_Wall_Node_kernel[2];
	Kernel Shear_kernels[3];
	Kernel Sort_kernel[3];
	Kernel Save_shear_kernel;
	Kernel Get_Umax_kernel;
	Kernel Save_Dists_kernel;
	Kernel Update_SS_kernel[2];
	Kernel Clump_Particle_kernels[2];

	Array1DP P;		//Particles
	Array2DNT Nod;	//temporary node array used to create NodI and NodC
	Array2DNI NodI;	//contains list of BL's near each node and a flag if
					//node contains wall (0 for no wall, 1 for bottom and 2 for top)
	Array2DNC NodC;	//Contains the coefficients used to calculate NodV
	Array2DNV NodV; //Contains velocities and temps used for interpolation at each 
	Array1DPP parP; //particle parameters
	Trparam trP;	//parameters used to re-release particles
	Array1DBL BL;	//Boundary Link info
	Array1Dv2i Bindicies_loc;
	
	Array1Dv2u RandList;//rng streams for device (one for each particle)
	Array1Dd Umax_val;
	Array2Di BLind_ind;	//Current index to save BL number when updateing NodI array
	Array2Du BL_dep;	//Number of particles deposited at a node since last update step
	Array1Dv4i Sind;		//Index of first node used in averaging to get shear at wall
	Array1Di BLindicies; //Index of BL array each index in SS arrays corresponds to	
	Array1Dv4d Weights;	//Weights applied to surrounding Shear stresses at nodes
	Array1Dv2i Shear_inds;	//locations of nodes where shear is calculated
	Array1Dv8d Shear_coeffs;	//coefficients used in calc. of shear at nodes
	Array1Dv2i Ploc;		//tracks the location of particles after sort
	Array1Dv2d PV;			//particle velocities
	Array1Dd Tau;			//Shear stresses calculated at boundary nodes
	Array1DP SortTmp;		//Temp array for particles during sort if there is an odd num of merges
	Array2Di Node_neigh;	//indicies of neighbors surrounding each node
	Array1Di Winds;			//indicies of wall nodes
	Array1Di Num_W_nodes;	///Number of nodes along the wall
	Array1DGL P_vbo;		//OpenGL object for particles
	Array2Dd SS_output;		//Array to store Shear Stresses for output
	Array1Ds Update_flag;	//indicates if particle has been updated or not (to avoid two
							//updates when it crosses into a new node
	BLbound Bounds;				//Indicies of BL's that are at start and end of tracer area
	Array1Di TR_indicies;		//List of indicies of NodI, NodC etc that are active for tracers
	Array1Dv2i Active_indicies;	//List of active nodes
	Array2Di Active_flag;	///indicates if a node in domain is active or not
	Array1Df Par_Color_List;
	Array1Di IO_inds;
	Array1Du IO_dists_save;
	Array1Di Reflect_info;
	Array1Di Reflect_inds;

	cl_int2 TrDomainSize;
	bool trSolverFlag, calcSSFlag, calcIOFlag, saveMacroStart;
	int Save_IO_Loc;
	cl_uint2 IO_inds_info;
	double Umean_Current;
	cl_int4 Ploc_inds;
	bool restartRunFlag;
	int Bnode_top_start;
	cl_uint Par_conc_inlet;
	int nActiveNodes;
	int nN, Nd;
	double X_release;
	double rmax;
	int Sort_timer, OpenGL_timer;	
	int localRange;
	cl_uint numMerges;
	int Shear_array_len;
	int numbl_bounds;
	int Num_wall_nodes, Num_wall_nodes_max;
	int Save_loc_SS;
	int numpargl;
	int Num_nodes_at_release;
	Array1Dd Tcrit_color;

//////////Simulation Parameters/////////////////
	int indRadiusSearch, maxBLPerNode;
	int parReleaseTime, stopDistX, timeBeforeReRelease;
	int timeBeforeStick, maxNumBLRoll, xReleasePos;
	int xStopPos, reduceDepStop1, reduceDepStop2, startThermoVel;
	double amtReduceDep, xMinVal, massFluxInlet;
	double cutoffRadius, foulSizeSwitchSoot, numEachPar, parVolMultiplier;





////////Particle Properties/////////////////////
	void save_particle_arrays();

	double WofA, Estar, WofA_s, Estar_s;
	Array1Dd D_p, D_p_real, M_p, V_par, D_dists, F_po, Kth_pars;
	Array2Dd Q_A_prime, Q_A, Tau_crit, R_d;
	double Tau_crit_max;
	int Mean_index, numStepsBtwSort;
	double sootNumConc, parThermalCond, mfpAir;
	double surfEnergySurf, poissonSurf, yModSurf;
	double surfEnergySoot, poissonSoot, yModSoot;
	double hamakerConst, depPorosity, wallCorrection;
	double liftCoeff;
////////////////////////////////////////////////	
	
	void setSourceDefines();
	void saveParams();
	void loadParams();
	void testRestartRun();
	void createKernels();
	void ini();
	void iniGL();
	void ini_trp();
	void ini_particle_properties();
	void ini_node();
	void ini_particles();
	void ini_shear_coeffs();
	int get_ind(cl_int2 ij);
	double weight_kernel_func(double Sval);
	void redistribute(int i);
	double rand1();
	void release_distribute();
	void ini_blinks();
	int get_type();
	void ini_rand();
	BOOL test_cross(cl_double2 vL0, cl_double2 vLd, cl_double2 vP0, cl_double2 vP1);
	void test_node(int i, int j, cl_double2 vL0, cl_double2 vL1, int blind, int bot_flag);
	void allocate_buffers();
	void setKernelArgs();
	void freeHostMem();
	void Release_Objects();
	void save2file();
	void Split_Node_arrays();
	void ini_sort();
#ifdef USE_OPENGL
	void sort_particles(cl_event *TR_prev_evt);
#else
	void sort_particles();
#endif
	
	void iniShear();
	void update_Shear_arrays();
	void update_trp();
	void Find_node_neighs();
	void Check_Neigh_Node(int i, int j, int *cur_ind, int nod_loc);
	void Re_Release_Par();
	void Update_Shear();
	void Update_Par_Rem_Args();
	void UpdateSS();
	void reset_SS_out();
	void saveDebug();
	void Save_SS(); 
	void initial_sort_particles();
	void Update_Transient_Coeffs();
	double get_offset(double xval);
	void save_box(double x1, double y1, double dx, double dy);
	void ini_particle_colors();
	void Clump_Particles();
	par AddPars(par p1, par p2);
	void sort_particles_for_clumping();
	void Update_IO_Dists();
	void Reset_IO_Dists();
	void saveRestartFiles();
	void save_restart_files_ini();    ///files that dont change after initialization (only need bin saved once)
	void call_update_wall_particles();
	void call_update_wall_particles(cl_event *fill_evt);
};

///////////////////////////////////////////*/
#endif // !defined(AFX_CLVARIABLESTR_H__INCLUDED_)