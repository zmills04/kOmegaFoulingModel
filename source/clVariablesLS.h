// clVariables.h: interface for the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLVARIABLESLS_H__INCLUDED_)
#define AFX_CLVARIABLESLS_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clProblem.h"

/////////////////////////////////////////////////////////////////////////
class clVariablesLS
{
public:

	// Func Pointer for calling loadParams
	std::function<void(void)> loadParamsPtr;
	
	clVariablesLS() : BL("lsbl"), Masses("lsmasses"), C("lsc"),
		C0("lsc0"), M("lbm"), M_o("lbm0"), s("lbms"), sf("lbmsf")
	{
		loadParamsPtr = std::bind(&clVariablesLS::loadParams, this);
	};
	virtual ~clVariablesLS()
	{
		//Release_Objects(); 
	};

	enum {find_ibb, find_dx_cur, find_dx};

	Array2Di BL;		// Contains Indicies of C that define wall boundaries
	Array1Dd Masses;	// Number of lattice sites in reduced domain containing 
						// Fluid (ind 0) and fouling layer (ind 1)
	Array1Dv2d C, C0;	//Location of boundary points (C is current location
						// and C0 is original location)
	Array2Dc M, M_o;	// Array defining the current state of lattice nodes with
						// fluid nodes = 0 and solid/fould nodes = 1. (M is current and
						// M_o is original)
	Array2Dc s, sf;		// s defines current state of nodes (1 for solid/fouling and 0 for fluid)
						// sf same but with fouling = 1 and 0 otherwise.
	Array2Dc M_red;		// Reduced form of M array

	//OpenGL Arrays
	Array1DGL LSt_vbo, LSb_vbo;   //Current Position of the Walls
	Array1DGL LSt0_vbo, LSb0_vbo; //Original Position of the Walls	
	Array1DGL LinesV, LinesH;	  //Gridlines (Uncomment OPENGL_GRIDLINES to use)
	

	double lsSpacing;
	int nN, nBL;		// Number of wall nodes and boundary links	
	bool saveMacroStart, restartRunFlag;
	//allocates buffers on device and writes values to them
	void allocateArrays();

	// Calls functions to initialize variables
	void ini();			

	//Initializes openGL variables
	void iniGL();

	// Fills M array
	void inifillMap();		

	// Fills M_o array
	void inifillMap0();	

	// Fill M_red array (and calls vlb.ini_transform_arrays)
	//void inifill_mapred();	
	
	//fills s and sf arrays
	void inifillS();		

	//Frees host memory that is no longer necessary after initialization
	//of device memory
	void freeHostMem();
	
	void loadParams();

	// Releases all memory associated with Arrays
	// called upon completion of simulation
	void releaseObjects();

	void saveParams();

	//Saves Output variables
	void save2file();
	
	//saves variables useful for debugging
	void saveDebug();

	void saveRestartFiles();
	
	
	
//////////Min Max functions (used throughout code) //////////////	
	
	//Returns int2 with floored values of the minimums of v0 and v1
	cl_int2 min2(cl_double2 v0, cl_double2 v1);
	//Returns int2 with ceiled values of the maximums of v0 and v1
	cl_int2 max2(cl_double2 v0, cl_double2 v1);


/////////////////// IBB variables ///////////////////////////////////

	Array1Dv2i ii0_array;	// Array containing locations of boundary nodes (within fluid)
	// one entry for each distribution direction that crosses boundary
	Array1Dv2i iis1_array;	// Array containing location of second node used in IBB calculations
	// each entry corresponds to each element in ii0_array
	Array1Di dir_array;		// distribution direction or corresponding node in ii0_array 
	// which crosses the boundary
	Array1Dd D_array;		// distance between node and boundary (normalized by distance traveled
	// between nodes for a given direction)
	Array2Di Bflag;			// Used to make sure that multiple Bouncebacks are not applied to the
	// same distribution function
	Array1Dv3i Bnodes;

	Array3Dd dXArr, dXArrCur;

	Array2Di Node_Numbers;

	Array1Ds IBB_flag_heat;
	Array1Di Node_loc_heat, Node_loc_heat_dist;
	Array1Dv8i Neighs_elbm_heat;
	Array1Dv2d dxdy_elbm_heat;
	Array2Dd ibb_coeff_heat, ibb_coeff_heat_temp;
	
	Array1Ds IBB_flag_fluid;
	Array2Dd ibb_coeff_fluid, ibb_coeff_fluid_temp;
	Array1Di Node_loc_fluid;
	Array1Dv8i Neighs_elbm_fluid;
	Array1Dv2d dxdy_elbm_fluid;

	int cur_el_elbm, fullsize_elbm_fluid, num_el_fluid, length_ibb_fluid;
	int num_el_heat, cur_len_heat;

	int cur_el;				// Index of current element to fill in the IBB arrays
	int num_el;				// Current number of entries in ibb arrays
	int length_ibb;			// Next multiple of WORKGROUPSIZE_IBB larger than num_el
	int nBnodes, curElBnodes, lengthBnodes, bot_ind_end;
	int fullsize_ibb_arrays, fullsize_Bnodes, fullsize_Bnodes_old;
	int imin_Bnode = X_RELEASE_POS - 10;
	int imax_Bnode = X_STOP_POS + 10;
	int LFNearest = LB_BOUNDARY_LINK_1 | LB_BOUNDARY_LINK_2 |
		LB_BOUNDARY_LINK_3 | LB_BOUNDARY_LINK_4;
///////////Functions for determining IBB variables//////////////////
	
	void ini_dX_arrays();

	void update_dXCur_array();
	
	//Initializes IBB arrays and calls update_IBB_arrays
	void ini_IBB_arrays();

	//Fills Arrays used to update IBB arrays in vlb class
	//Calls bcFindNodes
	void update_IBB_arrays();
	
	//Finds lattice sites surrounding Wall nodes 
	//Calls bcFindDirection
	void bcFindNodes(int bl, int findtype);

	//Finds which distributions will cross wall when streaming from
	//nodes found in bcFindNodes
	//Calls bcFindIntersection and if it returns true, calls bcSetBoundaryNode
	void bcFindDirection(int dir, int bl, cl_double2 vC0, cl_double2 vC1, 
		cl_double2 vCn, cl_int2 vCmin, cl_int2 vCmax, int findtype);

	//Determines if distribution will intersect with wall when streaming
	//Calls bcFindIntersectionLinePlane and testInside which must be true
	//for this function to return true
	BOOL bcFindIntersection(cl_int2 &vCcut0, cl_int2 &vCcut1, cl_double2 &vCcut,
		double &dist, cl_double2 vL0, cl_double2 vLd, cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);
	
	//Finds where intersection of distribution direction
	//and line crossing btw nodes of the BL is located (vC).
	//returns false if no point exists
	BOOL bcFindIntersectionLinePlane(cl_double2 &vC, double &dist, cl_double2 vL0, cl_double2 vLd,
		cl_double2 vP0, cl_double2 vP1);

	//Determines if intersection (vd) is located between v0 and v1
	//returns true if it is btw and false if not
	BOOL testInside(cl_double2 vd, cl_double2 v0, cl_double2 v1);

	//Fills values in IBB arrays (ii0,_array, iis1_array, dir_array and D_array)
	//With values used to create IBB arrays used in IBB_kernels.
	//Calls testperiodicindex to determine if distribution
	//crosses across periodic x boundary
	void bcSetBoundaryNode(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
		int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);
	
	void bcSetdXArr(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
		int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);

	void bcSetdXCurArr(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
		int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);


	//Determines if distribution crosses BL across periodic boundary
	BOOL testperiodicindex(cl_int2 i0, cl_int2 i1);

	void bcFindNodes_Heat(int bl);
	
	void bcFindDirection_Heat(int dir, int bl, cl_double2 vC0,
		cl_double2 vC1, cl_double2 vCn, cl_int2 vCmin, cl_int2 vCmax);
	
	void bcSetBoundaryNode_Heat(cl_int2 ii0, cl_int2 ii1, int dir,
		double dist, int bl, cl_double2 vCc, cl_double2 vC0,
		cl_double2 vC1, cl_double2 vCn);

	void create_IBB_arrays_heat();

	BOOL bcFindIntersection_Heat(cl_int2 &vCcut0, cl_int2 &vCcut1,
		cl_double2 &vCcut, double &dist, cl_double2 vL0, cl_double2 vLd,
		cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);

	// Empty functions to allow for templated function pointer wrapper
	// without having to use template specialization or SFINAE
	void createKernels() {}
	void setKernelArgs() {}

};

#endif // !defined(AFX_CLVARIABLESLS_H__INCLUDED_)