// clVariables.h: interface for the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLVARIABLESLS_H__INCLUDED_)
#define AFX_CLVARIABLESLS_H__INCLUDED_

#pragma once

#include "StdAfx.h"
#include "clProblem.h"

// TODO: Combine all arrays containing node information, to reduce to as few
//		as possible. This includes node info in vtr, vfd, vlb and vls.
class clVariablesLS
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

	clVariablesLS() : BL("lsbl"), Masses("lsmasses"), C0("lsc0"), C("lsc"),
		M("lbm"), ibbArr(WORKGROUPSIZE_IBB*2,"ibbArr", WORKGROUPSIZE_IBB, 1),
		ssArr(WORKGROUPSIZE_TR_SHEAR * 2, "ssArr", WORKGROUPSIZE_TR_SHEAR, 1),
		dXArr("dXArr"), dXArr0("dXArr0"), bFlag("bFlag")

	{
		loadParamsPtr = std::bind(&clVariablesLS::loadParams, this);
	};
	virtual ~clVariablesLS()
	{
		//Release_Objects(); 
	};

	enum { find_ibb, find_dx, find_dx0 };


	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	//////////////                                               ///////////////
	//////////////           KERNELS/REDUCTIONS/SOLVERS          ///////////////
	//////////////                                               ///////////////
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////


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




	Array1Dv2d C, C0;	//Location of boundary points (C is current location
						// and C0 is original location)

	Array2Dc M;	// Array defining the current state of lattice nodes
				// bit flags set for whether it is a solid node or fluid node
				// currently and at t = 0. (From this, foul node can be determined)


////////////////////////////////////////////////////////////////////////////	
//////////////                   Method Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////

	Array3Dd dXArr0;	// Spacing between nodes and wall, with dXarr being
	Array3Dd dXArr;		// the original spacing and dXarrCur being the current
						// spacing


	DynArray1Di ibbArr;	// Array containing location of distributions bouncing back 
						// (within fluid). 
						// Use decodeGlobalIdx(ii0_array[ind], &xval, &yval, &distval)
						// to get node location in x,y,dist coordinates.

	DynArray1Di ssArr;	// Array containing location of nodes at which shear is calculated 
						// Use vlb.NodeType to get directions pointing into wall 
						// (Shear stress calculation only considers N,S,E,W
						// directions crossing BL)

	Array1Dd Masses;	// Number of lattice sites containing 
						// Fluid (ind 0) and fouling layer (ind 1)

	Array2Di BL;		// Contains Indicies of C that define wall boundaries

	Array2Di bFlag;		// Used to make sure that multiple Bouncebacks are not applied to the
						// same distribution function by setting bits once a value is set
						// for a given direction.

	Array1Di ssArrInds; // index of first element in ssArr for a given x value
						// i.e. ssArrInds(2) is the first element in ssArr that is located at x = 2

	// These arrays were originally used for IBB implementation. New implementation
	// will just use ii0 array, which contains absolute position of bounced back distribution
	// and dXArrCur, which contains distances (trying to avoid reading from memory as much as possible)

	//Array1Dv2i Bnodes;	// when used previously, each element was 3d int vec with x and y being
						// the i and j value of a boundary node, and z was the location of node
						// info in ii0, iis1, dir, etc arrays 
						// should no longer be needed because ii0 will be the list of boundary nodes,
						// and it contains info about actual position while dXArr contains dist
						// information

	//Array1Dv2i iis1_array;	// Array containing location of second node used in IBB calculations
	// each entry corresponds to each element in ii0_array
	//Array1Di dir_array;		// distribution direction or corresponding node in ii0_array 
	// which crosses the boundary
	//Array1Dd D_array;		// distance between node and boundary (normalized by distance traveled
	// between nodes for a given direction)

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
	double lsSpacing;
	int nN, nBL;		// Number of wall nodes and boundary links	
	bool saveMacroStart, restartRunFlag;

	////////////////////////////////////////////////////////////////////////////	
	//////////////                Method Variables               ///////////////
	////////////////////////////////////////////////////////////////////////////	
	static const int LF[9]; // Flags used for setting bFlag array indicies
	int maxSSArrNodeX;	// max x position of boundary nodes tracked for ss calculations
	int minSSArrNodeX;	// min x position of boundary nodes tracked for ss calculations




	// Not sure about these variables (may not be necessary)
	//int imin_Bnode = X_RELEASE_POS - 10;
	//int imax_Bnode = X_STOP_POS + 10;
	//int LFNearest = LB_BOUNDARY_LINK_1 | LB_BOUNDARY_LINK_2 |
	//	LB_BOUNDARY_LINK_3 | LB_BOUNDARY_LINK_4;

		//int curElIBB;			// Index of current element to fill in ii0_array
	//int numElIBB;			// Current number of entries in ii0_array arrays
	//						// num_el ~= lengthIBB since array is buffered to allow
	//						// its size to increase without constant re-sizing
	//int lengthIBB;			// Next multiple of WORKGROUPSIZE_IBB larger than num_el
	//int nBnodes, curElBnodes, lengthBnodes, bot_ind_end;
	//int fullsize_ibb_arrays, fullsize_Bnodes, fullsize_Bnodes_old;


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
	void saveTimeData() {}

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
	void updateTimeData() {};


	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////	
	//////////////                                               ///////////////
	//////////////            CLASS SPECIFIC FUNCTIONS           ///////////////
	//////////////                                               ///////////////
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////


	//////////Min Max functions (used throughout code) //////////////	
	//Returns int2 with floored values of the minimums of v0 and v1
	cl_int2 min2(cl_double2 v0, cl_double2 v1);
	//Returns int2 with ceiled values of the maximums of v0 and v1
	cl_int2 max2(cl_double2 v0, cl_double2 v1);

	//Finds which distributions will cross wall when streaming from
	//nodes found in bcFindNodes
	//Calls bcFindIntersection and if it returns true, calls bcSetBoundaryNode
	void bcFindDirection(int dir, int bl, cl_double2 vC0, cl_double2 vC1,
		cl_double2 vCn, cl_int2 vCmin, cl_int2 vCmax, int findtype);

	//Determines if distribution will intersect with wall when streaming
	//Calls bcFindIntersectionLinePlane and testInside which must be true
	//for this function to return true
	bool bcFindIntersection(cl_int2& vCcut0, cl_int2& vCcut1, cl_double2& vCcut,
		double& dist, cl_double2 vL0, cl_double2 vLd, cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);

	//Finds where intersection of distribution direction
	//and line crossing btw nodes of the BL is located (vC).
	//returns false if no point exists
	bool bcFindIntersectionLinePlane(cl_double2& vC, double& dist, cl_double2 vL0, cl_double2 vLd,
		cl_double2 vP0, cl_double2 vP1);

	//Finds lattice sites surrounding Wall nodes 
	//Calls bcFindDirection
	void bcFindNodes(int bl, int findtype);

	//Fills values in IBB arrays (ii0,_array, iis1_array, dir_array and D_array)
	//With values used to create IBB arrays used in IBB_kernels.
	//Calls testperiodicindex to determine if distribution
	//crosses across periodic x boundary
	//void bcSetBoundaryNode(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
	//	int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);

	// Fills ibb array based on dXArr values (if not = 1.0 and a fluid node)
	// called by both iniIBBArray and updateIBBArray
	cl_bool2 fillIBBArray();

	//Determines if intersection (vd) is located between v0 and v1
	//returns true if it is btw and false if not
	bool testInside(cl_double2 vd, cl_double2 v0, cl_double2 v1);

	//Determines if distribution crosses a BL across periodic boundary
	bool testPeriodicIndex(cl_int2 i0, cl_int2 i1);


	////////////////////////////////////////////////////////////////////////////	
	//////////////            Initialization Functions           ///////////////
	////////////////////////////////////////////////////////////////////////////

	// Sets value in dXArr0 based on values passed to it.
	void bcSetdXArr0(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
		int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);

	// Compares M but flags in M array to count number of solid, fluid and fouling
	// layer nodes
	void iniCountSolid();

	// Fills dXArr and dXArr0
	void inidXArrays();

	// Fills M array
	void iniFillMap();

	// Fills M0 array
	void iniFillMap0();

	//Initializes IBB arrays and calls update_IBB_arrays
	//void iniIBBArray();


	////////////////////////////////////////////////////////////////////////////	
	//////////////              Updating Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	// Sets value in dXArr based on values passed to it.
	void bcSetdXArr(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
		int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vCn);

	// Function to update dXArrCur (either on host, or by calling kernels
	// on device)
	void updatedXArr();

	//Fills Arrays used to update IBB arrays
	void updateIBBArray();

	////////////////////////////////////////////////////////////////////////////	
	//////////////               Solving Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////	
	//////////////                Output Functions               ///////////////
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////	
	//////////////                Display Functions              ///////////////
	////////////////////////////////////////////////////////////////////////////

};

#endif // !defined(AFX_CLVARIABLESLS_H__INCLUDED_)