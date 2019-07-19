// clVariablesTR.h: interface for the clVariables class.
// Due to large number of methods and variables, the methods and
// variables associated with particle sorting, shear stress, particle
// properties, and displaying with opengl have been implemented
// as individual classes in separate source files. Instances of those
// classes have been included in this class. Each of these subclasses
// follow the same structure as the other clVariable classes.
//
// Note: nodes in vtr represent area enclosed by indicies 
//       [(i,j),(i+1,j),(i,j+1),(i+1,j+1)]
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////



// TODO: Edit comments to make sure that tracers is used to refer to
//		parcels being tracked in simulation domain, and particles refers
//		to the total number of particles being represented by the parcels.

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

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////    CONSTRUCTOR/DESTRUCTOR/ENUMS/FUNC_PTRS     ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	
	// Func Pointer for calling loadParams
	std::function<void(void)> loadParamsPtr;

	clVariablesTR() : ioDistsSave("ioDists"), blDep("blDep_temp"),
		reflectInfo("reflectInfo"), reflectInds("reflectInds"), 
		updateFlag("updateFlag"), PV("PV"), RandList("RandList"), 
		activeInds("activeInds", WORKGROUPSIZE_TR, 2),
		wallInds("wallInds", WORKGROUPSIZE_TR, 2)
	{
		loadParamsPtr = std::bind(&clVariablesTR::loadParams, this);
	};

	virtual ~clVariablesTR()
	{
	};

	// Sections of vtr split into own classes
	particleSort parSort;
	particleProperties parP;
	shearStress wallShear;
	particleDisplay glParticles;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////                                               ///////////////
//////////////           KERNELS/REDUCTIONS/SOLVERS          ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

	// Updates interpolation nodes with [0] temperature
	// and [1] velocity information (non wall-adjacent nodes)
	//Kernel trNodeKernel[2];

	// Updates interpolation nodes with [0] temperature
	// and [1] velocity information (wall-adjacent nodes)
	//Kernel trWallNodeKernel[2];

	// Updates particle locations in non-wall adjacent nodes 
	Kernel trUpdateParKernel;

	// Updates particle locations in wall adjacent nodes
	// performs calculations to deposit/reflect particles
	Kernel trWallParKernel[5];

	// Calculates particle distributions at inlet and outlet
	// and prepares it to be written to file
	Kernel saveDistsKernel;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
//////////////                                               ///////////////
//////////////            ARRAYS/PARTICLESTRUCTS             ///////////////
//////////////                                               ///////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	
////////////////////////////////////////////////////////////////////////////	
//////////////                   Data Arrays                 ///////////////
////////////////////////////////////////////////////////////////////////////
	
	// Boundary Link info
	BLinks BL;

	// Particles
	Par P;

////////////////////////////////////////////////////////////////////////////	
//////////////                   Method Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////

	// Indicies of active nodes (this is all nodes i.e. wall nodes + non-wall nodes)
	// Similar array is wallInds, which is a list of the wall node indicies 
	DynArray1Di activeInds;

	// indicates if a node in domain is active or not
	// size is trDomainSize
	// No longer needed as this info is stored in NodI.wallFlag
	// Array2Di activeFlag;

	// Number of particles deposited at a node since last update step
	// Will be added to vfl array storing total deposit (after some smoothing)
	// size is vls.nBL
	Array2Du blDep;

	// Current index to save BL number when updating NodI array
	// Used in legacy code. May no longer be needed, or its implementation
	// may need to be changed slightly
	//Array2Di blIndInd;

	// indicies of neighbors surrounding each node
	// Used in legacy code. May no longer be needed, or its implementation
	// may need to be changed slightly
	//Array2Di nodeNeigh;

	// Contains the coefficients used to calculate NodV
	// size - trDomainSize
	NodeC NodC;

	// contains list BL closest to center of trNode a flag if
	// node contains wall (0 for no wall, 1 for bottom and 2 for top)
	// size - trDomainSize
	NodeI NodI;

	// Contains velocities and temps used for interpolation at each node
	// size - trDomainSize
	// NodeV NodV;

	// temporary node array used to create NodI and NodC
	// Memory should be deallocated after use
	// Array2DNT *Nod;
	// NodeT* Nod;

	// Number of nodes along the wall
	//Array1Di numWallNodes;

	// particle velocities (info must be kept between updating of location and
	// deposit kernels
	// size - nN
	Array1Dv2d PV;

	// rng streams for device (one for each particle)
	Array1Dv2u RandList;

	// Stores information about particles reflecting
	// when updating particles near walls
	Array1Di reflectInfo;
	Array1Di reflectInds;

	// List of indicies of NodI, NodC etc that are active for tracers
	// Used in legacy code. May no longer be needed, or its implementation
	// may need to be changed slightly
	// Array1Di trIndicies;

	// indicates if particle has been updated or not (to avoid two
	// updates when it crosses into a new node)
	Array1Ds updateFlag;

	// indicies of wall nodes
	DynArray1Di wallInds;

////////////////////////////////////////////////////////////////////////////	
//////////////                   Output Arrays               ///////////////
////////////////////////////////////////////////////////////////////////////	

	// For saving inlet/outlet distributions
	//Array1Di ioInds;				// indicies of inlet and outlet nodes
									// to save dists in
									// shouldnt need this as the indicies 
									// can be easily calculated

	TimeData<cl_uint> ioDistsSave;	// distributions at inlet and outlet
									// to save as timedata

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

	bool restartRunFlag;
	bool trSolverFlag, saveMacroStart;
	bool calcIOFlag;		// Flag to save distribution at inlet outlet
							// as TimeData array

	int maxOutputLinesIO;	// number of lines in ioDistsSave

	int nN;					// Number of particles

	int indRadiusSearch;	// Radius of indicies to search when finding 
							// Shear Neighbor
	
	int timeBeforeReRelease;// Starting value for Par.timer (i.e. time a given particle
							// can remain in domain before it is re-released at inlet)
							// Should be a large number since this is not physical, and
							// mainly a way to make sure particles are eventually re-released
							// if they escape domain and its not caught, or they are unable to
							// escape from a recirculation region.

	int timeBeforeStick;	// starting value for Par.deptimer. Time a particle is temporarily
							// deposited on a BL and able to be removed before its considered
							// permanently deposited and re-released.

	int maxNumBLRoll;		// maximum number of BL particle can roll over as its looking for
							// local minimum to deposit on.
	
	double xReleaseVal;		// Actual x position particles are released from
	
	int xReleasePos;		// x index particles are released from

	double xMinVal;			// Minimum x value for particle before its re-released
							// Needs to be less than xReleaseVal

	double xStopVal;		// x position at which particles leave trdomain
							// and are re-released
	
	int xStopPos;			// max x index (calculated from ceil(cStopVal))							// and are re-released

	// X locations near inlet where amtRedDep is applied to limit the thickness of
	// deposit layer. This was necessary when performing laminar simulations to 
	// keep the deposit layer from forming too fast at inlet and crashing.
	// may not be as necessary with turb flow where shear stress will be higher
	// x < reduceDepStop1 - volume of deposited particle multiplied by amtRedDep
	// x > reduceDepStop1 and x < reduceDepStop2 - volume multiplied by interpolated
	//		value between amtReduceDep (at reduceDepStop1) and 1 (at reduceDepStop2)
	// x > reduceDepStop2 - volume of particle unchanged
	int reduceDepStop1, reduceDepStop2;	
	double amtReduceDep;	// must be between 0 and 1
	
	// Like the reduce deposit variable, this was another method to avoid rapid
	// layer formation at inlet causing simulation to become unstable. This defines
	// the x value at which thermophoretic velocity begins to be applied. (i.e. between
	// xReleasePos and startThermoVel no thermophoretic velocity is applied to particle)
	int startThermoVel;

	
	double cutoffRadius;	// radius used in search for boundary nodes near BL
							// when calculating interpolation weights for shear
							// calculation at BL
	
	double foulSizeSwitchSoot;	// thickness of deposit layer where surface properties
								// switch from metal to those of fouling layer
	
	double lsSpacing;		//starting length of BLinks

	// Legacy variables which may no longer be needed 
	// (Check and remove once done refactoring)

	int maxBLPerNode;		// Legacy, originally used to define size of BLind in NodI 
	int parReleaseTime;		// legacy, time to start releasing particles (probably unneeded
							// since the fluid and temp is solved before starting main time loop)

	int stopDistX;			// X location where particles stop and are re-released	
							// legacy code had stopDistX defined as 0 and not used 
							// anywhere so most likely unnecessary
	
	int massFluxInlet;		// not sure about this variable (define variable unused
							// in legacy code)

	int blSearchRad;		// Number of BLinks to look at on either side of NodI.BLind

////////////////////////////////////////////////////////////////////////////	
//////////////                Method Variables               ///////////////
////////////////////////////////////////////////////////////////////////////	

	// Indicies of BL's that are at start and end of tracer area
	BLbound Bounds;

	cl_int2 trDomainSize;	// size of tracer domain 
							// (x values between release and stop positions)
	int trDomainFullSizeX;	// padded size of trDomainSize.x

	cl_int2 trDomainXBounds;	// [min i index, max i index) 
							// last i index is actually (max i - 1) since in vtr nodes
							// span from i to i+1



//int Save_IO_Loc;			// no longer needed as TimeData class stores this info

	cl_uint2 IO_inds_info;		// used with calculation of IO distributions
								// possibly legacy

	double Umean_Current;		// current mean Ux velocity in domain
	int nActiveNodes;			// Number of active nodes in domain


	int Num_nodes_at_release;	// number of nodes in y direction at x location of release
								// possibly legacy

	double rmax;				//(double)RAND_MAX			

	// Number of wall nodes varies as walls foul, so arrays which
	// relate to the wall calculations need to be dynamic. These 
	// variables are used to track the actual number of nodes and
	// the full size of padded arrays. This is most likely no longer
	// needed if dynamicArrays are used.
	int Num_wall_nodes, Num_wall_nodes_max;


	int ioDistsArgIndex;


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

	// checks if neighboring node at (i,j) is an active node
	// adds it to Node_Neigh array if it does. This is a legacy code function
	// and either is no longer needed or needs significant refactoring
	//void checkNeighNode(int i, int j, int* cur_ind, int nod_loc);

	// fills NodI with information from BLinks
	void fillNodeI(int i, Array2Dd& dist2Center);



	// Converts [i,j] index to linear index
	int getInd(cl_int2 ij);
	
	// Gets y location of wall for a given x location
	double getOffset(double xval);

	// Generates the particle size represented by a tracer (according to
	// particle size distribution)
	int getType();
  	
	// Tests if vector between vL0 and vLd crosses line between vP0 and vP1
	bool testCross(cl_double2 vL0, cl_double2 vLd, cl_double2 vP0, cl_double2 vP1);

	// Tests if node contains a BL and sets flags in NodeI.wallFlag if so
	// needs updating to handle new way BL index is handled by NodeI (i.e. legacy
	// code stored all BL indicies in area defined by node, while new code only stores
	// index of one closest to center.
	void clVariablesTR::testNode(int i0, int j0, cl_double2 c0t, cl_double2 c1t,
		int blind, Array2Dd& dist2Center);

	// Generate random number between 0 and 1 from uniform distribution
	double rand1();

	// Calculates interpolation weight using basic kernel
	// TODO: implements various kernels, and allow for
	//		the kernel used to be selected with input parameter
	double weightKernelFunc(double Sval);


////////////////////////////////////////////////////////////////////////////	
//////////////            Initialization Functions           ///////////////
////////////////////////////////////////////////////////////////////////////
	
	// Initializes particles
	void iniParticles();

	// initializes BLinks and NodI data
	void iniBLandNodI();
	void iniBLinksBottom();
	void iniBLinksTop();
	void iniNodI();

	// Initializes random number array used in kernels
	void iniRand();

	// Two functions generate NodI and NodC structures of arrays

	void iniNode();			// generates temporary Nod structure of arrays
	void iniNodC(Array2Dd& distInfoX, Array2Dd& distInfoY,
		Array2Dd& distInfoX0, Array2Dd& distInfoY0, Array2Di& distInfoType);	// splits Nod into NodC and NodI

	// Initializes activeInds and wallInds arrays
	void iniIndexInfo();


	// Initializes arrays containing node neighbor arrays
	// Legacy code restructed arrays to remove any permanent wall nodes
	// which allowed wavy channel to be represented by straight channel.
	// This meant (i,j) and (i+1,j+1) were not necessarily neighbors, so
	// code needed to implement methods to find neighbors. Refactored code
	// does not do this restructuring, so may not need to use same techinques
	// as legacy code.
	// void iniNodeNeighs();

////////////////////////////////////////////////////////////////////////////	
//////////////              Updating Functions               ///////////////
////////////////////////////////////////////////////////////////////////////

	

////////////////////////////////////////////////////////////////////////////	
//////////////               Solving Functions               ///////////////
////////////////////////////////////////////////////////////////////////////
	
	// Calls kernels to update tracer locations
	void updateWallParticles();


////////////////////////////////////////////////////////////////////////////	
//////////////                Output Functions               ///////////////
////////////////////////////////////////////////////////////////////////////

	// Saves all data within box drawn in opengl window.
	// Very useful for debugging
	void saveBox(double x1, double y1, double dx, double dy);

	// files that dont change after initialization (only need bin saved once)
	// void saveRestartFilesIni();    

////////////////////////////////////////////////////////////////////////////	
//////////////                Display Functions              ///////////////
////////////////////////////////////////////////////////////////////////////





















	
};

///////////////////////////////////////////*/
#endif // !defined(AFX_CLVARIABLESTR_H__INCLUDED_)