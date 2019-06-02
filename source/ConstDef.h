// ConstDef.h : include file for constants
// (c) Alexander Alexeev, 2007
#if !defined(__CONSTDEF_H__INCLUDED_)
#define __CONSTDEF_H__INCLUDED_


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Switches and defines for simulation control/debugging/output control //
// These are not defined elsewhere (i.e. in parameter files) since they //
// must be known at compilation time                                    //
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

#define CHECK_KERNEL_STATUS
#define PRINT_ERROR_MESSAGES		
//#define DEBUG_TURBARR




#define OPENCL_VERSION_1_2
#define OPENMP_NUMBER_THREADS	4



#define CREATE_BIN_FILES
//#define RENAME_BIN_FILES
#define RENAME_DUMP_FILES

#define LOG_ERROR_IN_FILE

//#define DISABLE_TURBULENT_VISC		//Only turns off turbulent visc. in lattice boltzmann solver.


//#define INLET_OUTLET_BC
//#define MAINTAIN_FLOWRATE

//#define PRINT_OPENCL_INFO
#define PRINT_BUILD_INFO
//#define PROFILING

#define RESIZE_MULTIPLIER		1


////////AIR PROPERTIES (IN REAL UNITS)///////////////
#define PIPE_DIAMETER_REAL			2.e-3
#define VISC_AIR					4.1e-5
#define DENSITY_AIR					0.675
#define	THERMAL_CONDUCTIVITY_AIR	0.042
#define SOOT_NUMBER_CONCENTRATION	1.e10



//////////////////////////////////////////////////////////////////////////
///////////                 WORKGROUP SIZES                    ///////////
//////////////////////////////////////////////////////////////////////////

#define BLOCKSIZE					128

#define WORKGROUPSIZE_CSRMV			256
#define WORKGROUPSIZE_AXPBY			256


#define WORKGROUPSIZEX_LB			BLOCKSIZE
#define WORKGROUPSIZEY_LB			1

#define WORKGROUPSIZE_RED			256
#define WORKGROUPSIZE_IBB			256
#define WORKGROUPSIZE_INLET			16

#define WORKGROUPSIZEX_FD			16
#define WORKGROUPSIZEY_FD			16
#define WORKGROUPSIZE_NU			256

#define WORKGROUPSIZE_TR			64
#define WORKGROUPSIZE_TR_WALL		32
#define WORKGROUPSIZE_RERELEASE		64.
#define WORKGROUPSIZE_TR_GL			256
#define WORKGROUPSIZE_TR_SHEAR		64
#define WORKGROUPSIZE_SORT			256
#define WORKGROUPSIZE_PLOC			64
#define WORKGROUPSIZE_TR_DISTS		16
#define WORKGROUPSIZE_TR_WALL_REFLECT	16

#define WORKGROUPSIZE_FL_SHIFT		64
#define WORKGROUPSIZE_UPDATEM		64
#define WORKGROUPSIZE_UPDATEFD		64
#define WORKGROUPSIZE_UPDATEWALL	64
#define WORKGROUPSIZE_SHIFT_PAR		64
#define WORKGROUPSIZE_UPDATE_GL		256







//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Switches and defines for simulation control/debugging/output control //
// These are default values, that can be overridden in parameter file   //
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
#define DEVICE_ID						0



/////////////////////////////////////////////////////////////
//////                   Time Params                   //////
/////////////////////////////////////////////////////////////

#ifdef PROFILING
#define STOP_TIME					(1000)
#else
#define STOP_TIME					(200000000)
#endif

#define DUMP_STEP_START					(0)
#define DUMP_STEP_NUM					(500)
#define DUMP_STEPS_PER_SAVE				(1)
#define DUMP_STEP_START_SAVE			(0)
#define DUMP_STEPS_PER_DUMP_BIN			(1)

#define CURRENT_SAVE_STEP_NUM			(2000)
#define CURRENT_SAVE_STEP_NUM_NU		(2000)
#define CURRENT_SAVE_STEP_NUM_SS		(2000)
#define CURRENT_SAVE_STEP_NUM_IO		(2000)

#define CURRENT_SAVE_START_TIME		(0)
#define CURRENT_SAVE_START_TIME_NU	(0)
#define CURRENT_SAVE_START_TIME_SS	(0)
#define CURRENT_SAVE_START_TIME_IO	(0)

#define LB_STEPS					(1)
#define TR_STEPS_PER_LB				(5)
#define TR_STEPS_PER_LB_WALL		(3)
#define FD_STEPS_PER_LB				(1)

#define CLUMP_TIME					(100000)
/////////////////////////////////////////////////////////////
//////              Fluid Solver Params                //////
/////////////////////////////////////////////////////////////

// Fluid Initialization
#define PERTURB_VEL_FIELD           true
#define INI_TURB_PROFILE			true
#define RUN_LB_FD_FIRST             true
#define INI_POISEUILLE              true
#define FLUID_SAVE_MACROS_ON_START  false
#define USE_TURBULENCE_MODEL        true
// Initialization time (to obtain pseudo steady state)
// dump step is time between saves during this initialization
// LB ini steps is steps with just lb being solved before
// Steady state temperature is calculated
// LBFD is number of steps solving transient LB and temp.
#define TLBM_INI_DUMP_STEP			(200000)
#define NUM_LB_INI_STEPS			(4000000)
#define NUM_LBFD_INI_STEPS			(10000000)


// Perturb parameters
#define DUPLUS_MULT				0.025
#define EPSILON_MULT            250.
#define BETA_MULT               100.
#define ALPHA_MULT              55.75
#define SIGMA_MULT              0.00055



#define ROUGHNESS_FACTOR		9.

// Keep MU_NUMBER and FORCE_TERM_VAL as functions. (Do not set those directly here)
#define UX_INLET				(0.02)
#define RE_NUMBER				(10800.)
#define LB_TAU_NUMBER			(1.99)
#define MU_NUMBER				((1./LB_TAU_NUMBER - 0.5)/3.)
#define LB_RHO_NUMBER			1.

// Controls for iteratively finding dP necessary
// for a set flow rate
#define NUM_INTERVALS_PER_AVG	    500
#define TIME_BETWEEN_INTERVALS	    100
#define PAUSE_BETWEEN_ADJ		    100000
#define MAX_FLOWRATE_PERCENT_DIFF	0.02


#define KOMEGA_MAX_REL_TOL			1.0e-4
#define KOMEGA_MAX_ABS_TOL			1.0e-7
#define KOMEGA_MAX_ITERS			1000



/////////////////////////////////////////////////////////////
//////                Thermal Parameters               //////
/////////////////////////////////////////////////////////////
#define USE_THERMAL_SOLVER			true
#define CALC_NUSSELT				false
#define SOLVE_SS_TEMP               false


#define MAX_SS_FD_STEPS				(10000000)
#define THERMAL_SAVE_MACROS_ON_START            false

//Actual Temps (to re-dimensionalize)
#define LBT_ACTUAL_TEMP_MIN			(90.+273.)
#define LBT_ACTUAL_TEMP_MAX			(400.+273.)

//Wall Temperature
#define TFD_INI_TEMP		1.9
#define TFD_X_IN			2.0
#define TFD_WALL			1.0

#define PR_TURB_NUMBER			(0.9)
#define PR_NUMBER				(0.708)
#define ALPHA_NUMBER			(MU_NUMBER/PR_NUMBER)

#define NU_SKIP				10
#define THERMAL_MAX_REL_TOL			1.0e-4
#define THERMAL_MAX_ABS_TOL			1.0e-7
#define THERMAL_MAX_ITERS			1000



/////////////////////////////////////////////////////////////
//////               Particle Parameters               //////
/////////////////////////////////////////////////////////////
#define USE_PARTICLE_SOLVER		false
#define TR_SAVE_MACROS_ON_START	false
#define SAVE_WALL_SHEAR			true
#define X_RELEASE_POS			(490)
#define X_STOP_POS				(10500)
#define REDUCE_DEP_STOP1		(990)
#define REDUCE_DEP_STOP2		(1490)
#define AMT_REDUCE_DEP			0.1
#define START_THERMO_VEL		(500)
#define AVG_PAR_PER_RELEASE		(50)

#define X_MIN_VAL			(X_RELEASE_POS - 1.)


#define INDEX_RADIUS_SEARCH		30	//Radius of indicies to search when finding Shear Neighbors
#define MASS_FLUX_INLET			10
#define CUTOFF_RADIUS			5.
#define LIFT_COEFFICIENT		81.2
#define MAX_BL_PER_NODE			6
#define TRC_NUM_TRACERS			1048576//(4194304)//(1048576)
#define NUM_STEPS_BTW_SORT		(40)
#define PARTICLE_RELEASE_TIME	(0)
#define NUM_PAR_SIZES			(1)
#define MEAN_FREE_PATH_AIR		(68.0e-9)
#define FOUL_SIZE_SWITCH_SOOT2	0.5


#define STOP_DIST_X				0
#define K_coeff					0.4711

#define NUM_EACH_PAR			(50.) ////Min number of particles represented in a given Particle Size
#define PARTICLE_VOL_MULTIPLIER	(NUM_EACH_PAR/(1. - DEP_POROSITY))
#define TIME_BEFORE_RERELEASE	( 20000 / TR_STEPS_PER_LB_WALL )   //After the particle contacts surface this many times before fully depositing it is re-released
#define TIME_BEFORE_STICK		( 20000 / TR_STEPS_PER_LB_WALL )
#define MAX_NUM_BL_ROLL			(10)

/////////SOOT PROPERTIES (IN REAL UNITS)/////////////
#define SURF_ENERGY_SOOT			0.15
#define Y_MOD_SOOT					(35e9)
#define POISSON_SOOT				0.126
#define DEP_POROSITY				0.98
#define THERMAL_CONDUCTIVITY_FOUL	0.057	//W/(m*K)
#define DENSITY_SOOT				35.      //Kg/m^3
#define	HAMAKER_CONST				1.e-20
#define WALL_CORRECTION				1.7
#define THERMAL_CONDUCTIVITY_PARTICLE	5.	//W/(m*K)
#define HEAT_CAPACITY_SOOT			866.8	//J/(kg*K)
////////SURFACE PROPERTIES (IN REAL UNITS)///////////
#define SURF_ENERGY_SURF			1.37
#define Y_MOD_SURF					(210e9)
#define POISSON_SURF				0.29



/////////////////////////////////////////////////////////////
//////               OUTPUT PARAMETERS                 //////
/////////////////////////////////////////////////////////////

#define OUTPUT_MAX_LINES	128
#define OUTPUT_MAX_LINES_NU	10
#define OUTPUT_MAX_LINES_SS 10
#define OUTPUT_MAX_LINES_IO 100
#define REDUCE_RESULTS_SIZE	128

#define SAVE_SHEAR_LOC		"bothWalls"  //can be "bothWalls", "topWall" or "bottomWall"

#ifdef RENAME_DUMP_FILES
#define FLAG_DUMP	true
#else //RENAME_DUMP_FILES
#define FLAG_DUMP	false
#endif //RENAME_DUMP_FILES



#define OUTPUT_DIR				"results"
#define AVG_OUTPUT_FILE			"tout.txt"
#define UNIX_OUTPUT_FILE		"stdout.txt"
#define NUSSELT_OUTPUT_FILE		"Nusselt.txt"
#define STRESS_OUTPUT_FILE		"shear.txt"
#define IO_OUTPUT_FILE			"IO_dists.txt"
#define YAML_PARAM_FILE			"RunParam.yaml"

/////////////////////////////////////////////////////////////
//   OPENGL PARAMETERS  (vars located in particleDisplay)  //
/////////////////////////////////////////////////////////////
#define OPENGL_DISPLAY			false	 //turn openGL on or off
#define OPENGL_GRIDLINES		false    //Renders gridlines on window if defined
#define NUM_PAR_GL_DIV		1	//Reduces the number of particles plotted by increasing
#define POINT_SIZES			2.	//initial size of particles
#define LINE_SIZES			2.	//initial thickness of lines

#define screenWidth  1800		//width of opengl window
#define screenHeight 1000		//height of opengl window


/////////////////////////////////////////////////////////////
//////                 FL PARAMETERS                   //////
/////////////////////////////////////////////////////////////
#define USE_FL_SOLVER		false
#define FL_SAVE_ON_START	true

#define UPDATE_TIME				(8000)
#define NUM_INLET_OUTLET_NODES	30
#define LS_SPACING				1.
#define NUM_UPDATES_BTW_SMOOTHING	20
#define PERCENT_USED_IN_SMOOTHING	0.05
#define NEIGHS_PER_SIDE_SMOOTHING	3



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//              CONSTANTS (THESE SHOULD NEVER CHANGE)                   //
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

#define KO_C_MU_0_25		    (0.5477225575)


#define ARRAY_INCREASE_SIZE		100
#define NUMBER_DIMENSIONS		2
#define LB_NUMBER_COMPONENTS	3

#define LB_NUMBER_CONNECTIONS	9

#define LB_SOLID				0
#define LB_LIQUID1				1
#define LB_LIQUID2				2

#define LB_NO_BOUNDARY_LINK		0x0
#define LB_BOUNDARY_LINK_1		0x1
#define LB_BOUNDARY_LINK_2		0x2
#define LB_BOUNDARY_LINK_3		0x4
#define LB_BOUNDARY_LINK_4		0x8
#define LB_BOUNDARY_LINK_5		0x10
#define LB_BOUNDARY_LINK_6		0x20
#define LB_BOUNDARY_LINK_7		0x40
#define LB_BOUNDARY_LINK_8		0x80
#define LB_BOUNDARY_LINK_9		0x100
#define SHEAR_NODE_EXISTS		0x20000000

#define PI_NUMBER				3.1415926535897931
#define INDEX					12


#define POINT_VBO				1
#define LINE_VBO				0
#define CONST_COLOR_VBO			0
#define VAR_COLOR_VBO			1

#define NUM_ELEMENTS_IN_LBSIZE	2

#ifdef CREATE_BIN_FILES
#define SAVE_BIN_IN_LOAD
#endif

#ifdef SAVE_BIN_IN_LOAD
#define BIN_SAVE_DIR	"load" SLASH
#else //SAVE_BIN_IN_LOAD
#define BIN_SAVE_DIR	""
#endif //SAVE_BIN_IN_LOAD


#define C_DIR 0
#define E_DIR 1
#define W_DIR 2
#define N_DIR 3
#define S_DIR 4
#define NE_DIR 5
#define SW_DIR 6
#define SE_DIR 7
#define NW_DIR 8

#define C_REV C_DIR
#define E_REV W_DIR
#define W_REV E_DIR
#define N_REV S_DIR
#define S_REV N_DIR
#define NE_REV SW_DIR
#define NW_REV SE_DIR
#define SE_REV NW_DIR
#define SW_REV NE_DIR



//Type 1 Boundaries are IBB boundaries applied before collision step (q < 0.5)
//Type 2 Boundaries are IBB boundaries applied after collision step (q >= 0.5)
#define C_BOUND 0x8
#define E_BOUND 0x10
#define W_BOUND 0x20
#define N_BOUND 0x40
#define S_BOUND 0x80
#define NE_BOUND 0x400
#define SW_BOUND 0x800
#define SE_BOUND 0x1000
#define NW_BOUND 0x2000
#define E_BOUND_T1 0x4000
#define W_BOUND_T1 0x8000
#define N_BOUND_T1 0x10000
#define S_BOUND_T1 0x20000
#define NE_BOUND_T1 0x40000
#define SW_BOUND_T1 0x80000
#define SE_BOUND_T1 0x100000
#define NW_BOUND_T1 0x200000
#define E_BOUND_T2 0x400000
#define W_BOUND_T2 0x800000
#define N_BOUND_T2 0x1000000
#define S_BOUND_T2 0x2000000
#define NE_BOUND_T2 0x4000000
#define SW_BOUND_T2 0x8000000
#define SE_BOUND_T2 0x10000000
#define NW_BOUND_T2 0x20000000

#define SOLID_NODE 0x0
#define FLUID_NODE 0x1
#define FOULED_NODE 0x2
#define BOUNDARY_NODE 0x40000000
#define GHOST_NODE 0x4


// vls.M flag values 
#define M_SOLID_NODE	0b0001
#define M_FLUID_NODE	0b0010
#define M0_SOLID_NODE	0b0100
#define M0_FLUID_NODE	0b1000
#define M_FOUL_NODE		0b1001





#define OPTION_SAVE_MACRO_FIELDS			0x2
#define OPTION_LOCAL_UPDATE					0x8
#define OPTION_TIME_SAMPLE					0x20

#define REDUCE_TYPE_SUM				1
#define REDUCE_TYPE_MAX				2
#define REDUCE_TYPE_ABSMAX			3
#define REDUCE_TYPE_MIN				4




#ifdef DISABLE_TURBULENT_VISC
#ifdef PERTURB_VEL_FIELD
#undef PERTURB_VEL_FIELD
#endif
#ifndef USE_PARABOLIC_INI_VEL
#define USE_PARABOLIC_INI_VEL
#endif
#endif

#define TRP_TOP_LOC_IND		0
#define TRP_BOT_LOC_IND		1
#define TRP_UMAX_VAL_IND	2
#define TRP_BVAL_IND		3
#define TRP_OFFSET_Y_IND	4



#define ERROR_OCL_INITIALIZATION							-100
#define ERROR_BUFFER_ALLOCATION								-101
#define ERROR_MEMORY_ALLOCATION								-102
#define ERROR_MEMORY_COPY									-103
#define ERROR_MEMORY_OUT_OF_BOUNDS							-104
#define ERROR_OPENING_FILE									-105
#define ERROR_KERNEL_INITIALIZATION							-106
#define ERROR_CALLING_KERNEL								-107
#define ERROR_SETTING_KERNEL_ARG							-108
#define ERROR_IN_REDUCE_KERNEL								-109
#define ERROR_LOADING_PARAMS								-110
#define ERROR_INITIALIZING_KOMEGA							-111
#define ERROR_INITIALIZING_VLS								-112
#define ERROR_INITIALIZING_VTR								-113
#define ERROR_INITIALIZING_VFL								-114
#define ERROR_CREATING_BICGSTAB_SOLVER						-115
#define ERROR_IN_OPENGL_ARRAY								-116

#ifdef _DEBUG

#ifndef DEBUG_TURBARR
#define DEBUG_TURBARR
#endif

#ifndef PRINT_ERROR_MESSAGES
#define PRINT_ERROR_MESSAGES
#endif

#endif


#define CL_HPP_MINIMUM_OPENCL_VERSION 120


#ifdef OPENCL_VERSION_1_2
#define CL_HPP_TARGET_OPENCL_VERSION	120
#define CL_HPP_CL_1_2_DEFAULT_BUILD
#endif

#define getGlobalSizeMacro(actSize, wgSize)	\
	((int)ceil((double)actSize / (double)wgSize)* (int)wgSize)

#endif // !defined(__CONSTDEF_H__INCLUDED_)

///////////Possible places for optimization///////////////
