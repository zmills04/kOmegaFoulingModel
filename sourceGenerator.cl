#define FULLSIZEX		(2000)
#define FULLSIZEY		(200)
#define CHANNEL_LENGTH		(2000)
#define CHANNEL_HEIGHT		(204)
#define CHANNEL_LENGTH_FULL		(2048)
#define DIST_SIZE		(417792)
#define FULLSIZEXY		(400000)
#define FULLSIZEY_UT		(201)
#define DOMAIN_SIZE_Y		(204)
#define BLOCK_SIZE		(128)
#define WORKGROUPSIZE_RED		(256)
#define WORKGROUPSIZEX_LB		(128)
#define WORKGROUPSIZEY_LB		(1)
#define WORKGROUPSIZE_INLET		(16)
#define WORKGROUPSIZE_IBB		(256)
#define WORKGROUPSIZEX_FD		(16)
#define WORKGROUPSIZEY_FD		(16)
#define WORKGROUPSIZE_NU		(256)
#define WORKGROUPSIZE_RERELEASE		(64)
#define WORKGROUPSIZE_TR		(64)
#define WORKGROUPSIZE_TR_WALL		(32)
#define WORKGROUPSIZE_TR_GL		(256)
#define WORKGROUPSIZE_TR_SHEAR		(64)
#define WORKGROUPSIZE_SORT		(256)
#define WORKGROUPSIZE_PLOC		(64)
#define WORKGROUPSIZE_FL_SHIFT		(64)
#define WORKGROUPSIZE_FL_RAMP		(30)
#define WORKGROUPSIZE_UPDATEM		(64)
#define WORKGROUPSIZE_UPDATEFD		(64)
#define WORKGROUPSIZE_UPDATEWALL		(64)
#define WORKGROUPSIZE_UPDATE_GL		(256)
#define WORKGROUPSIZE_SHIFT_PAR		(64)
#define WORKGROUPSIZE_TR_DISTS		(16)
#define WORKGROUPSIZE_TR_WALL_REFLECT		(16)
#define OPENCL_VERSION_1_2
#define FTERM_VAL		(6.9444539848e-009)
#define INI_MASS		(      400000)
#define tau0		(  1.99445983)
#define UX_INLET_VAL		(        0.02)
#define RE_TURBULENT		(         180)
#define UTAU_VAL		(0.00083333390575)
#define YPLUS_WALL_NODE		(         0.9)
#define K_WALL_VALUE		(6.904534896e-008)
#define NUT_WALL_VALUE		(-0.00038129762506)
#define OMEGA_WALL_VALUE		(0.046340296752)
#define TIN		(           2)
#define TMIN		(         363)
#define TDIFF		(         310)
#define FULLSIZEXM1		(1999)
#define NU_MULTIPLIER		(         400)
#define LOCAL_SOLVE_SIZE_FD		(256)
#define PR_TURB_NUMBER		(         0.9)
#define PR_NUMBER		(           0)
#define TFD_X_IN_VAL		(           2)
#define START_THERMO_VEL		(500)
#define MAX_BL_PER_NODE		(6)
#define TRC_NUM_TRACERS		(1048576)
#define NUM_PAR_SIZES		(0)
#define FULLSIZEX_TR		(0)
#define FULLSIZEY_TR		(0)
#define MU_NUMBER		(0.00046296328097)
#define DTTR		(         0.2)
#define DTTR_WALL		(0.33333333333)
#define FOUL_SIZE_SWITCH_SOOT2		(         0.5)
#define KAIR		(8.9584647449e-009)
#define KSOOT		(1.2157916439e-008)
#define ALPHA_DEN_AIR		(1.3699991555e-005)
#define ALPHA_DEN_SOOT		(0.00057307116273)
#define ALPHA_TEST_AIR		(0.00065390293923)
#define ALPHA_TEST_SOOT		(2.1215369452e-005)
#define PAR_TIMER_START		(2222)
#define DEP_TIMER_START		(2222)
#define X_START_VAL_INT		(490)
#define X_MAX_VAL_INT		(10500)
#define X_START_VAL		(         490)
#define X_MIN_VAL		(     10499.5)
#define X_MAX_VAL		(10500)
#define Y_MIN_VAL		(         0.5)
#define Y_MAX_VAL		(       203.5)
#define XVALS_HEIGHT		(202)
#define NUM_PAR_GL_DIV		(1)
#define ALPHA_FLUID		(0.00065390293923)
#define ALPHA_FOUL		(2.1215369452e-005)
#define MAX_NUM_BL_ROLL		(10)
#define SORT_NUM_MERGES		(8)
#define PERCENT_USED_IN_SMOOTHING		(        0.05)
#define PERCENT_NOT_USED_IN_SMOOTHING		(        0.95)
#define NEIGHS_PER_SIDE_SMOOTHING		(3)
#ifndef OPENCL_VERSION_1_2
#pragma OPENCL EXTENSION cl_khr_subgroups : enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics : enable
#endif

#define FTYPE	double

#define MAX(A, B)		(((A) > (B)) ? (A) : (B))
#define MIN(A, B)		(((A) < (B)) ? (A) : (B))
#define EPSILON			1.11e-16
#define TRUE			1
#define FALSE			0
#define EQUALS2D(A, B)	((A.x == B.x) && (A.y == B.y))
#define NEQUALS2D(A, B)	((A.x != B.x) || (A.y != B.y))
#define CEPS			1.0e-5
#define RAND_MAX		4294967296.
#define MODFAST(A, B)	(((A) >= (B)) ? ((A) - (B)) : (((A) < 0) ? ((A) + (B)) : (A)))
#define CX_DEF_HI	(double4)(1.,-1.,1.,-1.)
#define CY_DEF_HI	(double4)(1.,-1.,-1.,1.)
#define CXX_DEF_HI	(double4)(1.,1.,1.,1.)
#define CXY_DEF	(double4)(1.,1.,-1.,-1.)

#define CX_DEF_HIF	(float4)(1.f,-1.f,1.f,-1.f)
#define CY_DEF_HIF	(float4)(1.f,-1.f,-1.f,1.f)
#define CXX_DEF_HIF	(float4)(1.f,1.f,1.f,1.f)
#define CXY_DEFF	(float4)(1.f,1.f,-1.f,-1.f)

#define PI_NUMBER				3.1415926535897931


enum{ MWC64X_A = 4294883355U };
enum{ MWC64X_M = 18446383549859758079UL };

#define WAVEFRONT_SIZE	64


#define CCxlo	(double4)(1., -1., 0., 0.)
#define CCxhi	(double4)(1., -1., 1., -1.)
#define	CCylo	(double4)(0., 0., 1., -1.)
#define CCyhi	(double4)(1., -1., -1., 1.)

#define CX8		(double8)(1., -1., 0., 0., 1., -1., 1., -1.)
#define CY8		(double8)(0., 0., 1., -1., 1., -1., -1., 1.)
#define WEIGH8	(double8)(1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.)

constant int CDirX[8] = { 1, -1, 0, 0, 1, -1, 1, -1 };
constant int CDirY[8] = { 0, 0, 1, -1, 1, -1, -1, 1 };
constant int RevDir[8] = { 1, 0, 3, 2, 5, 4, 7, 6 };

#ifdef USE_FLOATS
#define alpha_tolerance		1e-10f
#define entropy_tolerance	1e-7f
#define LOG49	(8.10930216216328769718e-01f)
#define LOG19	(2.19722457733621956422f)
#define LOG136	(3.58351893845610991463f)
#else
#define alpha_tolerance		1e-10
#define entropy_tolerance	1e-7
#define LOG49	(8.10930216216328769718e-01)
#define LOG19	(2.19722457733621956422)
#define LOG136	(3.58351893845610991463)
#endif


constant short IBB_Flag[8] = { 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80 };

constant double2 CXY_DOUBLE[8] = { 
	(double2)(1., 0.),
	(double2)(-1., 0.),
	(double2)(0., 1.),
	(double2)(0., -1.),
	(double2)(1., 1.),
	(double2)(-1., -1.),
	(double2)(1., -1.),
	(double2)(-1., 1.)};

#ifndef OPENCL_VERSION_1_2
void AtomicAdd(__global double *val, double delta)
{
	union
	{
		double f;
		ulong  i;
	} old;
	union
	{
		double f;
		ulong  i;
	} new;
	do
	{
		old.f = *val;
		new.f = old.f + delta;
	} while (atom_cmpxchg((volatile __global ulong *)val, old.i, new.i) != old.i);
}
#endif

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
#define BOUNDARY_NODE 0x3FFFFFF8
#define GHOST_NODE 0x4


#define OPTION_SAVE_MACRO_FIELDS			0x2
#define OPTION_LOCAL_UPDATE					0x8
#define OPTION_TIME_SAMPLE					0x20


#define C_NEIGH (0)
#define E_NEIGH (1)
#define W_NEIGH (-1)
#define N_NEIGH (CHANNEL_LENGTH_FULL)
#define S_NEIGH (-CHANNEL_LENGTH_FULL)
#define NE_NEIGH (CHANNEL_LENGTH_FULL+1)
#define NW_NEIGH (CHANNEL_LENGTH_FULL-1)
#define SE_NEIGH (-CHANNEL_LENGTH_FULL+1)
#define SW_NEIGH (-CHANNEL_LENGTH_FULL-1)


#define EVAL1(...) __VA_ARGS__

#define EMPTY()


#define CALL_ALL_DISTS(func)		EVAL1(func EMPTY() (C));\
									EVAL1(func EMPTY() (E));\
									EVAL1(func EMPTY() (W));\
									EVAL1(func EMPTY() (N));\
									EVAL1(func EMPTY() (S));\
									EVAL1(func EMPTY() (NE));\
									EVAL1(func EMPTY() (SW));\
									EVAL1(func EMPTY() (SE));\
									EVAL1(func EMPTY() (NW));\


#define CALL_ALL_DISTS_NO_C(func)	EVAL1(func EMPTY() (E));\
									EVAL1(func EMPTY() (W));\
									EVAL1(func EMPTY() (N));\
									EVAL1(func EMPTY() (S));\
									EVAL1(func EMPTY() (NE));\
									EVAL1(func EMPTY() (SW));\
									EVAL1(func EMPTY() (SE));\
									EVAL1(func EMPTY() (NW));\

#define CALL_WEST_DISTS(func)		EVAL1(func EMPTY() (W));\
									EVAL1(func EMPTY() (NW));\
									EVAL1(func EMPTY() (SW));\

#define CALL_EAST_DISTS(func)		EVAL1(func EMPTY() (E));\
									EVAL1(func EMPTY() (NE));\
									EVAL1(func EMPTY() (SE));\





#define KO_ALPHA1			(5./9.)
#define KO_ALPHA2			(0.44)
#define KO_BETA1			(3./40.)
#define KO_BETA2			(0.0828)
#define KO_BETA_STAR		(9./100.)
#define KO_SIGMA_K1			(0.85)
#define KO_SIGMA_K2			(1.)
#define KO_SIGMA_O1			(0.5)
#define KO_SIGMA_O2			(0.856)
#define KO_A1				(0.31)
#define KO_KAPPA			(0.41)
#define KO_C_MU_0_25		(0.5477225575)
#define KO_VON_KARMAN		(0.41)
#define KO_GAMMA1			(5./9.)
#define KO_GAMMA2			(0.44)


///All structures used in BD solver
//Uses attribute aligned and padding to ensure that
//	memory alignment matches alignment for structs on host
//structures defined in StdAfx.h for host side

typedef struct FD_coeff
{
	double4 cA;			///finite difference coefficients
	double4 cB;
	double4 cC;
	double Alpha[5];
	int Type;
	int indc;
	int inde;		//for outlet inde is i-1, indw is i-2
	int indw;
	int indn;
	int inds;
} FDcoeff __attribute__((aligned(64)));

typedef struct Particle
{
	double2 pos;		///position vector
	uint Num_rep;		///number of particles represented by this particle
	short type;			///type of particle (cooresponds to Param element
	short Dep_Flag; 	//-2 if waiting for re-release, -1 for not deposited, > 0 signifies the BL location it has deposited at.
	ushort Dep_timer;	//timer set to specified value once particle deposits and decrements at 0, particle is deposited
	ushort timer;		///time for use in re-releasing
	int loc;			//Node number particle is located within
} par __attribute__((aligned(32)));

typedef struct Particle_Param
{
	double Q_A_prime[2];	///Used in deposition/rebound calculation
	double Q_A[2];			///used in dep/reb
	double tau_crit[2];	///critical shear stress
	double Dp;			///particle diameter
	double Mp;			///particle mass
	double Kth;			///Thermophoretic Coeff
	double D_dist;		///% of particles distributed in bin
	double L_coeff;		///Lift coefficient (Lift force = L_coeff*Tau^1.5)
	double D_coeff;		///Drag coefficient (Drag force = D_coeff*Tau)
} Pparam __attribute__ ((aligned(32)));

typedef struct NodeCoeff
{
	double4 CoeffT00; //coefficients used to calculate NodeVar variables in update_Nodes kernel
	double4 CoeffT10;
	double4 CoeffT01;
	double4 CoeffT11;
	double4 CoeffU00;
	double4 CoeffU10;
	double4 CoeffU01;
	double4 CoeffU11;
	int4 neigh;		//Indicies of 4 lattice points enclosing square (U and T array are of size vlb.nX*(vlb.Channel_Height+1))
	char Pad[16];
} nodeC __attribute__ ((aligned(32)));

typedef struct NodeInfo
{
	short BLind[MAX_BL_PER_NODE];	//Inidices of three closest boundary links (-1 used to represent when less than three BL's near)
	short Wall_Flag;  //0 is no walls nearby, 1 if bottom wall nearby and 2 if top wall nearby, -1 if inactive
	char Pad[14 - 2 * MAX_BL_PER_NODE];
} nodeI __attribute__ ((aligned(16)));

typedef struct Neighbors
{
	int2 ii00;
	int2 ii10;
	int2 ii01;
	int2 ii11;
} Nact __attribute__((aligned(32)));

typedef struct NodeVar
{
	double4 Temps;
	double2 U00;
	double2 U10;
	double2 U01;
	double2 U11;
} nodeV __attribute__ ((aligned(32)));

typedef struct BL_bounds
{//Defines range of BL used in particle section of domain
	int MIN_BL_BOT;
	int	MAX_BL_BOT;
	int MIN_BL_TOP;
	int MAX_BL_TOP;
} BLbound __attribute__ ((aligned(16)));

typedef struct BL_Links
{
	double2 vP0;		//Location of left node
	double2 vP1;		//Location of right node
	double2 vTvec;	//tangential vector (points to the right)
	double2 vNvec;   //normal vector pointing into domain
	double Tau;		//Shear stress at location
	double blLen;	//length of BL
	int Node_loc;	//points to location BL is located in (or majority of BL when across two)
	int Color_ind;
	short P0ind;		//index of right node in vls.C array (maybe unnecessary)
	short P1ind;		//index of left node in vls.C array
	short dir;		//direction of shear stress
	short int_type;	//Designates if interface between particle and surface is soot-soot or soot-wall
} bLinks __attribute__ ((aligned(32)));

typedef struct Tr_Param
{//Used for re-releasing particles
	double Top_location;	//Y value of uppermost node
	double Bottom_location;	//Y value of lowermost node
	double umax_val;		//Max velocity at inlet
	double bval;			//spacing between upper and lower wall at particle inlet
	double offset_y;		//Location of wall at particle inlet
	double X_release;	//X location of release
	uint BL_rel_bot;	//index of bottom BL at particle inlet 
	uint BL_rel_top;	//index of top BL at particle inlet
	uint Uvals_start;	//location of starting point for Uvals used in re-distribution
	char Pad[4];		//padding
} Trparam __attribute__((aligned(64)));

typedef struct FoulInfo
{//used for updating fouling layer
	double4 WeightsL;	//Weights applied to BL deposits to left
	double4 WeightsR;	//Weights applied to BL deposits to right
	double2 vN;			//Normal vector of C node
	double disp;		//Current displacement distance	
	uint BL_ind;		//index of BL to left of node
	uint C_ind;			//index of wall node this struct corresponds to
} foulI __attribute__((aligned(32)));

typedef struct Ramp_info
{
	double Ybegin;
	double Coeff;
	uint IOind; //Index of LS nodes at beginning and end of TR area
	uint Cind; //direction traversed from IO_end when creating ramp
	char pad[8];

} rampI __attribute__((aligned(32)));

//Second Reduction step to obtain sum of x velocity and density
////This kernel is only used when FD temp solver is not being used
//__kernel __attribute__((reqd_work_group_size(1, 1, 1)))
//void LB_reduce_Ro_J_2(__global double *input,
//__global double *output,
//int savestep)
//{
//	int gid = get_global_id(0);
//
//	double temp = input[gid];
//	for (unsigned int s = 1; s < NUMBLOCKS; s++)
//	{
//		temp += input[s * 2 + gid];
//	}
//
//	output[savestep * 2 + gid] = temp;
//
//}

////Second Reduction step to obtain sum of x velocity and density
////when using FD solver. Also calculates bulk mean temp at inlet and outlet
//__kernel __attribute__((reqd_work_group_size(1, 1, 1)))
//void FD_reduce_Ro_J_T(__global double *input,
//__global double *output,
//__global double *T,
//__global double2 *J,
//int savestep)
//{
//	int gid = get_global_id(0);
//
//	if (gid < 2)
//	{
//		double temp = input[gid];
//		for (unsigned int s = 1; s < NUMBLOCKS; s++)
//		{
//			temp += input[s * 2 + gid];
//		}
//
//		output[savestep * 5 + gid] = temp;
//	}
//	else if (gid == 2)
//	{
//		double temp1 = 0., temp2 = 0.;
//		for (unsigned int s = 0; s < FULLSIZEY; s++)
//		{
//			int ind = s + X_START_VAL_INT * FULLSIZEY_UT;
//			temp1 += fabs(J[ind].x*T[ind]);
//			temp2 += fabs(J[ind].x);
//		}
//		output[savestep * 5 + gid] = temp1 / temp2;
//	}
//	else if (gid == 3)
//	{
//		double temp1 = 0., temp2 = 0.;
//		for (unsigned int s = 0; s < FULLSIZEY; s++)
//		{
//			int ind = s + X_MAX_VAL_INT * FULLSIZEY_UT;
//			temp1 += fabs(J[ind].x*T[ind]);
//			temp2 += fabs(J[ind].x);
//		}
//		output[savestep * 5 + gid] = temp1 / temp2;
//	}
//	else
//	{
//		double temp1 = 0., temp2 = 0.;
//		for (unsigned int s = 0; s < FULLSIZEY; s++)
//		{
//			int ind = s + (FULLSIZEX - 1) * FULLSIZEY_UT;
//			temp1 += fabs(J[ind].x*T[ind]);
//			temp2 += fabs(J[ind].x);
//		}
//		output[savestep * 5 + gid] = temp1 / temp2;
//	}
//}

//Calculates Nusselt along the boundaries
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
void FD_calc_Nu(__global bLinks *BL,
__global int *BLind,
__global double *Nu,
__global nodeV *NV,
__global double2 *J,
__global double *T,
int SaveLoc,
int NUMBLINKS)
{
	int ix = get_global_id(0);

	if (ix >= NUMBLINKS)
		return;

	int shiftval = SaveLoc * NUMBLINKS * 3;
	
	int ii = BLind[ix];
	
	bLinks BLtemp = BL[ii];
	double2 BLpos = BLtemp.vP0 + BLtemp.vTvec * BLtemp.blLen/2.;

	int2 Posi = convert_int2(BLpos);
	int loc = Posi.x*FULLSIZEY_TR + Posi.y;

	nodeV NVtemp = NV[loc];
	double4 TT = (NVtemp.Temps - TMIN) / TDIFF;
	
	double2 dX1 = BLpos - trunc(BLpos);
	double2 dX0 = 1. - dX1;
	
	double T0y = TT.x * dX0.y + TT.z * dX1.y;
	double T1y = TT.y * dX0.y + TT.w * dX1.y;

	double T0x = TT.x * dX0.x + TT.y * dX1.x;
	double T1x = TT.z * dX0.x + TT.w * dX1.x;

	double2 dTdx = (double2)(T1y - T0y, T1x - T0x);  ///-dTdx
	
	double dTdn = dot(dTdx, BLtemp.vNvec); 

	double Um = 0.;
	double Tm = 0.;
	for (int j = 0; j < FULLSIZEY; j++)
	{
		double Unorm = length(J[Posi.x*FULLSIZEY_UT + j]);
		Um += Unorm;
		Tm += Unorm * T[Posi.x*FULLSIZEY_UT + j];
	}
	Tm = Tm / Um;
	Nu[shiftval + ix * 3] = BLpos.x;
	Nu[shiftval + ix * 3 + 1] = BLpos.y;
	Nu[shiftval + ix * 3 + 2] = fabs(NU_MULTIPLIER * dTdn / Tm);
}

int bcFindIntersectionNusselt(double *dist, double2 vLd, double2 vPL, double2 vN)
{
	double den = dot(vN, vLd);
	if (den == 0.)
		return FALSE;

	*dist = dot(vN, vPL) / den;
	
	if ((*dist) >= 0. && (*dist) <= 1.)
	{
		return TRUE;
	}
	return FALSE;
}

double Get_dT(double T, double Tn, double dX_cur, double dX0)
{
	double Ts = (dX_cur == dX0) ? Tn : (KSOOT * Tn / (dX0 - dX_cur) + KAIR*T / dX_cur) / (KAIR / dX_cur + KSOOT / (dX0 - dX_cur));
	return (Ts - T) / dX_cur;
}


//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
//void FD_calc_Nu_Orig(__global double *T,
//__global double *Nu,
//__global double2 *J,
//__global double2 *LBdist_bot,
//__global double2 *LBdist_top,
//__global bLinks *BL,
//int SaveLoc,
//int Numpts,
//__global int2 *BLinds,
//int NodeX_offset,
//__global int2 *LBnodey,
//__global double *dX_cur,
//__global double *dX0,
//__global int *Stor)
//{
//	int ix = get_global_id(0);
//
//	if (ix >= Numpts)
//		return;
//
//	int shiftval = SaveLoc * Numpts * 5;
//
//	int ii = ix + NodeX_offset;
//
//	double Um = 0.;
//	double Tm = 0.;
//	for (int j = 0; j < FULLSIZEY; j++)
//	{
//		Um += J[ii*FULLSIZEY_UT + j].x;
//		Tm += J[ii*FULLSIZEY_UT + j].x * T[ii*FULLSIZEY_UT + j];
//	}
//
//	Tm /= Um;
//	int2 ibl = BLinds[ix];
//	double2 dTdn = 0.;
//
//	int jj = LBnodey[ix].x;
//	double2 Ts = 0.;
//	if (jj != -1)
//	{
//		double2 LBdist = LBdist_bot[ix];
//		int xdir = (LBdist.x <= 0.) ? 1 : -1;
//		int ydir = (LBdist.y <= 0.) ? 1 : -1;
//		int dXind = (LBdist.x <= 0.) ? 0 : 1;
//		int dYind = (LBdist.y <= 0.) ? 2 : 3;
//
//
//		int Cyred = Stor[ii*DOMAIN_SIZE_Y + jj];
//		int nXyred = Stor[(ii + xdir) * DOMAIN_SIZE_Y + jj];
//		int nYyred = Stor[(ii)* DOMAIN_SIZE_Y + jj + ydir];
//
//		dXind += 4 * (ii*FULLSIZEY + Cyred);
//		dYind += 4 * (ii*FULLSIZEY + Cyred);
//
//		double Tc = T[ii*FULLSIZEY_UT + Cyred];
//		double Tnx = (nXyred == -1) ? 0. : T[(ii + xdir)*FULLSIZEY_UT + nXyred];
//		double Tny = (nYyred == -1) ? 0. : T[(ii)*FULLSIZEY_UT + nYyred];
//
//		double2 dT = 0.;
//		if (LBdist.x != 0.)
//			dT.x = (double)xdir * Get_dT(Tc, Tnx, dX_cur[dXind], dX0[dXind]);
//		if (LBdist.y != 0.)
//			dT.y = (double)ydir * Get_dT(Tc, Tny, dX_cur[dYind], dX0[dYind]);
//
//		dTdn.x = dot(dT, BL[ibl.x].vNvec);
//		Ts.x = MAX(Tc - dTdn.x*length(LBdist), 0.);
//
//
//
//	}
//
//	jj = LBnodey[ix].y;
//	if (jj != -1)
//	{
//		double2 LBdist = LBdist_top[ix];
//		int xdir = (LBdist.x <= 0.) ? 1 : -1;
//		int ydir = (LBdist.y <= 0.) ? 1 : -1;
//		int dXind = (LBdist.x <= 0.) ? 0 : 1;
//		int dYind = (LBdist.y <= 0.) ? 2 : 3;
//
//
//		int Cyred = Stor[ii*DOMAIN_SIZE_Y + jj];
//		int nXyred = Stor[(ii + xdir) * DOMAIN_SIZE_Y + jj];
//		int nYyred = Stor[(ii)* DOMAIN_SIZE_Y + jj + ydir];
//
//		dXind += 4 * (ii*FULLSIZEY + Cyred);
//		dYind += 4 * (ii*FULLSIZEY + Cyred);
//
//		double Tc = T[ii*FULLSIZEY_UT + Cyred];
//		double Tnx = (nXyred == -1) ? 0. : T[(ii + xdir)*FULLSIZEY_UT + nXyred];
//		double Tny = (nYyred == -1) ? 0. : T[(ii)*FULLSIZEY_UT + nYyred];
//
//		double2 dT = 0.;
//		if (LBdist.x != 0.)
//			dT.x = (double)xdir * Get_dT(Tc, Tnx, dX_cur[dXind], dX0[dXind]);
//		if (LBdist.y != 0.)
//			dT.y = (double)ydir * Get_dT(Tc, Tny, dX_cur[dYind], dX0[dYind]);
//
//
//		dTdn.y = dot(dT, BL[ibl.y].vNvec);
//		Ts.y = MAX(Tc - dTdn.y*length(LBdist), 0.);
//
//	}
//
//	Nu[shiftval + ix * 5] = dTdn.x;
//	Nu[shiftval + ix * 5 + 1] = Ts.x;
//	Nu[shiftval + ix * 5 + 2] = dTdn.y;
//	Nu[shiftval + ix * 5 + 3] = Ts.y;
//	Nu[shiftval + ix * 5 + 4] = Tm;
//
//}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
void FD_calc_Nu_Orig(__global double *T,
__global double *Nu,
__global double2 *J,
__global double2 *LBdist_bot,
__global double2 *LBdist_top,
__global bLinks *BL,
int SaveLoc,
int Numpts,
__global int2 *BLinds,
int NodeX_offset,
__global int2 *LBnodey,
__global double *dX_cur,
__global double *dX0,
__global int *Stor)
{
	int ix = get_global_id(0);

	if (ix >= Numpts)
		return;

	int shiftval = SaveLoc * Numpts * 5;

	int ii = ix + NodeX_offset;

	double Um = 0.;
	double Tm = 0.;
	for (int j = 0; j < FULLSIZEY; j++)
	{
		Um += J[ii*FULLSIZEY_UT + j].x;
		Tm += J[ii*FULLSIZEY_UT + j].x * T[ii*FULLSIZEY_UT + j];
	}

	Tm /= Um;
	int2 ibl = BLinds[ix];
	double2 dTdn = 0.;
	int jj = LBnodey[ix].x;
	double2 Ts = 0.;

	if (jj != -1)
	{
		double2 LBdist = LBdist_bot[ix];
		int xdir = (LBdist.x <= 0.) ? 1 : -1;
		int ydir = (LBdist.y <= 0.) ? 1 : -1;
		int dXind = (LBdist.x <= 0.) ? 0 : 1;
		int dYind = (LBdist.y <= 0.) ? 2 : 3;


		int Cyred = Stor[ii*DOMAIN_SIZE_Y + jj];
		int nXyred = Stor[(ii + xdir) * DOMAIN_SIZE_Y + jj];
		int nYyred = Stor[(ii)* DOMAIN_SIZE_Y + jj + ydir];

		dXind += 4 * (ii*FULLSIZEY + Cyred);
		dYind += 4 * (ii*FULLSIZEY + Cyred);

		double Tc = T[ii*FULLSIZEY_UT + Cyred];
		double Tnx = (nXyred == -1) ? 0. : T[(ii + xdir)*FULLSIZEY_UT + nXyred];
		double Tny = (nYyred == -1) ? 0. : T[(ii)*FULLSIZEY_UT + nYyred];
		double2 dT_e = (double2)((double)xdir*(Tnx - Tc) / dX0[dXind], (double)ydir*(Tny - Tc) / dX0[dYind]);
		dTdn.x = dot(dT_e, BL[ibl.x].vNvec);
		Ts.x = MAX(Tc - dTdn.x*length(LBdist), 0.);
	}

	jj = LBnodey[ix].y;
	if (jj != -1)
	{
		double2 LBdist = LBdist_top[ix];
		int xdir = (LBdist.x <= 0.) ? 1 : -1;
		int ydir = (LBdist.y <= 0.) ? 1 : -1;
		int dXind = (LBdist.x <= 0.) ? 0 : 1;
		int dYind = (LBdist.y <= 0.) ? 2 : 3;


		int Cyred = Stor[ii*DOMAIN_SIZE_Y + jj];
		int nXyred = Stor[(ii + xdir) * DOMAIN_SIZE_Y + jj];
		int nYyred = Stor[(ii)* DOMAIN_SIZE_Y + jj + ydir];

		dXind += 4 * (ii*FULLSIZEY + Cyred);
		dYind += 4 * (ii*FULLSIZEY + Cyred);

		double Tc = T[ii*FULLSIZEY_UT + Cyred];
		double Tnx = (nXyred == -1) ? 0. : T[(ii + xdir)*FULLSIZEY_UT + nXyred];
		double Tny = (nYyred == -1) ? 0. : T[(ii)*FULLSIZEY_UT + nYyred];
		double2 dT_e = (double2)((double)xdir*(Tnx - Tc) / dX0[dXind], (double)ydir*(Tny - Tc) / dX0[dYind]);

		dTdn.y = dot(dT_e, BL[ibl.y].vNvec);
		Ts.y = MAX(Tc - dTdn.y*length(LBdist), 0.);

	}

	Nu[shiftval + ix * 5] = dTdn.x;
	Nu[shiftval + ix * 5 + 1] = Ts.x;
	Nu[shiftval + ix * 5 + 2] = dTdn.y;
	Nu[shiftval + ix * 5 + 3] = Ts.y;
	Nu[shiftval + ix * 5 + 4] = Tm;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
	void TR_shear_save_out(__global int *BLind,
	__global bLinks *BL_array,
	__global double *SSout,
	int NUMBLINKS, int save_loc)
{
	uint ii = get_global_id(0);
	if (ii >= NUMBLINKS)
		return;
	int jj = BLind[ii];
	bLinks BLtemp = BL_array[jj];
	double Tautemp = BLtemp.Tau * convert_double(BLtemp.dir);
	SSout[save_loc*NUMBLINKS + ii] = Tautemp;
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_DISTS, 1, 1)))
void TR_save_dists(__global int2 *Ploc_array,
__global par *P,
__global uint *IOinds,
__global int *Neighs,
__global uint *Output,
uint2 indicies, int save_loc, int num_sizes)
{
	int i = get_global_id(0);


	if (i >= indicies.y)
	{
		return;
	}

	int save_start = (i >= indicies.x) ? (save_loc*num_sizes*2 + num_sizes) : (save_loc*num_sizes*2);

	int zz = IOinds[i];
	int j = zz * 9;

	int ii = Neighs[j];

	int jj = j;

	while (jj < j + 9)
	{
		int Node_loc = Neighs[jj];
		int2 ploc = Ploc_array[Node_loc + 2];
		jj++;

		if (Node_loc == -1 || ploc.x == -1)
			continue;

		int kkt = ploc.x;
		while (kkt < ploc.y)
		{
			int kk = kkt;
			kkt++;
			par Pcur = P[kk];

			if (Pcur.loc != ii)
				continue;

			int ind_add = Pcur.type + save_start;
			atomic_add(&Output[ind_add], Pcur.Num_rep);
		}
	}
}


#define GET_GLOBAL_IDX(gx,gy)	(gx + gy*CHANNEL_LENGTH_FULL)


inline void compute_0th_moment(const FTYPE *Fi, FTYPE *out)
{
	*out = Fi[0] + Fi[E_DIR] + Fi[W_DIR] + Fi[N_DIR] + Fi[S_DIR] + Fi[NE_DIR] + Fi[SW_DIR] + Fi[SE_DIR] + Fi[NW_DIR];
}

inline void compute_1st_moment(const FTYPE *Fi, FTYPE *out, const FTYPE factor)
{
	out[0] = factor*(Fi[E_DIR] - Fi[W_DIR] +
		Fi[NE_DIR] - Fi[SW_DIR] + Fi[SE_DIR] - Fi[NW_DIR]);
	
	out[1] = factor*(Fi[N_DIR] - Fi[S_DIR] +
		Fi[NE_DIR] - Fi[SW_DIR] - Fi[SE_DIR] + Fi[NW_DIR]);
}

inline void compute_macro_quant(const FTYPE *fi, FTYPE *rho, FTYPE *v)
{
	compute_0th_moment(fi, rho);
	compute_1st_moment(fi, v, 1. / (*rho));
}

__kernel void Calc_Ro_J(__global double *__restrict__ rho_array,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ fvals,
	__global const int *__restrict__ Map)
{
	const int gx = get_global_id(0);			
	
	if (gx > CHANNEL_LENGTH) { return; }
	const int gy = get_global_id(1) + 1;
	int gi = GET_GLOBAL_IDX(gx, gy);
	const int type = Map[gi];					
	if ((type & FLUID_NODE) == 0) { return; }	

	double Fi[9];

#define	GET_DIST_MACRO(dir)				Fi[(dir ## _DIR)] = fvals[gi + DIST_SIZE * (dir ## _DIR)];

	CALL_ALL_DISTS(GET_DIST_MACRO);
#undef GET_DIST_MACRO

	double v0[2], rho;
	compute_macro_quant(Fi, &rho, v0);

	rho_array[gi] = rho;
	ivx[gi] = v0[0];
	ivy[gi] = v0[1];
}

inline void postCollisionBC(const FTYPE *Fi, const int orientation, const int gi, __global FTYPE *__restrict__ dist_out)
{
#define TEST_WRITE(dir)		if(orientation & (dir ## _BOUND)) { dist_out[gi + DIST_SIZE * (dir ## _REV)] = Fi[(dir ## _DIR)]; } 
	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE
}



__kernel void LB_collision_SRT_Fluid(__global double *__restrict__ rho_array,
	__global double *__restrict__ ovx,
	__global double *__restrict__ ovy,
	__global double *__restrict__ FA,
	__global double *__restrict__ FB,
	__global int *__restrict__ Map,
	__global double *__restrict__ dxvals,
	int options)
{

	const int gx = get_global_id(0);

	if (gx > CHANNEL_LENGTH) { return; }
	const int gy = get_global_id(1);
	const int gi = GET_GLOBAL_IDX(gx, gy);
	const int type = Map[gi];
	if ((type & FLUID_NODE) == 0) { return; }

	double Fi[9];

#define	GET_DIST_MACRO(dir)				Fi[(dir ## _DIR)] = FA[gi + DIST_SIZE * (dir ## _DIR)];

	CALL_ALL_DISTS(GET_DIST_MACRO);
#undef GET_DIST_MACRO

	
	// IBB for q < 0.5
#define TEST_WRITE(dir)		if(type & (dir ## _BOUND_T1)) { \
								double q2 = 2.*dxvals[gi + DIST_SIZE * ((dir ## _DIR) - 1)]; \
								Fi[(dir ## _REV)] = q2 * Fi[(dir ## _REV)] + (1. - q2) * Fi[(dir ## _DIR)]; \
							}

	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE

	
	double v0[2], rho;
	compute_macro_quant(Fi, &rho, v0);

	if ((options & OPTION_SAVE_MACRO_FIELDS))
	{
		rho_array[gi] = rho;
		ovx[gi] = v0[0] + FTERM_VAL / 2.;
		ovy[gi] = v0[1];
	}
	FTYPE Fval = v0[0] + FTERM_VAL;
	FTYPE xi_0 = rho*tau0;
	FTYPE xi_1 = (Fval)*(Fval)*rho;
	FTYPE xi_2 = (v0[0])*(v0[0]);
	FTYPE xi_3 = rho*xi_2;
	FTYPE xi_4 = rho*tau0*xi_2;
	FTYPE xi_5 = (v0[1])*(v0[1]);
	FTYPE xi_6 = rho*tau0*xi_5;
	FTYPE xi_7 = xi_1*xi_5;
	FTYPE xi_8 = xi_3*xi_5;
	FTYPE xi_9 = xi_4*xi_5;
	FTYPE xi_10 = Fval*rho;
	FTYPE xi_11 = (1. / 3.)*xi_10;
	FTYPE xi_12 = v0[0] * rho;
	FTYPE xi_13 = (1. / 3.)*xi_12;
	FTYPE xi_14 = tau0*xi_13;
	FTYPE xi_15 = Fval*rho*xi_5;
	FTYPE xi_16 = 0.5*xi_15;
	FTYPE xi_17 = v0[0] * rho*xi_5;
	FTYPE xi_18 = 0.5*xi_17;
	FTYPE xi_19 = tau0*xi_18;
	FTYPE xi_20 = (1.0 / 9.0)*xi_0;
	FTYPE xi_21 = -0.5*xi_7;
	FTYPE xi_22 = 0.5*xi_8;
	FTYPE xi_23 = -0.5*xi_9;
	FTYPE xi_24 = (1. / 3.)*xi_1 + xi_20 + xi_21 + xi_22 + xi_23 - (1. / 3.)*xi_3 + (1. / 3.)*xi_4 - (1.0 / 6.0)*xi_6;
	FTYPE xi_25 = v0[1] * rho*tau0;
	FTYPE xi_26 = (1. / 3.)*xi_25;
	FTYPE xi_27 = 0.5*v0[1];
	FTYPE xi_28 = xi_1*xi_27;
	FTYPE xi_29 = xi_27*xi_3;
	FTYPE xi_30 = 0.5*xi_2*xi_25;
	FTYPE xi_31 = -(1.0 / 6.0)*xi_1 + xi_20 + xi_21 + xi_22 + xi_23 + (1.0 / 6.0)*xi_3 - (1.0 / 6.0)*xi_4 + (1. / 3.)*xi_6;
	FTYPE xi_32 = (1.0 / 12.0)*xi_10;
	FTYPE xi_33 = (1.0 / 12.0)*xi_12;
	FTYPE xi_35 = 0.25*v0[1];
	FTYPE xi_36 = xi_10*xi_35;
	FTYPE xi_37 = xi_12*xi_35;
	FTYPE xi_39 = tau0*xi_33;
	FTYPE xi_40 = tau0*xi_37;
	FTYPE xi_41 = 0.25*xi_15;
	FTYPE xi_42 = 0.25*xi_17;
	FTYPE xi_44 = tau0*xi_42;
	FTYPE xi_45 = xi_0 / 36.;
	FTYPE xi_46 = (1.0 / 12.0)*xi_25;
	FTYPE xi_47 = (1.0 / 12.0)*xi_1;
	FTYPE xi_48 = -(1.0 / 12.0)*xi_3;
	FTYPE xi_49 = xi_1*xi_35;
	FTYPE xi_50 = xi_3*xi_35;
	FTYPE xi_51 = (1.0 / 12.0)*xi_4;
	FTYPE xi_52 = (1.0 / 12.0)*xi_6;
	FTYPE xi_53 = xi_35*xi_4;
	FTYPE xi_54 = 0.25*xi_7;
	FTYPE xi_55 = -0.25*xi_8;
	FTYPE xi_56 = 0.25*xi_9;
	FTYPE xi_57 = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 - xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
	Fi[0] += -Fi[0] * tau0 + (4.0 / 9.0)*xi_0 - (2. / 3.)*xi_1 + (2. / 3.)*xi_3 - (2. / 3.)*xi_4 - (2. / 3.)*xi_6 + xi_7 - xi_8 + xi_9;
	Fi[1] += -Fi[1] * tau0 + xi_11 - xi_13 + xi_14 - xi_16 + xi_18 - xi_19 + xi_24;
	Fi[2] += -Fi[2] * tau0 - xi_11 + xi_13 - xi_14 + xi_16 - xi_18 + xi_19 + xi_24;
	Fi[3] += -Fi[3] * tau0 + xi_26 - xi_28 + xi_29 - xi_30 + xi_31;
	Fi[4] += -Fi[4] * tau0 - xi_26 + xi_28 - xi_29 + xi_30 + xi_31;
	Fi[5] += -Fi[5] * tau0 + xi_32 - xi_33 + xi_36 - xi_37 + xi_39 + xi_40 + xi_41 - xi_42 + xi_44 + xi_57;
	Fi[6] += -Fi[6] * tau0 - xi_32 + xi_33 + xi_36 - xi_37 - xi_39 + xi_40 - xi_41 + xi_42 - xi_44 + xi_45 - xi_46 + xi_47 + xi_48 - xi_49 + xi_50 + xi_51 + xi_52 - xi_53 + xi_54 + xi_55 + xi_56;
	Fi[7] += -Fi[7] * tau0 + xi_32 - xi_33 - xi_36 + xi_37 + xi_39 - xi_40 + xi_41 - xi_42 + xi_44 + xi_45 - xi_46 + xi_47 + xi_48 - xi_49 + xi_50 + xi_51 + xi_52 - xi_53 + xi_54 + xi_55 + xi_56;
	Fi[8] += -Fi[8] * tau0 - xi_32 + xi_33 - xi_36 + xi_37 - xi_39 - xi_40 - xi_41 + xi_42 - xi_44 + xi_57;

	
#define TEST_WRITE(dir)		if(type & (dir ## _BOUND_T2)) { \
								double q2 = 0.5 / dxvals[gi + DIST_SIZE * ((dir ## _DIR) - 1)]; \
								Fi[(dir ## _DIR)] = q2 * Fi[(dir ## _DIR)] + (1. - q2) * Fi[(dir ## _REV)]; \
							}

	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE


	
	if (type & BOUNDARY_NODE)
	{
		postCollisionBC(Fi, type, gi, FB);
	}


	FB[gi] = Fi[0];

	int eneigh = (gx < CHANNEL_LENGTH - 1) ? 1 : -gx;
	int wneigh = (gx > 0) ? -1 : (CHANNEL_LENGTH - 1);
	
	if ((type & E_BOUND) == 0)
		FB[gi + eneigh + DIST_SIZE] = Fi[1];
	
	if ((type & W_BOUND) == 0)
		FB[gi + wneigh + DIST_SIZE*2] = Fi[2];
	
	if ((type & N_BOUND) == 0)
		FB[gi + N_NEIGH + DIST_SIZE*3] = Fi[3];
	
	if ((type & S_BOUND) == 0)
		FB[gi + S_NEIGH + DIST_SIZE * 4] = Fi[4];
		
	if ((type & NE_BOUND) == 0)
		FB[gi + N_NEIGH + eneigh + DIST_SIZE * 5] = Fi[5];
	
	if ((type & SW_BOUND) == 0)
		FB[gi + S_NEIGH + wneigh + DIST_SIZE * 6] = Fi[6];
	
	if ((type & SE_BOUND) == 0)
			FB[gi + S_NEIGH + eneigh + DIST_SIZE * 7] = Fi[7];
	
	if ((type & NW_BOUND) == 0)
		FB[gi + N_NEIGH + wneigh + DIST_SIZE * 8] = Fi[8];
	//int lx = get_local_id(0);

	//prop_fE[lx] = -1.0;

	//barrier(CLK_LOCAL_MEM_FENCE);
	//// Update the 0-th direction distribution
	//FB[gi] = Fi[C_DIR];
	//// Propagation in directions orthogonal to the X axis (global memory)
	//if ((type & N_BOUND) == 0)
	//	FB[gi + (DIST_SIZE * N_DIR + N_NEIGH)] = Fi[N_DIR];
	//if ((type & S_BOUND) == 0)
	//	FB[gi + (DIST_SIZE * S_DIR + S_NEIGH)] = Fi[S_DIR];

	//// E propagation in shared memory
	//if (lx < (BLOCK_SIZE - 1) && gx < (CHANNEL_LENGTH - 1))
	//{
	//	prop_fE[lx + 1] = Fi[E_DIR];
	//	prop_fNE[lx + 1] = Fi[NE_DIR];
	//	prop_fSE[lx + 1] = Fi[SE_DIR];
	//	// E propagation in global memory (at right block boundary)
	//}
	//else
	//{
	//	int xprop_loc = (gx == CHANNEL_LENGTH - 1) ? (0) : (gx + 1);
	//	int Eprop = GET_GLOBAL_IDX(xprop_loc, gy);

	//	FB[Eprop + DIST_SIZE] = Fi[E_DIR];
	//	if ((type & N_BOUND) == 0)
	//	{
	//		int NEprop = GET_GLOBAL_IDX(xprop_loc, gy + 1);
	//		FB[NEprop + DIST_SIZE * NE_DIR] = Fi[NE_DIR];
	//	}
	//	if ((type & S_BOUND) == 0)
	//	{
	//		int SEprop = GET_GLOBAL_IDX(xprop_loc, gy - 1);
	//		FB[SEprop + DIST_SIZE * SE_DIR] = Fi[SE_DIR];
	//	}
	//}

	//barrier(CLK_LOCAL_MEM_FENCE);

	//// Save locally propagated distributions into global memory.
	//// The leftmost thread is not updated in this block.
	//if (lx > 0 && gx < CHANNEL_LENGTH)
	//{
	//	if (prop_fE[lx] != -1.0)
	//	{
	//		FB[gi + (DIST_SIZE)*E_DIR] = prop_fE[lx];

	//		if ((type & N_BOUND) == 0)
	//			FB[gi + (DIST_SIZE * NE_DIR + N_NEIGH)] = prop_fNE[lx];
	//		if ((type & S_BOUND) == 0)
	//			FB[gi + (DIST_SIZE * SE_DIR + S_NEIGH)] = prop_fSE[lx];
	//	}
	//}

	//barrier(CLK_LOCAL_MEM_FENCE);

	//prop_fE[lx] = -1.0;
	//barrier(CLK_LOCAL_MEM_FENCE);

	//// W propagation in shared memory
	//// Note: propagation to ghost nodes is done directly in global memory as there
	//// are no threads running for the ghost nodes.
	//if (lx > 0 && gx > 0)
	//{
	//	prop_fW[lx - 1] = Fi[W_DIR];
	//	prop_fNW[lx - 1] = Fi[NW_DIR];
	//	prop_fSW[lx - 1] = Fi[SW_DIR];

	//	// W propagation in global memory (at left block boundary)
	//}
	//else
	//{
	//	int xprop_loc = (gx > 0) ? (gx - 1) : (CHANNEL_LENGTH - 1);
	//	int Wprop = GET_GLOBAL_IDX(xprop_loc, gy);

	//	FB[Wprop + DIST_SIZE * W_DIR] = Fi[W_DIR];
	//	if ((type & N_BOUND) == 0)
	//	{
	//		int NWprop = GET_GLOBAL_IDX(xprop_loc, gy + 1);
	//		FB[NWprop + DIST_SIZE * NW_DIR] = Fi[NW_DIR];
	//	}
	//	if ((type & S_BOUND) == 0)
	//	{
	//		int SWprop = GET_GLOBAL_IDX(xprop_loc, gy - 1);
	//		FB[SWprop + DIST_SIZE * SW_DIR] = Fi[SW_DIR];
	//	}
	//}

	//barrier(CLK_LOCAL_MEM_FENCE);
	//// The rightmost thread is not updated in this block.
	//if (lx < (BLOCK_SIZE - 1) && gx < CHANNEL_LENGTH-1)
	//{
	//	if (prop_fE[lx] != -1.0)
	//	{
	//		FB[gi + (DIST_SIZE * W_DIR)] = prop_fW[lx];
	//		if ((type & N_BOUND) == 0)
	//		{
	//			FB[gi + (DIST_SIZE * NW_DIR + N_NEIGH)] = prop_fNW[lx];
	//		}
	//		if ((type & S_BOUND) == 0)
	//		{
	//			FB[gi + (DIST_SIZE * SW_DIR + S_NEIGH)] = prop_fSW[lx];
	//		}
	//	}
	//}
}
#define DBGARR(dbgnum, dbgval)	DebugArrs[gi + dbgnum*DIST_SIZE] = dbgval

__kernel void LB_collision_SRT_Fluid_w_kOmega(__global double *__restrict__ rho_array,
	__global double *__restrict__ ovx,
	__global double *__restrict__ ovy,
	__global double *__restrict__ FA,
	__global double *__restrict__ FB,
	__global int *__restrict__ Map,
	__global double *__restrict__ dxvals,
	__global double *__restrict__ wallD,
	__global double *__restrict__ iNut,
	__global double *__restrict__ iK,
	__global double *__restrict__ iO,
	__global double *__restrict__ iSxy,
	int options)
{

	const int gx = get_global_id(0);

	if (gx > CHANNEL_LENGTH) { return; }
	const int gy = get_global_id(1);
	const int gi = GET_GLOBAL_IDX(gx, gy);
	const int type = Map[gi];
	if ((type & FLUID_NODE) == 0) { return; }



	double Fi[9];

#define	GET_DIST_MACRO(dir)				Fi[(dir ## _DIR)] = FA[gi + DIST_SIZE * (dir ## _DIR)];

	CALL_ALL_DISTS(GET_DIST_MACRO);
#undef GET_DIST_MACRO


	// IBB for q < 0.5
#define TEST_WRITE(dir)		if(type & (dir ## _BOUND_T1)) { \
								double q2 = 2.*dxvals[gi + DIST_SIZE * ((dir ## _DIR) - 1)]; \
								Fi[(dir ## _REV)] = q2 * Fi[(dir ## _REV)] + (1. - q2) * Fi[(dir ## _DIR)]; \
								}

	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE


	double v0[2], rho;
	compute_macro_quant(Fi, &rho, v0);

	double omegaval = iO[gi], kval = iK[gi], yval = wallD[gi];

	double F2arg = max(2.*sqrt(kval) / (KO_BETA_STAR*omegaval*yval), 500.*MU_NUMBER / (yval*yval*omegaval));
	double F2val = tanh(F2arg*F2arg); 

	double sxx_val = iSxy[gi], sxy_val = iSxy[gi + DIST_SIZE], syy_val = iSxy[gi + DIST_SIZE*2];
	F2val *= sqrt(2.*(sxx_val*sxx_val + 2.* sxy_val*sxy_val + syy_val*syy_val));

	//double nut_old = iNut[gi];
	double nut = KO_A1*kval / max(KO_A1*omegaval, F2val);
	//if ((type & (N_BOUND)) || (type & S_BOUND))
	//	nut = NUT_WALL_VALUE;
	iNut[gi] = nut;


#ifndef DISABLE_TURBULENT_VISC
	double tau_t = (1. / ((3.*(nut + MU_NUMBER)) + 0.5));
#else
	double tau_t = tau0;
#endif


	if ((options & OPTION_SAVE_MACRO_FIELDS))
	{
		rho_array[gi] = rho;
		ovx[gi] = v0[0] + FTERM_VAL / 2.;
		ovy[gi] = v0[1];
	}

	double Fval = v0[0] + FTERM_VAL;	
	double xi_0 = rho*tau_t;
	double xi_1 = (Fval)*(Fval);
	double xi_2 = rho*xi_1;
	double xi_3 = (v0[0])*(v0[0]);
	double xi_4 = rho*xi_3;
	double xi_5 = (2. / 3.)*rho*tau_t;
	double xi_6 = (v0[1])*(v0[1]);
	double xi_7 = xi_2*xi_6;
	double xi_8 = xi_4*xi_6;
	double xi_9 = v0[0] * rho;
	double xi_10 = (1. / 3.)*xi_9;
	double xi_11 = v0[0] * rho*xi_6;
	double xi_12 = 0.5*xi_11;
	double xi_13 = (1.0 / 9.0)*rho;
	double xi_14 = rho*xi_6;
	double xi_15 = -(1.0 / 6.0)*xi_14;
	double xi_16 = -0.5*xi_8;
	double xi_17 = xi_13 + xi_15 + xi_16 + (1. / 3.)*xi_4;
	double xi_18 = Fval*rho;
	double xi_19 = (1. / 3.)*xi_18;
	double xi_20 = Fval*rho*xi_6;
	double xi_21 = 0.5*xi_20;
	double xi_22 = -0.5*xi_7;
	double xi_23 = xi_13 + xi_15 + (1. / 3.)*xi_2 + xi_22;
	double xi_24 = v0[1] * rho;
	double xi_25 = (1. / 3.)*xi_24;
	double xi_26 = v0[1] * rho*xi_3;
	double xi_27 = 0.5*xi_26;
	double xi_28 = (1. / 3.)*xi_14;
	double xi_29 = xi_13 + xi_16 + xi_28 - (1.0 / 6.0)*xi_4;
	double xi_30 = v0[1] * rho*xi_1;
	double xi_31 = 0.5*xi_30;
	double xi_32 = xi_13 - (1.0 / 6.0)*xi_2 + xi_22 + xi_28;
	double xi_34 = (1.0 / 12.0)*xi_9;
	double xi_35 = 0.25*v0[1];
	double xi_36 = xi_35*xi_9;
	double xi_37 = 0.25*xi_11;
	double xi_38 = 0.0277777777777778*rho;
	double xi_39 = (1.0 / 12.0)*xi_24;
	double xi_40 = (1.0 / 12.0)*xi_4;
	double xi_41 = (1.0 / 12.0)*xi_14;
	double xi_42 = 0.25*xi_26;
	double xi_43 = 0.25*xi_8;
	double xi_44 = xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43;
	double xi_45 = (1.0 / 12.0)*xi_18;
	double xi_46 = xi_18*xi_35;
	double xi_47 = 0.25*xi_20;
	double xi_48 = (1.0 / 12.0)*xi_2;
	double xi_49 = 0.25*xi_30;
	double xi_50 = 0.25*xi_7;
	double xi_51 = xi_38 + xi_39 + xi_41 + xi_48 + xi_49 + xi_50;
	Fi[0] += -Fi[0] * tau_t + xi_0*xi_3*xi_6 + (4.0 / 9.0)*xi_0 - (2. / 3.)*xi_2 - xi_3*xi_5 + (2. / 3.)*xi_4 - xi_5*xi_6 + xi_7 - xi_8;

	double feq = xi_10 - xi_12 + xi_17;
	double sxx_new = Fi[1] - feq;
	Fi[1] += -Fi[1] * tau_t + feq*tau_t - feq + xi_19 - xi_21 + xi_23;


	feq = -xi_10 + xi_12 + xi_17;
	sxx_new += Fi[2] - feq;
	Fi[2] += -Fi[2] * tau_t + feq*tau_t - feq - xi_19 + xi_21 + xi_23;


	feq = xi_25 - xi_27 + xi_29;
	double syy_new = Fi[3] - feq;
	Fi[3] += -Fi[3] * tau_t + feq*tau_t - feq + xi_25 - xi_31 + xi_32;


	feq = -xi_25 + xi_27 + xi_29;
	syy_new += Fi[4] - feq;
	Fi[4] += -Fi[4] * tau_t + feq*tau_t - feq - xi_25 + xi_31 + xi_32;


	feq = xi_34 + xi_36 + xi_37 + xi_44;
	sxx_new += Fi[5] - feq;
	syy_new += Fi[5] - feq;
	double sxy_new = Fi[5] - feq;
	Fi[5] += -Fi[5] * tau_t + feq*tau_t - feq + xi_45 + xi_46 + xi_47 + xi_51;


	feq = -xi_34 + xi_36 - xi_37 + xi_38 - xi_39 + xi_40 + xi_41 - xi_42 + xi_43;
	sxx_new += Fi[6] - feq;
	syy_new += Fi[6] - feq;
	sxy_new += Fi[6] - feq;
	Fi[6] += -Fi[6] * tau_t + feq*tau_t - feq + xi_38 - xi_39 + xi_41 - xi_45 + xi_46 - xi_47 + xi_48 - xi_49 + xi_50;


	feq = xi_34 - xi_36 + xi_37 + xi_38 - xi_39 + xi_40 + xi_41 - xi_42 + xi_43;
	sxx_new += Fi[7] - feq;
	syy_new += Fi[7] - feq;
	sxy_new -= Fi[7] - feq;
	Fi[7] += -Fi[7] * tau_t + feq*tau_t - feq + xi_38 - xi_39 + xi_41 + xi_45 - xi_46 + xi_47 + xi_48 - xi_49 + xi_50;


	feq = -xi_34 - xi_36 - xi_37 + xi_44;
	sxx_new += Fi[8] - feq;
	syy_new += Fi[8] - feq;
	sxy_new -= Fi[8] - feq;
	Fi[8] += -Fi[8] * tau_t + feq*tau_t - feq - xi_45 - xi_46 - xi_47 + xi_51;

	iSxy[gi] = -1.5 * (tau0) * sxx_new;
	iSxy[gi + DIST_SIZE] = -1.5 * (tau0) * sxy_new;
	iSxy[gi + DIST_SIZE*2] = -1.5 * (tau0) * syy_new;



#define TEST_WRITE(dir)		if(type & (dir ## _BOUND_T2)) { \
								double q2 = 0.5 / dxvals[gi + DIST_SIZE * ((dir ## _DIR) - 1)]; \
								Fi[(dir ## _DIR)] = q2 * Fi[(dir ## _DIR)] + (1. - q2) * Fi[(dir ## _REV)]; \
								}

	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE



	if (type & BOUNDARY_NODE)
	{
		postCollisionBC(Fi, type, gi, FB);
	}


	FB[gi] = Fi[0];

	int eneigh = (gx < CHANNEL_LENGTH - 1) ? 1 : -gx;
	int wneigh = (gx > 0) ? -1 : (CHANNEL_LENGTH - 1);

	if ((type & E_BOUND) == 0)
		FB[gi + eneigh + DIST_SIZE] = Fi[1];

	if ((type & W_BOUND) == 0)
		FB[gi + wneigh + DIST_SIZE * 2] = Fi[2];

	if ((type & N_BOUND) == 0)
		FB[gi + N_NEIGH + DIST_SIZE * 3] = Fi[3];

	if ((type & S_BOUND) == 0)
		FB[gi + S_NEIGH + DIST_SIZE * 4] = Fi[4];

	if ((type & NE_BOUND) == 0)
		FB[gi + N_NEIGH + eneigh + DIST_SIZE * 5] = Fi[5];

	if ((type & SW_BOUND) == 0)
		FB[gi + S_NEIGH + wneigh + DIST_SIZE * 6] = Fi[6];

	if ((type & SE_BOUND) == 0)
		FB[gi + S_NEIGH + eneigh + DIST_SIZE * 7] = Fi[7];

	if ((type & NW_BOUND) == 0)
		FB[gi + N_NEIGH + wneigh + DIST_SIZE * 8] = Fi[8];

}


#undef DBGARR

////Interpolated bounce back kernel
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1)))
//void LB_IBB(__global int2 *ibb_loc,
//__global double2 *ibb_coeff,
//__global double *FB,
//int max_el)
//{
//	int i = get_global_id(0);
//
//	if (i >= max_el)
//		return;
//
//	double Fold = FB[ibb_loc[i].x];
//	double Fnew = ibb_coeff[i].x * Fold + ibb_coeff[i].y * FB[ibb_loc[i].y];
//	FB[ibb_loc[i].x] = Fnew;
//}


#define DBGARR(dbgnum, dbgval)	DebugArrs[gid + dbgnum*DIST_SIZE] = dbgval



__kernel void Update_kOmega_Diffusivities_Alphat(__global int *__restrict__ map,
	__global double *__restrict__ iK,
	__global double *__restrict__ iO,
	__global double *__restrict__ iNut,
	__global double *__restrict__ WallD,
	__global double *__restrict__ iDk,
	__global double *__restrict__ iDo,
	__global double *__restrict__ iF1,
	__global double *__restrict__ idkdxy,

	__global double *__restrict__ irho,
	__global double *__restrict__ iAlphat,
#ifndef DEBUG_ARRAYS
	__global double *__restrict__ dXcur)
#else
	__global double *__restrict__ dXcur,
	__global double *__restrict__ DebugArrs)
#endif
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;
	const int type = map[gid];
	if ((type & FLUID_NODE) == 0) { return; }
	
	double Xe_coeff = dXcur[gid], Xw_coeff = dXcur[gid + DIST_SIZE], Xc_coeff = dXcur[gid + DIST_SIZE * 2];
	double Yn_coeff = dXcur[gid + DIST_SIZE * 3], Ys_coeff = dXcur[gid + DIST_SIZE * 4], Yc_coeff = dXcur[gid + DIST_SIZE * 5];

	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double omegaval = iO[gid], kval = iK[gid], yval = WallD[gid];
	
	double dOdx = Xe_coeff*iO[gide] + Xw_coeff*iO[gidw] + Xc_coeff*omegaval;
	double dKdx = Xe_coeff*iK[gide] + Xw_coeff*iK[gidw] + Xc_coeff*kval;
	
	double dOdy = Yn_coeff*iO[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iO[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*omegaval;
	double dKdy = Yn_coeff*iK[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iK[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*kval;
	
	double dodx_mult = 2.*KO_SIGMA_O2 / omegaval;
	double dody_mult = dodx_mult*dKdy;
	dodx_mult *= dKdx;

	idkdxy[gid] = dodx_mult;
	idkdxy[gid+DIST_SIZE] = dody_mult;

	double CDkw = max(irho[gid] * (dodx_mult*dOdx + dody_mult*dOdy), 1.0e-10);
	
#ifdef DEBUG_ARRAYS
	DBGARR(0, dKdx);
	DBGARR(1, dKdy);
	DBGARR(2, dOdx);
	DBGARR(3, dOdy);
	DBGARR(4, CDkw);
#endif
	
	double F1val = min(max(sqrt(kval) / (KO_BETA_STAR*omegaval*yval), 500. * MU_NUMBER / (yval*yval*omegaval)), 4.*kval*KO_SIGMA_O2 / (CDkw*yval*yval));
	F1val = tanh(pown(F1val, 4));
	iF1[gid] = F1val;
	double nutval = iNut[gid];
	double sigma = KO_SIGMA_K1*F1val + KO_SIGMA_K2*(1. - F1val);
	iDk[gid] = MU_NUMBER + sigma*nutval;

	sigma = KO_SIGMA_O1*F1val + KO_SIGMA_O2*(1. - F1val);
	iDo[gid] = MU_NUMBER + sigma*nutval;

	iAlphat[gid] = MU_NUMBER / PR_NUMBER + nutval / PR_TURB_NUMBER;
}


__kernel void Update_kOmega_Diffusivities(__global int *__restrict__ map,
	__global double *__restrict__ iK,
	__global double *__restrict__ iO,
	__global double *__restrict__ iNut,
	__global double *__restrict__ WallD,
	__global double *__restrict__ iDk,
	__global double *__restrict__ iDo,
	__global double *__restrict__ iF1,
	__global double *__restrict__ idkdxy,

	__global double *__restrict__ irho,
#ifndef DEBUG_ARRAYS
	__global double *__restrict__ dXcur)
#else
	__global double *__restrict__ dXcur,
	__global double *__restrict__ DebugArrs)
#endif
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;
	const int type = map[gid];
	if ((type & FLUID_NODE) == 0) { return; }

	double Xe_coeff = dXcur[gid], Xw_coeff = dXcur[gid + DIST_SIZE], Xc_coeff = dXcur[gid + DIST_SIZE * 2];
	double Yn_coeff = dXcur[gid + DIST_SIZE * 3], Ys_coeff = dXcur[gid + DIST_SIZE * 4], Yc_coeff = dXcur[gid + DIST_SIZE * 5];

	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double omegaval = iO[gid], kval = iK[gid], yval = WallD[gid];

	double dOdx = Xe_coeff*iO[gide] + Xw_coeff*iO[gidw] + Xc_coeff*omegaval;
	double dKdx = Xe_coeff*iK[gide] + Xw_coeff*iK[gidw] + Xc_coeff*kval;

	double dOdy = Yn_coeff*iO[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iO[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*omegaval;
	double dKdy = Yn_coeff*iK[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iK[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*kval;

	double dodx_mult = 2.*KO_SIGMA_O2 / omegaval;
	double dody_mult = dodx_mult*dKdy;
	dodx_mult *= dKdx;

	idkdxy[gid] = dodx_mult;
	idkdxy[gid + DIST_SIZE] = dody_mult;

	double CDkw = max(irho[gid] * (dodx_mult*dOdx + dody_mult*dOdy), 1.0e-10);

#ifdef DEBUG_ARRAYS
	DBGARR(0, dKdx);
	DBGARR(1, dKdy);
	DBGARR(2, dOdx);
	DBGARR(3, dOdy);
	DBGARR(4, CDkw);
#endif

	double F1val = min(max(sqrt(kval) / (KO_BETA_STAR*omegaval*yval), 500. * MU_NUMBER / (yval*yval*omegaval)), 4.*kval*KO_SIGMA_O2 / (CDkw*yval*yval));
	F1val = tanh(pown(F1val, 4));
	iF1[gid] = F1val;
	double nutval = iNut[gid];
	double sigma = KO_SIGMA_K1*F1val + KO_SIGMA_K2*(1. - F1val);
	iDk[gid] = MU_NUMBER + sigma*nutval;

	sigma = KO_SIGMA_O1*F1val + KO_SIGMA_O2*(1. - F1val);
	iDo[gid] = MU_NUMBER + sigma*nutval;
}



__kernel void Update_kOmega_Coeffs_Implicit(__global int *__restrict__ IndArr,
	__global double *__restrict__ iK,
	__global double *__restrict__ iO,
	__global double *__restrict__ iNut,
	__global double *__restrict__ iDk,
	__global double *__restrict__ iDo,
	__global double *__restrict__ iF1,
	__global double *__restrict__ idkdo,
	__global double *__restrict__ iSxy,
	__global double *__restrict__ irho,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ kAmat,
	__global double *__restrict__ oAmat,
	__global double *__restrict__ kbvec,
	__global double *__restrict__ obvec,
	__global double *__restrict__ dXcur,
	__global double *__restrict__ WallD,
#ifndef DEBUG_ARRAYS
	double DTFD, double TurbInts, double TurbLScale_Inv)
#else
	double DTFD, double TurbInts, double TurbLScale_Inv,
	__global double *__restrict__ DebugArrs)
#endif
{

	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	int Cind = IndArr[gid];
	if (Cind < 0)
		return;

	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];


	double Xe_coeff = dXcur[gid], Xw_coeff = dXcur[gid + DIST_SIZE], Xc_coeff = dXcur[gid + DIST_SIZE * 2];
	double Yn_coeff = dXcur[gid + DIST_SIZE * 3], Ys_coeff = dXcur[gid + DIST_SIZE * 4], Yc_coeff = dXcur[gid + DIST_SIZE * 5];
	double Xe2_coeff = dXcur[gid + DIST_SIZE * 6], Xw2_coeff = dXcur[gid + DIST_SIZE * 7], Xc2_coeff = dXcur[gid + DIST_SIZE * 8];
	double Yn2_coeff = dXcur[gid + DIST_SIZE * 9], Ys2_coeff = dXcur[gid + DIST_SIZE * 10], Yc2_coeff = dXcur[gid + DIST_SIZE * 11];

	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double diffk_w = iDk[gidw], diffo_w = iDo[gidw]; ///will never be out of bounds because node (0,0) is always solid
	double diffk_e = iDk[gide], diffo_e = iDo[gide];
	double diffk_n = iDk[gid + CHANNEL_LENGTH_FULL], diffo_n = iDo[gid + CHANNEL_LENGTH_FULL];
	double diffk_s = iDk[gid - CHANNEL_LENGTH_FULL], diffo_s = iDo[gid - CHANNEL_LENGTH_FULL];

	double diffk = iDk[gid], diffo = iDo[gid];
	
	double diffk_dx = Xe_coeff * diffk_e + Xw_coeff * diffk_w + Xc_coeff * diffk;
	double diffo_dx = Xe_coeff * diffo_e + Xw_coeff * diffo_w + Xc_coeff * diffo;
	double diffk_dy = Yn_coeff * diffk_n + Ys_coeff * diffk_s + Yc_coeff * diffk;
	double diffo_dy = Yn_coeff * diffo_n + Ys_coeff * diffo_s + Yc_coeff * diffo;



#ifdef DEBUG_ARRAYS
	DBGARR(0, diffk_dx);
	DBGARR(1, diffk_dy);
	DBGARR(2, diffo_dx);
	DBGARR(3, diffo_dy);
#endif




	double Ux = ivx[gid], Uy = ivy[gid];

	double Jx_omega = Ux - diffo_dx - idkdo[gid];
	double Jy_omega = Uy - diffo_dy - idkdo[gid+DIST_SIZE];

	double Jx_k = Ux - diffk_dx;
	double Jy_k = Uy - diffk_dy;

	double omegaval = iO[gid], kval = iK[gid]; 
	double F1val = iF1[gid], nutval = iNut[gid];
	double Sxx = iSxy[gid], Sxy = iSxy[gid + DIST_SIZE], Syy = iSxy[gid + 2 * DIST_SIZE];
	double Pk = 2.*nutval*(Sxx*Sxx + Syy*Syy + Sxy*Sxy) - 2. / 3.*kval * (Sxx + Syy);
	Pk = min(Pk, 10.*KO_BETA_STAR*kval*omegaval);
	double gammaval = (F1val*KO_GAMMA1 + (1. - F1val)*KO_GAMMA2) / nutval;


	double Sc = -KO_BETA_STAR * omegaval;

#ifdef DEBUG_ARRAYS
	DBGARR(4, Jx_omega);
	DBGARR(5, Jy_omega);
	DBGARR(6, Jx_k);
	DBGARR(7, Jy_k);
	DBGARR(8, Pk);
	DBGARR(9, Sc);
#endif
	
	double Kc = 1. + DTFD * (Jx_k*Xc_coeff + Jy_k*Yc_coeff - diffk * (Xc2_coeff + Yc2_coeff)) - Sc;
	double Ke = DTFD * (Jx_k*Xe_coeff - diffk * Xe2_coeff);
	double Kw = DTFD * (Jx_k*Xw_coeff - diffk * Xw2_coeff);
	double Kn = DTFD * (Jy_k*Yn_coeff - diffk * Yn2_coeff);
	double Ks = DTFD * (Jy_k*Ys_coeff - diffk * Ys2_coeff);
	double Ksrc = kval + DTFD*Pk;
	
	double Wc = 1. + DTFD * (Jx_omega*Xc_coeff + Jy_omega*Yc_coeff - diffo * (Xc2_coeff + Yc2_coeff)) - Sc;
	double We = DTFD * (Jx_omega*Xe_coeff - diffo * Xe2_coeff);
	double Ww = DTFD * (Jx_omega*Xw_coeff - diffo * Xw2_coeff);
	double Wn = DTFD * (Jy_omega*Yn_coeff - diffo * Yn2_coeff);
	double Ws = DTFD * (Jy_omega*Ys_coeff - diffo * Ys2_coeff);
	double Wsrc = omegaval + DTFD*(gammaval * Pk);

#ifdef DEBUG_ARRAYS
	DBGARR(10, Kc);
	DBGARR(11, Ke);
	DBGARR(12, Kw);
	DBGARR(13, Kn);
	DBGARR(14, Ks);
	DBGARR(15, Ksrc);
	DBGARR(16, Wc);
	DBGARR(17, We);
	DBGARR(18, Ww);
	DBGARR(19, Wn);
	DBGARR(20, Ws);
	DBGARR(21, Wsrc);
#endif



	if (Sind >= 0 && Wind >= 0 && Eind >= 0 && Nind >= 0)
	{
		kAmat[Sind] = Ks;
		oAmat[Sind] = Ws;
		kAmat[Wind] = Kw;
		oAmat[Wind] = Ww;
		kAmat[Cind] = Kc;
		oAmat[Cind] = Wc;
		kAmat[Eind] = Ke;
		oAmat[Eind] = We;
		kAmat[Nind] = Kn;
		oAmat[Nind] = Wn;

		kbvec[gid] = Ksrc;
		obvec[gid] = Wsrc;
	}
	else
	{
		//if (Sind < 0)
		//{
		//	kAmat[Cind] = Kc + Ks;
		//	kAmat[Nind] = Kn;
		//	oAmat[Nind] = 0.;
		//}
		//else
		//{
		//	oAmat[Sind] = 0.;
		//	kAmat[Sind] = Ks;
		//	kAmat[Cind] = Kc + Kn;
		//}

		//
		//kAmat[Wind] = Kw;
		//oAmat[Wind] = 0.;
		//oAmat[Cind] = 1.;
		//kAmat[Eind] = Ke;
		//oAmat[Eind] = 0.;

		//if (Sind < 0)
		//{
		//	kAmat[Nind] = 1e-3;
		//	oAmat[Nind] = 1e-3;
		//}
		//else
		//{
		//	kAmat[Sind] = 1e-3;
		//	oAmat[Sind] = 1e-3;
		//}

		//kAmat[Wind] = 1e-3;
		//kAmat[Eind] = 1e-3;

		kAmat[Cind] = 1.;
		oAmat[Cind] = 1.;
		
		kbvec[gid] = K_WALL_VALUE;
		obvec[gid] = OMEGA_WALL_VALUE;
	}


}

#undef DBGARR
__kernel void Update_T_Coeffs(__global double *__restrict__ Temp,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ Amat,
	__global double *__restrict__ bvec,
	__global double *__restrict__ dXcur,
	__global int *__restrict__ IndArr,
	double alpha, double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	int Cind = IndArr[gid];
	if (Cind < 0)
		return;
	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];


	double Ux = ivx[gid], Uy = ivy[gid];
	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx = dx_e + dx_w;
	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy = dy_n + dy_s;

	double xi_0 = 1.0 / dx_e;
	double xi_1 = Ux*DTFD;
	double xi_2 = xi_0*xi_1;
	double xi_3 = 1.0 / dy_n;
	double xi_4 = Uy*DTFD;
	double xi_5 = xi_3*xi_4;
	double xi_6 = 1.0 / dx_w;
	double xi_7 = xi_1*xi_6;
	double xi_8 = 1.0 / dy_s;
	double xi_9 = xi_4*xi_8;
	double xi_10 = 2.0*alpha*DTFD;
	double xi_11 = 1.0 / dx;
	double xi_12 = 2.0*alpha*DTFD*xi_11;
	double xi_13 = 1.0 / dy;
	double xi_14 = 2.0*alpha*DTFD*xi_13;
	//Tc = -xi_0*xi_10*xi_6 - xi_10*xi_3*xi_8 + xi_2 + xi_5 - xi_7 - xi_9 + 1.0;
	//Te = -dx_w*xi_11*xi_2 + xi_0*xi_12;
	//Tw = dx_e*xi_11*xi_7 + xi_12*xi_6;
	//Tn = -dy_s*xi_13*xi_5 + xi_14*xi_3;
	//Ts = dy_n*xi_13*xi_9 + xi_14*xi_8;

	if (Sind >= 0)
		Amat[Sind] = dy_n*xi_13*xi_9 + xi_14*xi_8;

	if (Wind >= 0)
		Amat[Wind] = dx_e*xi_11*xi_7 + xi_12*xi_6;
	else if (i == 0)
		bvec[gid] = TFD_X_IN_VAL * (dx_e*xi_11*xi_7 + xi_12*xi_6);

	Amat[Cind] = -xi_0*xi_10*xi_6 - xi_10*xi_3*xi_8 + xi_2 + xi_5 - xi_7 - xi_9 + 1.0;

	if (Eind >= 0)
		Amat[Eind] = -dx_w*xi_11*xi_2 + xi_0*xi_12;
	else if (i == CHANNEL_LENGTH - 1)
		bvec[gid] = Temp[gid] * (-dx_w*xi_11*xi_2 + xi_0*xi_12);

	if (Nind >= 0)
		Amat[Nind] = -dy_s*xi_13*xi_5 + xi_14*xi_3;

}


//__kernel void Update_T_Coeffs_Implicit(__global double *__restrict__ Temp,
//	__global double *__restrict__ ivx,
//	__global double *__restrict__ ivy,
//	__global double *__restrict__ Amat,
//	__global double *__restrict__ bvec,
//	__global double *__restrict__ dXcur,
//	__global int *__restrict__ IndArr,
//	double alpha, double DTFD)
//{
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
//		return;
//
//	int gid = i + CHANNEL_LENGTH_FULL*j;
//
//	int Cind = IndArr[gid];
//	if (Cind < 0)
//		return;
//	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
//	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];
//
//
//	double Ux = ivx[gid], Uy = ivy[gid];
//	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx = dx_e + dx_w;
//	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy = dy_n + dy_s;
//
//	double xi_0 = 1.0 / dx_e;
//	double xi_1 = Ux*DTFD;
//	double xi_2 = xi_0*xi_1;
//	double xi_3 = 1.0 / dy_n;
//	double xi_4 = Uy*DTFD;
//	double xi_5 = xi_3*xi_4;
//	double xi_6 = 1.0 / dx_w;
//	double xi_7 = xi_1*xi_6;
//	double xi_8 = 1.0 / dy_s;
//	double xi_9 = xi_4*xi_8;
//	double xi_10 = 2.0*alpha*DTFD;
//	double xi_11 = 1.0 / dx;
//	double xi_12 = 2.0*alpha*DTFD*xi_11;
//	double xi_13 = 1.0 / dy;
//	double xi_14 = 2.0*alpha*DTFD*xi_13;
//	double Tc = xi_0*xi_10*xi_6 + xi_10*xi_3*xi_8 + xi_2 + xi_5 - xi_7 - xi_9 + 1.0;
//	double Te = dx_w*xi_11*xi_7 - xi_12*xi_6;
//	double Tw = -dx_e*xi_11*xi_2 - xi_0*xi_12;
//	double Tn = dy_s*xi_13*xi_9 - xi_14*xi_8;
//	double Ts = -dy_n*xi_13*xi_5 - xi_14*xi_3;
//	
//	double Tsrc = Temp[gid];
//	
//	if (Sind >= 0)
//		Amat[Sind] = Ts;
//
//	if (Wind >= 0)
//		Amat[Wind] = Tw;
//	else if (i == 0)
//		Tsrc -= TFD_X_IN_VAL * Tw;
//
//	Amat[Cind] = Tc;
//
//	if (Eind >= 0)
//		Amat[Eind] = Te;
//	else if (i == CHANNEL_LENGTH - 1)
//		Amat[Cind] += Te;
//
//	if (Nind >= 0)
//		Amat[Nind] = Tn;
//
//	bvec[gid] = Tsrc;
//
//	//double Tsrc = Temp[gid];
//
//	//if (Sind >= 0)
//	//	Amat[Sind] = DTFD*(-Uy / 2. - alpha);
//
//	//if (Wind >= 0)
//	//	Amat[Wind] = DTFD*(-Ux / 2. - alpha);
//	//else if (i == 0)
//	//	Tsrc += TFD_X_IN_VAL * DTFD*((Ux / 2. + alpha));
//
//	//Amat[Cind] = 1 + 4.*alpha;
//
//	//if (Eind >= 0)
//	//	Amat[Eind] = DTFD*(Ux / 2. - alpha);
//	//else if (i == CHANNEL_LENGTH - 1)
//	//	Amat[Cind] += DTFD*((Ux / 2. - alpha));
//
//	//if (Nind >= 0)
//	//	Amat[Nind] = DTFD*(Uy / 2. - alpha);
//
//	//bvec[gid] = Tsrc;
//
//
//}



__kernel void Update_T_Coeffs_Implicit(__global double *__restrict__ Temp,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ Amat,
	__global double *__restrict__ bvec,
	__global double *__restrict__ dXcur,
	__global int *__restrict__ IndArr,
	double alpha, double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	int Cind = IndArr[gid];
	if (Cind < 0)
		return;
	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];

	
	double dX_e = dXcur[gid], dX_w = dXcur[gid + DIST_SIZE], dX_c = dXcur[gid + DIST_SIZE * 2];
	double dY_n = dXcur[gid + DIST_SIZE * 3], dY_s = dXcur[gid + DIST_SIZE * 4], dY_c = dXcur[gid + DIST_SIZE * 5];
	double dX2_e = dXcur[gid + DIST_SIZE * 6], dX2_w = dXcur[gid + DIST_SIZE * 7], dX2_c = dXcur[gid + DIST_SIZE * 8];
	double dY2_n = dXcur[gid + DIST_SIZE * 9], dY2_s = dXcur[gid + DIST_SIZE * 10], dY2_c = dXcur[gid + DIST_SIZE * 11];

	double Ux = ivx[gid], Uy = ivy[gid];
	

	
	double xi_0 = Ux*DTFD;
	double xi_1 = Uy*DTFD;
	double xi_2 = alpha*DTFD;
	double Tc = -dX2_c*xi_2 + dX_c*xi_0 - dY2_c*xi_2 + dY_c*xi_1 + 1.0;
	double Te = -dX2_e*xi_2 + dX_e*xi_0;
	double Tw = -dX2_w*xi_2 + dX_w*xi_0;
	double Tn = -dY2_n*xi_2 + dY_n*xi_1;
	double Ts = -dY2_s*xi_2 + dY_s*xi_1;


	/// Testing w/ periodic domain and Twall = 1;
	double Tsrc = Temp[gid];

	if (Sind >= 0)
		Amat[Sind] = Ts;
	//else
	//	Tsrc -= Tw;

	if (Wind >= 0)
		Amat[Wind] = Tw;
	else if (i == 0)
		Tsrc -= TFD_X_IN_VAL * Tw;
	Amat[Cind] = Tc;

	if (Eind >= 0)
		Amat[Eind] = Te;
	else if (i == CHANNEL_LENGTH-1)
		Amat[Cind] += Te;

	if (Nind >= 0)
		Amat[Nind] = Tn;
	//else
	//	Tsrc -= Tn;

	bvec[gid] = Tsrc;








	//if (Sind >= 0)
	//	Amat[Sind] = Ts;

	//if (Wind >= 0)
	//	Amat[Wind] = Tw;
	//else if (i == 0)
	//	Tsrc -= TFD_X_IN_VAL * Tw;

	//Amat[Cind] = Tc;

	//if (Eind >= 0)
	//	Amat[Eind] = Te;
	//else if (i == CHANNEL_LENGTH - 1)
	//	Amat[Cind] += Te;

	//if (Nind >= 0)
	//	Amat[Nind] = Tn;

	//bvec[gid] = Tsrc;

}

__kernel void Update_T_Coeffs_Implicit_Turbulent(__global double *__restrict__ Temp,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ Amat,
	__global double *__restrict__ bvec,
	__global double *__restrict__ dXcur,
	__global int *__restrict__ IndArr,
	__global double *__restrict__ iAlphat,
	double alpha, double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	int Cind = IndArr[gid];
	if (Cind < 0)
		return;
	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];


	double dX_e = dXcur[gid], dX_w = dXcur[gid + DIST_SIZE], dX_c = dXcur[gid + DIST_SIZE * 2];
	double dY_n = dXcur[gid + DIST_SIZE * 3], dY_s = dXcur[gid + DIST_SIZE * 4], dY_c = dXcur[gid + DIST_SIZE * 5];
	double dX2_e = dXcur[gid + DIST_SIZE * 6], dX2_w = dXcur[gid + DIST_SIZE * 7], dX2_c = dXcur[gid + DIST_SIZE * 8];
	double dY2_n = dXcur[gid + DIST_SIZE * 9], dY2_s = dXcur[gid + DIST_SIZE * 10], dY2_c = dXcur[gid + DIST_SIZE * 11];

	double Ux = ivx[gid], Uy = ivy[gid], alphat = iAlphat[gid];
	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double dAdx = dX_e*iAlphat[gide] + dX_w*iAlphat[gidw] + dX_c*alphat;
	double dAdy = dY_n*iAlphat[gid+CHANNEL_LENGTH_FULL] + 
		dY_s*iAlphat[gide - CHANNEL_LENGTH_FULL] + dY_c*alphat;

	double Jx = DTFD*(Ux - dAdx);
	double Jy = DTFD*(Uy - dAdy);
	alphat *= DTFD;

	double Tc = 1. + (Jx*dX_c + Jy*dY_c - alphat * (dX2_c + dY2_c));
	double Te = (Jx*dX_e - alphat*dX2_e);
	double Tw = (Jx*dX_w - alphat*dX2_w);
	double Tn = (Jy*dY_n - alphat*dY2_n);
	double Ts = (Jy*dY_s - alphat*dY2_s);
	double Tsrc = Temp[gid];

	if (Sind >= 0)
		Amat[Sind] = Ts;

	if (Wind >= 0)
		Amat[Wind] = Tw;
	else if (i == 0)
		Tsrc -= TFD_X_IN_VAL * Tw;

	Amat[Cind] = Tc;

	if (Eind >= 0)
		Amat[Eind] = Te;
	else if (i == CHANNEL_LENGTH - 1)
		Amat[Cind] += Te;

	if (Nind >= 0)
		Amat[Nind] = Tn;

	bvec[gid] = Tsrc;
}





/*

__kernel void Test_BiCGStab(__global double *__restrict__ Temp,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ Amat,
	__global double *__restrict__ bvec,
	__global double *__restrict__ dXcur,
	__global int *__restrict__ IndArr,
	__global double *__restrict__ iAlphat,
	double alpha, double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	int Cind = IndArr[gid];
	if (Cind < 0)
		return;
	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];

	alpha *= DTFD;

	double Tc = 1. + 4 * alpha;
	double Te = -alpha;
	double Tw = -alpha;
	double Tn = -alpha;
	double Ts = -alpha;
	double Tsrc = Temp[gid];

	if (Sind >= 0)
		Amat[Sind] = Ts;

	if (Wind >= 0)
		Amat[Wind] = Tw;
	else if (i == 0)
		Tsrc -= TFD_X_IN_VAL * Tw;

	Amat[Cind] = Tc;

	if (Eind >= 0)
		Amat[Eind] = Te;

	if (Nind >= 0)
		Amat[Nind] = Tn;

	bvec[gid] = Tsrc;




}*/


