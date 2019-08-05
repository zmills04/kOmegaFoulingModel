#ifndef OPENCL_VERSION_1_2
#pragma OPENCL EXTENSION cl_khr_subgroups : enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics : enable
#endif

#define FTYPE	double

#define MAX(A, B)		(((A) > (B)) ? (A) : (B))
#define MIN(A, B)		(((A) < (B)) ? (A) : (B))
#define EPSILON			1.11e-16
#define true			1
#define false			0
#define EQUALS2D(A, B)	((A.x == B.x) && (A.y == B.y))
#define NEQUALS2D(A, B)	((A.x != B.x) || (A.y != B.y))
#define CEPS			1.0e-5
#define RAND_MAX		4294967296.
#define MODFAST(A, B)	(((A) >= (B)) ? ((A) - (B)) : (((A) < 0) ? ((A) + (B)) : (A)))
#define CX_DEF_HI	(double4)(1.,-1.,1.,-1.)
#define CY_DEF_HI	(double4)(1.,-1.,-1.,1.)
#define CXX_DEF_HI	(double4)(1.,1.,1.,1.)
#define CXY_DEF		(double4)(1.,1.,-1.,-1.)

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


//https://github.com/ddemidov/vexcl-experiments/blob/master/sort-by-key-atomic.cpp
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

bool AtomicMin(volatile __global double* __restrict__ source, const double newValDouble) {
	union {
		double f;
		ulong i;
	} newVal;
	union {
		double f;
		ulong i;
	} prevVal;
	do {
		prevVal.f = *source;
		newVal.f = min(prevVal.f, newValDouble);
	} while (atomic_cmpxchg((volatile __global ulong *)source, prevVal.i, newVal.i) != prevVal.i);


	// returns true if source == newValDouble (i.e. newValDouble was smaller and therefore
	// the new value in source)
	return (*source == newValDouble);

}

////////////////////////////////////////////////
////// Random number generator Functions  //////
////////////////////////////////////////////////

//! Return a 32-bit integer in the range [0..2^32)

uint2 MWC64X_NextUint(uint2 s, double* resd)
{
	uint res = s.x ^ s.y;
	uint X = s.x, C = s.y;

	uint Xn = MWC64X_A * X + C;
	uint carry = (uint)(Xn < C);				// The (Xn<C) will be zero or one for scalar
	uint Cn = mad_hi(MWC64X_A, X, carry);

	s.x = Xn;
	s.y = Cn;
	*resd = convert_double(res) / RAND_MAX;
	return s;
}

uint2 MWC64X_NextUint2(uint2 s, double2* resd)
{
	uint res = s.x ^ s.y;

	(*resd).x = convert_double(res) / RAND_MAX;
	uint X = s.x, C = s.y;

	uint Xn = MWC64X_A * X + C;
	uint carry = (uint)(Xn < C);				// The (Xn<C) will be zero or one for scalar
	uint Cn = mad_hi(MWC64X_A, X, carry);

	res = Xn ^ Cn;
	(*resd).y = convert_double(res) / RAND_MAX;
	Xn = MWC64X_A * X + C;
	carry = (uint)(Xn < C);				// The (Xn<C) will be zero or one for scalar
	Cn = mad_hi(MWC64X_A, X, carry);

	s.x = Xn;
	s.y = Cn;
	return s;
}


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



// Bit flags for vls.nType
//Type 1 Boundaries are IBB boundaries applied before collision step (q < 0.5)
//Type 2 Boundaries are IBB boundaries applied after collision step (q >= 0.5)

#ifndef IN_KERNEL_IBB
#define C_BOUND					0x0000 //just a placeholder. Not actually used
#define M_SOLID_NODE			0x0001
#define M_FLUID_NODE			0x0002
#define M0_SOLID_NODE			0x0004
#define M0_FLUID_NODE			0x0008
#define SOLID_BOUNDARY_NODE		0x0010
#define SOLID_BOUNDARY_NODE0	0x0020	// Solid boundary node in unfouled domain
#define NESW_BOUNDARY_NODE0		0x0040  // NESW_BOUNDARY_NODE in unfouled domain
#define BOUNDARY_NODE0			0x0080  // BOUNDARY_NODE in unfouled domain

#define E_BOUND					0x0100
#define W_BOUND					0x0200
#define N_BOUND					0x0400
#define S_BOUND					0x0800
#define NE_BOUND				0x1000
#define SW_BOUND				0x2000
#define SE_BOUND				0x4000
#define NW_BOUND				0x8000

#define M_FOUL_NODE				0x0009  // M0_FLUID_NODE | M_SOLID_NODE
#define NESW_BOUNDARY_NODE		0x0F00  // E_BOUND | W_BOUND | N_BOUND | S_BOUND
#define FD_BOUNDARY_NODE		0x0F10	//SOLID_BOUNDARY_NODE | NESW_BOUNDARY_NODE
#define BOUNDARY_NODE			0xFF00  // E_BOUND | W_BOUND | ... | NW_BOUND
#define FD_BOUNDARY_NODE0		0x0060	//SOLID_BOUNDARY_NODE0 | NESW_BOUNDARY_NODE0


#else
#define C_BOUND					0x00000000 //just a placeholder. Not actually used
#define M_SOLID_NODE			0x00000001
#define M_FLUID_NODE			0x00000002
#define M0_SOLID_NODE			0x00000004
#define M0_FLUID_NODE			0x00000008
#define SOLID_BOUNDARY_NODE		0x00000010	
#define SOLID_BOUNDARY_NODE0	0x00000020	// Solid boundary node in unfouled domain
#define NESW_BOUNDARY_NODE0		0x00000040  // NESW_BOUNDARY_NODE in unfouled domain
#define BOUNDARY_NODE0			0x00000080  // BOUNDARY_NODE in unfouled domain

#define E_BOUND					0x00000100
#define W_BOUND					0x00000200
#define N_BOUND					0x00000400
#define S_BOUND					0x00000800
#define NE_BOUND				0x00001000
#define SW_BOUND				0x00002000
#define SE_BOUND				0x00004000
#define NW_BOUND				0x00008000
#define E_BOUND_T1				0x00010000
#define W_BOUND_T1				0x00020000
#define N_BOUND_T1				0x00040000
#define S_BOUND_T1				0x00080000
#define NE_BOUND_T1				0x00100000
#define SW_BOUND_T1				0x00200000
#define SE_BOUND_T1				0x00400000
#define NW_BOUND_T1				0x00800000
#define E_BOUND_T2				0x01000000
#define W_BOUND_T2				0x02000000
#define N_BOUND_T2				0x04000000
#define S_BOUND_T2				0x08000000
#define NE_BOUND_T2				0x10000000
#define SW_BOUND_T2				0x20000000
#define SE_BOUND_T2				0x40000000
#define NW_BOUND_T2				0x80000000

#define M_FOUL_NODE				0x00000009  // M0_FLUID_NODE | M_SOLID_NODE
#define NESW_BOUNDARY_NODE		0x00000F00  // E_BOUND | W_BOUND | N_BOUND | S_BOUND
#define FD_BOUNDARY_NODE		0x00000F10	//SOLID_BOUNDARY_NODE | NESW_BOUNDARY_NODE
#define BOUNDARY_NODE			0xFFFFFF00
#define FD_BOUNDARY_NODE0		0x00000060	//SOLID_BOUNDARY_NODE0 | NESW_BOUNDARY_NODE0
#endif

// Resets node originally set as M0_FLUID_NODE to M0_SOLID_NODE
// These are mutually exclusive, so the M0_FLUID_NODE bit must be set
// to zero along with M0_SOLID_NODE bit being set to 1.
#define RESET_NODE_TO_M0_SOLID(nodeVal)		nodeVal &= !M0_FLUID_NODE;\
											nodeVal |= M0_SOLID_NODE;

// Converts node that is currently tagged as M_FLUID_NODE
// to a M_SOLID_NODE. These are mutually exclusive, so the
// M_FLUID_NODE bit must be set to zero along with 
// M_SOLID_NODE bit being set to 1.
#define RESET_NODE_TO_M_SOLID(nodeVal)		nodeVal &= !M_FLUID_NODE;\
											nodeVal |= M_SOLID_NODE;

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


#define WF_EMPTY			0x0000
#define WF_SOLID			0x0001
#define WF_FLUID			0x0002
#define WF_BOT_WALL			0x0004
#define WF_TOP_WALL			0x0008
#define WF_00_SOLID			0x0010
#define WF_10_SOLID			0x0020
#define WF_01_SOLID			0x0040
#define WF_11_SOLID			0x0080

#define WF_WALL				0x000C


#define WF_TEST_ALL_SOLID	0x00F0


// adding ibbRev[dir-1] to ibb index give index of the 
// reverse disribution 
__constant int ibbRev[8] = { DIST_SIZE, -DIST_SIZE, DIST_SIZE, -DIST_SIZE, DIST_SIZE, -DIST_SIZE, DIST_SIZE, -DIST_SIZE };

// Adding ibbNeigh to ibb index gives index of neighbor in reverse
// direction (for use with q >= 0.5). This may need to be set in 
// vlb.setSourceDefines as it requires calculations and im not sure if
// opencl can handle it.
__constant int ibbNeigh[8] = {	DIST_SIZE - 1, 
								-DIST_SIZE + 1,
								DIST_SIZE - CHANNEL_LENGTH_FULL, 
								-DIST_SIZE + CHANNEL_LENGTH_FULL,
								DIST_SIZE - (CHANNEL_LENGTH_FULL + 1),
								-DIST_SIZE + (CHANNEL_LENGTH_FULL + 1),
								DIST_SIZE - CHANNEL_LENGTH_FULL + 1, 
								-DIST_SIZE + CHANNEL_LENGTH_FULL - 1 }


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


#define GET_GLOBAL_IDX(gx,gy)	(gx + gy*CHANNEL_LENGTH_FULL)
#define GET_FULL_GLOBAL_IDX(gx,gy,gdir)	(gx + gy*CHANNEL_LENGTH_FULL + gdir*DIST_SIZE)

#define GET_TR_GLOBAL_IDX(gx,gy)	(gx + gy*FULLSIZEX_TR_PADDED)

void decodeFullGlobalIdx(const int gid, int* xval, int* yval, int* qval)
{
	*qval = gid / DIST_SIZE;
	int xyval = gid % DIST_SIZE;
	*xval = xyval % CHANNEL_LENGTH_FULL;
	*yval = xyval / CHANNEL_LENGTH_FULL;
}

void decodeGlobalIdx(const int gid, int* xval, int* yval)
{
	*xval = gid % CHANNEL_LENGTH_FULL;
	*yval = gid / CHANNEL_LENGTH_FULL;
}

double2 getLocFromGlobalIdx(const int gid)
{
	int2 ret = { { (gid % CHANNEL_LENGTH_FULL), (gid / CHANNEL_LENGTH_FULL) } };
	return convert_double2(ret);
}

int getTRPaddedPosFromLBPos(int gid)
{
	int xval, yval;
	decodeGlobalIdx(gid, &xval, &yval);
	return xval - TR_X_IND_START + FULLSIZEX_TR_PADDED * yval;
}

uint getRevDist(const int gid)
{
	// gid/DIST_SIZE = dist number
	// gid % DIST_SIZE = i + j*CHANNEL_LENGTH_FULL


	return gid % DIST_SIZE + RevDir[(gid / DIST_SIZE) + 1] *
		DIST_SIZE;
}

// BL_TOP_STOP and BL_BOT_STOP is index last element, so test for outside bounds
// is > BL_*_STOP not >= BL_*_STOP

// Convert index for vtr.BL.* arrays to vls.BL index
int convertTRBL2LSBL(int trbl_)
{
	return (trbl_ < NUM_BL_BOT) ? trbl_ + BL_BOT_START : trbl_ + BL_TOP_START;
}

// Convert index for vls.BL to index for vtr.BL.* arrays.
bool convertLSBL2TRBL(int lsbl_, int *trbl_)
{
	if (lsbl_ < BL_BOT_START || lsbl_ > BL_TOP_STOP || (lsbl_ > BL_BOT_STOP && lsbl_ < BL_BOT_STOP))
	{
		*trbl_ = -1;
		return false;
	}
	return (trbl_ < BL_BOT_STOP) ? trbl_ - BL_BOT_START : trbl_ - BL_TOP_START + NUM_BL_BOT;
}


int2 min2(double2 v0, double2 v1)
{
	if (v1.x < v0.x) v0.x = v1.x;
	if (v1.y < v0.y) v0.y = v1.y;
	return convert_int2(v0);
}

int2 max2(double2 v0, double2 v1)
{
	if (v1.x > v0.x) v0.x = v1.x;
	if (v1.y > v0.y) v0.y = v1.y;

	return convert_int2(ceil(v0));
}

int2 min4(double2 v0, double2 v1, double2 v2, double2 v3)
{
	if (v1.x < v0.x) v0.x = v1.x;
	if (v2.x < v0.x) v0.x = v2.x;
	if (v3.x < v0.x) v0.x = v3.x;
	if (v1.y < v0.y) v0.y = v1.y;
	if (v2.y < v0.y) v0.y = v2.y;
	if (v3.y < v0.y) v0.y = v3.y;

	return convert_int2(v0);
}

int2 max4(double2 v0, double2 v1, double2 v2, double2 v3)
{
	if (v1.x > v0.x) v0.x = v1.x;
	if (v2.x > v0.x) v0.x = v2.x;
	if (v3.x > v0.x) v0.x = v3.x;
	if (v1.y > v0.y) v0.y = v1.y;
	if (v2.y > v0.y) v0.y = v2.y;
	if (v3.y > v0.y) v0.y = v3.y;

	return convert_int2(ceil(v0));
}


