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


#define GET_GLOBAL_IDX(gx,gy)	(gx + gy*CHANNEL_LENGTH_FULL)

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


