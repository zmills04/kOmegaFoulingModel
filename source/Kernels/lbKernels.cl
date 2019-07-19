
////////////////////////////////////////////////////////////
////                                                    ////
//// Old Propagation Code, which has an error somewhere ////
////                                                    ////
////////////////////////////////////////////////////////////

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
	__global const NTYPE_TYPE *__restrict__ Map)
{
	const int gx = get_global_id(0);			
	
	if (gx > CHANNEL_LENGTH) { return; }
	const int gy = get_global_id(1) + 1;
	int gi = GET_GLOBAL_IDX(gx, gy);
	const NTYPE_TYPE type = Map[gi];
	if (type & M_SOLID_NODE) { return; }

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

inline void postCollisionBC(const FTYPE *Fi, const NTYPE_TYPE orientation, const int gi, __global FTYPE *__restrict__ dist_out)
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
	__global NTYPE_TYPE *__restrict__ Map,
#ifdef IN_KERNEL_IBB
	__global double *__restrict__ dxvals,
#endif
	int options)
{

	const int gx = get_global_id(0);

	if (gx > CHANNEL_LENGTH) { return; }
	const int gy = get_global_id(1);
	const int gi = GET_GLOBAL_IDX(gx, gy);
	const NTYPE_TYPE type = Map[gi];
	
	if (type & M_SOLID_NODE) { return; }

	double Fi[9];

#define	GET_DIST_MACRO(dir)				Fi[(dir ## _DIR)] = FA[gi + DIST_SIZE * (dir ## _DIR)];

	CALL_ALL_DISTS(GET_DIST_MACRO);
#undef GET_DIST_MACRO

#ifdef IN_KERNEL_IBB
	// IBB for q < 0.5
#define TEST_WRITE(dir)		if(type & (dir ## _BOUND_T1)) { \
								double q2 = 2.*dxvals[gi + DIST_SIZE * ((dir ## _DIR) - 1)]; \
								Fi[(dir ## _REV)] = q2 * Fi[(dir ## _REV)] + (1. - q2) * Fi[(dir ## _DIR)]; \
							}

	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE
#endif
	
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

#ifdef IN_KERNEL_IBB
#define TEST_WRITE(dir)		if(type & (dir ## _BOUND_T2)) { \
								double q2 = 0.5 / dxvals[gi + DIST_SIZE * ((dir ## _DIR) - 1)]; \
								Fi[(dir ## _DIR)] = q2 * Fi[(dir ## _DIR)] + (1. - q2) * Fi[(dir ## _REV)]; \
							}

	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE
#endif

	
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

}


#define DBGARR(dbgnum, dbgval)	DebugArrs[gi + dbgnum*DIST_SIZE] = dbgval

__kernel void LB_collision_SRT_Fluid_w_kOmega(__global double* __restrict__ rho_array,
	__global double* __restrict__ ovx,
	__global double* __restrict__ ovy,
	__global double* __restrict__ FA,
	__global double* __restrict__ FB,
	__global NTYPE_TYPE* __restrict__ Map,
#ifdef IN_KERNEL_IBB
	__global double* __restrict__ dxvals,
#endif
	__global double* __restrict__ wallD,
	__global double* __restrict__ iNut,
	__global double* __restrict__ iK,
	__global double* __restrict__ iO,
	__global double* __restrict__ iSxy,// strain rate tensor components
#ifdef USE_TURB_VEL_BC
	__global int* __restrict__ ssIndMap,
	__global double* __restrict__ wallShearArr,
#endif
	int options)
{

	const int gx = get_global_id(0);

	if (gx > CHANNEL_LENGTH) { return; }
	const int gy = get_global_id(1);
	const int gi = GET_GLOBAL_IDX(gx, gy);
	const NTYPE_TYPE type = Map[gi];
	if ((type & FLUID_NODE) == 0) { return; }



	double Fi[9];

#define	GET_DIST_MACRO(dir)				Fi[(dir ## _DIR)] = FA[gi + DIST_SIZE * (dir ## _DIR)];

	CALL_ALL_DISTS(GET_DIST_MACRO);
#undef GET_DIST_MACRO


#ifdef IN_KERNEL_IBB
	// IBB for q < 0.5
#define TEST_WRITE(dir)		if(type & (dir ## _BOUND_T1)) { \
								double q2 = 2.*dxvals[gi + DIST_SIZE * ((dir ## _DIR) - 1)]; \
								Fi[(dir ## _REV)] = q2 * Fi[(dir ## _REV)] + (1. - q2) * Fi[(dir ## _DIR)]; \
								}

	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE
#endif

	double v0[2], rho;
	compute_macro_quant(Fi, &rho, v0);

	// Read in and calculate necessary terms for calculation of 
	// turbulent viscosity
	double omegaval = iO[gi], kval = iK[gi], yval = wallD[gi];
	double F2arg = max(2.*sqrt(kval) / (KO_BETA_STAR*omegaval*yval), 500.*MU_NUMBER / (yval*yval*omegaval));
	double F2val = tanh(F2arg*F2arg); 


	// Multiply F2val by the magnitude of the strain rate S = sqrt(2*S_ij*S_ij)
	double sxx_val = iSxy[gi], sxy_val = iSxy[gi + DIST_SIZE], syy_val = iSxy[gi + DIST_SIZE*2];
	F2val *= sqrt(2. * (sxx_val * sxx_val + 2. * sxy_val * sxy_val + syy_val * syy_val));

	// Calculate Nu_turbulent
	double nut = KO_A1*kval / max(KO_A1*omegaval, F2val);
	
	// Use SST boundary condition for any node with N,E,S,W neighbors
	// on other side of boundary (i.e. these nodes have tau_wall already
	// calculated, which can allow for easy calculation of friction Velocity
	// Boundary condition for omega will be handled in the BiCGStab solver setup


#ifdef USE_TURB_VEL_BC
	if (type & NESW_BOUNDARY_NODE)
	{
		// get yplus value = sqrt(tau_wall/rho)*wall_distance/nu
		double yplus = yval*sqrt(fabs(wallShearArr[ssIndMap[gi]]) / rho) / MU_NUMBER;
		double velMag = sqrt(pown(v0[0], 2) + pown(v0[1], 2));
		double fricVelVisc = velMag / yplus;
		double fricVelLog = velMag / (log(yplus) / KO_KAPPA + 5.0);
		double fricVel = pow(pown(fricVelLog, 4), pown(fricVelVisc, 4));
		v0[0] *= fricVel / velMag;
		v0[1] *= fricVel / velMag;
	}
#endif


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

	//								   ____
	//								   \
	//                                  \	
	// Shear stress t_ab = -3*nu*tau *	/    C_ia*C_ib * (fneq_i)
	//								   /___i


	//								      ____
	//								      \
	//									   \	
	// Strain tensor S_ab = -1.5*tau/rho * /    C_ia*C_ib * (fneq_i) = t_ab/(2*rho*nu)
	//								      /___i

	// Storing strain rate Tensor components for use in kOmega and shear stress
	// kernels
	iSxy[gi] = -1.5 * (tau0) / rho * sxx_new;
	iSxy[gi + DIST_SIZE] = -1.5 * (tau0) / rho * sxy_new;
	iSxy[gi + DIST_SIZE*2] = -1.5 * (tau0) / rho * syy_new;


#ifdef IN_KERNEL_IBB
#define TEST_WRITE(dir)		if(type & (dir ## _BOUND_T2)) { \
								double q2 = 0.5 / dxvals[gi + DIST_SIZE * ((dir ## _DIR) - 1)]; \
								Fi[(dir ## _DIR)] = q2 * Fi[(dir ## _DIR)] + (1. - q2) * Fi[(dir ## _REV)]; \
								}

	CALL_ALL_DISTS_NO_C(TEST_WRITE);
#undef TEST_WRITE
#endif


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

// Interpolated bounce back kernel
// TODO: Test that implementation is correct
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1)))
void lbIBB(__global int *__restrict__ ibbArr,
	__global double *__restrict__ dXCurr,
	__global double *__restrict__ FB,
	int max_el)
{
	int i = get_global_id(0);

	if (i >= max_el)
		return;

	int loc = ibbArr[i];
	
	int dir = (loc / DIST_SIZE);
	int revDir = loc + ibbRev[dir-1];
	double twoQ = 2.*dXCurr[loc + (dir-1)*DIST_SIZE];

	if (dxVal <= 0.5)
	{
		FB[revDir] = twoQ * FB[revDir] + (1. - twoQ) * FB[loc];
	}
	else
	{
		revNeigh = loc + ibbNeigh[dir - 1];
		FB[revDir] = 1./twoQ * (FB[revDir] + (twoQ - 1.) * FB[revNeigh])
	}
}

