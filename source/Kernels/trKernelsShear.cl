// Called at beginning of TR solver. Compares shear at BL a particle is deposited on
// with critical shear for the particle size. Removes particle from wall if wall shear
// is sufficient, and places it back in the flow. Currently this is only performed on
// particles that deposited before the previous sort, and removed particles do not
// start to update their position until after the next sort

// TODO: Currently only tests particles that deposited before last sort
//		step. Should implement method to efficiently test all particles
//		even if they deposited after sort step.	

// TODO: Make sure indexing is done correctly, since there has been significant
//			changes to how variables have been passed.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RERELEASE, 1, 1)))
void TR_shear_removal(__global double2* __restrict__ lsc,
	__global ushort* __restrict__ BL_P01ind,
	__global double2* __restrict__ BL_vNvec,
	__global double* __restrict__ BL_blLen,
	__global double* __restrict__ BL_Tau,
	__global short* __restrict__ BL_Tau_Dir,
	__global short* __restrict__ BL_Surf_Type,
	__global double2* __restrict__ P_pos,
	__global uint* __restrict__ P_Num_rep,
	__global short* __restrict__ P_type,
	__global short* __restrict__ P_Dep_Flag,
	__global ushort* __restrict__ P_Dep_timer,
	__global int* __restrict__ P_loc,
	__global uint* BLdep,
	cl_int2 offsetAndMax) // offsetAndMax = {{ offset, maxel }}
{
	int i = get_global_id(0);

	if (i >= offsetAndMax.y)
		return;

	i += offsetAndMax.x;
	
	// skip if already deposited
	if (P_Dep_Flag < 0)
		return;

	int sizeInd = P_type[i];
	int depTimer = P_Dep_timer[i] - 1;
	int blInd = P_Dep_Flag[i];
	
	// if dep timer reaches zero, set particle to permanently deposited
	// and return.
	if (depTimer == 0)
	{
		int bldep_ind = blInd * NUM_PAR_SIZES + sizeInd;
		atomic_add(&BLdep[bldep_ind], P_Num_rep[i]);
		P_Dep_Flag[i] = -2;
		P_loc[i] = -2;
		return;
	}

	int blSurf = sizeInd*2 + BL_Surf_Type[blInd];
	double shear_surf = BL_Tau[blInd];
	
	// If insufficient shear, particle remains with decremented
	// dep_timer being written back to memory.
	if (fabs(shear_surf) <= paramTauCrit[blSurf])
	{
		P_Dep_timer[i] = depTimer;
		return;
	}

	// At sufficient shear, particle is re-entrained in flow
	double Fdrag = paramDcoeff[sizeInd];
	double Flift = paramLcoeff[sizeInd] * sqrt(shear_surf);

	// Currently normal vector.
	double2 vLvec = BL_vNvec[blInd];

	// tangent vector = shear_direction*[vN.y, -vN.x]
	double2 vDvec = Fdrag * double2(vLvec.y, -vLvec.x);
	
	// is actually ushort2, so index is blSurf*2
	ushort lscInd = BL_P01ind[blInd * 2];

	double blLen = BL_blLen[blInd];
	// C1 is point at center of BL slightly offset into fluid domain
	double2 C1 = lsc[lscInd] + double2(vLvec.y, -vLvec.x) * (blLen * 0.5) + vLvec * 0.001;

	// No longer need normal vector, so multiplied by lift magnitude to get lift vector
	vLvec *= Flift;

	// Resultant vector
	double2 vRvec = vLvec + vDvec;

	// displacement vector = F / mass * (dt*dt/2)
	double2 dC_rem = vRvec * shear_surf / paramMp[sizeInd] * DTTR_WALL * DTTR_WALL / 2.;
	
	double dC_mag = length(dC_rem);

	// make sure it is displaced sufficiently far from surface
	dC_mag = max(dC_mag, blLen / 8.);

	// C1 becomes new particle location
	C1 += dC_rem;

	// if particle is removed and goes outside of tracer range, it must be
	// caught here, or behavior is undefined
	P_Dep_Flag[i] = (C1.x < X_MIN_VAL || C1.x > X_MAX_VAL) ? (-2) : (-1);
	int2 Posi = convert_int2(C1);
	P_loc[i] = (C1.x < X_MIN_VAL || C1.x > X_MAX_VAL) ? (-2) : (pcur_loc);
	P_pos[i] = C1;
	P_timer[i] = PAR_TIMER_START;
}

//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////   Functions for calculating wall   /////////////////
/////////////////               shear                /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

// Calculates Fneq using central differences (used in calculation of wall shear)
void calcLBDists(__local double* Fvals, double* v0, double rho, int lid)
{
	double xi_3 = (v0[0]) * (v0[0]);
	double xi_4 = rho * xi_3;
	double xi_6 = (v0[1]) * (v0[1]);
	double xi_8 = xi_4 * xi_6;
	double xi_9 = v0[0] * rho;
	double xi_10 = (1. / 3.) * xi_9;
	double xi_11 = v0[0] * rho * xi_6;
	double xi_12 = 0.5 * xi_11;
	double xi_14 = rho * xi_6;
	double xi_17 = (1.0 / 9.0) * rho - ((1.0 / 6.0) * xi_14 + 0.5 * xi_8) + (1. / 3.) * xi_4;
	double xi_24 = v0[1] * rho;
	double xi_25 = (1. / 3.) * xi_24;
	double xi_34 = (1.0 / 12.0) * xi_9;
	double xi_36 = xi_9 * 0.25 * v0[1];
	double xi_37 = 0.25 * xi_11;
	double xi_38 = 0.0277777777777778 * rho;
	double xi_39 = (1.0 / 12.0) * xi_24;
	double xi_40 = (1.0 / 12.0) * xi_4;
	double xi_41 = (1.0 / 12.0) * xi_14;
	double xi_42 = 0.25 * v0[1] * rho * xi_3;
	double xi_43 = 0.25 * xi_8;
	double xi_44 = xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43;

	// calculation of equilibrium dists
	double8 feq;
	feq.s0 = xi_10 - xi_12 + xi_17;
	feq.s1 = -xi_10 + xi_12 + xi_17;
	feq.s2 = xi_25 - xi_27 + xi_29;
	feq.s3 = -xi_25 + xi_27 + xi_29;
	feq.s4 = xi_34 + xi_36 + xi_37 + xi_44;
	feq.f5 = -xi_34 + xi_36 - xi_37 + xi_38 - xi_39 + xi_40 + xi_41 - xi_42 + xi_43;
	feq.f6 = xi_34 - xi_36 + xi_37 + xi_38 - xi_39 + xi_40 + xi_41 - xi_42 + xi_43;
	feq.f7 = -xi_34 - xi_36 - xi_37 + xi_44;

	// Fneq = F - Feq
	Fvals[lid] -= feq;
}

// Calculates shear stress at LB nodes positioned along wall
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_shear_1(__global double* __restrict__ Ro_array,
	__global double* __restrict__ Ux_array,
	__global double* __restrict__ Uy_array,
	__global double8* __restrict__ CFSS,
	__global double* __restrict__ FA,
	__global double* __restrict__ Tau_array,
	__global int* __restrict__ Loc,
	__global double* __restrict__ nut,
	int TAU_ARRAY_SIZE,
	double dp,
	__local double8* Fvals)
{
	uint ii = get_global_id(0);
	if (ii >= (TAU_ARRAY_SIZE))
		return;

	uint ind = Loc[ii];
	int lid = get_local_id(0);

	Fvals[lid] = double8(FA[ind + DIST_SIZE], 
		FA[ind + 2 * DIST_SIZE],
		FA[ind + 3 * DIST_SIZE],
		FA[ind + 4 * DIST_SIZE],
		FA[ind + 5 * DIST_SIZE],
		FA[ind + 6 * DIST_SIZE],
		FA[ind + 7 * DIST_SIZE]);

	double U[2] = { Ux_array[ind], Uy_array[ind] };
	double Ro = Ro_array[indU];

	calcLBDists(Fvals, U, Ro, lid);
	
	//double visc_val = MU_NUMBER;
	double visc_val = nut[ind] + MU_NUMBER;
	visc_val *= (1. / ((3. * (visc_val)) + 0.5));

	double8 XCoeffs = CFSS[2 * ii] * visc_val;
	double8 YCoeffs = CFSS[2 * ii + 1] * visc_val;
	
	double2 Tau;
	
	Tau.x = dot(Fvals[lid].lo, XCoeffs.lo) + dot(Fvals[lid].hi, XCoeffs.hi);
	Tau.y = dot(Fvals[lid].lo, YCoeffs.lo) + dot(Fvals[lid].hi, YCoeffs.hi);
	
	double tau_val = length(Tau) / Ro;

	Tau_array[ii] = (Tau.x < 0) ? -tau_val : tau_val;
}


// This kernel has changes implemented, but the c++ code has not been updated
// to utilize the changes

#ifdef USE_OPENGL
// Calculates shear at each BL from average of values at surrounding nodes
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_shear_2(__global double4* __restrict__ Weights,
	__global double* __restrict__ tauArray,
	__global int4* __restrict__ sInd,
	__global double* __restrict__ blTau,
	__global short* __restrict__ blType,
	__global int* __restrict__ blColorInd,
	__global float* __restrict__ ls_color_b,
	__global float* __restrict__ ls_color_t,
	__global double* __restrict__ TCRIT_MAX2)
{
	uint ii = get_global_id(0);
	if (ii >= NUM_BL_TOTAL)
		return;

	double4 weight = Weights[ii];
	int4 i = sInd[ii];

	double4 Tautemp = (double4)(nodeTau[i.x], nodeTau[i.y], nodeTau[i.z], nodeTau[i.w]);
	double Tau_wall = dot(weight, Tautemp);
	blTau[ii] = Tau_wall;

	double TCRIT_MAX = TCRIT_MAX2[blType[ii]];
	
	if (Tau_wall < 0.)
	{
		fr = 0.f;
		if (Tau_wall < TCRIT_MAX / 2)
		{
			fb = convert_float((TCRIT_MAX - Tau_wall * 2.f) / TCRIT_MAX);
		}
		else if (Tau_wall < TCRIT_MAX)
		{
			double Tw2 = Tau_wall * 2. - TCRIT_MAX;
			fg -= convert_float((TCRIT_MAX - Tw2) / TCRIT_MAX);
		}
	}
	else
	{
		fb = 0.f;
		if (Tau_wall < TCRIT_MAX / 2)
		{
			fr = convert_float((TCRIT_MAX - Tau_wall * 2.) / TCRIT_MAX);
		}
		else if (Tau_wall < TCRIT_MAX)
		{
			double Tw2 = Tau_wall * 2. - TCRIT_MAX;
			fg -= convert_float((TCRIT_MAX - Tw2) / TCRIT_MAX);
		}
	}
	
	int LSind = blColorInd[ii];
	if (LSind < 0)
	{
		ls_color_t[-LSind * 3] = fr;
		ls_color_t[-LSind * 3 + 1] = fg;
		ls_color_t[-LSind * 3 + 2] = fb;
	}
	else
	{
		ls_color_b[LSind * 3] = fr;
		ls_color_b[LSind * 3 + 1] = fg;
		ls_color_b[LSind * 3 + 2] = fb;
	}
}

#else

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_shear_2(__global double4* __restrict__ Weights,
	__global double* __restrict__ tauArray,
	__global int4* __restrict__ sInd,
	__global double* __restrict__ blTau)
{
	uint ii = get_global_id(0);
	if (ii >= NUM_BL_TOTAL)
		return;

	double4 weight = Weights[ii];
	int4 i = sInd[ii];

	double4 Tautemp = (double4)(nodeTau[i.x], nodeTau[i.y], nodeTau[i.z], nodeTau[i.w]);
	blTau[ii] = dot(weight, Tautemp);
}

#endif

//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////    Functions for updating shear    /////////////////
/////////////////           coefficients             /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

// Updates the arrays storing coefficients/indicies used in the
// calculation of shear stress at nodes and boundary links
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_Find_Shear_Coeffs1(__global int* __restrict__ ssArr,//0
	__global int* __restrict__ nType,
	__global double* CFSS,
	int numBnodes)
{
	uint ii = get_global_id(0);
	if (ii >= numBnodes)
		return;

	int nodeInd = ssArr[ii];
	double2 nodeLoc = getLocFromGlobalIdx(nodeInd);
	
	double2 norm_vec = 0.;
	for (int j = max(ii - 50, 0); j < min(ii + 50, numBnodes); j++)
	{
		int neighInd = ssArr[j];
		double2 neighLoc = getLocFromGlobalIdx(neighInd);
		double2 dR = nodeLoc - neighLoc;
		if (length(dR) > SHEAR_CUTOFF_RADIUS)
			continue;
		
		int neighType = nType[neighInd];
		if (neighType & E_BOUND) { norm_vec.x += 1.; }
		if (neighType & W_BOUND) { norm_vec.x -= 1.; }
		if (neighType & N_BOUND) { norm_vec.y += 1.; }
		if (neighType & S_BOUND) { norm_vec.y -= 1.; }
	}

	norm_vec = normalize(norm_vec);

	int coeffs_ind = ii * 16;

	//double ssMultiplier = 3. * MU_NUMBER / (MU_NUMBER*3. + 0.5);

	for (int aa = 0; aa < 8; aa++)
	{
		double CCx = CXY_DOUBLE[aa].x, CCy = CXY_DOUBLE[aa].y;
		double sumk = (CCx * norm_vec.x + CCy * norm_vec.y);
		
		CFSS[coeffs_ind + aa] = 3. * sumk * (CCx - norm_vec.x * sumk);
		//CFSS[coeffs_ind + aa] = ssMultiplier * sumk * (CCx - norm_vec.x * sumk);

		CFSS[coeffs_ind + 8 + aa] = 3. * sumk * (CCy - norm_vec.y*sumk);
		//CFSS[coeffs_ind + 8 + aa] = ssMultiplier * sumk * (CCy - norm_vec.y*sumk);
	}
}

void reorganize_values(double4 * dRvals, int4 * index_neigh)
{
	if ((*dRvals).x >= (*dRvals).y)
		return;

	(*dRvals).xy = (*dRvals).yx;
	(*index_neigh).xy = (*index_neigh).yx;

	if ((*dRvals).y >= (*dRvals).z)
		return;

	(*dRvals).yz = (*dRvals).zy;
	(*index_neigh).yz = (*index_neigh).zy;

	if ((*dRvals).z >= (*dRvals).w)
		return;

	(*dRvals).zw = (*dRvals).wz;
	(*index_neigh).zw = (*index_neigh).wz;

}

// TODO: Compare various weight kernels for shear, can define yaml parameter
//		to control which weight kernel is used
double weight_kernel_func(double Sval)
{
	//if (Sval <= 1.)
	//{
	//	double tempval = 1. - Sval*Sval;
	//	return (35. / 32. * tempval * tempval * tempval);
	//}
	return 1. / sqrt(2. * PI_NUMBER) * exp(-0.5 * Sval * Sval);
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_Find_Shear_Coeffs2(__global double2* __restrict__ lsc,
	__global ushort* __restrict__ BL_P01ind,
	__global double2* __restrict__ BL_vNvec,
	__global double* __restrict__ BL_blLen,
	__global int* __restrict__ ssArr,
	__global int* __restrict__ ssArrInds,
	__global double4 * Weights,
	__global int4 * Sind)
{
	uint ii = get_global_id(0);
	if (ii >= NUM_BL_TOTAL)
		return;

	ushort lscInds = BL_P01ind[ii*2];
	double blLen = BL_blLen[ii]/2.;
	double2 vNvec = BL_vNvec[ii];
	// C1 is point at center of BL slightly offset into fluid domain
	double2 BLmiddle = lsc[lscInds] + (blLen) * double2(vNvec.y, -vNvec.x);

	int ival = convert_int(BLmiddle.x);

	int start_ind = ssArrInds[max(MIN_SS_ARR_NODE_X, ival - INDEX_RADIUS_SEARCH)];
	int stop_ind = ssArrInds[min(MAX_SS_ARR_NODE_X, ival + INDEX_RADIUS_SEARCH + 1)];

	int4 index_neigh = -1;
	double4 dRvals = 100.;

	for (int j = start_ind; j < stop_ind; j++)
	{
		int neighInd = ssArr[j];
		double2 neighLoc = getLocFromGlobalIdx(neighInd);
		double2 dR = BLmiddle - neighLoc;
		double dist = length(dR);
		if (dist < dRvals.x && dist < SHEAR_CUTOFF_RADIUS)
		{
			dRvals.x = dist;
			index_neigh.x = j;
			reorganize_values(&dRvals, &index_neigh);
		}
	}

	double4 weight_temp;
	weight_temp.x = weight_kernel_func(dRvals.x / SHEAR_CUTOFF_RADIUS);
	weight_temp.y = weight_kernel_func(dRvals.y / SHEAR_CUTOFF_RADIUS);
	weight_temp.z = weight_kernel_func(dRvals.z / SHEAR_CUTOFF_RADIUS);
	weight_temp.w = weight_kernel_func(dRvals.w / SHEAR_CUTOFF_RADIUS);

	if (index_neigh.x > index_neigh.y)
	{
		index_neigh.xy = index_neigh.yx;
		weight_temp.xy = weight_temp.yx;
	}

	if (index_neigh.z > index_neigh.w)
	{
		index_neigh.zw = index_neigh.wz;
		weight_temp.zw = weight_temp.wz;
	}

	if (index_neigh.x > index_neigh.z)
	{
		index_neigh.xz = index_neigh.zx;
		weight_temp.xz = weight_temp.zx;
	}

	if (index_neigh.y > index_neigh.w)
	{
		index_neigh.yw = index_neigh.wy;
		weight_temp.yw = weight_temp.wy;
	}

	if (index_neigh.y > index_neigh.z)
	{
		index_neigh.yz = index_neigh.zy;
		weight_temp.yz = weight_temp.zy;
	}

	Weights[ii] = weight_temp / dot(weight_temp, (double4)(1.));
	Sind[ii] = index_neigh;
}

#if defined(SAVE_FULL_SHEAR)

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_shear_save_out(__global double* __restrict__ blTau,
	__global double* __restrict__ SSout,
	int save_loc)
{
	uint ii = get_global_id(0);
	if (ii >= NUM_BL_TOTAL)
		return;
	SSout[ii + save_loc * NUM_BL_TOTAL] = blTau[ii];
}

#elif defined(SAVE_SHEAR_TOP)
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_shear_save_out(__global double* __restrict__ blTau,
	__global double* __restrict__ SSout,
	int save_loc)
{
	uint ii = get_global_id(0);
	if (ii >= NUM_BL_TOP)
		return;
	SSout[ii + save_loc * NUM_BL_TOP] = blTau[ii + NUM_BL_BOT];
}
#else //SAVE_SHEAR_BOT
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_shear_save_out(__global double* __restrict__ blTau,
	__global double* __restrict__ SSout,
	int save_loc)
{
	uint ii = get_global_id(0);
	if (ii >= NUM_BL_BOT)
		return;
	SSout[ii + save_loc * NUM_BL_BOT] = blTau[ii];
}
#endif