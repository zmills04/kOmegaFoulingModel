
uint find_SS_min(__global bLinks *BLlist, uint bl)
{
	int SSdir = BLlist[bl].dir;
	uint bl_next = bl + SSdir;
	double Taucur = BLlist[bl].Tau;
	double Taunext = BLlist[bl_next].Tau;

#ifdef MAX_NUM_BL_ROLL
	int count = 0;
	int max_roll = MAX_NUM_BL_ROLL;
	if (BLlist[bl].vP0.x <= X_START_VAL + 3.)
		max_roll = 3;
#endif

	while (1)
	{
#ifdef MAX_NUM_BL_ROLL
		if (Taucur > Taunext && (SSdir == BLlist[bl_next].dir) && count < max_roll)
		{
			bl = bl_next;
			Taucur = Taunext;
			bl_next += SSdir;
			Taunext = BLlist[bl_next].Tau;
			count++;
		}
		else
		{
			break;
		}
#else
		if (Taucur >= Taunext && (SSdir == BLlist[bl_next].dir))
		{
			bl = bl_next;
			Taucur = Taunext;
			bl_next += SSdir;
			Taunext = BLlist[bl_next].Tau;
		}
		else
		{
			break;
		}
#endif
	}
	return (bl);
}


//Called at beginning of TR solver.
// TODO: Currently only tests particles that deposited before last sort
//		step. Should implement method to efficiently test all particles
//		even if they deposited after sort step.	

// TODO: Make sure indexing is done correctly, since there has been significant
//			changes to how variables have been passed.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RERELEASE, 1, 1)))
void TR_shear_removal(__global double2* __restrict__ lsc,
	__global uint* __restrict__ BL_P01ind,
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
	__global double *__restrict__ PParam_tau_crit,
	__global double* __restrict__ PParam_D_coeff,
	__global double* __restrict__ PParam_L_coeff,
	__global double* __restrict__ PParam_Mp,
	__global uint* BLdep,
	uint offset,
	uint maxel)
{
	int i = get_global_id(0);

	if (i >= maxel)
		return;

	i += offset;

	
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
	if (fabs(shear_surf) <= PParam_tau_crit[blSurf])
	{
		P_Dep_timer[i] = depTimer;
		return;
	}

	// At sufficient shear, particle is re-entrained in flow
	double Fdrag = PParam_D_coeff[blSurf];
	double Flift = PParam_L_coeff[blSurf] * sqrt(shear_surf);

	// Currently normal vector.
	double2 vLvec = BL_vNvec[blInd];

	// tangent vector = shear_direction*[vN.y, -vN.x]
	double2 vDvec = Fdrag * double2(vLvec.y, -vLvec.x);
	
	// is actually uint2, so index is blSurf*2
	int lscInd = BL_P01ind[blInd * 2];

	double blLen = BL_blLen[blInd];
	// C1 is point at center of BL slightly offset into fluid domain
	double2 C1 = lsc[lscInd] + double2(vLvec.y, -vLvec.x) * (blLen * 0.5) + vLvec * 0.001;

	// No longer need normal vector, so multiplied by lift magnitude to get lift vector
	vLvec *= Flift;

	// Resultant vector
	double2 vRvec = vLvec + vDvec;

	// displacement vector = F / mass * (dt*dt/2)
	double2 dC_rem = vRvec * shear_surf / PParam_Mp[blSurf] * DTTR_WALL * DTTR_WALL / 2.;
	
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

// Calculates Fneq using central differences
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


//// This has not been updated for new implementations

//#ifdef USE_OPENGL
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
//void TR_shear_2(__global int* BLind,
//	__global double4* Weights,
//	__global double* Tau_array,
//	__global int4* Sind,
//	__global bLinks* BL_array,
//	int NUMBLINKS,
//	__global float* ls_color_b,
//	__global float* ls_color_t,
//	__global double* TCRIT_MAX2)
//{
//	uint ii = get_global_id(0);
//	if (ii >= NUMBLINKS)
//		return;
//	int jj = BLind[ii];
//	bLinks BLtemp = BL_array[jj];
//	int4 i = Sind[ii];
//	double4 weight = Weights[ii];
//	float fr = 1.f, fg = 1.f, fb = 1.f;
//	double4 Tautemp = (double4)(Tau_array[i.x], Tau_array[i.y], Tau_array[i.z], Tau_array[i.w]);
//	double Tau_wall = dot(weight, Tautemp);
//
//	double TCRIT_MAX = TCRIT_MAX2[BLtemp.int_type];
//
//	if (Tau_wall < 0.)
//	{
//		Tau_wall *= -1.;
//		BLtemp.dir = -1;
//		BLtemp.Tau = Tau_wall;
//		fr = 0.f;
//		if (Tau_wall < TCRIT_MAX / 2)
//		{
//			fb = convert_float((TCRIT_MAX - Tau_wall * 2.f) / TCRIT_MAX);
//		}
//		else if (Tau_wall < TCRIT_MAX)
//		{
//			double Tw2 = Tau_wall * 2. - TCRIT_MAX;
//			fg -= convert_float((TCRIT_MAX - Tw2) / TCRIT_MAX);
//		}
//	}
//	else
//	{
//		BLtemp.dir = 1;
//		BLtemp.Tau = Tau_wall;
//		fb = 0.f;
//		if (Tau_wall < TCRIT_MAX / 2)
//		{
//			fr = convert_float((TCRIT_MAX - Tau_wall * 2.) / TCRIT_MAX);
//		}
//		else if (Tau_wall < TCRIT_MAX)
//		{
//			double Tw2 = Tau_wall * 2. - TCRIT_MAX;
//			fg -= convert_float((TCRIT_MAX - Tw2) / TCRIT_MAX);
//		}
//	}
//
//	BL_array[jj] = BLtemp;
//
//	int LSind = BLtemp.Color_ind;
//	if (LSind < 0)
//	{
//		ls_color_t[-LSind * 3] = fr;
//		ls_color_t[-LSind * 3 + 1] = fg;
//		ls_color_t[-LSind * 3 + 2] = fb;
//	}
//	else
//	{
//		ls_color_b[LSind * 3] = fr;
//		ls_color_b[LSind * 3 + 1] = fg;
//		ls_color_b[LSind * 3 + 2] = fb;
//	}
//}
//#else

// Calculates shear at each BL from average of values at surrounding nodes
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1))) 
void TR_shear_2(__global int* __restrict__ BLind,
	__global double4* __restrict__ Weights,
	__global double* __restrict__ Tau_array,
	__global int4* __restrict__ Sind,
	__global double *__restrict__ BL_Tau,
	int NUMBLINKS)
{
	uint ii = get_global_id(0);
	if (ii >= NUMBLINKS)
		return;
	
	int jj = BLind[ii];
	double4 weight = Weights[ii];
	int4 i = Sind[ii];

	double4 Tautemp = (double4)(Tau_array[i.x], Tau_array[i.y], Tau_array[i.z], Tau_array[i.w]);
	BL_Tau[jj] = dot(weight, Tautemp);
}

//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////    Functions for updating shear    /////////////////
/////////////////           coefficients             /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
void TR_Find_Shear_Coeffs1(__global int3 * Bindicies_array,//0
	__global int* Stor,//1
	__global int* dir_array,//2
	__global int2 * ii0_array,//3
	__global double* CFSS,//4
	__global int2 * Loc,//5
	__global int* Bind_ind,//6
	int num_el,//7
	double Cutoff_Rad,//8
	int numBnodes,//9
	int Bindicies_top_start)//10
{
	uint ii = get_global_id(0);
	if (ii >= numBnodes)
		return;

	int3 Bindicies = Bindicies_array[ii];
	double2 ii0 = convert_double2(Bindicies.xy);
	int kk = (ii >= Bindicies_top_start) ? (Bindicies.x * 2 + 1) : (Bindicies.x * 2);

	atomic_xchg(&Bind_ind[kk], ii);

	Loc[ii] = (int2)(Bindicies.x * FULLSIZEY_UT, Bindicies.x * FULLSIZEY);
	Loc[ii] += Stor[Bindicies.x * DOMAIN_SIZE_Y + Bindicies.y];
	int ind = Bindicies.z;

	double2 norm_vec = 0.;
	int jmax = MIN(ind + 50, num_el);
	for (int j = MAX(ind - 50, 0); j < jmax; j++)
	{
		int dir = dir_array[j];
		if (dir > 4)
			continue;
		double2 ii1 = convert_double2(ii0_array[j]);
		double2 dR = ii0 - ii1;

		if (length(dR) <= Cutoff_Rad)
		{
			norm_vec += CXY_DOUBLE[RevDir[dir - 1]];
		}
	}

	norm_vec = normalize(norm_vec);

	int coeffs_ind = ii * 16;

	for (int aa = 0; aa < 8; aa++)
	{
		double CCx = CXY_DOUBLE[aa].x, CCy = CXY_DOUBLE[aa].y;
		double sumk = CCx * norm_vec.x * norm_vec.x;
		sumk += CCy * norm_vec.x * norm_vec.y;

		double sumj = CCx * norm_vec.x * (CCx - sumk);
		sumj += CCy * norm_vec.y * (CCx - sumk);

		//double taunumber = MU_NUMBER*3. + 0.5;

		//CFSS[coeffs_ind + aa] = sumj*MU_NUMBER*3. / taunumber;
		CFSS[coeffs_ind + aa] = sumj * 3.;
		double sumky = CCx * norm_vec.y * norm_vec.x;
		sumky += CCy * norm_vec.y * norm_vec.y;

		double sumjy = CCx * norm_vec.x * (CCy - sumky);
		sumjy += CCy * norm_vec.y * (CCy - sumky);

		CFSS[coeffs_ind + 8 + aa] = sumjy * 3.;
		//CFSS[coeffs_ind + 8 + aa] = sumjy*MU_NUMBER*3. / taunumber;
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
void TR_Find_Shear_Coeffs2(__global int* BLind,//0
	__global double4 * Weights,//1
	__global int4 * Sind,//2
	__global bLinks * BL,//3
	__global int3 * Bindicies,//4
	__global int2 * Bind_ind,//5
	int NUMBLINKS,//6
	double cutoff_radius,//7
	int BLinks_top_ind,//8
	int Bindicies_radius,//9
	int2 Bindicies_end)//10
{
	uint curind = get_global_id(0);
	if (curind >= NUMBLINKS)
		return;

	int i = BLind[curind];
	bLinks BLcur = BL[i];

	double2 BLmiddle = BLcur.vP0 + BLcur.vTvec * BLcur.blLen / 2.;
	int izz = convert_int(BLmiddle.x);
	int ind = Bind_ind[izz].x;
	while (ind == -1)
	{
		izz++;
		ind = Bind_ind[izz].x;
	}
	int start_ind = MAX(0, ind - Bindicies_radius);
	int stop_ind = MIN(Bindicies_end.x, ind + Bindicies_radius);

	//Bindicies_end.x is starting index of Bindicies located at top
	//Bindicies_end.y is total length of Bindicies

	if (i >= BLinks_top_ind)
	{
		ind = Bind_ind[izz].y;
		while (ind == -1)
		{
			izz--;
			ind = Bind_ind[izz].y;
		}
		start_ind = MAX(Bindicies_end.x, ind - Bindicies_radius);
		stop_ind = MIN(Bindicies_end.y, ind + Bindicies_radius);
	}

	int4 index_neigh = -1;
	double4 dRvals = 100.;

	int kk = start_ind;
	while (kk < stop_ind)
	{
		double2 ii0 = convert_double2(Bindicies[kk].xy);
		double2	dR = BLmiddle - ii0;
		double dist = length(dR);
		if (dist < dRvals.x)
		{
			dRvals.x = dist;
			index_neigh.x = kk;
			reorganize_values(&dRvals, &index_neigh);
		}

		kk++;
	}

	double4 weight_temp;
	weight_temp.x = weight_kernel_func(dRvals.x / cutoff_radius);
	weight_temp.y = weight_kernel_func(dRvals.y / cutoff_radius);
	weight_temp.z = weight_kernel_func(dRvals.z / cutoff_radius);
	weight_temp.w = weight_kernel_func(dRvals.w / cutoff_radius);

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

	Weights[curind] = weight_temp / dot(weight_temp, (double4)(1.));
	Sind[curind] = index_neigh;
}
