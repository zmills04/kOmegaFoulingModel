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



int bcFindIntersectionNusselt(double *dist, double2 vLd, double2 vPL, double2 vN)
{
	double den = dot(vN, vLd);
	if (den == 0.)
		return false;

	*dist = dot(vN, vPL) / den;
	
	if ((*dist) >= 0. && (*dist) <= 1.)
	{
		return true;
	}
	return false;
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


//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////    Functions for updating shear    /////////////////
/////////////////           coefficients             /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////
//
//// Re-organizes indicies of neighbor list and distances used to
//// calculate Nu value at wall from Nu value at neighboring nodes
//void reorganize_values_for_nu(double4* dRvals, int4* index_neigh)
//{
//	if ((*dRvals).x >= (*dRvals).y)
//		return;
//
//	(*dRvals).xy = (*dRvals).yx;
//	(*index_neigh).xy = (*index_neigh).yx;
//
//	if ((*dRvals).y >= (*dRvals).z)
//		return;
//
//	(*dRvals).yz = (*dRvals).zy;
//	(*index_neigh).yz = (*index_neigh).zy;
//
//	if ((*dRvals).z >= (*dRvals).w)
//		return;
//
//	(*dRvals).zw = (*dRvals).wz;
//	(*index_neigh).zw = (*index_neigh).wz;
//}
//
//// TODO: Compare various weight kernels for Nu, can define yaml parameter
////		to control which weight kernel is used
//double nu_weight_kernel_func(double Sval)
//{
//	return 1. / sqrt(2. * PI_NUMBER) * exp(-0.5 * Sval * Sval);
//}
//
//// Updates the arrays storing coefficients/indicies used in the
//// calculation of Nu at wall nodes
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
//void FD_Nu_Coeffs1(__global int2* __restrict__ nuArr,//same as ssArr
//	__global int* __restrict__ nType,
//	__global double* __restrict__ nuNormVecs,
//	int numBnodes)
//{
//	uint ii = get_global_id(0);
//	if (ii >= numBnodes)
//		return;
//
//	int nodeInd = ssArr[ii];
//	double2 nodeLoc = getLocFromGlobalIdx(nodeInd);
//
//	double2 norm_vec = 0.;
//	for (int j = max(ii - 50, 0); j < min(ii + 50, numBnodes); j++)
//	{
//		int neighInd = ssArr[j];
//		double2 neighLoc = getLocFromGlobalIdx(neighInd);
//		double2 dR = nodeLoc - neighLoc;
//		if (length(dR) > NU_CUTOFF_RADIUS)
//			continue;
//
//		int neighType = nType[neighInd];
//		if (neighType & E_BOUND) { norm_vec.x += 1.; }
//		if (neighType & W_BOUND) { norm_vec.x -= 1.; }
//		if (neighType & N_BOUND) { norm_vec.y += 1.; }
//		if (neighType & S_BOUND) { norm_vec.y -= 1.; }
//	}
//
//	nuNormVecs[ii] = normalize(norm_vec);
//}
//
//
//
//
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
//void FD_Nu_Coeffs2(__global double2 * __restrict__ lsc,
//	__global ushort * __restrict__ BL_P01ind,
//	__global double2 * __restrict__ BL_vNvec,
//	__global double* __restrict__ BL_blLen,
//	__global int* __restrict__ ssArr,
//	__global int* __restrict__ ssArrInds,
//	__global double4 * Weights,
//	__global int4 * Sind)
//{
//	uint ii = get_global_id(0);
//	if (ii >= NUM_BL_TOTAL)
//		return;
//
//	ushort lscInds = BL_P01ind[ii * 2];
//	double blLen = BL_blLen[ii] / 2.;
//	double2 vNvec = BL_vNvec[ii];
//	// C1 is point at center of BL slightly offset into fluid domain
//	double2 BLmiddle = lsc[lscInds] + (blLen)* double2(vNvec.y, -vNvec.x);
//
//	int ival = convert_int(BLmiddle.x);
//
//	int start_ind = ssArrInds[max(MIN_SS_ARR_NODE_X, ival - INDEX_RADIUS_SEARCH)];
//	int stop_ind = ssArrInds[min(MAX_SS_ARR_NODE_X, ival + INDEX_RADIUS_SEARCH + 1)];
//
//	int4 index_neigh = -1;
//	double4 dRvals = 100.;
//
//	for (int j = start_ind; j < stop_ind; j++)
//	{
//		int neighInd = ssArr[j];
//		double2 neighLoc = getLocFromGlobalIdx(neighInd);
//		double2 dR = BLmiddle - neighLoc;
//		double dist = length(dR);
//		if (dist < dRvals.x && dist < SHEAR_CUTOFF_RADIUS)
//		{
//			dRvals.x = dist;
//			index_neigh.x = j;
//			reorganize_values(&dRvals, &index_neigh);
//		}
//	}
//
//	double4 weight_temp;
//	weight_temp.x = weight_kernel_func(dRvals.x / SHEAR_CUTOFF_RADIUS);
//	weight_temp.y = weight_kernel_func(dRvals.y / SHEAR_CUTOFF_RADIUS);
//	weight_temp.z = weight_kernel_func(dRvals.z / SHEAR_CUTOFF_RADIUS);
//	weight_temp.w = weight_kernel_func(dRvals.w / SHEAR_CUTOFF_RADIUS);
//
//	if (index_neigh.x > index_neigh.y)
//	{
//		index_neigh.xy = index_neigh.yx;
//		weight_temp.xy = weight_temp.yx;
//	}
//
//	if (index_neigh.z > index_neigh.w)
//	{
//		index_neigh.zw = index_neigh.wz;
//		weight_temp.zw = weight_temp.wz;
//	}
//
//	if (index_neigh.x > index_neigh.z)
//	{
//		index_neigh.xz = index_neigh.zx;
//		weight_temp.xz = weight_temp.zx;
//	}
//
//	if (index_neigh.y > index_neigh.w)
//	{
//		index_neigh.yw = index_neigh.wy;
//		weight_temp.yw = weight_temp.wy;
//	}
//
//	if (index_neigh.y > index_neigh.z)
//	{
//		index_neigh.yz = index_neigh.zy;
//		weight_temp.yz = weight_temp.zy;
//	}
//
//	Weights[ii] = weight_temp / dot(weight_temp, (double4)(1.));
//	Sind[ii] = index_neigh;
//}
//
//#if defined(SAVE_FULL_SHEAR)
//
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
//void TR_shear_save_out(__global double* __restrict__ blTau,
//	__global double* __restrict__ SSout,
//	int save_loc)
//{
//	uint ii = get_global_id(0);
//	if (ii >= NUM_BL_TOTAL)
//		return;
//	SSout[ii + save_loc * NUM_BL_TOTAL] = blTau[ii];
//}
//
//#elif defined(SAVE_SHEAR_TOP)
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
//void TR_shear_save_out(__global double* __restrict__ blTau,
//	__global double* __restrict__ SSout,
//	int save_loc)
//{
//	uint ii = get_global_id(0);
//	if (ii >= NUM_BL_TOP)
//		return;
//	SSout[ii + save_loc * NUM_BL_TOP] = blTau[ii + NUM_BL_BOT];
//}
//#else //SAVE_SHEAR_BOT
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
//void TR_shear_save_out(__global double* __restrict__ blTau,
//	__global double* __restrict__ SSout,
//	int save_loc)
//{
//	uint ii = get_global_id(0);
//	if (ii >= NUM_BL_BOT)
//		return;
//	SSout[ii + save_loc * NUM_BL_BOT] = blTau[ii];
//}
//#endif






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

