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

