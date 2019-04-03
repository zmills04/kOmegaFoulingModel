//! Return a 32-bit integer in the range [0..2^32)
uint2 MWC64X_NextUint(uint2 s, double *resd)
{
	uint res = s.x ^ s.y;
	uint X = s.x, C = s.y;

	uint Xn = MWC64X_A*X + C;
	uint carry = (uint)(Xn<C);				// The (Xn<C) will be zero or one for scalar
	uint Cn = mad_hi(MWC64X_A, X, carry);

	s.x = Xn;
	s.y = Cn;
	*resd = convert_double(res)/RAND_MAX;
	return s;
}

uint2 MWC64X_NextUint2(uint2 s, double2 *resd)
{
	uint res = s.x ^ s.y;

	(*resd).x = convert_double(res) / RAND_MAX;
	uint X = s.x, C = s.y;

	uint Xn = MWC64X_A*X + C;
	uint carry = (uint)(Xn<C);				// The (Xn<C) will be zero or one for scalar
	uint Cn = mad_hi(MWC64X_A, X, carry);

	res = Xn ^ Cn;
	(*resd).y = convert_double(res) / RAND_MAX;
	Xn = MWC64X_A*X + C;
	carry = (uint)(Xn<C);				// The (Xn<C) will be zero or one for scalar
	Cn = mad_hi(MWC64X_A, X, carry);

	s.x = Xn;
	s.y = Cn;
	return s;
}

//#define USE_MAD

double2 get_Vadv(double2 U00, double2 U10, double2 U01, double2 U11, double2 dX1, double2 dX0)
{
#ifdef USE_MAD
	double2 U0y = mad(U00, dX0.y, U01 * dX1.y);
	double2 U1y = mad(U10, dX0.y, U11 * dX1.y);

	return mad(U0y, dX0.x, U1y * dX1.x);
#else
	double2 U0y = fma(U00, dX0.y, U01 * dX1.y);
	double2 U1y = fma(U10, dX0.y, U11 * dX1.y);

	return fma(U0y, dX0.x, U1y * dX1.x);
#endif
}

double2 get_Vth(double4 TT, double Kth, double2 dX1, double2 dX0)
{
#ifdef USE_MAD
	double T0y = mad(TT.x, dX0.y, TT.z * dX1.y);
	double T1y = mad(TT.y, dX0.y, TT.w * dX1.y);
	double Txy = mad(T0y, dX0.x, T1y * dX1.x);

	double T0x = mad(TT.x, dX0.x, TT.y * dX1.x);
	double T1x = mad(TT.z, dX0.x, TT.w * dX1.x);
#else
	double T0y = fma(TT.x, dX0.y, TT.z * dX1.y);
	double T1y = fma(TT.y, dX0.y, TT.w * dX1.y);
	double Txy = fma(T0y, dX0.x, T1y * dX1.x);

	double T0x = fma(TT.x, dX0.x, TT.y * dX1.x);
	double T1x = fma(TT.z, dX0.x, TT.w * dX1.x);
#endif
	double2 dTdx = (double2)(T0y - T1y, T0x - T1x);  ///-dTdx
	return dTdx * MU_NUMBER * Kth / Txy;
}

int bcFindIntersectionLinePlane(double2 *vC, double *dist, double2 vL0, double2 vLd, bLinks BLcur)
{
	double2 vN = BLcur.vNvec;
	double den = dot(vN, vLd);
	if (den == 0.)
		return FALSE;

	*dist = dot(vN, (BLcur.vP0 - vL0)) / den;
	*vC = vL0 + vLd * (*dist);

	if ((*dist) >= 0. && (*dist) <= 1.)
	{
		double vd0 = distance(*vC, BLcur.vP0), vd1 = distance(*vC, BLcur.vP1);
		if (fabs(BLcur.blLen - vd0 - vd1) < CEPS)
			return TRUE;
	}
	return FALSE;
}

int find_intersection(double2 *vCcut_use, double2 C1, double2 dC, short Bind[], __global bLinks *BLlist, int Sflag, int *bl_use)
{
	double dist_use = 100.;
	double dist;
	double2 vCcut;
	int Tflag = -1;
	for (int k = 0; k < MAX_BL_PER_NODE; k++)
	{
		int bl = Bind[k];
		if (bl == -1)
			continue;
		if (bcFindIntersectionLinePlane(&vCcut, &dist, C1, dC, BLlist[bl]) == TRUE)
		{
			if (dist < dist_use)
			{
				dist_use = dist;
				Tflag = (1)*Sflag;
				*bl_use = bl;
				*vCcut_use = vCcut;
			}
		}
	}
	return Tflag;
}


void reflect_across_bl(double2 vCcut, double2 vN, double2 *C1, double2 C2_C, double2 *dC)
{
	*C1 = vCcut + vN * 0.001;
	double2 R = (double2)(1.) - vN * vN * 2.;
	double xy = -2.*vN.x*vN.y;

	double2 PastBL = C2_C - (*C1);

	*dC = R * PastBL;
	*dC += xy * PastBL;
}


short find_SS_min(__global bLinks *BLlist, short bl)
{
	short SSdir = BLlist[bl].dir;
	short bl_next = bl + SSdir;
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


#ifdef OPENCL_VERSION_1_2
__kernel void find_umax(__global double2 *J_array, __global double *Umaxval, uint Jvals_start)
{
	Umaxval[0] = 0.;
	for (int i = 0; i < FULLSIZEY; i++)
	{
		int jj = i + Jvals_start;
		Umaxval[0] = MAX(Umaxval[0],J_array[jj].x);
	}
}
#else
__kernel void find_umax(__global double2 *J_array, __global double *Umaxval, uint Jvals_start)
{
	int i = get_global_id(0);
	int ii = i + Jvals_start;

	double Uval = J_array[ii].x;

	barrier(CLK_LOCAL_MEM_FENCE);

	Umaxval[0] = work_group_reduce_max(Uval);

}
#endif

double find_lineval(double Xloc, Trparam trP, __global double2 *J_array)
{
	if (Xloc < trP.Bottom_location)
	{
		return Xloc * J_array[trP.Uvals_start].x / trP.Bottom_location;
	}
	else if (Xloc < trP.Top_location)
	{
		double Xp = Xloc - trP.Bottom_location;
		double X1 = floor(Xp);
		int Xint = convert_int(X1) + trP.Uvals_start;
		double dX = Xp - X1;
		return J_array[Xint].x * (1. - dX) + dX * J_array[Xint + 1].x;
	}

	double Xp = Xloc - trP.Bottom_location;
	double X1 = floor(Xp);
	int Xint = convert_int(X1) + trP.Uvals_start;
	double dX = trP.Top_location - X1;
	double dX1 = trP.bval - Xp;
	return dX1 * J_array[Xint].x / dX;
}



int Get_Par_Type(double randval, __global Pparam *parP)
{
	int j = 0;

	while(j < NUM_PAR_SIZES - 1)
	{
		if (randval >= parP[j].D_dist && randval < parP[j + 1].D_dist)
			return j;
		j++;
	}
	return NUM_PAR_SIZES - 1;
}


//Called after sort.
//todo: include variable particle concentraion per TR particle to maintain constant mass flux
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RERELEASE, 1, 1))) 
void TR_release_par(__global par *P,
	__global uint2 *RandArray,
	Trparam trP,
	__global double2 *J_array,
	uint maxel,
	uint Conc_number,
	__global double *Umax_val,
	__global Pparam *parP)
{
	int i = get_global_id(0);

	if (i >= maxel)
		return;

	par Pcur = P[i];

	uint2 RandEl = RandArray[i];
	double Umv = Umax_val[0];
	double Bval = trP.bval;
	double2 randval;
	RandEl = MWC64X_NextUint2(RandEl, &randval);
	double Xt = randval.x * Bval;
	double Yt = randval.y * Umv;

	double Fxx = find_lineval(Xt, trP, J_array);
	
	while (Yt > Fxx)
	{
		RandEl = MWC64X_NextUint2(RandEl, &randval);
		Xt = randval.x * Bval;
		Yt = randval.y * Umv;
		Fxx = find_lineval(Xt, trP, J_array);
	}
	
	double randval_single;
	RandEl = MWC64X_NextUint(RandEl, &randval_single);
	RandArray[i] = RandEl;
	
	
	
	Pcur.type = Get_Par_Type(randval_single, parP);
	Pcur.timer = PAR_TIMER_START;
	Pcur.pos.y = Xt + trP.offset_y;
	Pcur.pos.x = trP.X_release;
	Pcur.Dep_Flag = -1;
	int2 Posi = convert_int2(Pcur.pos);
	Pcur.loc = Posi.x*FULLSIZEY_TR + Posi.y;
	Pcur.Num_rep = Conc_number;
	P[i] = Pcur;
}

//Called at beginning of TR solver.
//todo: currently is only updates shear removal for particles that were already deposited during last sort.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RERELEASE, 1, 1))) 
	void TR_shear_removal(__global bLinks *BLlist,
	__global par *P,
	__global Pparam *Pp,
	__global uint *BLdep,
	uint offset,
	uint maxel)
{
	int i = get_global_id(0);

	if (i >= maxel)
		return;

	i += offset;

	par Pcur = P[i];

	if (Pcur.Dep_Flag < 0)
		return;

	Pparam Parcur = Pp[Pcur.type];
	Pcur.Dep_timer--;
	if (Pcur.Dep_timer == 0)
	{
		int bldep_ind = Pcur.Dep_Flag*NUM_PAR_SIZES + Pcur.type;
		atomic_add(&BLdep[bldep_ind], Pcur.Num_rep);
		Pcur.Dep_Flag = -2;
		Pcur.loc = -2;
		P[i] = Pcur;
		return;
	}

	bLinks BLcur = BLlist[Pcur.Dep_Flag];
	double shear_surf = BLcur.Tau;

	if (shear_surf <= Parcur.tau_crit[BLcur.int_type])
	{
		P[i] = Pcur;
		return;
	}

	double Fdrag = Parcur.D_coeff;
	double Flift = Parcur.L_coeff * sqrt(shear_surf);
	
	double2 vLvec = Flift * BLcur.vNvec;
	double2 vDvec = (double)BLcur.dir * Fdrag * BLcur.vTvec;

	double2 vRvec = vLvec + vDvec;

	double2 C1 = BLcur.vP0 + BLcur.vTvec * BLcur.blLen / 2. + BLcur.vNvec * 0.001;

	double2 dC_rem = vRvec * shear_surf / Parcur.Mp * DTTR_WALL * DTTR_WALL / 2.;
	double dC_mag = length(dC_rem);
	dC_mag = (dC_mag > BLcur.blLen / 4.) ? (BLcur.blLen / 4.) : (dC_mag);
		
	Pcur.pos = C1 + dC_rem;
	Pcur.timer = PAR_TIMER_START;
	Pcur.Dep_Flag = (Pcur.pos.x < X_MIN_VAL || Pcur.pos.x > X_MAX_VAL) ? (-2) : (-1);
	int2 Posi = convert_int2(Pcur.pos);
	int pcur_loc = Posi.x*FULLSIZEY_TR + Posi.y;
	Pcur.loc = (Pcur.pos.x < X_MIN_VAL || Pcur.pos.x > X_MAX_VAL) ? (-2) : (pcur_loc);
	P[i] = Pcur;
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_update_nodes_along_wall_temp(__global double *T,
__global nodeC *nC,
__global nodeV *nV,
__global int *ActiveNodes,
__global int *Wall_inds,
int nWallNodes)
{
	int i = get_global_id(0);
	if (i >= nWallNodes)
	{
		return;
	}

	int j = Wall_inds[i];
	int ii = ActiveNodes[j];

	nodeC PN = nC[ii];

	double4 Tvals = (double4)(T[PN.neigh.x], T[PN.neigh.y], T[PN.neigh.z], T[PN.neigh.w]);

	nodeV nVtemp;
	nVtemp.Temps.x = dot(PN.CoeffT00, Tvals);
	nVtemp.Temps.y = dot(PN.CoeffT10, Tvals);
	nVtemp.Temps.z = dot(PN.CoeffT01, Tvals);
	nVtemp.Temps.w = dot(PN.CoeffT11, Tvals);

	nVtemp.Temps = nVtemp.Temps * TDIFF + TMIN;

	nV[ii] = nVtemp;
}

//Called after Shear removal kernel
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_update_nodes_along_wall_vel(__global double2 *U,
__global nodeC *nC,
__global nodeV *nV,
__global int *ActiveNodes,
__global int *Wall_inds,
int nWallNodes)
{
	int i = get_global_id(0);
	if (i >= nWallNodes)
	{
		return;
	}

	int j = Wall_inds[i];
	int ii = ActiveNodes[j];

	nodeC PN = nC[ii];
	nodeV nVtemp = nV[ii];

	double2 U00 = U[PN.neigh.x];
	double2 U10 = U[PN.neigh.y];
	double2 U01 = U[PN.neigh.z];
	double2 U11 = U[PN.neigh.w];
	double4 Uvels = (double4)(U00.x, U10.x, U01.x, U11.x);
	double4 Vvels = (double4)(U00.y, U10.y, U01.y, U11.y);

	nVtemp.U00 = (double2)(dot(PN.CoeffU00, Uvels), dot(PN.CoeffU00, Vvels));
	nVtemp.U10 = (double2)(dot(PN.CoeffU10, Uvels), dot(PN.CoeffU10, Vvels));
	nVtemp.U01 = (double2)(dot(PN.CoeffU01, Uvels), dot(PN.CoeffU01, Vvels));
	nVtemp.U11 = (double2)(dot(PN.CoeffU11, Uvels), dot(PN.CoeffU11, Vvels));

	nV[ii] = nVtemp;
}

//Called after Shear removal kernel
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR, 1, 1)))
void TR_update_nodes_temp(__global double *T,
__global nodeV *nV,
__global nodeI *nI,
__global nodeC *nC,
__global int *ActiveNodes,
int nActiveNodes)
{
	int i = get_global_id(0);
	if (i >= nActiveNodes)
	{
		return;
	}

	int ii = ActiveNodes[i];

	if (nI[ii].Wall_Flag > 0)
		return;
	
	int4 neigh = nC[ii].neigh;

	nV[ii].Temps = (double4)(T[neigh.x], T[neigh.y], T[neigh.z], T[neigh.w]) * TDIFF + TMIN;
}

//Called after Shear removal kernel
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR, 1, 1)))
void TR_update_nodes_vel(__global double2 *U,
__global nodeV *nV,
__global nodeI *nI,
__global nodeC *nC,
__global int *ActiveNodes,
int nActiveNodes)
{
	int i = get_global_id(0);
	if (i >= nActiveNodes)
	{
		return;
	}

	int ii = ActiveNodes[i];
	if (nI[ii].Wall_Flag > 0)
		return;

	int4 neigh = nC[ii].neigh;
	nodeV nVtemp = nV[ii];

	nVtemp.U00 = U[neigh.x];
	nVtemp.U10 = U[neigh.y];
	nVtemp.U01 = U[neigh.z];
	nVtemp.U11 = U[neigh.w];
	
	nV[ii] = nVtemp;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_update_par_along_wall(__global nodeV *nV,
__global int2 *Ploc_array,
__global par *P,
__global double2 *PV,
__global Pparam *Pp,
__global int *Neighs,
__global short *update_flag,
__global int *Wall_inds,
uint nWallNodes,
__global nodeI *nI,
__global bLinks *BLlist,
__global uint *index,
__global int *Out_array,
int neigh_ind)
{
	int i0 = get_global_id(0);
	if (i0 >= nWallNodes)
	{
		return;
	}

	int i = Wall_inds[i0];

	int jj = i * 9;
	int ii = Neighs[jj];

	__local nodeV nVtemp[WORKGROUPSIZE_TR_WALL];
	//__local par Pcur[WORKGROUPSIZE_TR_WALL];
	int lid = get_local_id(0);

	nVtemp[lid] = nV[ii];
	nodeI nItemp = nI[ii];

	jj += neigh_ind;

	int Node_loc = Neighs[jj];
	int2 ploc = Ploc_array[Node_loc + 2];

	if (Node_loc == -1 || ploc.x == -1)
		return;

	int kkt = ploc.x;

	while (kkt < ploc.y)
	{
		int kk = kkt;
		kkt++;
		par Pcur = P[kk];

		if (Pcur.loc != ii || update_flag[kk])
			continue;

		update_flag[kk] = 1;


		double2 dX1 = Pcur.pos - trunc(Pcur.pos);

		double2 Uadv = get_Vadv(nVtemp[lid].U00, nVtemp[lid].U10, nVtemp[lid].U01, nVtemp[lid].U11, dX1, (1. - dX1));
		Uadv += (Pcur.pos.x > START_THERMO_VEL) ? (get_Vth(nVtemp[lid].Temps, Pp[Pcur.type].Kth, dX1, (1. - dX1))) : (0.);

		double2 dC = (Uadv)* DTTR_WALL;

		double2 C1 = Pcur.pos;

		Pcur.pos = C1 + dC;

		int2 Posi = convert_int2(Pcur.pos);

		Pcur.loc = mad24(Posi.x, FULLSIZEY_TR, Posi.y);
		Pcur.Dep_Flag = -1;
		
		if (Pcur.pos.x >= X_MAX_VAL || Pcur.pos.x < X_MIN_VAL)
		{
			Pcur.loc = -2;
			Pcur.Dep_Flag = -2;
			P[kk] = Pcur;
			continue;
		}

		P[kk] = Pcur;

		int Sflag = (Pcur.pos.x <= X_START_VAL) ? (0) : (1);
		int bl;
		double2 vCcut;
	
		int Tflag = find_intersection(&vCcut, C1, dC, nItemp.BLind, BLlist, Sflag, &bl);

		if (Tflag == -1)
			continue;
		
		P[kk].timer -= 1;
		
		if (P[kk].timer <= 0)
		{
			P[kk].Dep_Flag = -2;
			P[kk].loc = -2;
			continue;
		}

		int ind_temp = atomic_inc(&index[Tflag]);

		int oset = TRC_NUM_TRACERS*Tflag;

		Out_array[oset + ind_temp * 2] = kk;
		Out_array[oset + ind_temp * 2 + 1] = bl;
		PV[oset + ind_temp * 2] = dC;
		PV[oset + ind_temp * 2 + 1] = vCcut;
	}
}



__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL_REFLECT, 1, 1)))
void TR_reflect_particles(__global bLinks *BLlist,
__global par *P,
__global double2 *PV,
__global int *Info,
int num_pars)
{
	int ii = get_global_id(0);

	if (ii >= num_pars)
		return;

	int lid = get_local_id(0);

	int par_ind = Info[ii * 2];
	par Pcur = P[par_ind];
	double2 dC = PV[ii * 2];
	double2 C1 = Pcur.pos - dC;

	if (Pcur.pos.x > X_MAX_VAL || Pcur.pos.x < (X_MIN_VAL) ||
		Pcur.pos.y < Y_MIN_VAL || Pcur.pos.y > Y_MAX_VAL)
	{
		P[par_ind].Dep_Flag = -2;
		P[par_ind].loc = -2;
		return;
	}

	reflect_across_bl(PV[ii * 2 + 1], BLlist[Info[ii * 2 + 1]].vNvec, &C1, Pcur.pos, &dC);
	Pcur.pos = C1 + dC;

	int2 Posi = convert_int2(Pcur.pos);
	Pcur.loc = (Pcur.pos.x < X_MIN_VAL) ? (-2) : (mad24(Posi.x, FULLSIZEY_TR, Posi.y));
	Pcur.Dep_Flag = (Pcur.pos.x < X_MIN_VAL) ? (-2) : (-1);
	P[par_ind] = Pcur;
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_update_par_contact_wall(__global bLinks *BLlist,
__global par *P,
__global double2 *PV,
__global Pparam *Pp,
__global uint2 *RandArray,
__global int *Info,
__global uint *index,
int offset,
int num_pars)
{
	int ii = get_global_id(0);

	if (ii >= num_pars)
		return;

	int par_ind = Info[offset + ii * 2];
	int bl = Info[offset + ii * 2 + 1];
	int lid = get_local_id(0);

	//__local par Pcur[WORKGROUPSIZE_TR_WALL];
	__local bLinks BLcur[WORKGROUPSIZE_TR_WALL];
	__local Pparam Parcur[WORKGROUPSIZE_TR_WALL];

	par Pcur = P[par_ind];

	if (Pcur.pos.x > X_MAX_VAL || Pcur.pos.x < (X_MIN_VAL) ||
		Pcur.pos.y < Y_MIN_VAL || Pcur.pos.y > Y_MAX_VAL)
	{
		P[par_ind].Dep_Flag = -2;
		P[par_ind].loc = -2;
		return;
	}

	Parcur[lid] = Pp[Pcur.type];
	
	BLcur[lid] = BLlist[bl];
	
	double shear_surf = BLcur[lid].Tau;
	
	short itype = BLcur[lid].int_type;

	int sh_flag = isgreaterequal(shear_surf, Parcur[lid].tau_crit[itype]);

	double2 Vimp_vec = PV[offset + ii * 2] / DTTR_WALL;

	double Vimp2 = dot(Vimp_vec, Vimp_vec);

	double Vreb2 = (sh_flag) ? ((shear_surf - Parcur[lid].tau_crit[itype]) / MU_NUMBER / 2. * Parcur[lid].Dp) :
		(Vimp2 - 2.*(Parcur[lid].Q_A_prime[itype] / Parcur[lid].Mp));

	double Sprob = (Vreb2 > 0.) ? (Parcur[lid].Q_A[itype] / (0.5 * Parcur[lid].Mp * Vreb2 - Parcur[lid].Q_A_prime[itype])) : (100.);
	
	double Rprob;
	RandArray[ii] = MWC64X_NextUint(RandArray[ii], &Rprob);

	Pcur.Dep_Flag = (Sprob < Rprob || sh_flag) ? (-1) : (bl);

	int ind_info = (Pcur.Dep_Flag > -1) ? 0 : 1;
	
	int ind_temp = atomic_inc(&index[ind_info]);

	Info[ind_temp * 2 + ind_info] = par_ind;



	if (Pcur.Dep_Flag == -1)
	{
		Vreb2 = sqrt(Vreb2) * DTTR_WALL;
		double2 dC = PV[offset + ii * 2];
		double2 C1 = Pcur.pos - dC;
		reflect_across_bl(PV[offset + ii * 2 + 1], BLcur[lid].vNvec, &C1, Pcur.pos, &dC);

		dC = normalize(dC)*Vreb2;
		Pcur.pos = C1 + dC;
		Pcur.Dep_Flag = bl;
		PV[ind_temp * 2] = dC;
		PV[ind_temp * 2 + 1] = C1;
	}

	P[par_ind] = Pcur;

}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_update_par_contact_wall2(__global bLinks *BLlist,
__global par *P,
__global double2 *PV,
__global int *Info,
__global uint *index,
int offset,
int num_pars)
{
	int ii = get_global_id(0);

	if (ii >= num_pars)
		return;

	int par_ind = Info[ii * 2 + 1];
	int lid = get_local_id(0);

	//__local par Pcur[WORKGROUPSIZE_TR_WALL];
	__local double2 dC[WORKGROUPSIZE_TR_WALL];
	__local double2 C1[WORKGROUPSIZE_TR_WALL];

	par Pcur = P[par_ind];
	dC[lid] = PV[ii * 2];
	C1[lid] = PV[ii * 2 + 1];

	if (Pcur.pos.x > X_MAX_VAL || Pcur.pos.x < (X_MIN_VAL) ||
		Pcur.pos.y < Y_MIN_VAL || Pcur.pos.y > Y_MAX_VAL)
	{
		P[par_ind].Dep_Flag = -2;
		P[par_ind].loc = -2;
		return;
	}

	double2 vCcut;
	double dist;
	int blneigh = (dC[lid].x >= 0) ? (1) : (-1);
	blneigh += Pcur.Dep_Flag;

	if (bcFindIntersectionLinePlane(&vCcut, &dist, C1[lid], dC[lid], BLlist[blneigh]))
	{
		int ind_temp = atomic_inc(&index[1]);

		Info[offset + ind_temp * 2] = par_ind;
		Info[offset + ind_temp * 2 + 1] = blneigh;
		PV[offset + ind_temp * 2] = dC[lid];
		PV[offset + ind_temp * 2 + 1] = vCcut;
	}

	int2 Posi = convert_int2(Pcur.pos);
	Pcur.loc = mad24(Posi.x, FULLSIZEY_TR, Posi.y);
	Pcur.Dep_Flag = -1;
	P[par_ind] = Pcur;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_deposit_particles_on_wall(__global bLinks *BLlist,
__global par *P,
__global int *Info,
int num_pars)
{
	int ii = get_global_id(0);

	if (ii >= num_pars)
		return;

	int par_ind = Info[ii * 2];
	int lid = get_local_id(0);

	//__local par Pcur[WORKGROUPSIZE_TR_WALL];

	/*Pcur[lid] = P[par_ind];

	Pcur[lid].Dep_Flag = find_SS_min(BLlist, Pcur[lid].Dep_Flag);
	Pcur[lid].Dep_timer = DEP_TIMER_START;
	bLinks BLdepcur = BLlist[Pcur[lid].Dep_Flag];
	Pcur[lid].pos = BLdepcur.vP0 + BLdepcur.vTvec*BLdepcur.blLen / 2.;
	Pcur[lid].loc = -1;
	P[par_ind] = Pcur[lid];*/

	par Pcur = P[par_ind];

	Pcur.Dep_Flag = find_SS_min(BLlist, Pcur.Dep_Flag);
	Pcur.Dep_timer = DEP_TIMER_START;
	bLinks BLdepcur = BLlist[Pcur.Dep_Flag];
	Pcur.pos = BLdepcur.vP0 + BLdepcur.vTvec*BLdepcur.blLen / 2.;
	Pcur.loc = -1;
	P[par_ind] = Pcur;
}



__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR, 1, 1)))
void TR_update_par_no_wall(__global nodeV *nV,
__global int2 *Ploc_array,
__global par *P,
__global Pparam *Pp,
__global int *Neighs,
__global short *update_flag,
__global nodeI *nI,
 uint nActiveNodes)//7
{
	int i = get_global_id(0);
	if (i >= nActiveNodes)
	{
		return;
	}
	int j = i * 9;
	int ii = Neighs[j];

	if (nI[ii].Wall_Flag > 0)
		return;

	nodeV nVtemp = nV[ii];
	
	double2 U00 = nVtemp.U00, U10 = nVtemp.U10, U01 = nVtemp.U01, U11 = nVtemp.U11;
	double4 Ttemp = nVtemp.Temps;

	int jj = j;

	while(jj < j+9)
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

			if (Pcur.loc != ii || update_flag[kk])
				continue;
			
			double2 dX1 = Pcur.pos - trunc(Pcur.pos);
			double2 dX0 = 1. - dX1;

			double2 Uadv = get_Vadv(U00, U10, U01, U11, dX1, dX0);
			Uadv += (Pcur.pos.x > START_THERMO_VEL) ? (get_Vth(Ttemp, Pp[Pcur.type].Kth, dX1, dX0)) : (0.);

			Pcur.pos += (Uadv)* DTTR;;

			int2 Posi = convert_int2(Pcur.pos);
			Pcur.loc = (Pcur.pos.x >= X_MAX_VAL || Pcur.pos.x < X_MIN_VAL) ? (-2) : (mad24(Posi.x,FULLSIZEY_TR,Posi.y));
			Pcur.Dep_Flag = (Pcur.pos.x >= X_MAX_VAL || Pcur.pos.x < X_MIN_VAL) ? (-2) : (-1);
			P[kk] = Pcur;
			update_flag[kk] = 1;
		}
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_GL, 1, 1))) 
void TR_opengl_par(__global par *P,
	__global float2 *gl_par,
	__global float *gl_color,
	__global float *gl_color_array)
{
	int i = get_global_id(0);


	if (i >= TRC_NUM_TRACERS)
		return;

	par Pcur = P[i];

	gl_par[i] = convert_float2(Pcur.pos);
	float Dep_mult = (Pcur.Dep_Flag > -2) ? (1.f) : (0.f);
	i *= 3;
	gl_color[i] = gl_color_array[Pcur.type * 3] * Dep_mult;
	gl_color[i+1] = gl_color_array[Pcur.type * 3+1] * Dep_mult;
	gl_color[i+2] = gl_color_array[Pcur.type * 3+2] * Dep_mult;
}

#ifdef USE_FLOATS
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1))) 
void TR_shear_1(__global double *Ro_array,
__global double2 *J_array,
__global double8 *CFSS,
#ifdef SOA_STORAGE
__global float *FA,
#else
__global float8 *FA,
#endif
__global double *Tau_array,
__global int2 *Loc,
int TAU_ARRAY_SIZE,
double dp,
__global float *alpha_array)
{
	uint ii = get_global_id(0);
	if (ii >= (TAU_ARRAY_SIZE))
		return;
	uint indU = Loc[ii].x;
	uint ind = Loc[ii].y;


#ifdef SOA_STORAGE
	double8 FAtemp;
	FAtemp.s0 = convert_double(FA[ind]);
	FAtemp.s1 = convert_double(FA[ind + FULLSIZEXY]);
	FAtemp.s2 = convert_double(FA[ind + 2 * FULLSIZEXY]);
	FAtemp.s3 = convert_double(FA[ind + 3 * FULLSIZEXY]);
	FAtemp.s4 = convert_double(FA[ind + 4 * FULLSIZEXY]);
	FAtemp.s5 = convert_double(FA[ind + 5 * FULLSIZEXY]);
	FAtemp.s6 = convert_double(FA[ind + 6 * FULLSIZEXY]);
	FAtemp.s7 = convert_double(FA[ind + 7 * FULLSIZEXY]);
#else
	double8 FAtemp = convert_double8(FA[ind]);
#endif
#ifdef USE_MRT_ELBM
	double relax_param = BETA_VAL;
#else
	double relax_param = alpha_array[ind] * BETA_VAL;
#endif
	double visc_val = 1. / 3. * (1. / relax_param - 0.5);
	double8 XCoeffs = CFSS[2 * ii] * relax_param*visc_val;
	double8 YCoeffs = CFSS[2 * ii + 1] * relax_param*visc_val;

	double2 U = J_array[indU];

	double Ro = Ro_array[indU];

	double2 op3u = sqrt(1. + 3.*U*U);
	double Psi = Ro * (2. - op3u.x)*(2. - op3u.y) / 9.;
	double2 B = (2. * U + op3u) / (1. - U);
	double2 Binv = 1. / B;

	FAtemp.s0 -= Psi*B.x;
	FAtemp.s1 -= Psi*Binv.x;
	FAtemp.s2 -= Psi*B.y;
	FAtemp.s3 -= Psi*Binv.y;
	FAtemp.s4 -= 0.25*Psi*B.x*B.y;
	FAtemp.s5 -= 0.25*Psi*Binv.x*Binv.y;
	FAtemp.s6 -= 0.25*Psi*B.x*Binv.y;
	FAtemp.s7 -= 0.25*Psi*Binv.x*B.y;
	double2 Tau;
	Tau.x = dot(FAtemp.lo, XCoeffs.lo) + dot(FAtemp.hi, XCoeffs.hi);
	Tau.y = dot(FAtemp.lo, YCoeffs.lo) + dot(FAtemp.hi, YCoeffs.hi);
	double Taumag = length(Tau) / Ro;
	if (U.x < 0.)
	{
		Tau_array[ii] = -Taumag;
	}
	else
	{
		Tau_array[ii] = Taumag;
	}
}

#else

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1))) 
void TR_shear_1(__global double *Ro_array,
__global double2 *J_array,
__global double8 *CFSS,
#ifdef SOA_STORAGE
__global double *FA,
#else
__global double8 *FA,
#endif
__global double *Tau_array,
__global int2 *Loc,
int TAU_ARRAY_SIZE,
double dp,
__global double *alpha_array)
{
	uint ii = get_global_id(0);
	if (ii >= (TAU_ARRAY_SIZE))
		return;
	uint indU = Loc[ii].x;
	uint ind = Loc[ii].y;


#ifdef SOA_STORAGE
	double8 FAtemp;
	FAtemp.s0 = (FA[ind]);
	FAtemp.s1 = (FA[ind + FULLSIZEXY]);
	FAtemp.s2 = (FA[ind + 2 * FULLSIZEXY]);
	FAtemp.s3 = (FA[ind + 3 * FULLSIZEXY]);
	FAtemp.s4 = (FA[ind + 4 * FULLSIZEXY]);
	FAtemp.s5 = (FA[ind + 5 * FULLSIZEXY]);
	FAtemp.s6 = (FA[ind + 6 * FULLSIZEXY]);
	FAtemp.s7 = (FA[ind + 7 * FULLSIZEXY]);
#else
	double8 FAtemp = FA[ind];
#endif
#ifdef USE_MRT_ELBM
	double relax_param = BETA_VAL;
#else
	double relax_param = alpha_array[ind] * BETA_VAL;
#endif
	double visc_val = 1. / 3. * (1. / relax_param - 0.5);
	double8 XCoeffs = CFSS[2 * ii] * relax_param*visc_val;
	double8 YCoeffs = CFSS[2 * ii + 1] * relax_param*visc_val;

	double2 U = J_array[indU];

	double Ro = Ro_array[indU];
	double2 op3u = sqrt(1. + 3.*U*U);
	double Psi = Ro * (2. - op3u.x)*(2. - op3u.y) / 9.;
	double2 B = (2. * U + op3u) / (1. - U);
	double2 Binv = 1. / B;

	FAtemp.s0 -= Psi*B.x;
	FAtemp.s1 -= Psi*Binv.x;
	FAtemp.s2 -= Psi*B.y;
	FAtemp.s3 -= Psi*Binv.y;
	FAtemp.s4 -= 0.25*Psi*B.x*B.y;
	FAtemp.s5 -= 0.25*Psi*Binv.x*Binv.y;
	FAtemp.s6 -= 0.25*Psi*B.x*Binv.y;
	FAtemp.s7 -= 0.25*Psi*Binv.x*B.y;
	double2 Tau;
	Tau.x = dot(FAtemp.lo, XCoeffs.lo) + dot(FAtemp.hi, XCoeffs.hi);
	Tau.y = dot(FAtemp.lo, YCoeffs.lo) + dot(FAtemp.hi, YCoeffs.hi);
	double Taumag = length(Tau) / Ro;
	if (U.x < 0.)
	{
		Tau_array[ii] = -Taumag;
	}
	else
	{
		Tau_array[ii] = Taumag;
	}
}
#endif

#ifdef USE_OPENGL
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
 void TR_shear_2(__global int *BLind,
	__global double4 *Weights,
	__global double *Tau_array,
	__global int4 *Sind,
	__global bLinks *BL_array,
	int NUMBLINKS,
	__global float *ls_color_b,
	__global float *ls_color_t,
	__global double *TCRIT_MAX2)
{
	uint ii = get_global_id(0);
	if (ii >= NUMBLINKS)
		return;
	int jj = BLind[ii];
	bLinks BLtemp = BL_array[jj];
	int4 i = Sind[ii];
	double4 weight = Weights[ii];
	float fr = 1.f, fg = 1.f, fb = 1.f;
	double4 Tautemp = (double4)(Tau_array[i.x], Tau_array[i.y], Tau_array[i.z], Tau_array[i.w]);
	double Tau_wall = dot(weight, Tautemp);

	double TCRIT_MAX = TCRIT_MAX2[BLtemp.int_type];

	if (Tau_wall < 0.)
	{
		Tau_wall *= -1.;
		BLtemp.dir = -1;
		BLtemp.Tau = Tau_wall;
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
		BLtemp.dir = 1;
		BLtemp.Tau = Tau_wall;
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

	BL_array[jj] = BLtemp;

	int LSind = BLtemp.Color_ind;
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
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1))) void TR_shear_2(__global int *BLind,
	__global double4 *Weights,
	__global double *Tau_array,
	__global int4 *Sind,
	__global bLinks *BL_array,
	int NUMBLINKS)
{
	uint ii = get_global_id(0);
	if (ii >= NUMBLINKS)
		return;
	int jj = BLind[ii];
	double4 weight = Weights[ii];
	int4 i = Sind[ii];
	bLinks BLtemp = BL_array[jj];
	double4 Tautemp = (double4)(Tau_array[i.x], Tau_array[i.y], Tau_array[i.z], Tau_array[i.w]);
	double Tau_wall = dot(weight, Tautemp);

	if (Tau_wall < 0.)
	{
		BLtemp.dir = -1;
		BLtemp.Tau = -Tau_wall;
	}
	else
	{
		BLtemp.dir = 1;
		BLtemp.Tau = Tau_wall;
	}

	BL_array[jj] = BLtemp;
}
#endif

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1))) 
	void TR_Find_Shear_Coeffs1(__global int3 *Bindicies_array,//0
	__global int *Stor,//1
	__global int *dir_array,//2
	__global int2 *ii0_array,//3
	__global double *CFSS,//4
	__global int2 *Loc,//5
	__global int *Bind_ind,//6
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
	int kk = (ii >= Bindicies_top_start) ? (Bindicies.x*2 + 1) : (Bindicies.x*2);

	atomic_xchg(&Bind_ind[kk], ii);

	Loc[ii] = (int2)(Bindicies.x * FULLSIZEY_UT, Bindicies.x * FULLSIZEY);
	Loc[ii] += Stor[Bindicies.x*DOMAIN_SIZE_Y+Bindicies.y];
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
	
	int coeffs_ind = ii*16;

	for (int aa = 0; aa < 8; aa++)
	{
		double CCx = CXY_DOUBLE[aa].x, CCy = CXY_DOUBLE[aa].y;
		double sumk = CCx*norm_vec.x*norm_vec.x;
		sumk += CCy*norm_vec.x*norm_vec.y;

		double sumj = CCx * norm_vec.x * (CCx - sumk);
		sumj += CCy * norm_vec.y * (CCx - sumk);

		//double taunumber = MU_NUMBER*3. + 0.5;

		//CFSS[coeffs_ind + aa] = sumj*MU_NUMBER*3. / taunumber;
		CFSS[coeffs_ind + aa] = sumj*3.;
		double sumky = CCx*norm_vec.y*norm_vec.x;
		sumky += CCy*norm_vec.y*norm_vec.y;

		double sumjy = CCx * norm_vec.x * (CCy - sumky);
		sumjy += CCy * norm_vec.y * (CCy - sumky);

		CFSS[coeffs_ind + 8 + aa] = sumjy*3.;
		//CFSS[coeffs_ind + 8 + aa] = sumjy*MU_NUMBER*3. / taunumber;
	}
}

void reorganize_values(double4 *dRvals, int4 *index_neigh)
{
	if((*dRvals).x >= (*dRvals).y)
		return;
	
	(*dRvals).xy = (*dRvals).yx;
	(*index_neigh).xy = (*index_neigh).yx;
	
	if((*dRvals).y >= (*dRvals).z)
		return;
		
	(*dRvals).yz = (*dRvals).zy;
	(*index_neigh).yz = (*index_neigh).zy;
	
	if((*dRvals).z >= (*dRvals).w)
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
	return 1. / sqrt(2.*PI_NUMBER) * exp(-0.5*Sval*Sval);
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1))) 
void TR_Find_Shear_Coeffs2(__global int *BLind,//0
	__global double4 *Weights,//1
	__global int4 *Sind,//2
	__global bLinks *BL,//3
	__global int3 *Bindicies,//4
	__global int2 *Bind_ind,//5
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
	
	double2 BLmiddle = BLcur.vP0 + BLcur.vTvec * BLcur.blLen/2.;
	int izz = convert_int(BLmiddle.x);
	int ind = Bind_ind[izz].x;
	while(ind == -1)
	{
		izz++;
		ind = Bind_ind[izz].x;
	}
	int start_ind = MAX(0, ind-Bindicies_radius);
	int stop_ind = MIN(Bindicies_end.x, ind+Bindicies_radius);
	
	//Bindicies_end.x is starting index of Bindicies located at top
	//Bindicies_end.y is total length of Bindicies
	
	if(i >= BLinks_top_ind)
	{
		ind = Bind_ind[izz].y;
		while(ind == -1)
		{
			izz--;
			ind = Bind_ind[izz].y;
		}
		start_ind = MAX(Bindicies_end.x, ind-Bindicies_radius);
		stop_ind = MIN(Bindicies_end.y, ind+Bindicies_radius);
	}
	
	int4 index_neigh = -1;
	double4 dRvals = 100.;
	
	int kk = start_ind;
	while(kk < stop_ind)
	{
		double2 ii0 = convert_double2(Bindicies[kk].xy);
		double2	dR = BLmiddle - ii0;
		double dist = length(dR);
		if(dist < dRvals.x)
		{
			dRvals.x = dist;
			index_neigh.x = kk;
			reorganize_values(&dRvals, &index_neigh);
		}
			
		kk++;
	}
	
	double4 weight_temp;
	weight_temp.x = weight_kernel_func(dRvals.x/cutoff_radius);
	weight_temp.y = weight_kernel_func(dRvals.y/cutoff_radius);
	weight_temp.z = weight_kernel_func(dRvals.z/cutoff_radius);
	weight_temp.w = weight_kernel_func(dRvals.w/cutoff_radius);

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
	
double get_offset(__global double2 *C0, double xval)
{
	int i = 0;
	while (TRUE)
	{
		if (C0[i].x < xval && C0[i + 1].x >= xval)
			break;
		i++;
	}
	return C0[i].y + (xval - C0[i].x) * (C0[i + 1].y - C0[i].y) / (C0[i + 1].x - C0[i].x);
}

__kernel void Rerelease_Clumps(__global par *P,
__global double2 *C0,
Trparam trP,
__global uint2 *RandArray,
int start_ind)
{
	int i = get_global_id(0) + start_ind;
	if (i >= TRC_NUM_TRACERS)
		return;
	
	par Ptemp;
	double2 Rprob;
	RandArray[i] = MWC64X_NextUint2(RandArray[i], &Rprob);
	Ptemp.pos.x = X_START_VAL + (0.75 + 0.25*Rprob.x)*(X_MAX_VAL - X_START_VAL);
	Ptemp.pos.y = get_offset(C0, Ptemp.pos.x) + 7.*trP.bval / 16. + Rprob.y*trP.bval / 8.;

	Ptemp.type = 0;
	Ptemp.timer = 1;
	Ptemp.Dep_Flag = -1;
	Ptemp.Num_rep = 0;
	Ptemp.Dep_timer = 1;
	int2 Posi = convert_int2(Ptemp.pos);
	Ptemp.loc = mad24(Posi.x, FULLSIZEY_TR, Posi.y);
	P[i] = Ptemp;
}


__kernel void Test_Bounds_Clumped_Particles(__global par *P,
	__global bLinks *BL,
	Trparam trP,
	int nBL,
	int num_pars,
	int start_ind)
{
	int i = get_global_id(0) + start_ind;

	if (i >= num_pars)
		return;
	par Pcur = P[i];
		
	double2 Parpos = Pcur.pos;

	if (Pcur.Dep_Flag == -2 || Pcur.Dep_Flag > -1)
		return;

	if (Parpos.x <= X_MIN_VAL || Parpos.x > X_MAX_VAL || isnan(Parpos.x))
	{
		P[i].Dep_Flag = -2;
		P[i].loc = -2;
		P[i].pos = (double2)(X_START_VAL, trP.offset_y + trP.bval / 2.);
		return;
	}

	int bl_ind_temp = 0;

	while (1)
	{
		if (BL[bl_ind_temp].vP0.x < Parpos.x && BL[bl_ind_temp].vP1.x >= Parpos.x)
			break;
		bl_ind_temp++;
	}

	int bl_bot_ind = bl_ind_temp;

	bl_ind_temp += (nBL / 2 - 10);

	while (1)
	{
		if (BL[bl_ind_temp].vP0.x < Parpos.x && BL[bl_ind_temp].vP1.x >= Parpos.x)
			break;

		bl_ind_temp++;
	}

	int bl_top_ind = bl_ind_temp;
		
	double2 V0 = BL[bl_bot_ind].vP0, V1 = BL[bl_bot_ind].vP1;
		
	double ybot_pos = V0.y + (Parpos.x - V0.x) * (V1.y - V0.y) / (V1.x - V0.x);

	V0 = BL[bl_top_ind].vP0, V1 = BL[bl_top_ind].vP1;
		
	double ytop_pos = V0.y + (Parpos.x - V0.x) * (V1.y - V0.y) / (V1.x - V0.x);

	if (Parpos.y <= ybot_pos || Parpos.y >= ytop_pos)
	{
		P[i].Dep_Flag = -2;
		P[i].loc = -2;
		P[i].pos = (double2)(X_START_VAL, trP.offset_y + trP.bval / 2.);
	}
}



