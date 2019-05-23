
//#define USE_MAD

double2 get_Vadv(double2 U00, double2 U10, double2 U01, double2 U11, double2 dX1, double2 dX0)
{
//	double2 U0y = fma(U00, dX0.y, U01 * dX1.y);
//	double2 U1y = fma(U10, dX0.y, U11 * dX1.y);
//
	return fma(fma(U00, dX0.y, U01 * dX1.y), dX0.x, fma(U10, dX0.y, U11 * dX1.y) * dX1.x);
}

double2 get_Vth(double4 TT, double Kth, double2 dX1, double2 dX0)
{
	double T0y = fma(TT.x, dX0.y, TT.z * dX1.y);
	double T1y = fma(TT.y, dX0.y, TT.w * dX1.y);
	double Txy = fma(T0y, dX0.x, T1y * dX1.x);

	double T0x = fma(TT.x, dX0.x, TT.y * dX1.x);
	double T1x = fma(TT.z, dX0.x, TT.w * dX1.x);

	double2 dTdx = (double2)(T0y - T1y, 
		fma(TT.x, dX0.x, TT.y * dX1.x) - 
		fma(TT.z, dX0.x, TT.w * dX1.x));  ///-dTdx
	return dTdx * MU_NUMBER * Kth / Txy;
}

int bcFindIntersectionLinePlane(double2 *vC, double *dist, double2 vL0, double2 vLd, bLinks BLcur)
{
	double2 vN = BLcur.vNvec;
	double den = dot(vN, vLd);
	if (den == 0.)
		return false;

	*dist = dot(vN, (BLcur.vP0 - vL0)) / den;
	*vC = vL0 + vLd * (*dist);

	if ((*dist) >= 0. && (*dist) <= 1.)
	{
		double vd0 = distance(*vC, BLcur.vP0), vd1 = distance(*vC, BLcur.vP1);
		if (fabs(BLcur.blLen - vd0 - vd1) < CEPS)
			return true;
	}
	return false;
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
		if (bcFindIntersectionLinePlane(&vCcut, &dist, C1, dC, BLlist[bl]) == true)
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




