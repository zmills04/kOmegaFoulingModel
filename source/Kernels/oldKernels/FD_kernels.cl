__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
void FD_update(__global double *T,
__global double4 *Aind,
__global double *AAA,
__global double *BBB,
__global double *AAA_base,
__global double *BBB_base,
__global double2 *A,
__global double2 *B,
__global double2 *C,
__global double2 *J_array,
double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= FULLSIZEX || j >= FULLSIZEY)
	{
		return;
	}
	int indUT = i * FULLSIZEY_UT + j;
	int indc = i * FULLSIZEY + j;

	double4 indM = Aind[indc];
	int indc5 = indc * 5;

	double2 U = J_array[indUT];
	
	double2 Avec = A[indc*2+1], Bvec = B[indc*2+1], Cvec = C[indc*2+1];

	double Btemp = (i > 0) ? (T[indUT]) : (T[indUT] - TIN * DTFD * ( U.x * Avec.x ));

	BBB[indc] = BBB_base[indc] + Btemp;

	AAA[indc5] = AAA_base[indc5] + DTFD*(U.x*Cvec.x + U.y*Cvec.y);

	AAA[indc5 + 1] = AAA_base[indc5 + 1] + indM.x * DTFD * U.x * Avec.x;

	AAA[indc5 + 2] = AAA_base[indc5 + 2] + indM.y * DTFD * U.x * Bvec.x;

	AAA[indc5 + 3] = AAA_base[indc5 + 3] + indM.z * DTFD * U.y * Avec.y;

	AAA[indc5 + 4] = AAA_base[indc5 + 4] + indM.w * DTFD * U.y * Bvec.y;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
void FD_update_base(__global double *Alpha,
__global double4 *Aind,
__global double *AAA_base,
__global double *BBB_base,
__global double4 *A,
__global double4 *B,
__global double4 *C,
double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= FULLSIZEX || j >= FULLSIZEY)
	{
		return;
	}

	int indc = i * FULLSIZEY + j;

	double4 indM = Aind[indc];
	int indc5 = indc * 5;

	double4 Avec = A[indc], Bvec = B[indc], Cvec = C[indc];

	BBB_base[indc] = (i > 0) ? (0.) : (TIN * DTFD * (Alpha[indc5 + 2] * Avec.x));

	AAA_base[indc5] = 1. - DTFD*Alpha[indc5] * (Cvec.x + Cvec.y);

	AAA_base[indc5 + 1] = -indM.x * DTFD * (Alpha[indc5 + 2] * Avec.x);

	AAA_base[indc5 + 2] = -indM.y * DTFD * (Alpha[indc5 + 1] * Bvec.x);

	AAA_base[indc5 + 3] = -indM.z * DTFD * (Alpha[indc5 + 4] * Avec.y);

	AAA_base[indc5 + 4] = -indM.w * DTFD * (Alpha[indc5 + 3] * Bvec.y);
}


//Updates the A and b matricies used to calculate steady state value of T
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
void FD_update_SS(__global double *T,
__global double *Alpha,
__global double4 *Aind,
__global double *AAA,
__global double *BBB,
__global double4 *A,
__global double4 *B,
__global double4 *C,
__global double2 *J_array,
double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= FULLSIZEX || j >= FULLSIZEY)
	{
		return;
	}
	int indUT = i * FULLSIZEY_UT + j;
	int indc = i * FULLSIZEY + j;
	double indw = Aind[indc].x, inde = Aind[indc].y;
	double inds = Aind[indc].z, indn = Aind[indc].w;
	int indc5 = indc * 5;
	double alpha_c = Alpha[indc5], alpha_e = Alpha[indc5 + 1], alpha_w = Alpha[indc5 + 2];
	double alpha_n = Alpha[indc5 + 3], alpha_s = Alpha[indc5 + 4];

	double2 U = J_array[indUT];

	double A2x = A[indc].x, B2x = B[indc].x, C2x = C[indc].x;
	double A2y = A[indc].y, B2y = B[indc].y, C2y = C[indc].y;
	double A1x = A[indc].z, B1x = B[indc].z, C1x = C[indc].z;
	double A1y = A[indc].w, B1y = B[indc].w, C1y = C[indc].w;

	AAA[indc5 + 4] = indn*(alpha_n*B2y - U.y*B1y);
	AAA[indc5 + 3] = inds*(alpha_s*A2y - U.y*A1y);


	if (i == 0)
	{
		BBB[indc] = -TIN * (alpha_w*A2x - U.x*A1x);
		AAA[indc5] = alpha_c*(C2x + C2y) - (U.x*C1x + U.y*C1y);
		AAA[indc5 + 1] = 0.;
		AAA[indc5 + 2] = alpha_e*B2x - U.x*B1x;

	}
	else if (i == FULLSIZEXM1)
	{
		BBB[indc] = 0.;
		AAA[indc5] = alpha_c*(C2x + C2y) - (U.x*C1x + U.y*C1y) + (alpha_e*B2x - U.x*B1x);
		AAA[indc5 + 1] = (alpha_w*A2x - U.x*A1x);
		AAA[indc5 + 2] = 0.;
	}
	else
	{
		BBB[indc] = 0.;
		AAA[indc5] = alpha_c*(C2x + C2y) - (U.x*C1x + U.y*C1y);
		AAA[indc5 + 1] = indw * (alpha_w*A2x - U.x*A1x);
		AAA[indc5 + 2] = inde * (alpha_e*B2x - U.x*B1x);
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
void FD_solve_no_res(__global double *T,
__global double *Tout,
__global double *AAA,
__global double *BBB,
__global int4 *Neighs)
{
	int ii = get_global_id(0);
	int jj = get_global_id(1);
	
	if (ii >= FULLSIZEX || jj >= FULLSIZEY)
	{
		return;
	}

	int indc = ii*FULLSIZEY + jj;
	int indc5 = indc * 5;
	int indUT = mad24(ii, FULLSIZEY_UT, jj);

	int4 Neigh = Neighs[indc];

	double Te = (Neigh.x > -1) ? (T[Neigh.x]) : (0.);
	double Tw = (Neigh.y > -1) ? (T[Neigh.y]) : (0.);
	double Tn = (Neigh.z > -1) ? (T[Neigh.z]) : (0.);
	double Ts = (Neigh.w > -1) ? (T[Neigh.w]) : (0.);

	double sum = AAA[indc5 + 1] * Tw + AAA[indc5 + 2] * Te + AAA[indc5 + 3] * Ts + AAA[indc5 + 4] * Tn;
	Tout[indUT] = (BBB[indc] - sum) / AAA[indc5];
}


#ifdef OPENCL_VERSION_1_2
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
void FD_solve(__global double *T,
__global double *Tout,
__global double *AAA,
__global double *BBB,
__global double *res_total,
__global int4 *Neighs)
{
	const int ii = get_global_id(0);
	const int jj = get_global_id(1);
	if (ii >= FULLSIZEX || jj >= FULLSIZEY)
	{
		return;
	}


	int indc = ii*FULLSIZEY + jj;
	int indc5 = indc * 5;
	int indUT = mad24(ii, FULLSIZEY_UT, jj);
	double oldval = T[indUT];
	int4 Neigh = Neighs[indc];

	double Te = (Neigh.x > -1) ? (T[Neigh.x]) : (0.);
	double Tw = (Neigh.y > -1) ? (T[Neigh.y]) : (0.);
	double Tn = (Neigh.z > -1) ? (T[Neigh.z]) : (0.);
	double Ts = (Neigh.w > -1) ? (T[Neigh.w]) : (0.);

	double sum = AAA[indc5 + 1] * Tw + AAA[indc5 + 2] * Te + AAA[indc5 + 3] * Ts + AAA[indc5 + 4] * Tn;

	double newval = (BBB[indc] - sum) / AAA[indc5];
	double resid = fabs(oldval - newval);
	Tout[indUT] = newval;
	res_total[indc] = resid;
}
#else
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
void FD_solve(__global double *T,
__global double *Tout,
__global double *AAA,
__global double *BBB,
__local volatile double *local_res,
__global double *res_total,
__global int4 *Neighs)
{
	const int ii = get_global_id(0);
	const int jj = get_global_id(1);

	if (ii >= FULLSIZEX || jj >= FULLSIZEY)
	{
		return;
	}


	const int lid = get_local_linear_id();
	const uint group_size = get_local_size(0)*get_local_size(1);

	int indc = ii*FULLSIZEY + jj;
	int indc5 = indc * 5;
	int indUT = mad24(ii, FULLSIZEY_UT, jj);
	double oldval = T[indUT];
	int4 Neigh = Neighs[indc];

	double Te = (Neigh.x > -1) ? (T[Neigh.x]) : (0.);
	double Tw = (Neigh.y > -1) ? (T[Neigh.y]) : (0.);
	double Tn = (Neigh.z > -1) ? (T[Neigh.z]) : (0.);
	double Ts = (Neigh.w > -1) ? (T[Neigh.w]) : (0.);

	double sum = AAA[indc5 + 1] * Tw + AAA[indc5 + 2] * Te + AAA[indc5 + 3] * Ts + AAA[indc5 + 4] * Tn;

	double newval = (BBB[indc] - sum) / AAA[indc5];
	double resid = fabs(oldval - newval);
	Tout[indUT] = newval;

	local_res[lid] = resid;
	barrier(CLK_LOCAL_MEM_FENCE);

	int i = group_size / 2;
	for (; i>WAVEFRONT_SIZE; i >>= 1)
	{
		if (lid < i)
			local_res[lid] = resid = resid + local_res[lid + i];
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (lid < i)
		resid = sub_group_reduce_add(resid + local_res[lid + i]);

	if (lid == 0)
		AtomicAdd(res_total, resid);
}
#endif
//Under/over relaxed iterative solver. Does not save residuals
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
	void FD_relaxed_solve_no_res(__global double *T,
	__global double *Tout,
	__global double *AAA,
	__global double *BBB,
	double relax_factor,
	__global int4 *Neighs)
{
	int ii = get_global_id(0);
	int jj = get_global_id(1);
	if (ii >= FULLSIZEX || jj >= FULLSIZEY)
	{
		return;
	}

	int indc = ii*FULLSIZEY + jj;
	int indc5 = indc * 5;
	int4 Neigh = Neighs[indc];

	double Te = (Neigh.x > -1) ? (T[Neigh.x]) : (0.);
	double Tw = (Neigh.y > -1) ? (T[Neigh.y]) : (0.);
	double Tn = (Neigh.z > -1) ? (T[Neigh.z]) : (0.);
	double Ts = (Neigh.w > -1) ? (T[Neigh.w]) : (0.);

	double sum = AAA[indc5 + 1] * Tw + AAA[indc5 + 2] * Te + AAA[indc5 + 3] * Ts + AAA[indc5 + 4] * Tn;
	
	double oldval = Tout[ii*FULLSIZEY_UT + jj];

	Tout[ii*FULLSIZEY_UT + jj] = (1. - relax_factor) * oldval + relax_factor * 
		(BBB[indc] - sum) / AAA[indc5];
}

//Under/over relaxed iterative solver. Saves residuals
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
	void FD_relaxed_solve(__global double *T,
	__global double *Tout,
	__global double *AAA,
	__global double *BBB,
	__global double *residuals,
	double relax_factor,
	__global int4 *Neighs)
{
	int ii = get_global_id(0);
	int jj = get_global_id(1);
	if (ii >= FULLSIZEX || jj >= FULLSIZEY)
	{
		return;
	}

	int indc = ii*FULLSIZEY + jj;
	int indc5 = indc * 5;
	int4 Neigh = Neighs[indc];

	double Te = (Neigh.x > -1) ? (T[Neigh.x]) : (0.);
	double Tw = (Neigh.y > -1) ? (T[Neigh.y]) : (0.);
	double Tn = (Neigh.z > -1) ? (T[Neigh.z]) : (0.);
	double Ts = (Neigh.w > -1) ? (T[Neigh.w]) : (0.);

	double sum = AAA[indc5 + 1] * Tw + AAA[indc5 + 2] * Te + AAA[indc5 + 3] * Ts + AAA[indc5 + 4] * Tn;

	double oldval = T[ii*FULLSIZEY_UT + jj];

	double newval = (1. - relax_factor) * oldval + relax_factor * (BBB[indc] - sum) / AAA[indc5];
	Tout[ii*FULLSIZEY_UT + jj] = newval;

	residuals[indc] = fabs(oldval - newval);
}

