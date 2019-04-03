//Reduction kernel used to sum array
//Reduces entire array to smaller array with num elements = NUMBLOCKS
//can be called recursively until reduced to single element or followed by
//reduce_generic_2 to sum NUMBLOCK_elements into single element
__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
void reduce_generic(__global double *input,
	__global double *output,
	__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = input[stride] + input[stride + 1];

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid] = sdata[0];
	}

}


__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
void reduce_generic_absmax(__global double *input,
__global double *output,
__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = max(fabs(input[stride]), fabs(input[stride + 1]));

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = max(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid] = sdata[0];
	}

}

__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
void reduce_generic_max(__global double *input,
__global double *output,
__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = max(input[stride],input[stride + 1]);

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = max(sdata[tid],sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid] = sdata[0];
	}

}

__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
void reduce_generic_min(__global double *input,
__global double *output,
__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = min(input[stride], input[stride + 1]);

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = min(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid] = sdata[0];
	}

}


__kernel void reduce_generic_max_2_updated(__global double *input,
	__global double *output,
	__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = max(input[stride],input[stride + 1]);

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = max(sdata[tid],sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[0] = sdata[0];
	}
}

__kernel void reduce_generic_absmax_2_updated(__global double *input,
	__global double *output,
	__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = max(fabs(input[stride]), fabs(input[stride + 1]));

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = max(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[0] = sdata[0];
	}
}

// Need to Test
__kernel void reduce_generic_min_2_updated(__global double *input,
	__global double *output,
	__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = min(input[stride], input[stride + 1]);

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = min(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[0] = sdata[0];
	}
}

// Need to Test
__kernel void reduce_generic_max_2_updated_multi(__global double *input,
	__global double *output,
	__local double *sdata,
	int saveloc)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = max(input[stride],input[stride + 1]);

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = max(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[saveloc] = sdata[0];
	}
}

// Need to Test
__kernel void reduce_generic_min_2_updated_multi(__global double *input,
	__global double *output,
	__local double *sdata,
	int saveloc)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = min(input[stride], input[stride + 1]);

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = min(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[saveloc] = sdata[0];
	}
}

// Need to Test
__kernel void reduce_generic_absmax_2_updated_multi(__global double *input,
	__global double *output,
	__local double *sdata,
	int saveloc)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = max(fabs(input[stride]), fabs(input[stride + 1]));

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = max(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[saveloc] = sdata[0];
	}
}

// Need to Test
__kernel void reduce_generic_2_updated(__global double *input,
	__global double *output,
	__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = input[stride] + input[stride + 1];

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[0] = sdata[0];
	}
}


// Need to Test
__kernel void reduce_generic_2_updated_multi(__global double *input,
	__global double *output,
	__local double *sdata,
	int saveloc)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = input[stride] + input[stride + 1];

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[saveloc] = sdata[0];
	}
}



__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
void reduce_generic_vec_x(__global double2 *input,
__global double *output,
__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = input[stride].x + input[stride + 1].x;

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid] = sdata[0];
	}

}



//Reduces output of reduce1 to single value using single work-item
__kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void reduce_generic_2(__global double *input,
__global double *output)
{
	double temp = input[0];
	for (unsigned int s = 1; s < NUMBLOCKS; s++)
	{
		temp += input[s];
	}
	output[0] = (temp);
}

//First reduction step to reduce both Ro and J arrays
//Used to generate data for output arrays.
__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
void LB_reduce_Ro_J(__global double *inputRo,
__global double2 *inputJ,
__global double *output,
__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;

	sdata[tid * 2] = inputRo[stride] + inputRo[(stride + 1)];
	sdata[tid * 2 + 1] = inputJ[stride].x + inputJ[(stride + 1)].x;

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid * 2] += sdata[(tid + s) * 2];
			sdata[tid * 2 + 1] += sdata[(tid + s) * 2 + 1];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid * 2] = sdata[0];
		output[bid * 2 + 1] = sdata[1];
	}
}


#ifndef OPENCL_VERSION_1_2
__kernel  __attribute__((reqd_work_group_size(WORKGROUPSIZE_RED, 1, 1)))
void reductionSubgrp(__global double *input,
__global double *output,
__local volatile double *localBuffer,
int skip, int start,
__global double *Mass)
{
	const uint id = get_global_id(0);
	const uint lid = get_local_id(0);
	const uint group_size = get_local_size(0);

	// initialize shared memory contents
	double res = input[id*skip + start] / Mass[0];
	localBuffer[lid] = res;
	barrier(CLK_LOCAL_MEM_FENCE);

	// local memory reduction
	int i = group_size / 2;
	for (; i>WAVEFRONT_SIZE; i >>= 1) {
		if (lid < i)
			localBuffer[lid] = res = res + localBuffer[lid + i];
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// subgroup reduction (introduced in OpenCL 2.0)
	if (lid < i)
		res = sub_group_reduce_add(res + localBuffer[lid + i]);

	// atomic reduce in global memory
	if (lid == 0)
		AtomicAdd(output, res);
}
#endif

//Second Reduction step to obtain sum of x velocity and density
//This kernel is only used when FD temp solver is not being used
__kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void LB_reduce_Ro_J_2(__global double *input,
__global double *output,
int savestep)
{
	int gid = get_global_id(0);

	double temp = input[gid];
	for (unsigned int s = 1; s < NUMBLOCKS; s++)
	{
		temp += input[s * 2 + gid];
	}

	output[savestep * 2 + gid] = temp;

}

//Second Reduction step to obtain sum of x velocity and density
//when using FD solver. Also calculates bulk mean temp at inlet and outlet
__kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void FD_reduce_Ro_J_T(__global double *input,
__global double *output,
__global double *T,
__global double2 *J,
int savestep)
{
	int gid = get_global_id(0);

	if (gid < 2)
	{
		double temp = input[gid];
		for (unsigned int s = 1; s < NUMBLOCKS; s++)
		{
			temp += input[s * 2 + gid];
		}

		output[savestep * 5 + gid] = temp;
	}
	else if (gid == 2)
	{
		double temp1 = 0., temp2 = 0.;
		for (unsigned int s = 0; s < FULLSIZEY; s++)
		{
			int ind = s + X_START_VAL_INT * FULLSIZEY_UT;
			temp1 += fabs(J[ind].x*T[ind]);
			temp2 += fabs(J[ind].x);
		}
		output[savestep * 5 + gid] = temp1 / temp2;
	}
	else if (gid == 3)
	{
		double temp1 = 0., temp2 = 0.;
		for (unsigned int s = 0; s < FULLSIZEY; s++)
		{
			int ind = s + X_MAX_VAL_INT * FULLSIZEY_UT;
			temp1 += fabs(J[ind].x*T[ind]);
			temp2 += fabs(J[ind].x);
		}
		output[savestep * 5 + gid] = temp1 / temp2;
	}
	else
	{
		double temp1 = 0., temp2 = 0.;
		for (unsigned int s = 0; s < FULLSIZEY; s++)
		{
			int ind = s + (FULLSIZEX - 1) * FULLSIZEY_UT;
			temp1 += fabs(J[ind].x*T[ind]);
			temp2 += fabs(J[ind].x);
		}
		output[savestep * 5 + gid] = temp1 / temp2;
	}
}

//Final kernel called in recursive reduce method for calculating
//mass correction factor
__kernel void LB_recursive_reduce_Ro(__global double *input,
	__global double *output,
	__local double *sdata,
	__global double *Mass0)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = input[stride] + input[stride + 1];

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[0] = (Mass0[0] - sdata[0]) / Mass0[0] / 81.;
		//output[0] = sdata[0];
	}

}

//Calculates Nusselt along the boundaries
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
void FD_calc_Nu(__global bLinks *BL,
__global int *BLind,
__global double *Nu,
__global nodeV *NV,
__global double2 *J,
__global double *T,
int SaveLoc,
int NUMBLINKS)
{
	int ix = get_global_id(0);

	if (ix >= NUMBLINKS)
		return;

	int shiftval = SaveLoc * NUMBLINKS * 3;
	
	int ii = BLind[ix];
	
	bLinks BLtemp = BL[ii];
	double2 BLpos = BLtemp.vP0 + BLtemp.vTvec * BLtemp.blLen/2.;

	int2 Posi = convert_int2(BLpos);
	int loc = Posi.x*FULLSIZEY_TR + Posi.y;

	nodeV NVtemp = NV[loc];
	double4 TT = (NVtemp.Temps - TMIN) / TDIFF;
	
	double2 dX1 = BLpos - trunc(BLpos);
	double2 dX0 = 1. - dX1;
	
	double T0y = TT.x * dX0.y + TT.z * dX1.y;
	double T1y = TT.y * dX0.y + TT.w * dX1.y;

	double T0x = TT.x * dX0.x + TT.y * dX1.x;
	double T1x = TT.z * dX0.x + TT.w * dX1.x;

	double2 dTdx = (double2)(T1y - T0y, T1x - T0x);  ///-dTdx
	
	double dTdn = dot(dTdx, BLtemp.vNvec); 

	double Um = 0.;
	double Tm = 0.;
	for (int j = 0; j < FULLSIZEY; j++)
	{
		double Unorm = length(J[Posi.x*FULLSIZEY_UT + j]);
		Um += Unorm;
		Tm += Unorm * T[Posi.x*FULLSIZEY_UT + j];
	}
	Tm = Tm / Um;
	Nu[shiftval + ix * 3] = BLpos.x;
	Nu[shiftval + ix * 3 + 1] = BLpos.y;
	Nu[shiftval + ix * 3 + 2] = fabs(NU_MULTIPLIER * dTdn / Tm);
}

int bcFindIntersectionNusselt(double *dist, double2 vLd, double2 vPL, double2 vN)
{
	double den = dot(vN, vLd);
	if (den == 0.)
		return FALSE;

	*dist = dot(vN, vPL) / den;
	
	if ((*dist) >= 0. && (*dist) <= 1.)
	{
		return TRUE;
	}
	return FALSE;
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



//Sums absolute values of array
__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
void reduce_generic_abs(__global double *input,
	__global double *output,
	__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = fabs(input[stride]) + fabs(input[stride + 1]);

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid] = sdata[0];
	}

}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_SHEAR, 1, 1)))
	void TR_shear_save_out(__global int *BLind,
	__global bLinks *BL_array,
	__global double *SSout,
	int NUMBLINKS, int save_loc)
{
	uint ii = get_global_id(0);
	if (ii >= NUMBLINKS)
		return;
	int jj = BLind[ii];
	bLinks BLtemp = BL_array[jj];
	double Tautemp = BLtemp.Tau * convert_double(BLtemp.dir);
	SSout[save_loc*NUMBLINKS + ii] = Tautemp;
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

//__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
//void reduce_J_only1(__global double2 *input,
//__global double *output,
//__local double *sdata)
//{
//	// load shared mem
//	unsigned int tid = get_local_id(0);
//	unsigned int bid = get_group_id(0);
//	unsigned int gid = get_global_id(0);
//
//	unsigned int localSize = get_local_size(0);
//	unsigned int stride = gid * 2;
//
//
//	sdata[tid] = input[stride].x + input[(stride + 1)].x;
//
//	barrier(CLK_LOCAL_MEM_FENCE);
//	// do reduction in shared mem
//	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
//	{
//		if (tid < s)
//		{
//			sdata[tid] += sdata[(tid + s)];
//		}
//		barrier(CLK_LOCAL_MEM_FENCE);
//	}
//
//	// write result for this block to global mem
//	if (tid == 0)
//	{
//		output[bid] = sdata[0];
//	}
//
//}
//
//
//__kernel __attribute__((work_group_size_hint(WORKGROUPSIZE_RED, 1, 1)))
//void reduce_J_Final(__global double *input,
//__global double *output,
//__global double *Mass)
//{
//	// load shared mem
//	unsigned int tid = get_local_id(0);
//
//	double Jsum = input[tid];
//	Jsum = work_group_reduce_add(Jsum);
//	if (tid == 0)
//	{
//		output[0] += Jsum / Mass[0];
//	}
//}
