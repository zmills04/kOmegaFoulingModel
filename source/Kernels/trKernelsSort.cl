// Kernels associated with particle sorting, which include merge
// sorting functions, particle release kernel, and kernels 
// associated with clumping method.
// First set of functions taken from AMD's Bolt library (altered
// for use with particle data). Remaining kernels implemented
// by me.
//
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
/////////////////                                  ///////////////////
///////////////// Merge sort methods adapted from  ///////////////////
/////////////////     AMD's bolt library           ///////////////////
/////////////////                                  ///////////////////
//////////////////////////////////////////////////////////////////////
int lowerBoundBinarylocal(local int* __restrict__ data_loc, int left, 
	int right, int searchVal_loc)
{
	int firstIndex = left;
	int lastIndex = right;

	while (firstIndex < lastIndex)
	{
		int midIndex = (firstIndex + lastIndex) / 2;

		if (data_loc[midIndex] < searchVal_loc)
		{
			firstIndex = midIndex + 1;
		}
		else
		{
			lastIndex = midIndex;
		}
	}
	return firstIndex;
}

int upperBoundBinarylocal(local int* __restrict__ data_loc, int left, 
	int right, int searchVal_loc)
{
	int upperBound = lowerBoundBinarylocal(data_loc, left, right, searchVal_loc);

	if (upperBound != right)
	{
		int mid = 0;
		int upperValue_loc = data_loc[upperBound];
		while ((upperValue_loc == searchVal_loc) && (upperBound < right))
		{
			mid = (upperBound + right) / 2;
			if (data_loc[mid] == searchVal_loc)
			{
				upperBound = mid + 1;
			}
			else
			{
				right = mid;
				upperBound++;
			}
			upperValue = data_loc[upperBound];
		}
	}
	return upperBound;
}

int lowerBoundLinear(global int* __restrict__ data_loc, int left, int right, int searchVal_loc)
{
	int firstIndex = left;
	int lastIndex = right;

	while (firstIndex < lastIndex)
	{
		if (data_loc[firstIndex] < searchVal_loc)
		{
			firstIndex = firstIndex + 1;
		}
		else
		{
			break;
		}
	}

	return firstIndex;
}

int lowerBoundBinary(__global int* __restrict__ source_ptr_loc, int left, int right, int searchVal_loc)
{
	int firstIndex = left;
	int lastIndex = right;

	while (firstIndex < lastIndex)
	{
		int midIndex = (firstIndex + lastIndex) / 2;

		if (source_ptr_loc[midIndex] < searchVal_loc)
		{
			firstIndex = midIndex + 1;
		}
		else
		{
			lastIndex = midIndex;
		}
	}
	return firstIndex;
}

int upperBoundBinary(__global int* __restrict__ source_ptr_loc, int left, int right, int searchVal_loc)
{
	int upperBound = lowerBoundBinary(source_ptr_loc, left, right, searchVal_loc);

	if (upperBound != right)
	{
		int mid = 0;
		int upperValue = source_ptr_loc[upperBound];
		while ((searchVal_loc == upperValue_loc) && (upperBound < right))
		{
			mid = (upperBound + right) / 2;

			if (source_ptr_loc[mid] == searchVal_loc)
			{
				upperBound = mid + 1;
			}
			else
			{
				right = mid;
				upperBound++;
			}
			upperValue = source_ptr_loc[upperBound];
		}
	}
	return upperBound;
}


// This version sorts loc array only and tracks original location
// Sorting of remaining arrays performed in Sort_update_loc
kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
void Sort_merge_global(__global int* __restrict__ source_loc,
	__global int* __restrict__ source_ptr_orig, 
	__global int* __restrict__ result_loc,
	__global int* __restrict__ result_ptr_orig,
	const int srcLogicalBlockSize)
{
	size_t globalID = get_global_id(0);// * get_global_size(0) + get_global_id(0);
	size_t groupID = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
	size_t localID = get_local_id(0);// +get_local_size(0) + get_local_id(0);
	size_t wgSize = get_local_size(0);// * get_local_size(1);


	if (globalID >= TRC_NUM_TRACERS)
		return;
	int srcBlockNum = globalID / srcLogicalBlockSize;
	int srcBlockIndex = globalID % srcLogicalBlockSize;


	int dstLogicalBlockSize = srcLogicalBlockSize << 1;
	int leftBlockIndex = globalID & ~(dstLogicalBlockSize - 1);

	leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
	leftBlockIndex = (((leftBlockIndex) < (TRC_NUM_TRACERS)) ?
		(leftBlockIndex) : (TRC_NUM_TRACERS));
	int rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (TRC_NUM_TRACERS)) ?
		(leftBlockIndex + srcLogicalBlockSize) : (TRC_NUM_TRACERS));

	int insertionIndex = 0;

	int search_val = source_loc[globalID];

	if ((srcBlockNum & 0x1) == 0)
	{
		insertionIndex = lowerBoundBinary(source_loc, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
	}
	else
	{
		insertionIndex = upperBoundBinary(source_loc, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
	}

	int dstIndex = srcBlockIndex + insertionIndex;
	dstIndex += (srcBlockNum / 2) * dstLogicalBlockSize;

	result_loc[dstIndex] = source_loc[globalID];
	result_loc_orig[dstIndex] = source_loc_orig[globalID];
}

kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
void Sort_merge_local(__global int* __restrict__ source_loc,
	__global int* __restrict__ source_loc_orig,
	__local int* __restrict__ lds_loc,
	__local int* __restrict__ lds_loc_orig,
	__local int* __restrict__ lds2_loc,
	__local int* __restrict__ lds2_loc_orig)
{
	size_t gloId = get_global_id(0);// *get_global_size(0) + get_global_id(0);
	size_t groId = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
	size_t locId = get_local_id(0);// +get_local_size(0) + get_local_id(0);
	size_t wgSize = get_local_size(0);// *get_local_size(1);


	if (gloId < TRC_NUM_TRACERS)
	{
		lds_loc[locId] = source_loc[gloId];
		lds_loc_orig[locId] = gloId;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	int end = wgSize;
	if ((groId + 1) * (wgSize) >= TRC_NUM_TRACERS)
		end = TRC_NUM_TRACERS - (groId * wgSize);

	int numMerges = SORT_NUM_MERGES;
	int pass;
	for (pass = 1; pass <= numMerges; ++pass)
	{
		int srcLogicalBlockSize = 1 << (pass - 1);
		if (gloId < TRC_NUM_TRACERS)
		{
			int srcBlockNum = (locId) / srcLogicalBlockSize;
			int srcBlockIndex = (locId) % srcLogicalBlockSize;

			int dstLogicalBlockSize = srcLogicalBlockSize << 1;
			int leftBlockIndex = (locId) & ~(dstLogicalBlockSize - 1);

			leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
			leftBlockIndex = (((leftBlockIndex) < (end)) ? (leftBlockIndex) : (end));
			int rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (end)) ? (leftBlockIndex + srcLogicalBlockSize) : (end));

			int insertionIndex = 0;
			if (pass % 2 != 0)
			{
				if ((srcBlockNum & 0x1) == 0)
				{
					insertionIndex = lowerBoundBinarylocal(lds_loc, leftBlockIndex, rightBlockIndex, lds_loc[locId]) - leftBlockIndex;
				}
				else
				{
					insertionIndex = upperBoundBinarylocal(lds_loc, leftBlockIndex, rightBlockIndex, lds_loc[locId]) - leftBlockIndex;
				}
			}
			else
			{
				if ((srcBlockNum & 0x1) == 0)
				{
					insertionIndex = lowerBoundBinarylocal(lds2_loc, leftBlockIndex, rightBlockIndex, lds2_loc[locId]) - leftBlockIndex;
				}
				else
				{
					insertionIndex = upperBoundBinarylocal(lds2_loc, leftBlockIndex, rightBlockIndex, lds2_loc[locId]) - leftBlockIndex;
				}
			}
			int dstIndex = srcBlockIndex + insertionIndex;
			dstIndex += (srcBlockNum / 2) * dstLogicalBlockSize;
			if (pass % 2 != 0)
			{
				lds2_loc[dstIndex] = lds_loc[locId];
				lds2_loc_orig[dstIndex] = lds_loc_orig[locId];
			}
			else
			{
				lds_loc[dstIndex] = lds2_loc[locId];
				lds_loc_orig[dstIndex] = lds2_loc_orig[locId];
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (gloId < TRC_NUM_TRACERS)
	{
		source_loc[gloId] = lds_loc[locId];
		source_loc_orig[gloId] = lds_loc_orig[locId];
	}
}


//////////////////////////////////////////////////////////////////////
/////////////////                                  ///////////////////
///////////////// Function sorting arrays based on ///////////////////
/////////////////     merge results and setting    ///////////////////
/////////////////         Ploc_array values        ///////////////////
/////////////////                                  ///////////////////
//////////////////////////////////////////////////////////////////////


// TODO: see if ploc_array can be replaced by int elements and use next index 
// as stop point (currently uses int2 with start and stop index).

// Creates Ploc array, which specifies the starting index in particle arrays
// for a given cell index.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_PLOC, 1, 1)))
void Sort_update_loc(global double2* __restrict__ source_pos, 
	global uint* __restrict__ source_Num_rep,
	global short* __restrict__ source_type,
	global short* __restrict__ source_Dep_Flag,
	global ushort* __restrict__ source_Dep_timer,
	global ushort* __restrict__ source_timer,
	global double2* __restrict__ result_pos,
	global uint* __restrict__ result_Num_rep,
	global short* __restrict__ result_type,
	global short* __restrict__ result_Dep_Flag,
	global ushort* __restrict__ result_Dep_timer,
	global ushort* __restrict__ result_timer,
	__global int* __restrict__ P_loc,
#ifdef ODD_NUM_MERGES
	__global int* __restrict__ result_loc,
#endif
	__global int* __restrict__ reOrderInfo,
	__global int* Ploc_array)
{
	int i = get_global_id(0);
	if (i >= TRC_NUM_TRACERS)
		return;

	// Reorder particle arrays
	int oldLoc = reOrderInfo[i];
	result_pos[i] = source_pos[oldLoc]; 
	result_Num_rep[i] = source_Num_rep[oldLoc];
	result_type[i] = source_type[oldLoc];
	result_Dep_Flag[i] = source_Dep_Flag[oldLoc];
	result_Dep_timer[i] = source_Dep_timer[oldLoc];
	result_timer[i] = source_timer[oldLoc];

	// Fill Ploc_array
	int ploc = P_loc[i];

#ifdef ODD_NUM_MERGES
	result_loc[i] = ploc;
#endif

	if (i == TRC_NUM_TRACERS - 1)
	{
		Ploc_array[2 * (ploc + 2) + 1] = i + 1;
		return;
	}

	int next_loc = P_loc[i + 1];

	if (i == 0)
	{
		Ploc_array[2 * (ploc + 2)] = 0;
	}

	if (ploc != next_loc)
	{
		Ploc_array[2 * (ploc + 2) + 1] = i + 1;
		Ploc_array[2 * (next_loc + 2)] = i + 1;
		for (int k = ploc + 1; k < next_loc; k++)
		{
			Ploc_array[2 * (k + 2)] = -1;
		}
	}
}

//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
///////////////// Functions for releasing particles  /////////////////
/////////////////       after fully depositing       /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

int Get_Par_Type(double randval, __global double* __restrict__ parP_D_dist)
{
	int j = 0;

	while (j < NUM_PAR_SIZES - 1)
	{
		if (randval >= parP_D_dist[j] && randval < parP_D_dist[j+1])
			return j;
		j++;
	}
	return NUM_PAR_SIZES - 1;
}


// Returns velocity value at y = Xloc. Used to randomly generate particles
// in distribution corresponding to velocity distribution.
double find_lineval(double Xloc, __global double* trP, __global double* U_array)
{
	if (Xloc < trP[TRP_BOT_LOC_IND])
	{
		return Xloc * U_array[UVALS_START_INDEX] / trP[TRP_BOT_LOC_IND];
	}
	else if (Xloc < trP[TRP_TOP_LOC_IND])
	{
		double Xp = Xloc - trP[TRP_BOT_LOC_IND];
		double X1 = floor(Xp);
		int Xint = convert_int(X1) + UVALS_START_INDEX;
		double dX = Xp - X1;
		return U_array[Xint] * (1. - dX) + dX * U_array[Xint + 1];
	}

	double Xp = Xloc - trP[TRP_BOT_LOC_IND];
	double X1 = floor(Xp);
	int Xint = convert_int(X1) + UVALS_START_INDEX;
	double dX = trP[TRP_TOP_LOC_IND] - X1;
	double dX1 = trP[TRP_BVAL_IND] - Xp;
	return dX1 * U_array[Xint] / dX;
}



// Called after sort.
// Rereleases particles according to velocity distribution.
// Done by randomly generating location for particle in box w/ Height
// equal to height of channel and width equal to maximum velocity.
// Particle is kept if its location falls within velocity distribution plotted
// within box and re-generated if it falls outside.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RERELEASE, 1, 1)))
void TR_release_par(global double2* __restrict__ source_pos, 
	__global uint* __restrict__ source_Num_rep,
	__global short* __restrict__ source_type,
	__global short* __restrict__ source_Dep_Flag,
	__global ushort* __restrict__ source_timer,
	__global int* __restrict__ source_loc,
	__global uint2* RandArray,
	__global double* __restrict__ trP,
	__global double* U_array,
	__global double* __restrict__ parP_D_dist,
	uint maxel,
	uint Conc_number)
{
	int i = get_global_id(0);

	if (i >= maxel)
		return;

	if (source_timer[i] > 0)
	{
		source_timer[i] -= 1;
		return;
	}


	uint2 RandEl = RandArray[i];
	double Umv = trP[TRP_UMAX_IND];
	double Bval = trP[TRP_BVAL_IND];
	double2 randval;
	RandEl = MWC64X_NextUint2(RandEl, &randval);
	double Xt = randval.x * Bval;
	double Yt = randval.y * Umv;

	double Fxx = find_lineval(Xt, trP, U_array);

	while (Yt > Fxx)
	{
		RandEl = MWC64X_NextUint2(RandEl, &randval);
		Xt = randval.x * Bval;
		Yt = randval.y * Umv;
		Fxx = find_lineval(Xt, trP, U_array);
	}

	double randval_single;
	RandEl = MWC64X_NextUint(RandEl, &randval_single);
	RandArray[i] = RandEl;

	source_type[i] = Get_Par_Type(randval_single, parP_D_dist);
	source_timer[i] = PAR_TIMER_START;
	source_pos[i].y = Xt + trP[TRP_OFFSET_Y_IND];
	source_pos[i].x = TRP_X_RELEASE;
	source_Dep_Flag[i] = -1;
	int2 Posi = convert_int2(source_pos[i]);
	source_loc[i] = Posi.x + FULLSIZEX_TR * Posi.y;
	source_Num_rep[i] = Conc_number;

}

//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////  Functions for clumping particles  /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

// Not very efficient, because it is only called a few hundred times
// per simulation. Function finds BL above and below particle and determines if it
// is located between them (and therefore is in domain).
__kernel void Test_Bounds_Clumped_Particles(__global double2* __restrict__ P_pos,
	__global int* __restrict__ P_Dep_Flag,
	__global int* __restrict__ P_loc,
	__global double* __restrict__ trP,
	__global double2* __restrict__ lsc,
	int num_pars,
	int start_ind)
{
	int i = get_global_id(0) + start_ind;

	if (i >= num_pars)
		return;

	double2 Parpos = P_pos[i];

	int depFl = P_Dep_Flag[i];

	// return if particle is deposited
	if (depFl == -2 || depFl > -1)
		return;

	// Mainly to return particles that have escaped domain back to beginning
	if (Parpos.x <= X_MIN_VAL || Parpos.x > X_MAX_VAL || isnan(Parpos.x))
	{
		P_Dep_Flag[i] = -2;
		P_loc[i] = -2;
		P_pos[i] = (double2)(X_START_VAL, trP[TRP_OFFSET_Y_IND] + trP[TRP_BVAL_IND] / 2.);
		return;
	}

	int lsc_ind_temp = LSC_REL_BOT;

	while (1)
	{
		if (lsc[lsc_ind_temp].x < Parpos.x && lsc[lsc_ind_temp+1].x >= Parpos.x)
			break;
		lsc_ind_temp++;
	}

	int lsc_bot_ind = lsc_ind_temp;

	lsc_ind_temp += (LSC_NN / 2 - 20);

	while (1)
	{
		if (lsc[lsc_ind_temp+1].x < Parpos.x && lsc[lsc_ind_temp].x >= Parpos.x)
			break;

		lsc_ind_temp++;
	}

	int lsc_top_ind = lsc_ind_temp;

	double2 V0 = lsc[lsc_bot_ind], V1 = lsc[lsc_bot_ind+1];

	double ybot_pos = V0.y + (Parpos.x - V0.x) * (V1.y - V0.y) / (V1.x - V0.x);

	V0 = lsc[lsc_top_ind+1], V1 = lsc[lsc_top_ind];

	double ytop_pos = V0.y + (Parpos.x - V0.x) * (V1.y - V0.y) / (V1.x - V0.x);

	// if the particle is outside of either of these boundary layers,
	// it is re-located by reflecting across BL in y direction only.
	if (Parpos.y <= ybot_pos || Parpos.y >= ytop_pos)
	{
		// Reflect particle across BL
		Parpos.y += (Parpos.y <= ybot_pos) ? 
			(2. * (ybot_pos - Parpos.y)) : (2. * (ytop_pos - Parpos.y));
		P_pos[i] = Parpos;

		// Update loc array with new position as it may have changed
		int2 Posi = convert_int2(Parpos);
		P_loc[i] = Posi.x + FULLSIZEX_TR * Posi.y;
	}
}

//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////  Find Umax Kernel                  /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

#ifdef OPENCL_VERSION_1_2
__kernel void find_umax(__global double* input,
	__global double* trP,
	__local double* sdata)
{
	const int tid = get_global_id(0); // = local_id(0)
	const int localSize = get_local_size(0);
	const int stride = gid * 2;

	double input1 = (stride < FULLSIZEY) ? input[UVALS_START_INDEX + stride] : 0.;
	double input2 = (stride+1 < FULLSIZEY) ? input[UVALS_START_INDEX + stride+1] : 0.;
	
	sdata[tid] = max(input1, input2);

	barrier(CLK_LOCAL_MEM_FENCE);
	for (int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = max(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0) 
	{
		trP[TRP_UMAX_VAL_IND] = sdata[0];
	}
}
#else

__kernel void find_umax(__global double* U_array,
	__global double* trP)
{
	int i = get_global_id(0);

	double Uval = (i < FULLSIZEY) ? U_array[UVALS_START_INDEX + i] : 0.;
	
	barrier(CLK_LOCAL_MEM_FENCE);

	trP[TRP_UMAX_VAL_IND] = work_group_reduce_max(Uval);
}

#endif





//// These macros are used to avoid clutter in the functions.
//
//#define GLOBAL_SOURCE	global double2* __restrict__ source_pos,\
//						global uint* __restrict__ source_Num_rep,\
//						global short* __restrict__ source_type,\
//						global short* __restrict__ source_Dep_Flag,\
//						global ushort* __restrict__ source_Dep_timer,\
//						global ushort* __restrict__ source_timer,\
//						global int* __restrict__ source_loc
//
//#define GLOBAL_RESULT	global double2* __restrict__ result_pos,\
//						global uint* __restrict__ result_Num_rep,\
//						global short* __restrict__ result_type,\
//						global short* __restrict__ result_Dep_Flag,\
//						global ushort* __restrict__ result_Dep_timer,\
//						global ushort* __restrict__ result_timer,\
//						global int* __restrict__ result_loc
//
//
//// Probably an elegant way to only create a single full macro which has
//// variable prefixes passed to it, but not worth effort at this point
//// For these result_*[result_ind] = source_*[source_ind]
//#define SET_SOURCE_TO_RESULT(result_ind, source_ind) \
//						result_pos[result_ind] = source_pos[source_ind];\
//						result_Num_rep[result_ind] = source_Num_rep[source_ind];\
//						result_type[result_ind] = source_type[source_ind];\
//						result_Dep_Flag[result_ind] = source_Dep_Flag[source_ind];\
//						result_Dep_timer[result_ind] = source_Dep_timer[source_ind];\
//						result_timer[result_ind] = source_timer[source_ind];\
//						result_loc[result_ind] = source_loc[source_ind];
//
//#define SET_LDS_TO_SOURCE(result_ind, source_ind) \
//						source_pos[result_ind] = lds_pos[source_ind];\
//						source_Num_rep[result_ind] = lds_Num_rep[source_ind];\
//						source_type[result_ind] = lds_type[source_ind];\
//						source_Dep_Flag[result_ind] = lds_Dep_Flag[source_ind];\
//						source_Dep_timer[result_ind] = lds_Dep_timer[source_ind];\
//						source_timer[result_ind] = lds_timer[source_ind];\
//						source_loc[result_ind] = lds_loc[source_ind];
//
//#define SET_SOURCE_TO_LDS(result_ind, source_ind) \
//						lds_pos[result_ind] = source_pos[source_ind];\
//						lds_Num_rep[result_ind] = source_Num_rep[source_ind];\
//						lds_type[result_ind] = source_type[source_ind];\
//						lds_Dep_Flag[result_ind] = source_Dep_Flag[source_ind];\
//						lds_Dep_timer[result_ind] = source_Dep_timer[source_ind];\
//						lds_timer[result_ind] = source_timer[source_ind];\
//						lds_loc[result_ind] = source_loc[source_ind];
//
//#define SET_LDS_TO_LDS2(result_ind, source_ind) \
//						lds2_pos[result_ind] = lds_pos[source_ind];\
//						lds2_Num_rep[result_ind] = lds_Num_rep[source_ind];\
//						lds2_type[result_ind] = lds_type[source_ind];\
//						lds2_Dep_Flag[result_ind] = lds_Dep_Flag[source_ind];\
//						lds2_Dep_timer[result_ind] = lds_Dep_timer[source_ind];\
//						lds2_timer[result_ind] = lds_timer[source_ind];\
//						lds2_loc[result_ind] = lds_loc[source_ind];
//
//#define SET_LDS2_TO_LDS(result_ind, source_ind) \
//						lds_pos[result_ind] = lds2_pos[source_ind];\
//						lds_Num_rep[result_ind] = lds2_Num_rep[source_ind];\
//						lds_type[result_ind] = lds2_type[source_ind];\
//						lds_Dep_Flag[result_ind] = lds2_Dep_Flag[source_ind];\
//						lds_Dep_timer[result_ind] = lds2_Dep_timer[source_ind];\
//						lds_timer[result_ind] = lds2_timer[source_ind];\
//						lds_loc[result_ind] = lds2_loc[source_ind];


// Implementation sorting entire particle data arrays (not debugged yet)
//kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
//void Sort_merge_global(GLOBAL_SOURCE,
//	GLOBAL_RESULT,
//	const int srcLogicalBlockSize)
//{
//	size_t globalID = get_global_id(0);// * get_global_size(0) + get_global_id(0);
//	size_t groupID = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
//	size_t localID = get_local_id(0);// +get_local_size(0) + get_local_id(0);
//	size_t wgSize = get_local_size(0);// * get_local_size(1);
//
//
//	if (globalID >= TRC_NUM_TRACERS)
//		return;
//	int srcBlockNum = globalID / srcLogicalBlockSize;
//	int srcBlockIndex = globalID % srcLogicalBlockSize;
//
//
//	int dstLogicalBlockSize = srcLogicalBlockSize << 1;
//	int leftBlockIndex = globalID & ~(dstLogicalBlockSize - 1);
//
//	leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
//	leftBlockIndex = (((leftBlockIndex) < (TRC_NUM_TRACERS)) ?
//		(leftBlockIndex) : (TRC_NUM_TRACERS));
//	int rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (TRC_NUM_TRACERS)) ?
//		(leftBlockIndex + srcLogicalBlockSize) : (TRC_NUM_TRACERS));
//
//	int insertionIndex = 0;
//	
//	int search_val = source_loc[globalID];
//	
//	if ((srcBlockNum & 0x1) == 0)
//	{
//		insertionIndex = lowerBoundBinary(source_loc, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
//	}
//	else
//	{
//		insertionIndex = upperBoundBinary(source_loc, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
//	}
//
//	int dstIndex = srcBlockIndex + insertionIndex;
//	dstIndex += (srcBlockNum / 2) * dstLogicalBlockSize;
//
//	SET_SOURCE_TO_RESULT(dstIndex, globalID);
//}
//
//kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
//void Sort_merge_local(GLOBAL_SOURCE,
//	local double2* __restrict__ lds_pos,
//	local uint* __restrict__ lds_Num_rep,
//	local short* __restrict__ lds_type,
//	local short* __restrict__ lds_Dep_Flag,
//	local ushort* __restrict__ lds_Dep_timer,
//	local ushort* __restrict__ lds_timer,
//	local int* __restrict__ lds_loc,
//	local double2* __restrict__ lds2_pos,
//	local uint* __restrict__ lds2_Num_rep,
//	local short* __restrict__ lds2_type,
//	local short* __restrict__ lds2_Dep_Flag,
//	local ushort* __restrict__ lds2_Dep_timer,
//	local ushort* __restrict__ lds2_timer,
//	local int* __restrict__ lds2_loc)
//{
//	size_t gloId = get_global_id(0);// *get_global_size(0) + get_global_id(0);
//	size_t groId = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
//	size_t locId = get_local_id(0);// +get_local_size(0) + get_local_id(0);
//	size_t wgSize = get_local_size(0);// *get_local_size(1);
//
//
//	if (gloId < TRC_NUM_TRACERS)
//	{
//		SET_SOURCE_TO_LDS(locID, gloID);
//	}
//	barrier(CLK_LOCAL_MEM_FENCE);
//	int end = wgSize;
//	if ((groId + 1) * (wgSize) >= TRC_NUM_TRACERS)
//		end = TRC_NUM_TRACERS - (groId * wgSize);
//
//	int numMerges = SORT_NUM_MERGES;
//	int pass;
//	for (pass = 1; pass <= numMerges; ++pass)
//	{
//		int srcLogicalBlockSize = 1 << (pass - 1);
//		if (gloId < TRC_NUM_TRACERS)
//		{
//			int srcBlockNum = (locId) / srcLogicalBlockSize;
//			int srcBlockIndex = (locId) % srcLogicalBlockSize;
//
//			int dstLogicalBlockSize = srcLogicalBlockSize << 1;
//			int leftBlockIndex = (locId) & ~(dstLogicalBlockSize - 1);
//
//			leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
//			leftBlockIndex = (((leftBlockIndex) < (end)) ? (leftBlockIndex) : (end));
//			int rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (end)) ? (leftBlockIndex + srcLogicalBlockSize) : (end));
//
//			int insertionIndex = 0;
//			if (pass % 2 != 0)
//			{
//				if ((srcBlockNum & 0x1) == 0)
//				{
//					insertionIndex = lowerBoundBinarylocal(lds_loc, leftBlockIndex, rightBlockIndex, lds_loc[locId]) - leftBlockIndex;
//				}
//				else
//				{
//					insertionIndex = upperBoundBinarylocal(lds_loc, leftBlockIndex, rightBlockIndex, lds_loc[locId]) - leftBlockIndex;
//				}
//			}
//			else
//			{
//				if ((srcBlockNum & 0x1) == 0)
//				{
//					insertionIndex = lowerBoundBinarylocal(lds2_loc, leftBlockIndex, rightBlockIndex, lds2_loc[locId]) - leftBlockIndex;
//				}
//				else
//				{
//					insertionIndex = upperBoundBinarylocal(lds2_loc, leftBlockIndex, rightBlockIndex, lds2_loc[locId]) - leftBlockIndex;
//				}
//			}
//			int dstIndex = srcBlockIndex + insertionIndex;
//			dstIndex += (srcBlockNum / 2) * dstLogicalBlockSize;
//			if (pass % 2 != 0)
//			{
//				SET_LDS_TO_LDS2(dstIndex, locId);
//			}
//			else
//			{
//				SET_LDS2_TO_LDS(dstIndex, locId);
//			}
//		}
//		barrier(CLK_LOCAL_MEM_FENCE);
//	}
//	if (gloId < TRC_NUM_TRACERS)
//	{
//		SET_LDS_TO_SOURCE(gloId, locId);
//	}
//}
//
//
//// Creates Ploc array, which specifies the starting index in particle arrays
//// for a given cell index.
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_PLOC, 1, 1)))
//void Sort_update_loc(__global int* __restrict__ P_loc,
//	__global int* Ploc_array)
//{
//	int i = get_global_id(0);
//	if (i >= TRC_NUM_TRACERS)
//		return;
//	int ploc = P_loc[i];
//
//	if (i == TRC_NUM_TRACERS - 1)
//	{
//		Ploc_array[2 * (ploc + 2) + 1] = i + 1;
//		return;
//	}
//
//	int next_loc = P_loc[i + 1];
//
//	if (i == 0)
//	{
//		Ploc_array[2 * (ploc + 2)] = 0;
//	}
//
//	if (ploc != next_loc)
//	{
//		Ploc_array[2 * (ploc + 2) + 1] = i + 1;
//		Ploc_array[2 * (next_loc + 2)] = i + 1;
//		for (int k = ploc + 1; k < next_loc; k++)
//		{
//			Ploc_array[2 * (k + 2)] = -1;
//		}
//	}
//}


// Working Version used to sort Par struct array (for reference when debugging new
// implementation)

//uint lowerBoundBinarylocal(local par* data, uint left, uint right, par searchVal)
//{
//	uint firstIndex = left;
//	uint lastIndex = right;
//
//	while (firstIndex < lastIndex)
//	{
//		uint midIndex = (firstIndex + lastIndex) / 2;
//		par midValue = data[midIndex];
//
//		if (midValue.loc < searchVal.loc)
//		{
//			firstIndex = midIndex + 1;
//		}
//		else
//		{
//			lastIndex = midIndex;
//		}
//	}
//	return firstIndex;
//}
//
//uint upperBoundBinarylocal(local par* data, uint left, uint right, par searchVal)
//{
//	uint upperBound = lowerBoundBinarylocal(data, left, right, searchVal);
//
//	if (upperBound != right)
//	{
//		uint mid = 0;
//		par upperValue = data[upperBound];
//		while ((upperValue.loc == searchVal.loc) && (upperBound < right))
//		{
//			mid = (upperBound + right) / 2;
//			par midValue = data[mid];
//			if (midValue.loc == searchVal.loc)
//			{
//				upperBound = mid + 1;
//			}
//			else
//			{
//				right = mid;
//				upperBound++;
//			}
//			upperValue = data[upperBound];
//		}
//	}
//	return upperBound;
//}
//
//uint lowerBoundLinear(global par* data, uint left, uint right, par searchVal)
//{
//	uint firstIndex = left;
//	uint lastIndex = right;
//
//	while (firstIndex < lastIndex)
//	{
//		par dataVal = data[firstIndex];
//
//		if (dataVal.loc < searchVal.loc)
//		{
//			firstIndex = firstIndex + 1;
//		}
//		else
//		{
//			break;
//		}
//	}
//
//	return firstIndex;
//}
//
//uint lowerBoundBinary(__global par *source_ptr, uint left, uint right, par searchVal)
//{
//	uint firstIndex = left;
//	uint lastIndex = right;
//
//	while (firstIndex < lastIndex)
//	{
//		uint midIndex = (firstIndex + lastIndex) / 2;
//		par midValue = source_ptr[midIndex];
//
//		if (midValue.loc < searchVal.loc)
//		{
//			firstIndex = midIndex + 1;
//		}
//		else
//		{
//			lastIndex = midIndex;
//		}
//	}
//	return firstIndex;
//}
//
//uint upperBoundBinary(__global par *source_ptr, uint left, uint right, par searchVal)
//{
//	uint upperBound = lowerBoundBinary(source_ptr, left, right, searchVal);
//
//	if (upperBound != right)
//	{
//		uint mid = 0;
//		par upperValue = source_ptr[upperBound];
//		while ((searchVal.loc == upperValue.loc) && (upperBound < right))
//		{
//			mid = (upperBound + right) / 2;
//			par midValue = source_ptr[mid];
//			if (midValue.loc == searchVal.loc)
//			{
//				upperBound = mid + 1;
//			}
//			else
//			{
//				right = mid;
//				upperBound++;
//			}
//			upperValue = source_ptr[upperBound];
//		}
//	}
//	return upperBound;
//}
//
//
//kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
//	void Sort_merge_global(global par* source_ptr,
//	global par* result_ptr,
//	const uint srcVecSize,
//	const uint srcLogicalBlockSize)
//{
//	size_t globalID = get_global_id(0);// * get_global_size(0) + get_global_id(0);
//	size_t groupID = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
//	size_t localID = get_local_id(0);// +get_local_size(0) + get_local_id(0);
//	size_t wgSize = get_local_size(0);// * get_local_size(1);
//
//
//	if (globalID >= srcVecSize)
//		return;
//	uint srcBlockNum = globalID / srcLogicalBlockSize;
//	uint srcBlockIndex = globalID % srcLogicalBlockSize;
//
//
//	uint dstLogicalBlockSize = srcLogicalBlockSize << 1;
//	uint leftBlockIndex = globalID & ~(dstLogicalBlockSize - 1);
//
//	leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
//	leftBlockIndex = (((leftBlockIndex) < (srcVecSize)) ? (leftBlockIndex) : (srcVecSize));
//	uint rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (srcVecSize)) ? (leftBlockIndex + srcLogicalBlockSize) : (srcVecSize));
//
//	uint insertionIndex = 0;
//	par search_val = source_ptr[globalID];
//	if ((srcBlockNum & 0x1) == 0)
//	{
//		insertionIndex = lowerBoundBinary(source_ptr, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
//	}
//	else
//	{
//		insertionIndex = upperBoundBinary(source_ptr, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
//	}
//
//	uint dstBlockIndex = srcBlockIndex + insertionIndex;
//	uint dstBlockNum = srcBlockNum / 2;
//
//	result_ptr[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = source_ptr[globalID];
//
//
//}
//
//kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
//	void Sort_merge_local(global par* data_ptr,
//	const uint vecSize,
//	local par* lds,
//	local par* lds2)
//{
//	size_t gloId = get_global_id(0);// *get_global_size(0) + get_global_id(0);
//	size_t groId = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
//	size_t locId = get_local_id(0);// +get_local_size(0) + get_local_id(0);
//	size_t wgSize = get_local_size(0);// *get_local_size(1);
//
//
//	if (gloId < vecSize)
//	{
//		lds[locId] = data_ptr[gloId];
//	}
//	barrier(CLK_LOCAL_MEM_FENCE);
//	uint end = wgSize;
//	if ((groId + 1)*(wgSize) >= vecSize)
//		end = vecSize - (groId*wgSize);
//
//	uint numMerges = SORT_NUM_MERGES;
//	uint pass;
//	for (pass = 1; pass <= numMerges; ++pass)
//	{
//		uint srcLogicalBlockSize = 1 << (pass - 1);
//		if (gloId < vecSize)
//		{
//			uint srcBlockNum = (locId) / srcLogicalBlockSize;
//			uint srcBlockIndex = (locId) % srcLogicalBlockSize;
//
//			uint dstLogicalBlockSize = srcLogicalBlockSize << 1;
//			uint leftBlockIndex = (locId)& ~(dstLogicalBlockSize - 1);
//
//			leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
//			leftBlockIndex = (((leftBlockIndex) < (end)) ? (leftBlockIndex) : (end));
//			uint rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (end)) ? (leftBlockIndex + srcLogicalBlockSize) : (end));
//
//			uint insertionIndex = 0;
//			if (pass % 2 != 0)
//			{
//				if ((srcBlockNum & 0x1) == 0)
//				{
//					insertionIndex = lowerBoundBinarylocal(lds, leftBlockIndex, rightBlockIndex, lds[locId]) - leftBlockIndex;
//				}
//				else
//				{
//					insertionIndex = upperBoundBinarylocal(lds, leftBlockIndex, rightBlockIndex, lds[locId]) - leftBlockIndex;
//				}
//			}
//			else
//			{
//				if ((srcBlockNum & 0x1) == 0)
//				{
//					insertionIndex = lowerBoundBinarylocal(lds2, leftBlockIndex, rightBlockIndex, lds2[locId]) - leftBlockIndex;
//				}
//				else
//				{
//					insertionIndex = upperBoundBinarylocal(lds2, leftBlockIndex, rightBlockIndex, lds2[locId]) - leftBlockIndex;
//				}
//			}
//			uint dstBlockIndex = srcBlockIndex + insertionIndex;
//			uint dstBlockNum = srcBlockNum / 2;
//			if (pass % 2 != 0)
//				lds2[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = lds[locId];
//			else
//				lds[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = lds2[locId];
//		}
//		barrier(CLK_LOCAL_MEM_FENCE);
//	}
//	if (gloId < vecSize)
//	{
//		data_ptr[gloId] = lds[locId];
//	}
//}
//
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_PLOC, 1, 1)))
//	void Sort_update_loc(__global par *P,
//	__global int *Ploc_array)
//{
//	int i = get_global_id(0);
//	if (i >= TRC_NUM_TRACERS)
//		return;
//	int ploc = P[i].loc;
//	
//	if (i == TRC_NUM_TRACERS - 1)
//	{
//		Ploc_array[2*(ploc + 2)+1] = i+1;
//		return;
//	}
//	
//	int next_loc = P[i + 1].loc;
//
//	if(i == 0)
//	{
//		Ploc_array[2*(ploc + 2)] = 0;
//	}
//	
//	if (ploc != next_loc)
//	{
//		Ploc_array[2*(ploc + 2)+1] = i + 1;
//		Ploc_array[2*(next_loc + 2)] = i + 1;
//		for (int k = ploc + 1; k < next_loc; k++)
//		{
//			Ploc_array[2*(k + 2)] = -1;
//		}
//	}
//}
//













