//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////     These functions are used by    /////////////////
/////////////////   multiple kernels in this file.   /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////


// Calculates advective velocity of particles
double2 get_Vadv(double2 U00, double2 U10, double2 U01, double2 U11, double2 dX1, double2 dX0)
{
	//	double2 U0y = fma(U00, dX0.y, U01 * dX1.y);
	//	double2 U1y = fma(U10, dX0.y, U11 * dX1.y);
	//
	return fma(fma(U00, dX0.y, U01 * dX1.y), dX0.x, fma(U10, dX0.y, U11 * dX1.y) * dX1.x);
}

// calculates thermophoretic velocity of particles
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


//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////   Functions for updating tracers   /////////////////
/////////////////             along walls.           /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////


// Returns corrected velocites and temperatures for 
// interpolating velocity and temperature of tracers
// near walls.
void trGetWallNodeTempVel(__global double* __restrict__ UxVals,
	__global double* __restrict__ UyVals,
	__global double* __restrict__ Tvals,
	__global double4* __restrict__ nodeCU00,
	__global double4* __restrict__ nodeCU01,
	__global double4* __restrict__ nodeCU10,
	__global double4* __restrict__ nodeCU11,
	__global double4* __restrict__ nodeCT00,
	__global double4* __restrict__ nodeCT01,
	__global double4* __restrict__ nodeCT10,
	__global double4* __restrict__ nodeCT11,
	__local double2* UV00,
	__local double2* UV10,
	__local double2* UV01,
	__local double2* UV11,
	__local double4* TNode,
	int i, int j, int lid)
{
	// x index in full domain
	int iLB = i + TR_X_IND_START;

	// linear index in full domain
	int iLBLin = GET_GLOBAL_IDX(iLB, j);

	// linear index in tr domain
	int iTRLin = GET_TR_GLOBAL_IDX(i, j);

	double4 Uvel = (double4)(UxVals[iLBLin],
		UxVals[iLBLin + 1],
		UxVals[iTRLin + CHANNEL_LENGTH_FULL],
		UxVals[iTRLin + CHANNEL_LENGTH_FULL + 1]);

	double4 Vvel = (double4)(UyVals[iLBLin],
		UyVals[iLBLin + 1],
		UyVals[iTRLin + CHANNEL_LENGTH_FULL],
		UyVals[iTRLin + CHANNEL_LENGTH_FULL + 1]);

	double4 Tval = (double4)(Tvals[iLBLin],
		Tvals[iLBLin + 1],
		Tvals[iTRLin + CHANNEL_LENGTH_FULL],
		Tvals[iTRLin + CHANNEL_LENGTH_FULL + 1]);


	UV00[lid] = (double2)(dot(CoeffU00[iTRLin], Uvel),
		dot(CoeffU00[iTRLin], Vvel));
	UV10[lid] = (double2)(dot(CoeffU10[iTRLin], Uvel),
		dot(CoeffU10[iTRLin], Vvel));
	UV01[lid] = (double2)(dot(CoeffU01[iTRLin], Uvel),
		dot(CoeffU01[iTRLin], Vvel));
	UV11[lid] = (double2)(dot(CoeffU11[iTRLin], Uvel),
		dot(CoeffU11[iTRLin], Vvel));

	TNode[lid] = (double4)(dot(CoeffT00, Tval), dot(CoeffT10, Tval),
		dot(CoeffT01, Tval), dot(CoeffT11, Tval));

	TNode[lid] *= TDIFF;

	TNode[lid] += TMIN;
}

// finds intersection of two lines (in this case, used to find 
// intersection of tracer trajectory and BL). Name comes from legacy
// code which was 3D.
bool bcFindIntersectionLinePlane(double2* vC, double* dist, double2* vL0,
	double2* vLd, double2* vP0, double2* vP1, double2* vN, double blLen)
{
	double den = dot(*vN, *vLd);
	if (den == 0.)
		return false;

	*dist = dot(*vN, (*vP0 - *vL0)) / den;
	*vC = (*vL0) + (*vLd) * (*dist);

	if ((*dist) >= 0. && (*dist) <= 1.)
	{
		double vd0 = distance(*vC, *vP0), vd1 = distance(*vC, *vP1);
		if (fabs(blLen - vd0 - vd1) < CEPS)
			return true;
	}
	return false;
}


//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////    Kernels for updating tracers    /////////////////
/////////////////             along walls.           /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

// Finds BL index that particle crosses when updating
// tracers near walls. Returns -1 if doesnt cross BL, 
// 0 if crosses when located behind release location
// and 1 if crosses when in front of release locaion
int find_intersection(__global double2* __restrict__ Cvals,
	__global ushort2* __restrict__ P01ind,
	__global double2* __restrict__ vNvals,
	__global double* __restrict__ blLen,
	double2* vCcut_use, double2* C1, double2* dC,
	int Sflag, int wallType, int* bl_use)
{
	double dist_use = 100.;
	double dist;
	double2 vCcut;
	int blSearchMin, blSearchMax;
	if (wallType == 1)
	{
		blSearchMin = max(*bl_use - BL_SEARCH_RAD, MIN_BL_BOT);
		blSearchMax = max(*bl_use + 1 + BL_SEARCH_RAD, MAX_BL_BOT);
	}
	else
	{
		blSearchMin = max(*bl_use - BL_SEARCH_RAD, MIN_BL_TOP);
		blSearchMax = max(*bl_use + 1 + BL_SEARCH_RAD, MAX_BL_TOP);
	}



	int Tflag = -1;
	for (int k = blSearchMin; k < blSearchMax; k++)
	{
		ushort2 p01 = P01ind[k];

		if (bcFindIntersectionLinePlane(&vCcut, &dist, C1,
			dC, &Cvals[p01.x], &Cvals[p01.y], &vNvals[k],
			blLen[k]) == true)
		{
			if (dist < dist_use)
			{
				dist_use = dist;
				Tflag = (1) * Sflag;
				*bl_use = k;
				*vCcut_use = vCcut;
			}
		}
	}
	return Tflag;
}

// Reflect tracer across the point vCcut on BL defined by vN
void reflect_across_bl(double2 vCcut, double2 vN, double2* C1, double2 C2_C, double2* dC)
{
	*C1 = vCcut + vN * 0.001;
	double2 R = (double2)(1.) - vN * vN * 2.;
	double xy = -2. * vN.x * vN.y;

	double2 PastBL = C2_C - (*C1);

	*dC = R * PastBL;
	*dC += xy * PastBL;
}

// Finds the BL that has possesses the local minimum shear magnitude
// for particles to deposit on when rolling after contacting surface.

//TODO: figure out if this needs to handle a tracer rolling outside of 
//		trDomain range
short find_SS_min(__global double *__restrict__ blTau, short bl)
{
	double Taucur = blTau[bl];
	short SSdir = (Taucur < 0.) ? -1 : 1;
	short bl_next = bl + SSdir;
	double Taunext = blTau[bl_next];

#ifdef MAX_NUM_BL_ROLL
	int count = 0;
	int max_roll = MAX_NUM_BL_ROLL;
	//if (BLlist[bl].vP0.x <= X_START_VAL + 3.)
	//	max_roll = 3;
#endif

	while (1)
	{
#ifdef MAX_NUM_BL_ROLL
		if ((fabs(Taucur) >= fabs(Taunext)) && (Taucur * Taunext >= 0.) && count < max_roll)
		{
			bl = bl_next;
			Taucur = Taunext;
			bl_next += SSdir;
			Taunext = blTau[bl_next];
			count++;
		}
		else
		{
			break;
		}
#else
		if ((fabs(Taucur) > fabs(Taunext)) && (Taucur * Taunext >= 0.))
		{
			bl = bl_next;
			Taucur = Taunext;
			bl_next += SSdir;
			Taunext = blTau[bl_next];
		}
		else
		{
			break;
		}
#endif
	}
	return (bl);
}


// Updates locations of tracers which are located in nodes 
// close to wall. Checks for particles which have been sorted into
// the surrounding nodes to see if they have entered the node being solved at.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_update_par_along_wall(
	// NodC arrays
	__global double4* __restrict__ nodeCU00,
	__global double4* __restrict__ nodeCU10,
	__global double4* __restrict__ nodeCU01,
	__global double4* __restrict__ nodeCU11,
	__global double4* __restrict__ nodeCT00,
	__global double4* __restrict__ nodeCT10,
	__global double4* __restrict__ nodeCT01,
	__global double4* __restrict__ nodeCT11,
	// Par arrays
	__global int* __restrict__ parLoc,
	__global double2* __restrict__ parPos,
	__global short* __restrict__ parType,
	__global int* __restrict__ parDepFlag,
	__global uint* __restrict__ parTimer,
	// BL arrays
	__global ushort2* __restrict__ P01ind,
	__global double2* __restrict__ vNvals,
	__global double* __restrict__ blLen,
	// NodI arrays
	__global short* __restrict__ wallFlag,
	__global int* __restrict__ blIndInfo,
	// Remaining variables not part of trStructs
	__global int* __restrict__ wallInds, // array containing indices of wall nodes
	__global double* __restrict__ UxVals, //macroscopic vel
	__global double* __restrict__ UyVals, //macroscopic vel
	__global double* __restrict__ TVals, //macroscopic temp
	__global double2* __restrict__ Cvals, // vls.C array
	__global int2* Ploc_array, // sort info for particles
	__global short* update_flag, // used to ensure particles only updated once
	__global int* Out_array, // index info of particles crossing bl
	__global double2* PV, // displacement info of particles crossing bl
	__global uint* index, // tracks current index of PV/Out_array to save to
	uint nWallNodes) // number of wall nodes (changes as wall fouls)
{
	int i0 = get_global_id(0);
	if (i0 >= nWallNodes)
	{
		return;
	}

	// wall locations spread across domain, so
	// node location is obtained from wallInds array
	int iLin = wallInds[i0];
	int i, j;
	decodeGlobalIdx(iLin, &i, &j);


	int lid = get_local_id(0);

	__local UV00[WORKGROUPSIZE_TR_WALL];
	__local UV10[WORKGROUPSIZE_TR_WALL];
	__local UV01[WORKGROUPSIZE_TR_WALL];
	__local UV11[WORKGROUPSIZE_TR_WALL];
	__local nodeT[WORKGROUPSIZE_TR_WALL];

	// get corrected wall velocity and temperature values
	// for this trNode
	trGetWallNodeTempVel(UxVals, UyVals, Tvals, nodeCU00,
		nodeCU01, nodeCU10, nodeCU11, nodeCT00, nodeCT01,
		nodeCT10, nodeCT11, UV00, UV10, UV01, UV11,
		nodeT, i, j, lid);


	// This seems a bit convoluted and would be easier to
	// implement with 2 for loops and a while loop, but this
	// would lead to idle cores as each loop would require
	// all cores to finish inner loops before incrementing.
	// This should allow for a all cores to complete their
	// calculations without remaining idle waiting for other
	// cores to complete a loop. THis is at the expense of extra
	// integer calculations to calculate the next index.
	int imin = min(i - 1, 0), imax = (i + 2, FULLSIZEX_TR);
	int jmin = min(j - 1, 0), jmax = (j + 2, FULLSIZEY_TR);
	int iNeigh = imin, jNeigh = jmin;
	int neighLin = GET_TR_GLOBAL_IDX(iNeigh, jNeigh);
	int2 ploc = Ploc_array[neighLin + 2];
	while (true)
	{
		// if no particle in node, ploc.x = ploc.y = -1
		// if reached end of particles in node ploc.x = ploc.y
		if (wallFlag[neighLin] == -1 || ploc.x == ploc.y)
		{
			iNeigh++;
			if (iNeigh == imax)
			{
				jNeigh++;
				if (jNeigh == jmax)
					break;
				iNeigh = imin;
			}
			neighLin = GET_TR_GLOBAL_IDX(iNeigh, jNeigh);
			ploc = Ploc_array[neighLin + 2];

			// need to make sure that new node has particles,
			// and if not continue to catch if statment above
			// and increment iNeigh
			if (ploc.x == ploc.y)
				continue;
		}
		// ploc.x needs to be updated here since there are multiple
		// continue statements, and we dont want to have to 
		// increment it right before each of these. Instead kk is 
		// used as the current index and ploc.x will be incremented now
		int kk = ploc.x;
		ploc.x++;

		if (parLoc[kk] != ii || updateFlag[kk])
			continue;

		// set flag to indicate that this particle was updated just in
		// case it end up in another node
		updateFlag[kk] = 1;

		// Check to see if parTimer has reached zero, 
		// if so, deposit particle
		parTimer[kk] -= 1;
		if (parTimer[kk] <= 0)
		{
			parDepFlag[kk] = -2;
			parLoc[kk] = -2;
			continue;
		}

		double2 pos = parPos[kk];

		double2 dX1 = pos - trunc(pos);

		double2 Uadv = get_Vadv(UV00[lid], UV10[lid], UV01[lid], UV11[lid], dX1, (1. - dX1));
		Uadv += (pos.x > START_THERMO_VEL) ? (get_Vth(nodeT[lid], paramKth[parType[kk]], dX1, (1. - dX1))) : (0.);

		double2 dC = (Uadv)* DTTR_WALL;

		double2 C1 = pos;

		pos += dC;

		int2 Posi = convert_int2(pos);

		parLoc[kk] = Posi.x + FULLSIZEX_TR * Posi.y;

		// If particle is outside of tr domain, it is considered deposited, and 
		// needs to be re-released.
		if (pos.x >= X_MAX_VAL || pos.x < X_MIN_VAL)
		{
			parLoc[kk] = -2;
			parDepFlag[kk] = (Ppos.x >= X_MAX_VAL) ? -3 : -2;
			continue;
		}
		parPos[kk] = pos;

		// Remaining code checks to see if particle crosses any Blinks,
		// and stores necessary information in arrays for subsequent
		// kernels to test for reflection/deposition

		// if particle is between X_MIN_VAL and X_START_VAL,
		// it will be reflected in TR_reflect_particles (code always
		// reflects and does not let particles deposit at this location)
		// Sflag indicates whether or not the particle is located here
		int Sflag = (pos.x <= X_START_VAL) ? (0) : (1);
		int bl = blIndInfo[neighLin];
		double2 vCcut;

		// Tests for crossing BLink
		int Tflag = find_intersection(Cvals, P01ind, vNvals, blLen,
			&vCcut, &C1, &dC, Sflag, wallFlag[neighLin], &bl);

		// if didnt cross any BLs, continue.
		if (Tflag == -1)
			continue;

		// index (vtr.reflectInds) stores number of particles that have
		// crossed the BLink before X_START_VAL (index[0]) and after 
		// (index[1])
		int ind_temp = atomic_inc(&index[Tflag]);


		// Assuming that less than 1/8 of particles are
		// crossing BL behind starting x value. The cpu code
		// will check this is not violated. 


		// offset is TRC_NUM_TRACERS/4 because two variables are
		// stored in each array.
		int oset = REFLECT_INFO_OFFSET * Tflag;

		// Store information necessary for reflecting particles
		// and or depositing on surface
		Out_array[oset + ind_temp * 2] = kk;
		Out_array[oset + ind_temp * 2 + 1] = bl;
		PV[oset + ind_temp * 2] = dC;
		PV[oset + ind_temp * 2 + 1] = vCcut;

	}
}


// Reflects particles without allowing for possible deposition along BL (used
// with particles that are located behind release location)


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL_REFLECT, 1, 1)))
void TR_reflect_particles(
	// Par arrays
	__global int* __restrict__ parLoc,
	__global double2* __restrict__ parPos,
	__global int* __restrict__ parDepFlag,
	// BL arrays
	__global double2* __restrict__ vNvals,
	// Remaining variables not part of trStructs
	__global double2* __restrict__ PV,
	__global int* __restrict__ Info,
	int nReflect)
{
	int ii = get_global_id(0);

	if (ii >= nReflect)
		return;

	// Get index of particle to reflect
	int par_ind = Info[ii * 2];

	// get current particle position
	double2 pPos = parPos[par_ind];

	// get displacement calculated previously
	double2 dC = PV[ii * 2];
	
	// previous location
	double2 C1 = pPos - dC;

	// get reflected distance.
	reflect_across_bl(PV[ii * 2 + 1], vNvals[Info[ii * 2 + 1]], &C1, pPos, &dC);
	pPos = C1 + dC;

	// make sure that particle is still in tracer domain
	if (pPos.x >= X_MAX_VAL || pPos.x < (X_MIN_VAL) ||
		pPos.y < Y_MIN_VAL || pPos.y >= Y_MAX_VAL)
	{
		parDepFlag[par_ind] = (Ppos.x >= X_MAX_VAL) ? -3 : -2;
		parLoc[par_ind] = -2;
		return;
	}


	int2 Posi = convert_int2(pPos);

	// save info
	parLoc[par_ind] = Posi.x + FULLSIZEX_TR * Posi.y;
	parDepFlag[par_ind] = -1; // just in case this value hasnt been set yet
	parPos[par_ind] = pPos; 
}

//TODO: test whether it is faster to set number of particles using
//		set kernel argument vs using buffer;
//TODO:	test if using constant global address space is faster for
//		particle Properties arrays

// Updates tracers that have crossed a wall between release location
// and stop location. If insufficient forces on particles, they will
// deposit on BL, otherwise they will reflect.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_update_par_contact_wall(
	// Par arrays
	__global double2* __restrict__ parPos,
	__global short* __restrict__ parType,
	__global int* __restrict__ parDepFlag,
	// BL arrays
	__global double2* __restrict__ vNvals,
	__global double* __restrict__ blTau,
	__global short* __restrict__ blType,
	// Remaining variables not part of trStructs
	__global double2 * PV,
	__global uint2 * RandArray,
	__global int* Info,
	__global uint * index,
	int nReflect)
{
	int ii = get_global_id(0);

	if (ii >= nReflect)
		return;

	// Get index info stored by previous kernel
	int par_ind = Info[REFLECT_INFO_OFFSET + ii * 2];
	int bl = Info[REFLECT_INFO_OFFSET + ii * 2 + 1];
	
	
	double2 pPos = parPos[par_ind];

	// check to make sure that the particle hasnt been
	// displaced outside of x or y bounds

	//TODO: Check to see if this is still necessary.
	//		This may have been implemented when dealing
	//		with bugs in original laminar implementation
	//		and maybe not be necessary.
	if (pPos.x >= X_MAX_VAL || pPos.x < (X_MIN_VAL) ||
		pPos.y < Y_MIN_VAL || pPos.y > Y_MAX_VAL)
	{
		parDepFlag[par_ind] = (Ppos.x >= X_MAX_VAL) ? -3 : -2;
		parLoc[par_ind] = -2;
		return;
	}

	// get indices for PParam arrays
	short distInd = parType[par_ind];	// surface agnostic params
	short itype = 2 * distInd + blType[bl]; // params dependent on surface

	// Current shear at BL
	double shear_surf = fabs(blTau[bl]);

	// Test shear to see if greater than critical for particle size and surface type
	int sh_flag = isgreaterequal(shear_surf, paramTauCrit[itype]);

	// Impact velocity
	double2 Vimp_vec = PV[REFLECT_INFO_OFFSET + ii * 2] / DTTR_WALL;

	// Vimp^2
	double Vimp2 = dot(Vimp_vec, Vimp_vec);

	// square of rebound velocity
	double Vreb2 = (sh_flag) ? ((shear_surf - paramTauCrit[itype]) / MU_NUMBER / 2. * paramDp[distInd]) :
		(Vimp2 - 2. * (paramQaPrime[itype] / paramMp[distInd]));

	// Calculate probabilty for deposit based on velocities, particle size and surface properties
	double Sprob = (Vreb2 > 0.) ? (paramQa[itype] / (0.5 * paramMp[distInd] * Vreb2 - paramQaPrime[itype])) : (100.);

	// Generate uniform random number and test if particle deposited
	// If shear stress is greater than critical value, or if random number
	// is greater than Sprob particle doesnt deposit and depFlagTemp = -1,
	// otherwise it deposits and depFlagTemp = bl index.
	double Rprob;
	RandArray[ii] = MWC64X_NextUint(RandArray[ii], &Rprob);
	int depFlagTemp = (Sprob < Rprob || sh_flag) ? (-1) : (bl);

	// If particle has deposited, increment element 2, else increment element 0
	int ind_info = (depFlagTemp > -1) ? 2 : 0;
	int ind_temp = atomic_inc(&index[ind_info]);

	// Save information in Info array to use in subsequent kernel calls
	// using first half of Info array to store info, since we dont want to
	// overwrite any information in second half of array.
	// First 1/4 of Info now stores particle index of deposited particles to
	// be updated by trWallParKernel[4] and 1/4-1/2 stores index of particles
	// to re-test for intersection with neighboring BL).
	// Note: only storing par_ind in array now, so no need to multiply ind_temp by 2
	Info[ind_info * REFLECT_INFO_OFFSET / 4 + ind_temp] = par_ind;


	if (depFlagTemp == -1)
	{
		// Reflect across BL with Vreb providing displacement distance
		// of rebounded particle
		Vreb2 = sqrt(Vreb2) * DTTR_WALL;
		double2 dC = PV[REFLECT_INFO_OFFSET + ii * 2];
		double2 C1 = Pcur.pos - dC;
		reflect_across_bl(PV[REFLECT_INFO_OFFSET + ii * 2 + 1], vNvals[par_ind], &C1, pPos, &dC);


		// Information still needs to be stored to test rebounded particle didnt cross
		// another boundary link
		dC = normalize(dC) * Vreb2;
		pPos = C1 + dC;

		// Dont need to store this info for deposited particles, so there
		// is no need to use the offset.
		PV[ind_temp * 2] = dC;
		PV[ind_temp * 2 + 1] = C1;

		// bl particle reflected off of is stored in parDepFlag,
		// but this will be corrected in next kernel called (it 
		// should actually be -1 if it reflected off of BL since that
		// indicates it hasnt deposited yet)
		depFlagTemp = bl; 
	}

	parPos[par_ind] = pPos;
	parDepFlag[par_ind] = depFlagTemp;

}

// Tests neighboring BL of BL which particle reflected off of to see if it 
// crossed after reflecting off previous BL.  
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_update_par_contact_wall2(
	// Par arrays
	__global int* __restrict__ parLoc,
	__global double2* __restrict__ parPos,
	__global int* __restrict__ parDepFlag,
	// BL arrays
	__global double2* __restrict__ vNvals,
	__global ushort2* __restrict__ P01ind,
	__global double* __restrict__ blLen,
	// Remaining variables not part of trStructs
	__global double2 *Cvals,
	__global double2 * PV,
	__global int* Info,
	__global uint * index,
	int num_pars)
{
	int ii = get_global_id(0);

	if (ii >= num_pars)
		return;

	// Get index of particle from last kernel
	int par_ind = Info[ii];
	int lid = get_local_id(0);

	__local double2 dC[WORKGROUPSIZE_TR_WALL];
	__local double2 C1[WORKGROUPSIZE_TR_WALL];

	// Get positions and displacements
	double2 pPos = parPos[par_ind];
	dC[lid] = PV[ii * 2];
	C1[lid] = PV[ii * 2 + 1];


	// Make sure particle didnt leave domain
	if (pPos.x >= X_MAX_VAL || pPos.x < (X_MIN_VAL) ||
		pPos.y < Y_MIN_VAL || pPos.y > Y_MAX_VAL)
	{
		parDepFlag[par_ind] = (Ppos.x >= X_MAX_VAL) ? -3 : -2;
		parLoc[par_ind] = -2;
		return;
	}

	double2 vCcut;
	double dist;

	int bl = parDepFlag[par_ind];

	// get Bl neighboring BL that particle reflected off of
	int blneigh = (dC[lid].x >= 0) ? (1) : (-1);
	blneigh += bl;

	ushort2 p01 = P01ind[bl];

	// test to see if particle crossed neighboring bl
	if (bcFindIntersectionLinePlane(&vCcut, &dist, C1[lid], 
		dC[lid], &Cvals[p01.x], &Cvals[p01.y], &vNvals[bl],	blLen[bl])
	{
		// if crossed neighboring BL, save necessary info for an additional iteration
		int ind_temp = atomic_inc(&index[1]);

		// TR_update_par_contact_wall will read from 2nd half of these arrays
		Info[REFLECT_INFO_OFFSET + ind_temp * 2] = par_ind;
		Info[REFLECT_INFO_OFFSET + ind_temp * 2 + 1] = blneigh;
		PV[REFLECT_INFO_OFFSET + ind_temp * 2] = dC[lid];
		PV[REFLECT_INFO_OFFSET + ind_temp * 2 + 1] = vCcut;
	}

	int2 Posi = convert_int2(pPos);
	parLoc[par_ind] = Posi.x + FULLSIZEX_TR * Posi.y;
	parDepFlag[par_ind] = -1;
	parPos[par_ind] = pPos;
}


// Finds BL with local minimum in shear stress near where particles deposited on wall
// and sets appropriate values of par struct.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void TR_deposit_particles_on_wall(
	// Par Arrays
	__global int* __restrict__ parLoc,
	__global double2* __restrict__ parPos,
	__global int* __restrict__ parDepFlag,
	__global ushort* __restrict__ parDepTimer,
	// BL Arrays
	__global ushort2* __restrict__ P01ind,
	__global double *__restrict__ blTau,
	__global double2* __restrict__ Cvals,
	__global int* __restrict__ Info,
	int num_pars)
{
	int ii = get_global_id(0);

	if (ii >= num_pars)
		return;

	// information was stored in 1/4-1/2 of par_ind range,
	// so we need to offset it when reading
	int par_ind = Info[REFLECT_INFO_OFFSET/4 + ii];

	int lid = get_local_id(0);
		
	int pDepFlag = parDepFlag[par_ind];

	pDepFlag = find_SS_min(blTau, pDepFlag);
	
	ushort2 p01 = P01ind[pDepFlag];
	double2 C1val = Cvals[p01.x];
	double2 vTvec = Cvals[p01.y] - C1val;
	
	parPos[par_ind] = C1val + normalize(vTvec) * length(vTvec) / 2.;
	parLoc[par_ind] = -1;
	parDepTimer[par_ind] = DEP_TIMER_START;
	parDepFlag[par_ind] = pDepFlag;
}







//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////   Functions for updating tracers   /////////////////
/////////////////          away from walls.          /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////



// Returns velocity and temperature for 
// all nodes for interpolating velocity and temperature of tracers
// within center of domain (away from walls).
void trGetNodeTempVel(__global double* __restrict__ UxVals,
	__global double* __restrict__ UyVals,
	__global double* __restrict__ Tvals,
	double2* UV00, double2* UV10, double2* UV01,
	double2* UV11, double4* TNode,
	int i, int j)
{
	// x index in full domain
	int iLB = i + TR_X_IND_START;

	// linear index in full domain
	int iLBLin = GET_GLOBAL_IDX(iLB, j);

	*UV00 = (double2)(UxVals[iLBLin], UyVals[iLBLin]);
	*UV10 = (double2)(UxVals[iLBLin + 1], UyVals[iLBLin + 1]);
	*UV01 = (double2)(UxVals[iLBLin + CHANNEL_LENGTH_FULL],
					  UyVals[iLBLin + CHANNEL_LENGTH_FULL]);
	*UV11 = (double2)(UxVals[iLBLin + CHANNEL_LENGTH_FULL + 1],
					  UyVals[iLBLin + CHANNEL_LENGTH_FULL + 1]);

	*TNode = (double4)(Tvals[iLBLin],
		Tvals[iLBLin + 1],
		Tvals[iTRLin + CHANNEL_LENGTH_FULL],
		Tvals[iTRLin + CHANNEL_LENGTH_FULL + 1]);

	*TNode *= TDIFF;

	*TNode += TMIN;
}



//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////    Kernel for updating tracers     /////////////////
/////////////////          away from walls.          /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////


//Kernel to update positions of particles located away from walls

// TODO: test whether its faster to use array of nodes to update, or use entire 
//		trDomain and return on solid and near-wall nodes
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR, 1, 1)))
void TR_update_par_no_wall(
	// Par arrays
	__global int* __restrict__ parLoc,
	__global double2* __restrict__ parPos,
	__global short* __restrict__ parType,
	__global int* __restrict__ parDepFlag,
	__global uint* __restrict__ parTimer,
	// NodI arrays
	__global short* __restrict__ wallFlag,
	// Remaining variables not part of trStructs
	__global double* __restrict__ UxVals, //macroscopic vel
	__global double* __restrict__ UyVals, //macroscopic vel
	__global double* __restrict__ TVals, //macroscopic temp
	__global double2* __restrict__ Cvals, // vls.C array
	__global int2 *__restrict__ Ploc_array,
	__global short *__restrict__ update_flag,
	__global int* __restrict__ activeNodes,
	uint nActiveNodes)
{
	int i = get_global_id(0);
	if (i >= nActiveNodes)
	{
		return;
	}

	// get index of nodes aways from wall
	int iLin = activeNodes[i0];

	//activeNodes contains wall nodes as well, which are handled by wall 
	//kernels, so we can skip those nodes
	if (wallFlag[neighLin] != 0)
		return;

	int i, j;
	decodeGlobalIdx(iLin, &i, &j);


	int lid = get_local_id(0);

	__local UV00[WORKGROUPSIZE_TR_WALL];
	__local UV10[WORKGROUPSIZE_TR_WALL];
	__local UV01[WORKGROUPSIZE_TR_WALL];
	__local UV11[WORKGROUPSIZE_TR_WALL];
	__local nodeT[WORKGROUPSIZE_TR_WALL];

	// get wall velocities and temperatures for interpolation
	trGetNodeTempVel(UxVals, UyVals, Tvals, uble2 * UV00,
		UV10, UV01, UV11, TNode, i, j);

	// This seems a bit convoluted and would be easier to
	// implement with 2 for loops and a while loop, but this
	// would lead to idle cores as each loop would require
	// all cores to finish inner loops before incrementing.
	// This should allow for a all cores to complete their
	// calculations without remaining idle waiting for other
	// cores to complete a loop. THis is at the expense of extra
	// integer calculations to calculate the next index.
	int imin = min(i - 1, 0), imax = (i + 2, FULLSIZEX_TR);
	int jmin = min(j - 1, 0), jmax = (j + 2, FULLSIZEY_TR);
	int iNeigh = imin, jNeigh = jmin;
	int neighLin = GET_TR_GLOBAL_IDX(iNeigh, jNeigh);
	int2 ploc = Ploc_array[neighLin + 2];

	while (true)
	{
		// if no particle in node, ploc.x = ploc.y = -1
		// if reached end of particles in node ploc.x = ploc.y
		if (wallFlag[neighLin] == -1 || ploc.x == ploc.y)
		{
			iNeigh++;
			if (iNeigh == imax)
			{
				jNeigh++;
				if (jNeigh == jmax)
					break;
				iNeigh = imin;
			}
			neighLin = GET_TR_GLOBAL_IDX(iNeigh, jNeigh);
			ploc = Ploc_array[neighLin + 2];

			// need to make sure that new node has particles,
			// and if not continue to catch if statment above
			// and increment iNeigh
			if (ploc.x == ploc.y)
				continue;
		}
		// ploc.x needs to be updated here since there are multiple
		// continue statements, and we dont want to have to 
		// increment it right before each of these. Instead kk is 
		// used as the current index and ploc.x will be incremented now
		int kk = ploc.x;
		ploc.x++;

		if (parLoc[kk] != ii || updateFlag[kk])
			continue;

		// set flag to indicate that this particle was updated just in
		// case it end up in another node
		updateFlag[kk] = 1;

		// Check to see if parTimer has reached zero, 
		// if so, deposit particle
		parTimer[kk] -= TIMER_DECREMENT;
		if (parTimer[kk] <= 0)
		{
			parDepFlag[kk] = -2;
			parLoc[kk] = -2;
			continue;
		}

		double2 pos = parPos[kk];

		double2 dX1 = pos - trunc(pos);
		double2 dX0 = 1. - dX1;
		double2 Uadv = get_Vadv(UV00[lid], UV10[lid], UV01[lid], UV11[lid], dX1, dX0);
		Uadv += (pos.x > START_THERMO_VEL) ? (get_Vth(nodeT[lid], paramKth[parType[kk]], dX1, dX0)) : (0.);

		double2 dC = (Uadv)* DTTR;

		double2 C1 = pos;

		pos += dC;
		parPos[kk] = pos;

		if (pos.x >= X_MAX_VAL || pos.x < (X_MIN_VAL))
		{
			parDepFlag[kk] = (pos.x >= X_MAX_VAL) ? -3 : -2;
			parLoc[kk] = -2;
			continue;
		}
		int2 Posi = convert_int2(Pcur.pos);
		parLoc[kk] = Posi.x + FULLSIZEX_TR * Posi.y;
	}
}





