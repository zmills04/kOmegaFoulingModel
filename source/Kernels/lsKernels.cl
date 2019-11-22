// This cl file contains kernels for updating LS arrays during the fouling
// layer update step (along with vtr.BL arrays). 

// LS arrays and where they are updated:
// C - updated in kernel found in flKernels.cl
// C0 - not updated
// nType - Nodes which fall in the updated fouling layer are 
//		   updated in updateM kernel implemented in this file.
//			NESW_BOUNDARY_NODE information is updated along with
//			vlb.ibb arrays in kernel in lbKernels.cl
// nTypePrev - updated each time from nType after updateM is called and
//				and before updateBoundaryNodes is called
// dXArr0 - not updated
// dXArr - updated in updateBLdXArr 
// ibbArr - dynamic array, so updated on host in updateBoundaryArrays
// ssArr - dynamic array, so updated on host in updateBoundaryArrays
// ssArrIndMap - updated along with ssArr  
// Masses - updated by reduce kernels
// BL - not updated
// bFlag - not updated
// ssArrInds - updated with ssArr
// lsMap - not updated

// calculates cross product of two vectors
double Cross2(double2 v1, double2 v2)
{//absolute value of cross product
	return fabs(v1.x * v2.y - v1.y * v2.x);
}

// tests if point L0 falls with triangle defined by P0,P1,P2
// by comparing area of P0,P1,P2 triangle with sum of areas
// of three triangles formed by P0,P1,P2 and L0
int test_inside(double2 P0, double2 P1, double2 P2, double2 L0)
{
	double2 v10 = P1 - P0, v20 = P2 - P0, v21 = P2 - P1, vd0 = L0 - P0, vd1 = L0 - P1;
	double stot = Cross2(v10, v20), s0 = Cross2(v10, vd0), s1 = Cross2(vd0, v20), s2 = Cross2(v21, vd1);

	if (fabs(stot - s0 - s1 - s2) < CEPS)
		return 1;
	return 0;
}


// Determines if a node falls within two triangles formed by C0[i], C0[i+1]
// C[i] and C[i+1], and if so, sets the node as a fouling layer node.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEM, 1, 1)))
void updateNType(__global NTYPE_TYPE * __restrict__ Map,
	__global double2 * __restrict__ C,
	__global double2 * __restrict__ C0,
	__global double* __restrict__ Ux_array,
	__global double* __restrict__ Uy_array,
	__global double* __restrict__ Ro_array,
	__global double* __restrict__ kAmat,
	__global double* __restrict__ kbmat,
	__global double* __restrict__ kxmat,
	__global double* __restrict__ oAmat,
	__global double* __restrict__ obmat,
	__global double* __restrict__ oxmat,
	__global double* __restrict__ nutArray,
	__global int* __restrict__ IndArr)
{
	int i = get_global_id(0);

	if (i >= LSC_UPDATE_M_NN || i == LSC_UPDATE_M_SKIP)
		return;

	double2 P0 = C0[i], P1 = C0[i + 1], P2 = C[i], P3 = C[i + 1];

	bool neq02 = EQUALS2D(P0, P2);
	bool neq13 = EQUALS2D(P1, P3);

	// if neq02 and neq13 == true, no fouling layer,
	// so skip
	if (neq02 && neq13)
		return;

	int2 Cmax = max4(P0, P1, P2, P3), Cmin = min4(P0, P1, P2, P3);

	if (Cmax.y >= CHANNEL_HEIGHT)
		Cmax.y = CHANNEL_HEIGHT - 1;
	if (Cmin.y >= CHANNEL_HEIGHT)
		return;

	if (Cmax.y < 0)
		return;
	if (Cmin.y < 0)
		Cmin.y = 0;

	if (Cmax.x >= CHANNEL_LENGTH)
		Cmax.x = CHANNEL_LENGTH - 1;
	if (Cmin.x >= CHANNEL_LENGTH)
		return;

	if (Cmax.x < 0)
		return;
	if (Cmin.x < 0)
		Cmin.x = 0;


	// Iterates through all surrounding nodes
	for (int ii = Cmin.x; ii <= Cmax.x; ii++)
	{
		for (int jj = Cmin.y; jj <= Cmax.y; jj++)
		{
			int gid = ii + CHANNEL_LENGTH_FULL * jj;
			NTYPE_TYPE type = Map[gid];

			// If already a solid node (either wall or fouling layer),
			// nothing needs to be changed
			if (type & M_SOLID_NODE)
				continue;

			// If P0 == P2, no fouling layer, or if test_inside(P0,P1,P2) == 0
			// node does not fall within layer, so skip
			// same for P1 == P3 and test_inside(P3,P2,P1)
			if ((!neq13 && test_inside(P3, P2, P1, (double2)(ii, jj)) ) ||
				(!neq02 && test_inside(P0, P1, P2, (double2)(ii, jj)) ))
			{
				// Change from fluid node to solid node (will still have M0_FLUID_NODE
				// flag, so it can be tested to see if it is a fouled node
				type ^= SWITCH_FLUID_TO_SOLID;
				//type |= M_SOLID_NODE;


				// Set A matrix C location to 1. and rest of terms
				// along with element in b vector to 0. so it always
				// solves as 0.
				int KOind = IndArr[gid + DIST_SIZE * 4];
				kAmat[KOind] = 0.;
				oAmat[KOind] = 0.;

				KOind = IndArr[gid + DIST_SIZE * 2];
				kAmat[KOind] = 0.;
				oAmat[KOind] = 0.;

				KOind = IndArr[gid];
				kAmat[KOind] = 1.;
				oAmat[KOind] = 1.;

				KOind = IndArr[gid + DIST_SIZE];
				kAmat[KOind] = 0.;
				oAmat[KOind] = 0.;

				KOind = IndArr[gid + DIST_SIZE * 3];
				kAmat[KOind] = 0.;
				oAmat[KOind] = 0.;

				kbmat[gid] = 0.;
				obmat[gid] = 0.;

				// Setting macroscopic arrays to
				// appropriate values
				kxmat[gid] = 0.;
				oxmat[gid] = 0.;
				nutArray[gid] = 0.;
				Ux_array[gid] = 0.;
				Uy_array[gid] = 0.;
				Ro_array[gid] = 1.;
			}
			Map[gid] = type;
		}
	}
}




// Finds the intersection between a node and a line through two points defined by a BL
// fluidBN is boundary node on fluid side and solidBN is inside solid. vCcut is actual location
// of the intersection, vL0 is location of LB node being tested, vLd is unit vector pointing
// toward BL, vC0 and vC1 are wall nodes
int bcFindIntersectionLS(double* dist, int *dir, double2 vL0, double2 *vLd,
	double2 vC0, double2 vC1, double vP10Length, double2 vN)
{
	// Calculate distance between vL0 and intersection
	double den = dot(vN, *vLd);
	
	// if perpendicular, will not cross so return 0 (false)
	if (den == 0.)
		return 0;

	// vN is pointing away from BL, vLd is pointing away from node,
	// so if den < 0, vLd is pointing toward BL and we should continue
	// to test that distribution direction
	// if den > 0, vLd is pointing away from BL, so we should test opposite
	// distribution direction.
	if (den > 0.)
	{
		*dir = RevDir[*dir];
		*vLd *= -1.0;
		den *= -1.0;
	}

	*dist = dot(vN, (vC0 - vL0)) / den;

	// Looking for node closest to wall, so dist > 1
	// will have a node closer to wall in that direction.
	// Testing to make sure dist > 0. just in case
	if ((*dist) < 0. || (*dist) > 1.)
		return 0;

	// Location of intersection between two lines
	double2 vCcut = vL0 + (*vLd) * (*dist);

	// Test to see if intersection falls between vC0 and vC1
	// if falls outside, return 0.
	double2 vd0 = vCcut - vC0, vd1 = vCcut - vC1;
	if (fabs(vP10Length - length(vd0) - length(vd1)) >= CEPS)
		return 0;

	// Limiting the minimum value of dX to avoid instabilities that
	// it might lead to.
	//(*dist) = max((*dist), CEPS);

	return 1;
}

//// This kernel also updates vtr.BL arrays even though it is grouped with LS kernels
// Updates vtr.BL arrays, vls.nType and ibb arrays by iterating through BL
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEM, 1, 1)))
void updateBoundaryNodes(__global ushort2* __restrict__ BL_P01ind,
	__global double2* __restrict__ BL_vNvec,
	__global double* __restrict__ BL_blLen,
	__global int* __restrict__ BL_nodeLoc,
	__global short* __restrict__ BL_intType,
	__global double2* C,
	__global double2* Cnot,
	__global NTYPE_TYPE* Map,
	__global double* dXArr
#ifndef IN_KERNEL_IBB
	, __global int* __restrict__ ibbArr,
	__global double* __restrict__ ibbDistArr,
	__global int* ibbIndCount,
	int ibbArrSize
#endif
)
{
	int i = get_global_id(0);
	if (i >= FULL_BL_ARRAY_SIZE)
		return;

	ushort2 p01Ind = BL_P01ind[i];

	double2 vC0 = C[p01Ind.x], vC1 = C[p01Ind.y];
	double2 vC0not = Cnot[p01Ind.x], vC1not = Cnot[p01Ind.y];
	
	double2 vT = vC1 - vC0;
	double2 vTnot = vC1not - vC0not;


	double2 CenterDist = (vC0 + vT / 2.) - (vC0not + vTnot / 2.);
	BL_intType[i] = (length(CenterDist) >= FOUL_SIZE_SWITCH_SOOT2) ? ((ushort)1) : ((ushort)0);

	double blLen = length(vT);
	BL_blLen[i] = blLen;
	vT /= blLen;
	double2 vCn = (double2)(-vT.y, vT.x);
	BL_vNvec[i] = vCn;

	double2 centpos = vC0 + vT * blLen / 2.;
	int2 blind = convert_int2(centpos);

	// Not sure about this. Cannot figure out it this is something new I was
	// trying to implemement, or if it was from legacy code. Keeping here, but
	// commented just in case it is needed and just has a bug in implementation.
	//BL_blInd[i] = (blind.x - TR_X_IND_START) + (FULLSIZEX_TR_PADDED)* blind.y;

	BL_nodeLoc[i] = (blind.x >= TR_X_IND_START && blind.x < TR_X_IND_STOP) ?
		((blind.x - TR_X_IND_START) + (FULLSIZEX_TR_PADDED)* blind.y) : -1;

	int2 vCmin = min2(vC0, vC1), vCmax = max2(vC0, vC1);
	//vCmin -= 1;
	//vCmax += 1;
	if (vCmin.x < 0) vCmin.x = 0;
	if (vCmin.x >= CHANNEL_LENGTH) return;
	if (vCmax.x < 0) return;
	if (vCmax.x >= CHANNEL_LENGTH) vCmax.x = CHANNEL_LENGTH-1;

	if (vCmin.y < 0) vCmin.y = 0;
	if (vCmin.y >= CHANNEL_HEIGHT) return;
	if (vCmax.y < 0) return;
	if (vCmax.y >= CHANNEL_HEIGHT) vCmax.y = CHANNEL_HEIGHT-1;


	for (int ii = vCmin.x; ii <= vCmax.x; ii++)
	{
		for (int jj = vCmin.y; jj <= vCmax.y; jj++)
		{
			int iiF_1D = ii + CHANNEL_LENGTH_FULL * jj;
			NTYPE_TYPE ntype = Map[iiF_1D];

			if (ntype & M_SOLID_NODE)
				continue;

			double dist;

			int dir = 0;
			double2 Cdir = (double2)(1., 0.);
			double2 vL0 = convert_double2((int2)(ii, jj));
			
			if (bcFindIntersectionLS(&dist, &dir, vL0, &Cdir,
				vC0, vC1, blLen, vCn))
			{
				// With simple straight inlet and outlet that doesmt 
				// have any fouling, will never have E or W direction
				// of first or last node cross BL, therefore no
				// need to worry about indexing out of bounds when
				// getting neighboring index
				int iiS_1D = iiF_1D + CDirX[dir];

				//////////// Set required info for IBB.//////////
#ifndef IN_KERNEL_IBB
				// get index of next element to add to array
				int ibbind = atomic_add(&ibbIndCount[0], 1);
				ibbind = min(ibbind, ibbArrSize-1);

				ibbDistArr[ibbind] = dist;
				ibbArr[ibbind] = iiF_1D + DIST_SIZE * (dir + 1);
#endif
				//////////// Set required info for dXArr //////////

				// x         x
				//
				//		/\
				// x   /  \  x
				//    /    \
				// If the wall forms a shape similar to above, where 
				// a BL falls in between two nodes, this will be accounted
				// for in the ibb method, but ignored in the finite difference
				// methods.
				NTYPE_TYPE ntype2 = Map[iiS_1D];
				if (ntype2 & M_SOLID_NODE)
				{
					// Set corresponding element in dXArr
					dXArr[iiF_1D + dir * DIST_SIZE] = max(dist, MINIMUM_DX_DIST);

					// if node is M0_FLUID_NODE, it is now a fouling
					// layer node and the distance needs to be tracked
					if (ntype2 & M0_FLUID_NODE)
						dXArr[iiS_1D + RevDir[dir] * DIST_SIZE] = max(1. - dist, MINIMUM_DX_DIST);
#ifdef GPU_DEBUG
					else
						dXArr[iiS_1D + RevDir[dir] * DIST_SIZE] = 0.0;
#endif
				}
			}

			// Below is the same as above but for the 3 other sets of
			// directions
			dir = 2;
			Cdir = (double2)(0., 1.);

			if (bcFindIntersectionLS(&dist, &dir, vL0, &Cdir,
				vC0, vC1, blLen, vCn))
			{
				// no need to worry about y being periodic
				int iiS_1D = ii + CHANNEL_LENGTH_FULL * (jj + CDirY[dir]);

				//////////// Set required info for IBB.//////////
#ifndef IN_KERNEL_IBB
				// get index of next element to add to array
				int ibbind = atomic_add(&ibbIndCount[0], 1);
				ibbind = min(ibbind, ibbArrSize-1);

				ibbDistArr[ibbind] = dist;
				ibbArr[ibbind] = iiF_1D + DIST_SIZE * (dir + 1);
#endif
				//////////// Set required info for dXArr //////////
				NTYPE_TYPE ntype2 = Map[iiS_1D];
				if (ntype2 & M_SOLID_NODE)
				{
					// Set corresponding element in dXArr
					dXArr[iiF_1D + dir * DIST_SIZE] = max(dist, MINIMUM_DX_DIST);

					// if node is M0_FLUID_NODE, it is now a fouling
					// layer node and the distance needs to be tracked
					if (ntype2 & M0_FLUID_NODE)
						dXArr[iiS_1D + RevDir[dir]*DIST_SIZE] = max(1. - dist, MINIMUM_DX_DIST);
#ifdef GPU_DEBUG
					else
						dXArr[iiS_1D + RevDir[dir] * DIST_SIZE] = 0.0;
#endif
				}
			}
					   			 
#ifndef IN_KERNEL_IBB
			dir = 4;
			Cdir = (double2)(1., 1.);

			if (bcFindIntersectionLS(&dist, &dir, vL0, &Cdir,
				vC0, vC1, blLen, vCn))
			{
				//////////// Set required info for IBB.//////////

				// get index of next element to add to array
				int ibbind = atomic_add(&ibbIndCount[0], 1);
				ibbind = min(ibbind, ibbArrSize-1);

				ibbDistArr[ibbind] = dist;
				ibbArr[ibbind] = iiF_1D + DIST_SIZE * (dir + 1);

			}

			dir = 6;
			Cdir = (double2)(1., -1.);

			if (bcFindIntersectionLS(&dist, &dir, vL0, &Cdir,
				vC0, vC1, blLen, vCn))
			{

				//////////// Set required info for IBB.//////////

				// get index of next element to add to array
				int ibbind = atomic_add(&ibbIndCount[0], 1);
				ibbind = min(ibbind, ibbArrSize-1);

				ibbDistArr[ibbind] = dist;
				ibbArr[ibbind] = iiF_1D + DIST_SIZE * (dir + 1);
			}
#endif
		}
	}
}


// if the ibb arrays need to be reallocated, they will also need to be updated again
// so instead of calling updateBoundaryNodes again which does other things,
// this can be called which will only update ibb arrays
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEFD, 1, 1)))
void updateIBBOnly(__global double2* C,
	__global NTYPE_TYPE* Map,
	__global int* __restrict__ ibbArr,
	__global double* __restrict__ ibbDistArr,
	__global int* ibbIndCount,
	int ibbArrSize)
{
	int i = get_global_id(0);

	if (i >= LSC_UPDATE_M_NN || i == LSC_UPDATE_M_SKIP)
		return;

	double2 vC0 = C[i], vC1 = C[i + 1];

	double2 vT = vC1 - vC0;

	double blLen = length(vT);
	vT /= blLen;
	double2 vCn = (double2)(-vT.y, vT.x);


	int2 vCmin = min2(vC0, vC1), vCmax = max2(vC0, vC1);
	//vCmin -= 1;
	//vCmax += 1;

	// Note, simple bounce back is implemented at x = 0 and x = CHANNEL_LENGTH-1
	// so these locations are skipped in the search. (simple bounce back is automatically
	// implemented in collision kernel, so this is what LB code defaults to when the ibb
	// information arrays do not include a distribution direction)

	if (vCmin.x < 1) vCmin.x = 1;
	if (vCmin.x >= (CHANNEL_LENGTH-1)) return;
	if (vCmax.x < 0) return;
	if (vCmax.x >= (CHANNEL_LENGTH-1)) vCmax.x = CHANNEL_LENGTH-2;

	if (vCmin.y < 0) vCmin.y = 0;
	if (vCmin.y >= CHANNEL_HEIGHT) return;
	if (vCmax.y < 0) return;
	if (vCmax.y >= CHANNEL_HEIGHT) vCmax.y = CHANNEL_HEIGHT - 1;


	for (int ii = vCmin.x; ii <= vCmax.x; ii++)
	{
		for (int jj = vCmin.y; jj <= vCmax.y; jj++)
		{
			int ii0_1D = ii + CHANNEL_LENGTH_FULL * jj;
			NTYPE_TYPE ntype = Map[ii0_1D];

			if (ntype & M_SOLID_NODE)
				continue;

			double2 vL0 = convert_double2((int2)(ii, jj));

			for(int dir = 0; dir < 8; dir+=2)
			{

				double2 Cdir = CXY_DOUBLE[dir];
				double dist;
				
				if (bcFindIntersectionLS(&dist, &dir, vL0, &Cdir,
					vC0, vC1, blLen, vCn))
				{
					//////////// Set required info for IBB.//////////

					// get index of next element to add to array
					int ibbind = atomic_add(&ibbIndCount[0], 1);
				
					ibbind = min(ibbind, ibbArrSize-1);

					ibbDistArr[ibbind] = dist;
					ibbArr[ibbind] = ii0_1D + DIST_SIZE * (dir + 1);
				}
			}
		}
	}
}





__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEM, 1, 1)))
void updateNodeBoundaryInfo(__global NTYPE_TYPE* Map)
{
	int i = get_global_id(0);
	if (i >= CHANNEL_LENGTH)
		return;

	int j = get_global_id(1);

	int gid = GET_GLOBAL_IDX(i, j);

	NTYPE_TYPE ntype = Map[gid] & RESET_NTYPE_NODE;
	NTYPE_TYPE ntypeMod = 0x0000;
	NTYPE_TYPE nCompare = ((ntype & 0x0003) ^ 0x0003);


	// test each direction
	// east
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j * CHANNEL_LENGTH_FULL);
	if ((Map[gide] & nCompare))
		ntypeMod |= (E_BOUND);

	// West
	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j * CHANNEL_LENGTH_FULL);
	if ((Map[gidw] & nCompare))
		ntypeMod |= (W_BOUND);

	// North
	if (j < (CHANNEL_HEIGHT - 1))
	{
		// North
		if ((Map[gid + CHANNEL_LENGTH_FULL] & nCompare))
			ntypeMod |= (N_BOUND);
		// NE
		if ((Map[gide + CHANNEL_LENGTH_FULL] & nCompare))
			ntypeMod |= (NE_BOUND);
		// NW		
		if ((Map[gidw + CHANNEL_LENGTH_FULL] & nCompare))
			ntypeMod |= (NW_BOUND);
	}

	// South
	if (j > 0)
	{
		// South
		if ((Map[gid - CHANNEL_LENGTH_FULL] & nCompare))
			ntypeMod |= (S_BOUND);
		// SE
		if ((Map[gide - CHANNEL_LENGTH_FULL] & nCompare))
			ntypeMod |= (SE_BOUND);
		// SW		
		if ((Map[gidw - CHANNEL_LENGTH_FULL] & nCompare))
			ntypeMod |= (SW_BOUND);
	}

	// set all bits which are re-calculated here to 0
	// i.e. set SOLID_BOUNDARY_NODE, and all direction bits to 0
	

	// if node is a solid, and any neighs in E,W,N,S direction are not solid,
	// then this node is a solid boundary node.
	// !ntypeMod will set any boundary direction bit to 1 if node
	// in that direction is fluid, & 0x0F00 will have a non-zero bit only
	// if N,E,S,W has a non-zero, so will only be true if neighbor is a fluid.
	if (ntype & M_SOLID_NODE)
	{
		ntypeMod = ((ntypeMod & 0x0F00)) ? SOLID_BOUNDARY_NODE : 0x0000;
	}

	Map[gid] = ntype | ntypeMod;
}




//// Finds the intersection between a node and a line through two points defined by a BL
//// fluidBN is boundary node on fluid side and solidBN is inside solid. vCcut is actual location
//// of the intersection, vL0 is location of LB node being tested, vLd is unit vector pointing
//// toward BL, vC0 and vC1 are wall nodes
//int bcFindIntersectionLS(int2* fluidBN, int2* solidBN, double* dist,
//	double2 vL0, double2 vLd, double2 vC0, double2 vC1)
//{
//
//	// vector between two points BL points
//	double2 vP10 = vC1 - vC0;
//
//	// normal vector
//	double2 vN = (double2)(-vP10.y, vP10.x);
//
//	// Calculate distance between vL0 and intersection
//	double den = dot(vN, vLd);
//	if (den == 0.)
//		return 0;
//	*dist = dot(vN, (vC0 - vL0)) / den;
//
//	// Location of intersection between two lines
//	double2 vCcut = vL0 + vLd * (*dist);
//
//	// Test to see if intersection falls between vC0 and vC1
//	// if falls outside, return 0.
//	double2 vd0 = vCcut - vC0, vd1 = vCcut - vC1;
//	if (fabs(length(vP10) - length(vd0) - length(vd1)) >= CEPS)
//		return 0;
//
//	// should never occur, but placed her just in case.
//	if ((*dist) < 0.)
//		return 0;
//
//	// number of lattice spacings between vL0 and first node beyond
//	// BL
//	int n_cut = (int)ceil((*dist));
//
//	// 
//	*solidBN = convert_int2(vL0) + convert_int2(vLd) * n_cut;
//	*fluidBN = *solidBN - convert_int2(vLd);
//
//	// Limiting the minimum value of dX to avoid instabilities that
//	// it might lead to.
//	(*dist) = max((*dist), CEPS);
//
//	return 1;
//}


//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEFD, 1, 1)))
//void updateBLdXArr(__global ushort* __restrict__ BL_P01ind,
//	__global double2* __restrict__ BL_vNvec,
//	__global double2* __restrict__ BL_blLen,
//	__global int* __restrict__ BL_nodeLoc,
//	__global short* __restrict__ BL_intType,
//	__global double2* C,
//	__global double2* Cnot,
//	__global NTYPE_TYPE* Map,
//	__global double* dXArr)
//{
//	int i = get_global_id(0);
//	if (i >= NUM_BL_TOTAL)
//		return;
//
//	ushort2 p01Ind = BL_P01ind[i];
//
//	double2 vC0 = C[p01Ind.x], vC1 = C[p01Ind.x];
//	double2 vC0not = Cnot[p01Ind.x], vC1not = Cnot[p01Ind.x];
//
//	double2 vT = vC1 - vC0;
//	double2 vTnot = vC1not - vC0not;
//
//	double2 CenterDist = vC0 + vT / 2. - vC0not - vTnot / 2.;
//	BL_intType[i] = (length(CenterDist) >= FOUL_SIZE_SWITCH_SOOT2) ? ((ushort)1) : ((ushort)0);
//
//	blLen = length(vT);
//	vT /= blLen;
//	double2 vCn = double2(vT.y, -vT.x);
//	BL_vNvec[i] = vCn;
//
//	double2 centpos = vP0 + vTvec * blLen / 2.;
//	int2 blind = convert_int2(centpos);
//	BL_blInd[i] = (blind.x - TR_X_IND_START) + (FULLSIZEX_TR_PADDED)* blind.y;
//
//	BL_nodeLoc[i] = (blind.x >= TR_X_IND_START && blind.x < TR_X_IND_STOP) ?
//		((blind.x - TR_X_IND_START) + (FULLSIZEX_TR_PADDED)* blind.y) : -1;
//
//	int2 vCmin = min2(vC0, vC1), vCmax = max2(vC0, vC1);
//
//	if (vCmin.x < 0) vCmin.x = 0;
//	if (vCmin.x >= CHANNEL_LENGTH) return;
//	if (vCmax.x < 0) return;
//	if (vCmax.x >= CHANNEL_LENGTH) vCmax.x = CHANNEL_LENGTH;
//
//	if (vCmin.y < 0) vCmin.y = 0;
//	if (vCmin.y >= CHANNEL_HEIGHT) return;
//	if (vCmax.y < 0) return;
//	if (vCmax.y >= CHANNEL_HEIGHT) vCmax.y = CHANNEL_HEIGHT - 1;
//
//	int dir = 0;
//	double dist;
//	double2 Cdir = (double2)(1., 0.);
//	double Cndir = dot(vCn, Cdir);
//
//	int2 vCcut0, vCcut1;
//	int kk, j;
//
//	// only considering NESW directions, so splitting it up into EW and NS tests
//	if (Cndir != 0.)
//	{
//		if (Cndir > 0.)
//		{
//			dir = 1;
//			Cdir.x = -1.;
//		}
//
//		// Instead of testing all nodes in x and y directions, just iterating
//		// through y direction and testing if Cdir crosses BL regardless of distance
//		// away. The closest node will then be determined from the distance between
//		// the node being tested from and the BL.
//		kk = (dir == 0) ? kk = vCmin.x : vCmax.x;
//
//		for (j = vCmin.y; j <= vCmax.y; j++)
//		{
//			// tests if Cdir crosses BL, and finds BN on either side of
//			// BL.
//			if (bcFindIntersectionLS(&vCcut0, &vCcut1, &dist,
//				convert_double2((int2)(kk, j)), Cdir, vC0, vC1))
//			{
//				int ii0_1D = vCcut0.x + CHANNEL_LENGTH_FULL * vCcut0.y;
//				int ii1_1D = vCcut1.x + CHANNEL_LENGTH_FULL * vCcut1.y;
//				// Should be unnecessary, but testing just in case
//				if (Map[ii0_1D] & M_SOLID_NODE)
//					continue;
//
//				NTYPE_TYPE ntype = Map[ii1_1D];
//				if (ntype & M_FLUID_NODE)
//					continue;
//
//				// Need to reset appropriate bits corresponding to the type of
//				// boundaries as Map was reset to Map0
//				Map[ii1_1D] |= SOLID_BOUNDARY_NODE;
//
//				Map[ii0_1D] |= (dir == 0) ? E_BOUND : W_BOUND;
//
//				// Set corresponding element to distance
//				dX_cur[ii0_1D + dir * DIST_SIZE] = dist;
//
//				// if node is M0_FLUID_NODE, it is now a fouling
//				// layer node and the distance needs to be tracked
//				if (ntype & M0_FLUID_NODE)
//					dX_cur[ii1_1D + RevDir[dir]] = 1. - dist;
//			}
//		}
//	}
//
//	dir = 2;
//	Cdir = (double2)(0., 1.);
//	Cndir = dot(vCn, Cdir);
//	if (Cndir != 0.)
//	{
//		if (Cndir > 0.)
//		{
//			dir = 3;
//			Cdir = (double2)(0., -1.);
//		}
//
//		j = (dir == 2) ? kk = vCmin.y : vCmax.y;
//
//		for (kk = vCmin.x; kk <= vCmax.x; kk++)
//		{
//			// tests if Cdir crosses BL, and finds BN on either side of
//			// BL.
//			if (bcFindIntersectionLS(&vCcut0, &vCcut1, &dist,
//				convert_double2((int2)(kk, j)), Cdir, vC0, vC1))
//			{
//				int ii0_1D = vCcut0.x + CHANNEL_LENGTH_FULL * vCcut0.y;
//				int ii1_1D = vCcut1.x + CHANNEL_LENGTH_FULL * vCcut1.y;
//				// Should be unnecessary, but testing just in case
//				if (Map[ii0_1D] & M_SOLID_NODE)
//					continue;
//
//				NTYPE_TYPE ntype = Map[ii1_1D];
//				if (ntype & M_FLUID_NODE)
//					continue;
//
//				// Need to reset appropriate bits corresponding to the type of
//				// boundaries as Map was reset to Map0
//				Map[ii1_1D] = ntype | SOLID_BOUNDARY_NODE;
//
//				Map[ii0_1D] |= (dir == 2) ? N_BOUND : S_BOUND;
//
//				// Set corresponding element to distance
//				dX_cur[ii0_1D + dir * DIST_SIZE] = dist;
//
//				// if node is M0_FLUID_NODE, it is now a fouling
//				// layer node and the distance needs to be tracked
//				if (ntype & M0_FLUID_NODE)
//					dX_cur[ii1_1D + RevDir[dir]] = 1. - dist;
//			}
//		}
//	}
//}

