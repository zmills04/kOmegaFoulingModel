__kernel void Update_T_Coeffs_Implicit(__global double *__restrict__ Temp,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ Amat,
	__global double *__restrict__ bvec,
	__global double *__restrict__ dXcur,
	__global int *__restrict__ IndArr,
	double alpha, double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	int Cind = IndArr[gid];
	if (Cind < 0)
		return;
	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];

	
	double dX_e = dXcur[gid], dX_w = dXcur[gid + DIST_SIZE], dX_c = dXcur[gid + DIST_SIZE * 2];
	double dY_n = dXcur[gid + DIST_SIZE * 3], dY_s = dXcur[gid + DIST_SIZE * 4], dY_c = dXcur[gid + DIST_SIZE * 5];
	double dX2_e = dXcur[gid + DIST_SIZE * 6], dX2_w = dXcur[gid + DIST_SIZE * 7], dX2_c = dXcur[gid + DIST_SIZE * 8];
	double dY2_n = dXcur[gid + DIST_SIZE * 9], dY2_s = dXcur[gid + DIST_SIZE * 10], dY2_c = dXcur[gid + DIST_SIZE * 11];

	double Ux = ivx[gid], Uy = ivy[gid];
	

	
	double xi_0 = Ux*DTFD;
	double xi_1 = Uy*DTFD;
	double xi_2 = alpha*DTFD;
	double Tc = -dX2_c*xi_2 + dX_c*xi_0 - dY2_c*xi_2 + dY_c*xi_1 + 1.0;
	double Te = -dX2_e*xi_2 + dX_e*xi_0;
	double Tw = -dX2_w*xi_2 + dX_w*xi_0;
	double Tn = -dY2_n*xi_2 + dY_n*xi_1;
	double Ts = -dY2_s*xi_2 + dY_s*xi_1;


	/// Testing w/ periodic domain and Twall = 1;
	double Tsrc = Temp[gid];

	if (Sind >= 0)
		Amat[Sind] = Ts;
	//else
	//	Tsrc -= Tw;

	if (Wind >= 0)
		Amat[Wind] = Tw;
	else if (i == 0)
		Tsrc -= TFD_X_IN_VAL * Tw;
	Amat[Cind] = Tc;

	if (Eind >= 0)
		Amat[Eind] = Te;
	else if (i == CHANNEL_LENGTH-1)
		Amat[Cind] += Te;

	if (Nind >= 0)
		Amat[Nind] = Tn;
	//else
	//	Tsrc -= Tn;

	bvec[gid] = Tsrc;








	//if (Sind >= 0)
	//	Amat[Sind] = Ts;

	//if (Wind >= 0)
	//	Amat[Wind] = Tw;
	//else if (i == 0)
	//	Tsrc -= TFD_X_IN_VAL * Tw;

	//Amat[Cind] = Tc;

	//if (Eind >= 0)
	//	Amat[Eind] = Te;
	//else if (i == CHANNEL_LENGTH - 1)
	//	Amat[Cind] += Te;

	//if (Nind >= 0)
	//	Amat[Nind] = Tn;

	//bvec[gid] = Tsrc;

}

__kernel void Update_T_Coeffs_Implicit_Turbulent(__global double *__restrict__ Temp,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ Amat,
	__global double *__restrict__ bvec,
	__global double *__restrict__ dXcur,
	__global int *__restrict__ IndArr,
	__global double *__restrict__ iAlphat,
	double alpha, double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	int Cind = IndArr[gid];
	if (Cind < 0)
		return;
	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];


	double dX_e = dXcur[gid], dX_w = dXcur[gid + DIST_SIZE], dX_c = dXcur[gid + DIST_SIZE * 2];
	double dY_n = dXcur[gid + DIST_SIZE * 3], dY_s = dXcur[gid + DIST_SIZE * 4], dY_c = dXcur[gid + DIST_SIZE * 5];
	double dX2_e = dXcur[gid + DIST_SIZE * 6], dX2_w = dXcur[gid + DIST_SIZE * 7], dX2_c = dXcur[gid + DIST_SIZE * 8];
	double dY2_n = dXcur[gid + DIST_SIZE * 9], dY2_s = dXcur[gid + DIST_SIZE * 10], dY2_c = dXcur[gid + DIST_SIZE * 11];

	double Ux = ivx[gid], Uy = ivy[gid], alphat = iAlphat[gid];
	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double dAdx = dX_e*iAlphat[gide] + dX_w*iAlphat[gidw] + dX_c*alphat;
	double dAdy = dY_n*iAlphat[gid+CHANNEL_LENGTH_FULL] + 
		dY_s*iAlphat[gide - CHANNEL_LENGTH_FULL] + dY_c*alphat;

	double Jx = DTFD*(Ux - dAdx);
	double Jy = DTFD*(Uy - dAdy);
	alphat *= DTFD;

	double Tc = 1. + (Jx*dX_c + Jy*dY_c - alphat * (dX2_c + dY2_c));
	double Te = (Jx*dX_e - alphat*dX2_e);
	double Tw = (Jx*dX_w - alphat*dX2_w);
	double Tn = (Jy*dY_n - alphat*dY2_n);
	double Ts = (Jy*dY_s - alphat*dY2_s);
	double Tsrc = Temp[gid];

	if (Sind >= 0)
		Amat[Sind] = Ts;

	if (Wind >= 0)
		Amat[Wind] = Tw;
	else if (i == 0)
		Tsrc -= TFD_X_IN_VAL * Tw;

	Amat[Cind] = Tc;

	if (Eind >= 0)
		Amat[Eind] = Te;
	else if (i == CHANNEL_LENGTH - 1)
		Amat[Cind] += Te;

	if (Nind >= 0)
		Amat[Nind] = Tn;

	bvec[gid] = Tsrc;
}

// Since only used one, this kernel will work for both turb and non turb (which is what turbFlag is for)
__kernel void steady_State_T_Coeffs(__global double *__restrict__ Temp,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ Amat,
	__global double *__restrict__ bvec,
	__global double *__restrict__ dXcur,
	__global int *__restrict__ IndArr,
	__global double *__restrict__ iAlphat,
	double alpha, int turbFlag)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	int Cind = IndArr[gid];
	if (Cind < 0)
		return;
	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];


	double dX_e = dXcur[gid], dX_w = dXcur[gid + DIST_SIZE], dX_c = dXcur[gid + DIST_SIZE * 2];
	double dY_n = dXcur[gid + DIST_SIZE * 3], dY_s = dXcur[gid + DIST_SIZE * 4], dY_c = dXcur[gid + DIST_SIZE * 5];
	double dX2_e = dXcur[gid + DIST_SIZE * 6], dX2_w = dXcur[gid + DIST_SIZE * 7], dX2_c = dXcur[gid + DIST_SIZE * 8];
	double dY2_n = dXcur[gid + DIST_SIZE * 9], dY2_s = dXcur[gid + DIST_SIZE * 10], dY2_c = dXcur[gid + DIST_SIZE * 11];

	double Ux = ivx[gid], Uy = ivy[gid], alphat = iAlphat[gid];
	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double dAdy = 0., dAdx = 0.;

	if (turbFlag)
	{
		dAdx = dX_e*iAlphat[gide] + dX_w*iAlphat[gidw] + dX_c*alphat;
		dAdy = dY_n*iAlphat[gid + CHANNEL_LENGTH_FULL] +
			dY_s*iAlphat[gide - CHANNEL_LENGTH_FULL] + dY_c*alphat;
	}
	double Jx = (Ux - dAdx);
	double Jy = (Uy - dAdy);

	double Tc = (Jx*dX_c + Jy*dY_c - alphat * (dX2_c + dY2_c));
	double Te = (Jx*dX_e - alphat*dX2_e);
	double Tw = (Jx*dX_w - alphat*dX2_w);
	double Tn = (Jy*dY_n - alphat*dY2_n);
	double Ts = (Jy*dY_s - alphat*dY2_s);
	double Tsrc = 0.;

	if (Sind >= 0)
		Amat[Sind] = Ts;

	if (Wind >= 0)
		Amat[Wind] = Tw;
	else if (i == 0)
		Tsrc -= TFD_X_IN_VAL * Tw;

	Amat[Cind] = Tc;

	if (Eind >= 0)
		Amat[Eind] = Te;
	else if (i == CHANNEL_LENGTH - 1)
		Amat[Cind] += Te;

	if (Nind >= 0)
		Amat[Nind] = Tn;

	bvec[gid] = Tsrc;
}


//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////      Functions for updating Nu     /////////////////
/////////////////      calculation coefficients      /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

// returns true if vector pointing from node at vL0 in direction perp.
// to BL intersects BL, otherwise false. vC is location vector intersects
// BL, vP0 and vP1 are points defining BL, and dist
// is magnitude of vector between vC and vL0.
bool nuFindIntersectionNormal(double2* vC, double* dist, double2 vL0,
	double2 vP0, double2 vP1)
{
	double2 vP10 = vP1 - vP0;
	double2 vLd = vL0 - vP0;
	double2 vNm = double2(-vP10.y, vP10.x);
	double2 vN = normalize(vNm);
	*dist = dot(vLd, vN);
	double2 vNdist = vN * (*dist);

	*vC = vL0 - vNdist;

	double2 vd0 = *vC - vP0, vd1 = *vC - vP1;

	if (fabs(length(vP10) - length(vd0) - length(vd1)) < EPSILON)
		return true;
	return false;
}


// Finds node which has shortest vector pointing to it from BL in direction
// normal to BL. 
void nuFindNearestNode(__global double2* __restrict__ Cvals,
	__global int* __restrict__ nType,
	__global int* nuYInds,
	__global double2* __restrict__ nuDistVec,
	__global double* __restrict__ nuDist,
	const ushort2 BLind, const int shiftInd)
{
	double2 C0 = Cvals[BLind.x], C1 = vls.C[BLind.y];
	int2 C0i = convert_int2_rtz(fmin(C0, C1));
	int2 C1i = convert_int2_rtp(fmax(C0, C1));
	C0i -= 1;
	C1i += 1;
	C0i = max(C0i, int2(TR_X_IND_START, 0));
	C1i = min(C1i, int2(TR_X_IND_STOP - 1, CHANNEL_HEIGHT - 1));
	if ((C0i.x >= TR_X_IND_STOP) || C1i.x < TR_X_IND_START)
		return;
	
	double dist;
	double2 vCcut, vN;

	for (int ii = C0i.x; ii <= C1i.x; ii++)
	{
		for (int jj = C0i.y; jj <= C1i.y; jj++)
		{
			if (nType(GET_GLOBAL_IDX(ii,jj)) & M_SOLID_NODE)
				continue;

			double2 L0 = double2((double)ii, (double)jj);
			if(nuFindIntersectionNormal(&vCcut, &dist, L0, C0, C1))
			{
				if(AtomicMin(&nuDist[ii - shiftInd], dist)
				{
					nuDistVec[ii - shiftInd] = vCcut;
					nuYInds[ii - shiftInd] = jj;
				}
			}
		}
	}
}


// Updates arrays used to calculate Nu values along wall.
// This is called right before calculating Nu every time, because
// the calculation of Nu is performed less frequently than the 
// updating of the wall location.

// TODO: Set up calculation of Nu to occur more frequently 
//		and average these values in between each save point.
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
void updateNuCoeff1(__global double2* __restrict__ Cvals,
	__global ushort2* __restrict__ BLinds,
	__global int* __restrict__ nType,
	__global int* __restrict__ nuYInds,
	__global double2* __restrict__ nuDistVec,
	__global double* __restrict__ nuDist)
{

	int i = get_global_id(0);

	if (i >= NU_BL_PER_SIDE_NU)
		return;

	nuFindNearestNode(Cvals, nType, nuYInds, nuDistVec,
		nuDist, BLind[i + MIN_BL_BOT], TR_X_IND_START);

	nuFindNearestNode(Cvals, nType, nuYInds, nuDistVec,
		nuDist, BLind[i + MIN_BL_TOP], TR_X_IND_START - NU_NUM_NODES);

}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
void updateNuCoeff2(__global int* __restrict__ nuYInds,
	__global double2* __restrict__ nuDistVec)
{
	int i = get_global_id(0);

	if (i >= NU_NUM_NODES*2)
		return;

	int jj = nuYInds[i];
	int nodeType = 0x0;
	if (jj != -1)
	{
		if (nuDistVec[i].x < 0.)
		{
			nodeType |= 0x1;
		}
		if (nuDistVec[i].y < 0.)
		{
			nodeType |= 0x2;
		}
	}
	nuYInds[i] = (nuYInds(i) << 2) | nodeType;
}


// Calculates bulk mean temperature
double getTbm(__global double* __restict__ T_array,
	__global double* __restict__ Ux_array,
	__global double* __restict__ Uy_array,
	int ival)
{
	double Um = 0.;
	double Tm = 0.;
	for (int j = 0; j < FULLSIZEY; j++)
	{
		int linInd = ival + j * FULLSIZE_Y;
		double Uval = sqrt(Ux_array[linInd] * Ux_array[linInd] + Uy_array[linInd] * Uy_array[linInd]);
		Um += Uval;
		Tm += Uval * T_array[linInd];
	}
	return Tm / Um;
}


// TODO: make sure that there is not any sign errors in calculation
//		 of dT or Ts.
double getTandDT(__global double* __restrict__ Tvals,
	__global double* __restrict__ dXvals,
	__global double* __restrict__ dX0vals,
	double2 vNorm, int linInd, int nType, double* dTval)
{
	double Tnode = Tvals[linInd];
	double Tdiff, dx, dy;
	if (nType & 0x1)
	{
		Tdiff = Tvals[linInd + 1] - Tnode; // east neigh
		dx = dXvals[linInd * 8]; // dxe
		dx0 = dX0vals[linInd * 8];
	}
	else
	{
		Tneigh = Tnode - Tvals[linInd - 1]; // west neigh
		dx = dXvals[linInd * 8 + 1]; // dxw
		dx0 = dX0vals[linInd * 8 + 1];
	}

	double2 dTdx;
	double cCoeff = KSOOT * dx / (KAIR * (dx0 - dx));
	dTdx.x = ((dx0 - dx) > CEPS) ?
		(KSOOT * (Tdiff) / (KAIR * (dx0 - dx) * (1. + cCoeff))) :
		((Tdiff) / dx);

	if (nType & 0x2)
	{
		Tdiff = Tvals[linInd + CHANNEL_LENGTH_FULL] - Tnode; // north neigh
		dx = dXvals[linInd * 8 + 3]; // dyn
		dx0 = dX0vals[linInd * 8 + 3];
	}
	else
	{
		Tdiff = Tnode - Tvals[linInd - CHANNEL_LENGTH_FULL]; // south neigh
		dx = dXvals[linInd * 8 + 4]; // dys
		dx0 = dX0vals[linInd * 8 + 4];
	}


	dTdx.y = KSOOT * Tdiff / (KAIR * (dx0 - dx) + KSOOT * dx);

	*dTval = dot(dTdx, normalize(vNorm));

	Tnode -= (*dTval) * length(vNorm);


	return Tnode;
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
void FD_calc_Nu(__global double* __restrict__ T_array,
	__global double* __restrict__ Nu_array,
	__global double* __restrict__ Ux_array,
	__global double* __restrict__ Uy_array,
	__global int* __restrict__ nuYInds,
	__global double2 * __restrict__ nuDist,
	__global double* __restrict__ dXvals,
	__global double* __restrict__ dX0vals,
	int timeIndex)
{
	int ix = get_global_id(0);

	if (ix >= NU_NUM_NODES)
		return;

	int ii = ix + TR_X_IND_START;

	double Tm = getTbm(T_array, Ux_array, Uy_array, ii);

	int jj = nuYInds[ix];

	int nType = (jj & 0x3);
	jj = jj >> 2;

	double2 dTdn = 0.;
	double2 Ts = 0.;

	if (jj != -1)
	{
		int linInd = ii + jj * CHANNEL_LENGTH_FULL;
		double2 vNorm = nuDist[ix];

		Ts.x = getTandDT(T_array, dXvals, dX0vals,
			vNorm, linInd, nType, &(dTdn.x));
	}

	jj = nuYInds[ix + FULLSIZEX_TR];

	nType = (jj & 0x3);
	jj = jj >> 2;

	if (jj != -1)
	{
		int linInd = ii + jj * CHANNEL_LENGTH_FULL;
		double2 vNorm = nuDist[ix + FULLSIZEX_TR];

		Ts.y = getTandDT(T_array, dXvals, dX0vals,
			vNorm, linInd, nType, &(dTdn.y));
	}

	timeIndex *= NU_NUM_NODES_FULLSIZE;
	
	Nu_Array[shiftval + ix * 5]		= dTdn.x;
	Nu_array[shiftval + ix * 5 + 1] = Ts.x;
	Nu_array[shiftval + ix * 5 + 2] = dTdn.y;
	Nu_array[shiftval + ix * 5 + 3] = Ts.y;
	Nu_array[shiftval + ix * 5 + 4] = Tm;
}

// Calculates Bulk Mean at beginning of each period
// global size = (number of periods+1, Height of Domain)
// local size = (1, Height of Domain)
#ifdef OPENCL_VERSION_1_2
__kernel void getNuMean(__global double* __restrict__ Uxvals,
	__global double* __restrict__ Uyvals,
	__global double* __restrict__ Tvals,
	__global double* __restrict__ Nuvals,
	const int saveNum,
	__local double* UTdata,
	__local double* Udata)
{
	const int tid = get_local_id(1);
	const int gid = get_global_id(0);
	const int xloc = WAVY_SECTION_START + gid * WAVY_PERIOD_LENGTH;
	const int localSize = get_local_size(1);
	const int stride = tid * 2;

	// top row of domain will always have solid nodes, so no need to worry if 
	double Uxy1 = 0., Uxy2 = 0., T1 = 0., T2 = 0.;
	if (stride < CHANNEL_HEIGHT)
	{
		int fullLoc = xloc + stride * CHANNEL_LENGTH_FULL;
		Uxy1 = sqrt(pown(Uxvals[fullLoc], 2) + pown(Uyvals[fullLoc], 2));
		T1 = Tvals[fullLoc];
	}
	if ((stride + 1) < CHANNEL_HEIGHT)
	{
		int fullLoc = xloc + (stride + 1) * CHANNEL_LENGTH_FULL;
		Uxy2 = sqrt(pown(Uxvals[fullLoc], 2) + pown(Uyvals[fullLoc], 2));
		T2 = Tvals[fullLoc];
	}

	UTdata[tid] = Uxy1 * T1 + Uxy2 * T2;
	Udata[tid] = Uxy1 + Uxy2;

	barrier(CLK_LOCAL_MEM_FENCE);
	for (int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			UTdata[tid] += UTdata[tid + s];
			Udata[tid] += Udata[tid + s];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		Nuvals[2 * gid + saveNum * NUM_NU_LOCS * 2] = UTdata[0];
		Nuvals[2 * gid + 1 + saveNum * NUM_NU_LOCS * 2] = Udata[0];
	}
}
#else

__kernel void getNuMean(__global double* __restrict__ Uxvals,
	__global double* __restrict__ Uyvals,
	__global double* __restrict__ Tvals,
	__global double* __restrict__ Nuvals,
	const int saveNum)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	const int xloc = WAVY_SECTION_START + i * WAVY_PERIOD_LENGTH;
	int fullLoc = xloc + j * CHANNEL_LENGTH_FULL;

	double Uval = (j < CHANNEL_HEIGHT) ?
		sqrt(pown(Uxvals[fullLoc], 2) + pown(Uyvals[fullLoc], 2)) : 0.;
	double Tval = (j < CHANNEL_HEIGHT) ? Tvals[fullLoc] : 0.;

	Tval *= Uval;

	barrier(CLK_LOCAL_MEM_FENCE);

	Nuvals[2 * i + saveNum * NUM_NU_LOCS * 2] = work_group_reduce_max(Uval);
	Nuvals[2 * i + 1 + saveNum * NUM_NU_LOCS * 2] = work_group_reduce_max(Tval);
}

#endif



//__kernel void Update_T_Coeffs(__global double* __restrict__ Temp,
//	__global double *__restrict__ ivx,
//	__global double *__restrict__ ivy,
//	__global double *__restrict__ Amat,
//	__global double *__restrict__ bvec,
//	__global double *__restrict__ dXcur,
//	__global int *__restrict__ IndArr,
//	double alpha, double DTFD)
//{
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
//		return;
//
//	int gid = i + CHANNEL_LENGTH_FULL*j;
//
//	int Cind = IndArr[gid];
//	if (Cind < 0)
//		return;
//	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
//	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];
//
//
//	double Ux = ivx[gid], Uy = ivy[gid];
//	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx = dx_e + dx_w;
//	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy = dy_n + dy_s;
//
//	double xi_0 = 1.0 / dx_e;
//	double xi_1 = Ux*DTFD;
//	double xi_2 = xi_0*xi_1;
//	double xi_3 = 1.0 / dy_n;
//	double xi_4 = Uy*DTFD;
//	double xi_5 = xi_3*xi_4;
//	double xi_6 = 1.0 / dx_w;
//	double xi_7 = xi_1*xi_6;
//	double xi_8 = 1.0 / dy_s;
//	double xi_9 = xi_4*xi_8;
//	double xi_10 = 2.0*alpha*DTFD;
//	double xi_11 = 1.0 / dx;
//	double xi_12 = 2.0*alpha*DTFD*xi_11;
//	double xi_13 = 1.0 / dy;
//	double xi_14 = 2.0*alpha*DTFD*xi_13;
//	//Tc = -xi_0*xi_10*xi_6 - xi_10*xi_3*xi_8 + xi_2 + xi_5 - xi_7 - xi_9 + 1.0;
//	//Te = -dx_w*xi_11*xi_2 + xi_0*xi_12;
//	//Tw = dx_e*xi_11*xi_7 + xi_12*xi_6;
//	//Tn = -dy_s*xi_13*xi_5 + xi_14*xi_3;
//	//Ts = dy_n*xi_13*xi_9 + xi_14*xi_8;
//
//	if (Sind >= 0)
//		Amat[Sind] = dy_n*xi_13*xi_9 + xi_14*xi_8;
//
//	if (Wind >= 0)
//		Amat[Wind] = dx_e*xi_11*xi_7 + xi_12*xi_6;
//	else if (i == 0)
//		bvec[gid] = TFD_X_IN_VAL * (dx_e*xi_11*xi_7 + xi_12*xi_6);
//
//	Amat[Cind] = -xi_0*xi_10*xi_6 - xi_10*xi_3*xi_8 + xi_2 + xi_5 - xi_7 - xi_9 + 1.0;
//
//	if (Eind >= 0)
//		Amat[Eind] = -dx_w*xi_11*xi_2 + xi_0*xi_12;
//	else if (i == CHANNEL_LENGTH - 1)
//		bvec[gid] = Temp[gid] * (-dx_w*xi_11*xi_2 + xi_0*xi_12);
//
//	if (Nind >= 0)
//		Amat[Nind] = -dy_s*xi_13*xi_5 + xi_14*xi_3;
//
//}


//__kernel void Update_T_Coeffs_Implicit(__global double *__restrict__ Temp,
//	__global double *__restrict__ ivx,
//	__global double *__restrict__ ivy,
//	__global double *__restrict__ Amat,
//	__global double *__restrict__ bvec,
//	__global double *__restrict__ dXcur,
//	__global int *__restrict__ IndArr,
//	double alpha, double DTFD)
//{
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
//		return;
//
//	int gid = i + CHANNEL_LENGTH_FULL*j;
//
//	int Cind = IndArr[gid];
//	if (Cind < 0)
//		return;
//	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
//	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];
//
//
//	double Ux = ivx[gid], Uy = ivy[gid];
//	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx = dx_e + dx_w;
//	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy = dy_n + dy_s;
//
//	double xi_0 = 1.0 / dx_e;
//	double xi_1 = Ux*DTFD;
//	double xi_2 = xi_0*xi_1;
//	double xi_3 = 1.0 / dy_n;
//	double xi_4 = Uy*DTFD;
//	double xi_5 = xi_3*xi_4;
//	double xi_6 = 1.0 / dx_w;
//	double xi_7 = xi_1*xi_6;
//	double xi_8 = 1.0 / dy_s;
//	double xi_9 = xi_4*xi_8;
//	double xi_10 = 2.0*alpha*DTFD;
//	double xi_11 = 1.0 / dx;
//	double xi_12 = 2.0*alpha*DTFD*xi_11;
//	double xi_13 = 1.0 / dy;
//	double xi_14 = 2.0*alpha*DTFD*xi_13;
//	double Tc = xi_0*xi_10*xi_6 + xi_10*xi_3*xi_8 + xi_2 + xi_5 - xi_7 - xi_9 + 1.0;
//	double Te = dx_w*xi_11*xi_7 - xi_12*xi_6;
//	double Tw = -dx_e*xi_11*xi_2 - xi_0*xi_12;
//	double Tn = dy_s*xi_13*xi_9 - xi_14*xi_8;
//	double Ts = -dy_n*xi_13*xi_5 - xi_14*xi_3;
//	
//	double Tsrc = Temp[gid];
//	
//	if (Sind >= 0)
//		Amat[Sind] = Ts;
//
//	if (Wind >= 0)
//		Amat[Wind] = Tw;
//	else if (i == 0)
//		Tsrc -= TFD_X_IN_VAL * Tw;
//
//	Amat[Cind] = Tc;
//
//	if (Eind >= 0)
//		Amat[Eind] = Te;
//	else if (i == CHANNEL_LENGTH - 1)
//		Amat[Cind] += Te;
//
//	if (Nind >= 0)
//		Amat[Nind] = Tn;
//
//	bvec[gid] = Tsrc;
//
//	//double Tsrc = Temp[gid];
//
//	//if (Sind >= 0)
//	//	Amat[Sind] = DTFD*(-Uy / 2. - alpha);
//
//	//if (Wind >= 0)
//	//	Amat[Wind] = DTFD*(-Ux / 2. - alpha);
//	//else if (i == 0)
//	//	Tsrc += TFD_X_IN_VAL * DTFD*((Ux / 2. + alpha));
//
//	//Amat[Cind] = 1 + 4.*alpha;
//
//	//if (Eind >= 0)
//	//	Amat[Eind] = DTFD*(Ux / 2. - alpha);
//	//else if (i == CHANNEL_LENGTH - 1)
//	//	Amat[Cind] += DTFD*((Ux / 2. - alpha));
//
//	//if (Nind >= 0)
//	//	Amat[Nind] = DTFD*(Uy / 2. - alpha);
//
//	//bvec[gid] = Tsrc;
//
//
//}


