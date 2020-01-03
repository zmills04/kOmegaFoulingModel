// Update_T_coeffs_Implicit now handles both turbulent and non-turbulent flows

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
//	double dX_e = dXcur[gid], dX_w = dXcur[gid + DIST_SIZE], dX_c = dXcur[gid + DIST_SIZE * 2];
//	double dY_n = dXcur[gid + DIST_SIZE * 3], dY_s = dXcur[gid + DIST_SIZE * 4], dY_c = dXcur[gid + DIST_SIZE * 5];
//	double dX2_e = dXcur[gid + DIST_SIZE * 6], dX2_w = dXcur[gid + DIST_SIZE * 7], dX2_c = dXcur[gid + DIST_SIZE * 8];
//	double dY2_n = dXcur[gid + DIST_SIZE * 9], dY2_s = dXcur[gid + DIST_SIZE * 10], dY2_c = dXcur[gid + DIST_SIZE * 11];
//
//	double Ux = ivx[gid], Uy = ivy[gid];
//	
//
//	
//	double xi_0 = Ux*DTFD;
//	double xi_1 = Uy*DTFD;
//	double xi_2 = alpha*DTFD;
//	double Tc = -dX2_c*xi_2 + dX_c*xi_0 - dY2_c*xi_2 + dY_c*xi_1 + 1.0;
//	double Te = -dX2_e*xi_2 + dX_e*xi_0;
//	double Tw = -dX2_w*xi_2 + dX_w*xi_0;
//	double Tn = -dY2_n*xi_2 + dY_n*xi_1;
//	double Ts = -dY2_s*xi_2 + dY_s*xi_1;
//
//
//	/// Testing w/ periodic domain and Twall = 1;
//	double Tsrc = Temp[gid];
//
//	if (Sind >= 0)
//		Amat[Sind] = Ts;
//	//else
//	//	Tsrc -= Tw;
//
//	if (Wind >= 0)
//		Amat[Wind] = Tw;
//	else if (i == 0)
//		Tsrc -= TFD_X_IN_VAL * Tw;
//	Amat[Cind] = Tc;
//
//	if (Eind >= 0)
//		Amat[Eind] = Te;
//	else if (i == CHANNEL_LENGTH-1)
//		Amat[Cind] += Te;
//
//	if (Nind >= 0)
//		Amat[Nind] = Tn;
//	//else
//	//	Tsrc -= Tn;
//
//	bvec[gid] = Tsrc;
//
//
//
//
//
//
//
//
//	//if (Sind >= 0)
//	//	Amat[Sind] = Ts;
//
//	//if (Wind >= 0)
//	//	Amat[Wind] = Tw;
//	//else if (i == 0)
//	//	Tsrc -= TFD_X_IN_VAL * Tw;
//
//	//Amat[Cind] = Tc;
//
//	//if (Eind >= 0)
//	//	Amat[Eind] = Te;
//	//else if (i == CHANNEL_LENGTH - 1)
//	//	Amat[Cind] += Te;
//
//	//if (Nind >= 0)
//	//	Amat[Nind] = Tn;
//
//	//bvec[gid] = Tsrc;
//
//}

// TODO: check to see how much penalty is associated with using dXArr values to calculate
//		d(nut)/dx values rather than dX0Arr, and how significant its effect on results is.
//		Implement if not significant performance penalty and has non-negligible effect.

// TODO: Check to see if setting the dXe, dXc and dXw coefficients to zero for outlet nodes
//		provides better results than adding coefficient for east node to center node coefficient
//		and using 0 for east coefficient in A matrix.

// TODO: check to see if it is faster to not use boundary elements in sparse matrix for Temp calculation
	


// Left off at creating kernels for calculating dAlpha and dRhoCp (see below old kernel for calculating
// steady state coefficients. Mistake made in implementing source term for CHT, by having it only be applied
// when k-omega is used. This should have its own set of ifdef's which will allow it to be turned on and off
// for testing (using ifdefs in kernels and flag in parameter file wont add overhead since it is only implemented
// in kernels outside of initialization of variables and ifdefs will avoid need for if statements everywhere)


// TODO: Set all diagonal elements to 1.0 at beginning of simulation
//		Since preconditioning with diagonal preconditioner, these elements will
//		always be 1.
__kernel void Update_T_Coeffs_Implicit(__global int* __restrict__ IndArr,//0
	__global NTYPE_TYPE* __restrict__ Map,//1
	__global double* __restrict__ dXcur,//2
	__global double* __restrict__ Amat,//3
	__global double* __restrict__ bvec,//4
	__global double *__restrict__ Temp,//5
	__global double* __restrict__ dAlpha,//6
	__global double *__restrict__ ivx,//7
	__global double *__restrict__ ivy,//8
#ifdef USING_KOMEGA_SOLVER
	__global double* __restrict__ iNut,//9
#endif
#ifdef USING_CHT_SOURCE_CORRECTION
	__global double* __restrict__ dRhoCp,
#endif
	double DTFD)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;

	NTYPE_TYPE ntype = Map[gid];

	// return if inside wall of channel (not fouling layer)
	// solid boundary nodes do not need coefficients updated,
	// so no need to handle these nodes here.
	if (ntype & M0_SOLID_NODE)
		return;

	//double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx_c = dx_e+dx_w;
	//double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy_c = dy_n+dy_s;

	double dx_e = max(dXcur[gid], 0.1), dx_w = max(dXcur[gid + DIST_SIZE], 0.1), dx_c = max(dx_e + dx_w, 0.1);
	double dy_n = max(dXcur[gid + DIST_SIZE * 2], 0.1), dy_s = max(dXcur[gid + DIST_SIZE * 3], 0.1), dy_c = max(dy_n + dy_s, 0.1);


	double dXe = dx_w / (dx_e * dx_c), dXw = -dx_e / (dx_w * dx_c), dXc = (dx_e - dx_w) / (dx_e * dx_w);
	double dYn = dy_s / (dy_n * dy_c), dYs = -dy_n / (dy_s * dy_c), dYc = (dy_n - dy_s) / (dy_n * dy_s);

	double dX2e = 2. / (dx_e * dx_c), dX2w = 2. / (dx_w * dx_c), dX2c = -2. / (dx_e * dx_w);
	double dY2n = 2. / (dy_n * dy_c), dY2s = 2. / (dy_s * dy_c), dY2c = -2. / (dy_n * dy_s);

	if (fabs(dXc) > 10.)
	{
		dXe = 1.0 / dx_c;
		dXw = -1.0 / dx_c;
		dXc = 0.0;		
	}

	if (fabs(dYn) > 10.)
	{
		dYn = 1.0 / dy_c;
		dYs = -1.0 / dy_c;
		dYc = 0.0;
	}

	

	double alphaVal = (ntype & M_SOLID_NODE) ? LB_ALPHA_FOUL : LB_ALPHA_FLUID;

#ifdef USING_KOMEGA_SOLVER
	double alphatVal = iNut[gid] * PR_TURB_NUMBER_INV;
#endif

//#ifdef USING_CHT_SOURCE_CORRECTION
//	double Jx = ivx[gid];
//	double Jy = ivy[gid];
//
//	double Dsource = Jx * dRhoCp[gid * 2] + Jy * dRhoCp[gid * 2 + 1];
//
//	#ifdef USING_KOMEGA_SOLVER
//		int gidw = ((i > 0) ? -1 : (CHANNEL_LENGTH - 1)) + gid;
//		int gide = ((i < CHANNEL_LENGTH - 1) ? 1 : (1 - CHANNEL_LENGTH_FULL)) + gid;
//		Jx -= dAlpha[gid * 2] + PR_TURB_NUMBER_INV * (dXe * iNut[gide] + dXw * iNut[gidw] + dXc * alphatVal);
//		Jy -= dAlpha[gid * 2+1] + PR_TURB_NUMBER_INV * (dYn * iNut[gid+CHANNEL_LENGTH_FULL] + dYs * iNut[gid-CHANNEL_LENGTH_FULL] + dYc * alphatVal);
//		alphaVal += alphatVal;
//	#else
//		Jx -= dAlpha[gid * 2];
//		Jy -= dAlpha[gid * 2 + 1];
//	#endif
//
//	double Tsrc = Temp[gid] * (1. + max(Dsource, 0.));
//	double Ccoeff = 1. + DTFD*(Jx*dXc + Jy*dYc - alphaVal*(dX2c+dY2c) - min(Dsource,0.));
//#else
//	#ifdef USING_KOMEGA_SOLVER
//		int gidw = ((i > 0) ? -1 : (CHANNEL_LENGTH - 1)) + gid;
//		int gide = ((i < CHANNEL_LENGTH - 1) ? 1 : (1 - CHANNEL_LENGTH_FULL)) + gid;
//		double Jx = ivx[gid] - (dAlpha[gid * 2] + PR_TURB_NUMBER_INV * (dXe * iNut[gide] + dXw * iNut[gidw] + dXc * alphatVal));
//		double Jy = ivy[gid] - (dAlpha[gid * 2 + 1] + PR_TURB_NUMBER_INV * (dYn * iNut[gid + CHANNEL_LENGTH_FULL] + 
//			dYs * iNut[gid - CHANNEL_LENGTH_FULL] + dYc * alphatVal));
//		alphaVal += alphatVal;
//	#else
//		double Jx = ivx[gid] - dAlpha[gid * 2];
//		double Jy = ivy[gid] - dAlpha[gid * 2 + 1];
//	#endif
//
//	double Tsrc = Temp[gid];
//	double Ccoeff = 1. + DTFD * (Jx * dXc + Jy * dYc - alphaVal * (dX2c + dY2c));
//#endif




	int gidw = ((i > 0) ? -1 : (CHANNEL_LENGTH - 1)) + gid;
	int gide = ((i < CHANNEL_LENGTH - 1) ? 1 : (1 - CHANNEL_LENGTH_FULL)) + gid;
	double Jx = ivx[gid] - (PR_TURB_NUMBER_INV * (dXe * iNut[gide] + dXw * iNut[gidw] + dXc * alphatVal));
	double Jy = ivy[gid] - (PR_TURB_NUMBER_INV * (dYn * iNut[gid + CHANNEL_LENGTH_FULL] +
		dYs * iNut[gid - CHANNEL_LENGTH_FULL] + dYc * alphatVal));
	alphaVal += alphatVal;

	double Tsrc = min(max(0.0, Temp[gid]), TFD_X_IN_VAL);
	double Ccoeff = 1. + DTFD * (Jx * dXc + Jy * dYc - alphaVal * (dX2c + dY2c));





	// No BCs to handle for south coefficient, so calculate and set element in A
	// directly
	double Scoeff = DTFD * (Jy * dYs - alphaVal * dY2s);
	
	// W coefficient is 0 for inlet nodes, so Wcoeff must be multiplied
	// by the inlet temperature and moved to the LHS of equation
	double Wcoeff = DTFD * (Jx * dXw - alphaVal * dX2w);
	Tsrc -= (i == 0) ? (TFD_X_IN_VAL * Wcoeff) : 0.;
	
	// E coefficient is 0 for outlet nodes and Ecoeff is added to
	// Ccoeff to essentially set Temp of east node to temp of outlet node
	double Ecoeff = DTFD * (Jx * dXe - alphaVal * dX2e);
	Ccoeff += (i == CHANNEL_LENGTH-1) ? (Ecoeff) : 0.;



	// now that Ccoeff will not change, setting Amat values with all values divided by Ccoeff
	// Ccoeff inverted to multiply other coefficients by rather than dividing it
	Ccoeff = 1. / Ccoeff;

	Amat[IndArr[gid + DIST_SIZE * 4]] = Scoeff * Ccoeff;

	if (i > 0)
		Amat[IndArr[gid + DIST_SIZE * 2]] = Wcoeff * Ccoeff;

	if (i < CHANNEL_LENGTH - 1)
		Amat[IndArr[gid + DIST_SIZE]] = Ecoeff * Ccoeff;

	// Set C coefficient after checking to see if Ecoeff is added to it
	Amat[IndArr[gid]] = 1.;

	// North coefficient Same situation as S coefficient, so set directly
	Amat[IndArr[gid + DIST_SIZE * 3]] = DTFD*(Jy * dYn - alphaVal * dY2n) * Ccoeff;

	// Set corresponding element in b vector.
	bvec[gid] = Tsrc * Ccoeff;
}


__kernel void Steady_State_T_Coeffs(__global int* __restrict__ IndArr,
	__global NTYPE_TYPE* __restrict__ Map,
	__global double* __restrict__ dXcur,
	__global double* __restrict__ Amat,
	__global double* __restrict__ bvec,
	__global double* __restrict__ Temp,
	__global double* __restrict__ dAlpha,
#ifdef USING_KOMEGA_SOLVER
	__global double* __restrict__ iNut,
#endif
#ifdef USING_CHT_SOURCE_CORRECTION
	__global double* __restrict__ dRhoCp,
#endif
	__global double* __restrict__ ivx,
	__global double* __restrict__ ivy)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH)
		return;

	int gid = i + CHANNEL_LENGTH_FULL * j;

	NTYPE_TYPE ntype = Map[gid];

	// return if inside wall of channel (not fouling layer)
	// solid boundary nodes do not need coefficients updated,
	// so no need to handle these nodes here.
	if (ntype & M0_SOLID_NODE)
		return;


	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx_c = dx_e + dx_w;
	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy_c = dy_n + dy_s;

	double dXe = dx_w / (dx_e * dx_c), dXw = -dx_e / (dx_w * dx_c), dXc = (dx_e - dx_w) / (dx_e * dx_w);
	double dYn = dy_s / (dy_n * dy_c), dYs = -dy_n / (dy_s * dy_c), dYc = (dy_n - dy_s) / (dy_n * dy_s);

	double dX2e = 2. / (dx_e * dx_c), dX2w = 2. / (dx_w * dx_c), dX2c = -2. / (dx_e * dx_w);
	double dY2n = 2. / (dy_n * dy_c), dY2s = 2. / (dy_s * dy_c), dY2c = -2. / (dy_n * dy_s);

	double alphaVal = LB_ALPHA_FLUID;

	double Jx = ivx[gid];
	double Jy = ivy[gid];

	double Tsrc = 0.;
	double Ccoeff = (Jx * dXc + Jy * dYc - alphaVal * (dX2c + dY2c));

	// No BCs to handle for south coefficient, so calculate and set element in A
	// directly
	Amat[IndArr[gid + DIST_SIZE * 4]] = (Jy * dYs - alphaVal * dY2s);

	// W coefficient is 0 for inlet nodes, so Wcoeff must be multiplied
	// by the inlet temperature and moved to the LHS of equation
	double Wcoeff = (Jx * dXw - alphaVal * dX2w);
	if (i > 0)
		Amat[IndArr[gid + DIST_SIZE * 2]] = Wcoeff;
	else
		Tsrc -= TFD_X_IN_VAL * Wcoeff;

	// E coefficient is 0 for outlet nodes and Ecoeff is added to
	// Ccoeff to essentially set Temp of east node to temp of outlet node
	double Ecoeff = (Jx * dXe - alphaVal * dX2e);
	if (i < CHANNEL_LENGTH - 1)
		Amat[IndArr[gid + DIST_SIZE]] = Ecoeff;
	else
		Ccoeff += Ecoeff;

	// Set C coefficient after checking to see if Ecoeff is added to it
	Amat[IndArr[gid]] = Ccoeff;

	// North coefficient Same situation as S coefficient, so set directly
	Amat[IndArr[gid + DIST_SIZE * 3]] = (Jy * dYn - alphaVal * dY2n);

	// Set corresponding element in b vector.
	bvec[gid] = Tsrc;
}



//__kernel void Steady_State_T_Coeffs(__global int* __restrict__ IndArr,
//	__global NTYPE_TYPE* __restrict__ Map,
//	__global double* __restrict__ dXcur,
//	__global double* __restrict__ Amat,
//	__global double* __restrict__ bvec,
//	__global double* __restrict__ Temp,
//	__global double* __restrict__ dAlpha,
//#ifdef USING_KOMEGA_SOLVER
//	__global double* __restrict__ iNut,
//#endif
//#ifdef USING_CHT_SOURCE_CORRECTION
//	__global double* __restrict__ dRhoCp,
//#endif
//	__global double* __restrict__ ivx,
//	__global double* __restrict__ ivy)
//{
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//	if (i >= CHANNEL_LENGTH)
//		return;
//
//	int gid = i + CHANNEL_LENGTH_FULL * j;
//
//	NTYPE_TYPE ntype = Map[gid];
//
//	// return if inside wall of channel (not fouling layer)
//	// solid boundary nodes do not need coefficients updated,
//	// so no need to handle these nodes here.
//	if (ntype & M0_SOLID_NODE)
//		return;
//
//	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx_c = dx_e + dx_w;
//	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy_c = dy_n + dy_s;
//
//	double dXe = dx_w / (dx_e * dx_c), dXw = -dx_e / (dx_w * dx_c), dXc = (dx_e - dx_w) / (dx_e * dx_w);
//	double dYn = dy_s / (dy_n * dy_c), dYs = -dy_n / (dy_s * dy_c), dYc = (dy_n - dy_s) / (dy_n * dy_s);
//
//	double dX2e = 2. / (dx_e * dx_c), dX2w = 2. / (dx_w * dx_c), dX2c = -2. / (dx_e * dx_w);
//	double dY2n = 2. / (dy_n * dy_c), dY2s = 2. / (dy_s * dy_c), dY2c = -2. / (dy_n * dy_s);
//
//	int i_w = (i > 0) ? i - 1 : (CHANNEL_LENGTH - 1);
//	int i_e = (i >= CHANNEL_LENGTH) ? 0 : i + 1;
//	int gidw = i_w + CHANNEL_LENGTH_FULL * j;
//	int gide = i_e + CHANNEL_LENGTH_FULL * j;
//
//
//
//	double alphaVal = LB_ALPHA_FLUID;
//
//
//#ifdef USING_KOMEGA_SOLVER
//	double alphatVal = iNut[gid] * PR_TURB_NUMBER_INV;
//#endif
//
//
//
//#ifdef USING_CHT_SOURCE_CORRECTION
//	double Jx = ivx[gid];
//	double Jy = ivy[gid];
//
//	double Dsource = Jx * dRhoCp[gid * 2] + Jy * dRhoCp[gid * 2 + 1];
//
//#ifdef USING_KOMEGA_SOLVER
//	Jx -= dAlpha[gid * 2] + PR_TURB_NUMBER_INV * (dXe * iNut[gide] + dXw * iNut[gidw]) + dXc * alphatVal;
//	Jy -= dAlpha[gid * 2 + 1] + PR_TURB_NUMBER_INV * (dYn * iNut[gid + CHANNEL_LENGTH_FULL] + dYs * iNut[gid - CHANNEL_LENGTH_FULL]) + dYs * alphatVal;
//	alphaVal += alphatVal;
//#else
//	Jx -= dAlpha[gid * 2];
//	Jy -= dAlpha[gid * 2 + 1];
//#endif
//
//	double Tsrc = Temp[gid] * (max(Dsource, 0.));
//	double Ccoeff = (Jx * dXc + Jy * dYc - alphaVal * (dX2c + dY2c) - min(Dsource, 0.));
//
//#else
//
//#ifdef USING_KOMEGA_SOLVER
//	double Jx = ivx[gid] - (dAlpha[gid * 2] + PR_TURB_NUMBER_INV * (dXe * iNut[gide] + dXw * iNut[gidw]) + dXc * alphatVal);
//	double Jy = ivy[gid] - (dAlpha[gid * 2 + 1] + PR_TURB_NUMBER_INV * (dYn * iNut[gid + CHANNEL_LENGTH_FULL] +
//		dYs * iNut[gid - CHANNEL_LENGTH_FULL]) + dYs * alphatVal);
//	alphaVal += alphatVal;
//#else
//	double Jx = ivx[gid] - dAlpha[gid * 2];
//	double Jy = ivy[gid] - dAlpha[gid * 2 + 1];
//#endif
//
//	double Tsrc = 0.;
//	double Ccoeff = (Jx * dXc + Jy * dYc - alphaVal * (dX2c + dY2c));
//#endif
//
//	// No BCs to handle for south coefficient, so calculate and set element in A
//	// directly
//	Amat[IndArr[gid + DIST_SIZE * 4]] = (Jy * dYs - alphaVal * dY2s);
//
//	// W coefficient is 0 for inlet nodes, so Wcoeff must be multiplied
//	// by the inlet temperature and moved to the LHS of equation
//	double Wcoeff = (Jx * dXw - alphaVal * dX2w);
//	if (i > 0)
//		Amat[IndArr[gid + DIST_SIZE * 2]] = Wcoeff;
//	else
//		Tsrc -= TFD_X_IN_VAL * Wcoeff;
//
//	// E coefficient is 0 for outlet nodes and Ecoeff is added to
//	// Ccoeff to essentially set Temp of east node to temp of outlet node
//	double Ecoeff = (Jx * dXe - alphaVal * dX2e);
//	if (i < CHANNEL_LENGTH - 1)
//		Amat[IndArr[gid + DIST_SIZE]] = Ecoeff;
//	else
//		Ccoeff += Ecoeff;
//
//	// Set C coefficient after checking to see if Ecoeff is added to it
//	Amat[IndArr[gid]] = Ccoeff;
//
//	// North coefficient Same situation as S coefficient, so set directly
//	Amat[IndArr[gid + DIST_SIZE * 3]] = (Jy * dYn - alphaVal * dY2n);
//
//	// Set corresponding element in b vector.
//	bvec[gid] = Tsrc;
//}



//Old version of steady state coefficient kernel

//// Since only used one, this kernel will work for both turb and non turb (which is what turbFlag is for)
//__kernel void steady_State_T_Coeffs(__global double *__restrict__ Temp,
//	__global double *__restrict__ ivx,
//	__global double *__restrict__ ivy,
//	__global double *__restrict__ Amat,
//	__global double *__restrict__ bvec,
//	__global double *__restrict__ dXcur,
//	__global int *__restrict__ IndArr,
//	__global double *__restrict__ iAlphat,
//	double alpha, int turbFlag)
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
//	double dX_e = dXcur[gid], dX_w = dXcur[gid + DIST_SIZE], dX_c = dXcur[gid + DIST_SIZE * 2];
//	double dY_n = dXcur[gid + DIST_SIZE * 3], dY_s = dXcur[gid + DIST_SIZE * 4], dY_c = dXcur[gid + DIST_SIZE * 5];
//	double dX2_e = dXcur[gid + DIST_SIZE * 6], dX2_w = dXcur[gid + DIST_SIZE * 7], dX2_c = dXcur[gid + DIST_SIZE * 8];
//	double dY2_n = dXcur[gid + DIST_SIZE * 9], dY2_s = dXcur[gid + DIST_SIZE * 10], dY2_c = dXcur[gid + DIST_SIZE * 11];
//
//	double Ux = ivx[gid], Uy = ivy[gid], alphat = iAlphat[gid];
//	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
//	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);
//
//	double dAdy = 0., dAdx = 0.;
//
//	if (turbFlag)
//	{
//		dAdx = dX_e*iAlphat[gide] + dX_w*iAlphat[gidw] + dX_c*alphat;
//		dAdy = dY_n*iAlphat[gid + CHANNEL_LENGTH_FULL] +
//			dY_s*iAlphat[gide - CHANNEL_LENGTH_FULL] + dY_c*alphat;
//	}
//	double Jx = (Ux - dAdx);
//	double Jy = (Uy - dAdy);
//
//	double Tc = (Jx*dX_c + Jy*dY_c - alphat * (dX2_c + dY2_c));
//	double Te = (Jx*dX_e - alphat*dX2_e);
//	double Tw = (Jx*dX_w - alphat*dX2_w);
//	double Tn = (Jy*dY_n - alphat*dY2_n);
//	double Ts = (Jy*dY_s - alphat*dY2_s);
//	double Tsrc = 0.;
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
//}


//////////////////////////////////////////////////////////////////////
/////////////////                                    /////////////////
/////////////////    Functions for updating Alpha    /////////////////
/////////////////    Rho and Cp derivative arrays    /////////////////
/////////////////                                    /////////////////
//////////////////////////////////////////////////////////////////////

// del initially = dXcurrent, divided by dx0 to get dXcur/dX0
void calcAlphaRhoCpBtwNodes(double del, const double dx0,
	const NTYPE_TYPE ntype1,
	const NTYPE_TYPE ntype2,
	double* alphaVal, double* rhoCpVal)
{
	double k_c = (ntype1 & M_FLUID_NODE) ? K_AIR_LB : K_SOOT_LB;
	double cp_c = (ntype1 & M_FLUID_NODE) ? CP_AIR_LB : CP_SOOT_LB;
	double ro_c = (ntype1 & M_FLUID_NODE) ? RHO_AIR_LB : RHO_SOOT_LB;

	double k_n = (ntype2 & M_FLUID_NODE) ? K_AIR_LB : K_SOOT_LB;
	double cp_n = (ntype2 & M_FLUID_NODE) ? CP_AIR_LB : CP_SOOT_LB;
	double ro_n = (ntype2 & M_FLUID_NODE) ? RHO_AIR_LB : RHO_SOOT_LB;

	// if very thin region of fouling layer, just use values based on node type
	if (fabs(del - dx0) < CEPS)
	{
		*rhoCpVal = ro_c * cp_c;
		*alphaVal = k_c / *rhoCpVal;
		return;
	}

	del /= dx0;

	double coeffA = k_n * ro_c * cp_c;
	double coeffB = k_n * ro_c * cp_n + k_n * ro_n * cp_c + k_c * ro_c * cp_c;
	double coeffC = k_n * ro_n * cp_n + k_c * ro_c * cp_n + k_c * ro_n * cp_c;
	double coeffD = k_c * ro_n * cp_n;

	// calculate alpha based on harmonic mean interpolation
	*alphaVal = k_c * k_n / (coeffA * del * del * del + coeffB * del * del * (1.0 - del) +
		coeffC * del * (1.0 - del) * (1.0 - del) + coeffD * (1.0 - del) * (1.0 - del) * (1.0 - del));

	// calculate alpha based on linear interp of rho*cp
	*rhoCpVal = ro_c * cp_c * del + ro_n * cp_n * (1. - del);
}


// del initially = dXcurrent, divided by dx0 to get dXcur/dX0
void calcAlphaBtwNodes(double del, const double dx0,
	const NTYPE_TYPE ntype1,
	const NTYPE_TYPE ntype2,
	double* alphaVal)
{
	double k_c = (ntype1 & M_FLUID_NODE) ? K_AIR_LB : K_SOOT_LB;
	double cp_c = (ntype1 & M_FLUID_NODE) ? CP_AIR_LB : CP_SOOT_LB;
	double ro_c = (ntype1 & M_FLUID_NODE) ? RHO_AIR_LB : RHO_SOOT_LB;

	double k_n = (ntype2 & M_FLUID_NODE) ? K_AIR_LB : K_SOOT_LB;
	double cp_n = (ntype2 & M_FLUID_NODE) ? CP_AIR_LB : CP_SOOT_LB;
	double ro_n = (ntype2 & M_FLUID_NODE) ? RHO_AIR_LB : RHO_SOOT_LB;

	if (fabs(del - dx0) < CEPS)
	{
		*alphaVal = k_c / ro_c * cp_c;
		return;
	}

	del /= dx0;
	
	double coeffA = k_n * ro_c * cp_c;
	double coeffB = k_n * ro_c * cp_n + k_n * ro_n * cp_c + k_c * ro_c * cp_c;
	double coeffC = k_n * ro_n * cp_n + k_c * ro_c * cp_n + k_c * ro_n * cp_c;
	double coeffD = k_c * ro_n * cp_n;

	// calculate alpha based on harmonic mean interpolation
	*alphaVal = k_c * k_n / (coeffA * del * del * del + coeffB * del * del * (1.0 - del) +
		coeffC * del * (1.0 - del) * (1.0 - del) + coeffD * (1.0 - del) * (1.0 - del) * (1.0 - del));

}


// Updates the derivative arrays for Alpha and rhoCp
__kernel void updateDerivativeArrays(__global NTYPE_TYPE* __restrict__ Map,
	__global double* __restrict__ dXArr,
	__global double* __restrict__ dX0Arr,
#ifdef USING_CHT_SOURCE_CORRECTION
	__global double* __restrict__ rhoCpDer,
#endif
	__global double* __restrict__ alphaDer)
{

	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH)
		return;
	
	int gid = GET_GLOBAL_IDX(i, j);

	NTYPE_TYPE ntype = Map[gid];
	
	if (ntype & M0_SOLID_NODE)
		return;

	int gide = ((i < CHANNEL_LENGTH - 1) ? 1 : 1 - CHANNEL_LENGTH) + gid;
	int gidw = ((i > 0) ? -1 : (CHANNEL_LENGTH - 1)) + gid;

	double dxe = dX0Arr[gid], dxw = dX0Arr[gid + DIST_SIZE];
	double dyn = dX0Arr[gid + DIST_SIZE * 2], dys = dX0Arr[gid + DIST_SIZE * 3];

	double alpha_e, alpha_w, alpha_n, alpha_s;

#ifdef USING_CHT_SOURCE_CORRECTION
	double rhoCp_e, rhoCp_w, rhoCp_n, rhoCp_s;

	calcAlphaRhoCpBtwNodes(dXArr[gid], dxe, ntype, Map[gide], &alpha_e, &rhoCp_e);
	
	calcAlphaRhoCpBtwNodes(dXArr[gid + DIST_SIZE], dxw, ntype, Map[gidw], &alpha_w, &rhoCp_w);
	
	calcAlphaRhoCpBtwNodes(dXArr[gid + DIST_SIZE * 2], dyn, ntype,
		Map[gid + CHANNEL_LENGTH_FULL], &alpha_n, &rhoCp_n);
	
	calcAlphaRhoCpBtwNodes(dXArr[gid + DIST_SIZE * 3], dys, ntype,
		Map[gid - CHANNEL_LENGTH_FULL], &alpha_s, &rhoCp_s);

	double roCp_c = 0.25 * (rhoCp_e + rhoCp_w + rhoCp_n + rhoCp_s);

	rhoCpDer[gid*2] = (2. * dxw * rhoCp_e / (dxe * (dxe + dxw)) - 2. * dxe * rhoCp_w / (dxw * (dxe + dxw))) / roCp_c +
		2. * (dxe - dxw) / (dxe * dxw);

	rhoCpDer[gid*2+1] = (2. * dys * rhoCp_n / (dyn * (dyn + dys)) - 2. * dyn * rhoCp_s / (dys * (dyn + dys))) / roCp_c +
		2. * (dyn - dys) / (dyn * dys);

#endif
	calcAlphaBtwNodes(dXArr[gid], dxe, ntype, Map[gide], &alpha_e);

	calcAlphaBtwNodes(dXArr[gid + DIST_SIZE], dxw, ntype, Map[gidw], &alpha_w);

	calcAlphaBtwNodes(dXArr[gid + DIST_SIZE * 2], dyn, ntype,
		Map[gid + CHANNEL_LENGTH_FULL], &alpha_n);

	calcAlphaBtwNodes(dXArr[gid + DIST_SIZE * 3], dys, ntype,
		Map[gid - CHANNEL_LENGTH_FULL], &alpha_s);

	double alpha_c = 0.25 * (alpha_e + alpha_w + alpha_n + alpha_s);

	alphaDer[gid*2] = 2. * dxw * alpha_e / (dxe * (dxe + dxw)) - 2. * dxe * alpha_w / (dxw * (dxe + dxw)) +
		2. * (dxe - dxw) * alpha_c / (dxe * dxw);

	alphaDer[gid*2+1] = 2. * dys * alpha_n / (dyn * (dyn + dys)) - 2. * dyn * alpha_s / (dys * (dyn + dys)) +
		2. * (dyn - dys) * alpha_c / (dyn * dys);

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
	double2 vNm = (double2)(-vP10.y, vP10.x);
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
	__global NTYPE_TYPE* __restrict__ nType,
	__global int* nuYInds,
	__global double2* __restrict__ nuDistVec,
	__global double* __restrict__ nuDist,
	const ushort2 BLind, const int shiftInd)
{
	double2 C0 = Cvals[BLind.x], C1 = Cvals[BLind.y];
	int2 C0i = convert_int2_rtz(fmin(C0, C1));
	int2 C1i = convert_int2_rtp(fmax(C0, C1));
	C0i -= 1;
	C1i += 1;

	C0i = max(C0i, (int2)(TR_X_IND_START, 0));
	C1i = min(C1i, (int2)(TR_X_IND_STOP - 1, CHANNEL_HEIGHT - 1));
	if ((C0i.x >= TR_X_IND_STOP) || C1i.x < TR_X_IND_START)
		return;

	double dist;
	double2 vCcut;

	for (int ii = C0i.x; ii <= C1i.x; ii++)
	{
		for (int jj = C0i.y; jj <= C1i.y; jj++)
		{
			if (nType[GET_GLOBAL_IDX(ii,jj)] & M_SOLID_NODE)
				continue;

			double2 L0 = (double2)((double)ii, (double)jj);
			if(nuFindIntersectionNormal(&vCcut, &dist, L0, C0, C1))
			{
				if(AtomicMin(&nuDist[ii - shiftInd], dist))
				{
					//nuDist[ii - shiftInd] = dist; // this is for testing and should not stay
					nuDistVec[ii - shiftInd] = L0 - vCcut;
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
	__global NTYPE_TYPE* __restrict__ nType,
	__global int* __restrict__ nuYInds,
	__global double2* __restrict__ nuDistVec,
	__global double* __restrict__ nuDist)
{

	int i = get_global_id(0);

	if (i >= NU_NUM_NODES)
		return;

	nuFindNearestNode(Cvals, nType, nuYInds, nuDistVec,
		nuDist, BLinds[i + MIN_BL_BOT], TR_X_IND_START);

	nuFindNearestNode(Cvals, nType, nuYInds, nuDistVec,
		nuDist, BLinds[i + MIN_BL_TOP], TR_X_IND_START - NU_NUM_NODES);

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
	nuYInds[i] = (nuYInds[i] << 2) | nodeType;
}


// Calculates bulk mean temperature
double getTbm(__global double* __restrict__ T_array,
	__global double* __restrict__ Ux_array,
	__global double* __restrict__ Uy_array,
	int ival)
{
	double Um = 0.;
	double Tm = 0.;
	for (int j = 0; j < FULLSIZEY; j++)
	{
		int linInd = ival + j * CHANNEL_HEIGHT;
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

	int neighInd = (nType & 0x1) ? 1 : -1;
	int dirInd = (nType & 0x1) ? 0 : DIST_SIZE;
	dirInd += linInd;

	double Tdiff = Tnode - Tvals[linInd + neighInd];
	double dx = dXvals[dirInd];
	double dx0 = dX0vals[dirInd];


	double2 dTdx;
	dTdx.x = KSOOT * Tdiff / (KAIR * (dx0 - dx) + KSOOT * dx);

	neighInd = (nType & 0x2) ? CHANNEL_LENGTH_FULL : -CHANNEL_LENGTH_FULL;
	dirInd = (nType & 0x2) ? (DIST_SIZE * 2) : (DIST_SIZE * 3);
	dirInd += linInd;


	Tdiff = Tnode - Tvals[linInd + neighInd];
	dx = dXvals[dirInd];
	dx0 = dX0vals[dirInd];


	dTdx.y = KSOOT * Tdiff / (KAIR * (dx0 - dx) + KSOOT * dx);

	*dTval = dot(dTdx, normalize(fabs(vNorm)));

	Tnode -= (*dTval) * length(vNorm);


	return Tnode;
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_NU, 1, 1)))
void FD_calc_Nu(__global double* __restrict__ T_array,
	__global double* __restrict__ Nu_array,
	__global double* __restrict__ Ux_array,
	__global double* __restrict__ Uy_array,
	__global int* __restrict__ nuYInds,
	__global double2* __restrict__ nuDist,
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

	double dTdn[2] = { 0. };
	double2 Ts = 0.;

	timeIndex *= NU_NUM_NODES_FULLSIZE;
	Nu_array[timeIndex + ix * 5] = Tm;


	if (jj != -1)
	{
		int linInd = ii + jj * CHANNEL_LENGTH_FULL;
		double2 vNorm = nuDist[ix];

		Ts.x = getTandDT(T_array, dXvals, dX0vals,
			vNorm, linInd, nType, &dTdn[0]);
	}

	Nu_array[timeIndex + ix * 5 + 1] = (jj > -1) ? dTdn[0] : -1000.;
	Nu_array[timeIndex + ix * 5 + 2] = (jj > -1) ? Ts.x : -1000.;

	jj = nuYInds[ix + NU_NUM_NODES];

	nType = (jj & 0x3);
	jj = jj >> 2;

	if (jj != -1)
	{
		int linInd = ii + jj * CHANNEL_LENGTH_FULL;
		double2 vNorm = nuDist[ix + NU_NUM_NODES];

		Ts.y = getTandDT(T_array, dXvals, dX0vals,
			vNorm, linInd, nType, &dTdn[1]);
	}

	Nu_array[timeIndex + ix * 5 + 3] = (jj > -1) ? dTdn[1] : -1000.;
	Nu_array[timeIndex + ix * 5 + 4] = (jj > -1) ? Ts.y : -1000.;

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
		//Nuvals[0] = UTdata[0];
		//Nuvals[1] = Udata[0];
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


