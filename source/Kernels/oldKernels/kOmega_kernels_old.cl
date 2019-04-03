#define DBGARR(dbgnum, dbgval)	DebugArrs[gid + dbgnum*DIST_SIZE] = dbgval



__kernel void Update_kOmega_Diffusivities(__global int *__restrict__ map,
	__global double *__restrict__ iK,
	__global double *__restrict__ iO,
	__global double *__restrict__ iNut,
	__global double *__restrict__ WallD,
	__global double *__restrict__ iDk,
	__global double *__restrict__ iDo,
	__global double *__restrict__ iF1,
	__global double *__restrict__ idkdo,
	__global double *__restrict__ irho,
#ifndef DEBUG_ARRAYS
	__global double *__restrict__ dXcur)
#else
	__global double *__restrict__ dXcur,
	__global double *__restrict__ DebugArrs)
#endif
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i >= CHANNEL_LENGTH || j >= CHANNEL_HEIGHT)
		return;

	int gid = i + CHANNEL_LENGTH_FULL*j;
	const int type = map[gid];
	if ((type & FLUID_NODE) == 0) { return; }





	iF1[gid] = 0;

	iDk[gid] = ALPHA_FLUID;
	iDo[gid] = ALPHA_FLUID;


	return;
}
//
//
//
//
//
//
//
//
//
//
//
//
//
//	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx = dx_e + dx_w;
//	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy = dy_n + dy_s;
//
//	double Xe_coeff = dx_w / (dx_e*dx), Xw_coeff = -dx_e / (dx_w*dx), Xc_coeff = (dx_e - dx_w) / (dx_e*dx_w);
//	double Yn_coeff = dy_s / (dy_n*dy), Ys_coeff = -dy_n / (dy_s*dy), Yc_coeff = (dy_n - dy_s) / (dy_n*dy_s);
//
//	if (dx_w < 1.)
//	{
//		Xe_coeff = 1. / dx_e;
//		Xc_coeff = -Xe_coeff;
//		Xw_coeff = 0.;
//	}
//	else if (dx_e < 1.)
//	{
//		Xc_coeff = 1. / dx_w;
//		Xw_coeff = -Xc_coeff;
//		Xe_coeff = 0.;
//	}
//
//	if (dy_s < 1.)
//	{
//		Yn_coeff = 1. / dy_n;
//		Yc_coeff = -Yn_coeff;
//		Ys_coeff = 0.;
//	}
//	else if (dy_n < 1.)
//	{
//		Yc_coeff = 1. / dy_s;
//		Ys_coeff = -Yc_coeff;
//		Yn_coeff = 0.;
//	}
//
//#ifdef DEBUG_ARRAYS
//	DBGARR(0, Xe_coeff);
//	DBGARR(1, Xw_coeff);
//	DBGARR(2, Xc_coeff);
//	DBGARR(3, Yn_coeff);
//	DBGARR(4, Ys_coeff);
//	DBGARR(5, Yc_coeff);
//#endif
//
//
//	double omegaval = iO[gid], kval = iK[gid], yval = WallD[gid];
//	
//	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
//	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);
//	
//	double dOdx = Xe_coeff*iO[gide] + Xw_coeff*iO[gidw] + Xc_coeff*omegaval;
//	double dKdx = Xe_coeff*iK[gide] + Xw_coeff*iK[gidw] + Xc_coeff*kval;
//
//	double dOdy = Yn_coeff*iO[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iO[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*omegaval;
//	double dKdy = Yn_coeff*iK[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iK[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*kval;
//
//	double dko_mult = 2.*KO_SIGMA_O2 / omegaval * (dOdx*dKdx + dOdy*dKdy);
//	double CDkw = max(dko_mult * irho[gid], 1.0e-10);
//	idkdo[gid] = dko_mult;
//
//#ifdef DEBUG_ARRAYS
//	DBGARR(6, dKdx);
//	DBGARR(7, dKdy);
//	DBGARR(8, dOdx);
//	DBGARR(9, dOdy);
//	DBGARR(10, CDkw);
//#endif
//	
//	double F1val = min(max(sqrt(kval) / (KO_BETA_STAR*omegaval*yval), 500. * MU_NUMBER / (yval*yval*omegaval)), 4.*kval*KO_SIGMA_O2 / (CDkw*yval*yval));
//	F1val = tanh(pown(F1val, 4));
//	iF1[gid] = F1val;
//
//	double sigma = KO_SIGMA_K1*F1val + KO_SIGMA_K2*(1. - F1val);
//	iDk[gid] = MU_NUMBER + sigma*iNut[gid];
//
//	sigma = KO_SIGMA_O1*F1val + KO_SIGMA_O2*(1. - F1val);
//	iDo[gid] = MU_NUMBER + sigma*iNut[gid];
//}

#undef DBGARR



double get_omega_wall(double WallD)
{
	double Utau = 180.*MU_NUMBER / 100.;
	double yplus = WallD * Utau / MU_NUMBER;
	double omega_vis = 6.*MU_NUMBER / (0.075*yplus*yplus), omega_log = Utau / (0.3*yplus*KO_VON_KARMAN);
	return sqrt(omega_vis*omega_vis + omega_log*omega_log);


}



__kernel void Update_kOmega_Coeffs_Explicit(__global int *__restrict__ IndArr,
	__global double *__restrict__ iK,
	__global double *__restrict__ iO,
	__global double *__restrict__ iNut,
	__global double *__restrict__ iDk,
	__global double *__restrict__ iDo,
	__global double *__restrict__ iF1,
	__global double *__restrict__ idkdo,
	__global double *__restrict__ iSxy,
	__global double *__restrict__ irho,
	__global double *__restrict__ ivx,
	__global double *__restrict__ ivy,
	__global double *__restrict__ kAmat,
	__global double *__restrict__ oAmat,
	__global double *__restrict__ kbvec,
	__global double *__restrict__ obvec,
	__global double *__restrict__ dXcur,
	__global double *__restrict__ WallD,
#ifndef DEBUG_ARRAYS
	double DTFD, double TurbInts, double TurbLScale_Inv)
#else
	double DTFD, double TurbInts, double TurbLScale_Inv,
	__global double *__restrict__ DebugArrs)
#endif
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


	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx = dx_e + dx_w;
	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy = dy_n + dy_s;

	double Xe_coeff = dx_w / (dx_e*dx), Xw_coeff = -dx_e / (dx_w*dx), Xc_coeff = (dx_e - dx_w) / (dx_e*dx_w);
	double Yn_coeff = dy_s / (dy_n*dy), Ys_coeff = -dy_n / (dy_s*dy), Yc_coeff = (dy_n - dy_s) / (dy_n*dy_s);

	if (dx_w < 1.)
	{
		Xe_coeff = 1. / dx_e;
		Xc_coeff = -Xe_coeff;
		Xw_coeff = 0.;
	}
	else if (dx_e < 1.)
	{
		Xc_coeff = 1. / dx_w;
		Xw_coeff = -Xc_coeff;
		Xe_coeff = 0.;
	}

	if (dy_s < 1.)
	{
		Yn_coeff = 1. / dy_n;
		Yc_coeff = -Yn_coeff;
		Ys_coeff = 0.;
	}
	else if (dy_n < 1.)
	{
		Yc_coeff = 1. / dy_s;
		Ys_coeff = -Yc_coeff;
		Yn_coeff = 0.;
	}

	double Xe2_coeff = 2. / (dx_e*dx), Xw2_coeff = 2. / (dx_w*dx), Xc2_coeff = -2. / (dx_e*dx_w);
	double Yn2_coeff = 2. / (dy_n*dy), Ys2_coeff = 2. / (dy_s*dy), Yc2_coeff = -2. / (dy_n*dy_s);

	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double Ux = ivx[gid], Uy = ivy[gid];

	double diffk_w = iDk[gidw], diffo_w = iDo[gidw]; ///will never be out of bounds because node (0,0) is always solid
	double diffk_e = iDk[gide], diffo_e = iDo[gide];
	double diffk_n = iDk[gid + CHANNEL_LENGTH_FULL], diffo_n = iDo[gid + CHANNEL_LENGTH_FULL];
	double diffk_s = iDk[gid - CHANNEL_LENGTH_FULL], diffo_s = iDo[gid - CHANNEL_LENGTH_FULL];

	double diffk = iDk[gid], diffo = iDo[gid];
	double diffk_dx, diffk_dy, diffo_dx, diffo_dy;

	diffk_dx = Xe_coeff * diffk_e + Xw_coeff * diffk_w + Xc_coeff * diffk;
	diffo_dx = Xe_coeff * diffo_e + Xw_coeff * diffo_w + Xc_coeff * diffo;
	diffk_dy = Yn_coeff * diffk_n + Ys_coeff * diffk_s + Yc_coeff * diffk;
	diffo_dy = Yn_coeff * diffo_n + Ys_coeff * diffo_s + Yc_coeff * diffo;


	double Kc = -DTFD*(Ux*Xc_coeff + Uy*Yc_coeff - ((Xc2_coeff + Yc2_coeff)*diffk + Xc_coeff*diffk_dx + Yc_coeff*diffk_dy)) + 1.0;
	double Ke = -DTFD*(Ux*Xe_coeff - Xe2_coeff*diffk - Xe_coeff*diffk_dx);
	double Kw = -DTFD*(Ux*Xw_coeff - Xw2_coeff*diffk - Xw_coeff*diffk_dx);
	double Kn = -DTFD*(Uy*Yn_coeff - Yn2_coeff*diffk - Yn_coeff*diffk_dy);
	double Ks = -DTFD*(Uy*Ys_coeff - Ys2_coeff*diffk - Ys_coeff*diffk_dy);

	double Wc = -DTFD*(Ux*Xc_coeff + Uy*Yc_coeff - ((Xc2_coeff + Yc2_coeff)*diffo + Xc_coeff*diffo_dx + Yc_coeff*diffo_dy)) + 1.0;
	double We = -DTFD*(Ux*Xe_coeff - Xe2_coeff*diffo - Xe_coeff*diffo_dx);
	double Ww = -DTFD*(Ux*Xw_coeff - Xw2_coeff*diffo - Xw_coeff*diffo_dx);
	double Wn = -DTFD*(Uy*Yn_coeff - Yn2_coeff*diffo - Yn_coeff*diffo_dy);
	double Ws = -DTFD*(Uy*Ys_coeff - Ys2_coeff*diffo - Ys_coeff*diffo_dy);

	double Ksrc = 0.;
	double Wsrc = 0.;

	if (Sind >= 0)
	{
		kAmat[Sind] = Ks;
		oAmat[Sind] = Ws;
	}

	if (Wind >= 0)
	{
		kAmat[Wind] = Kw;
		oAmat[Wind] = Ww;
	}

	kAmat[Cind] = Kc;
	oAmat[Cind] = Wc;

	if (Eind >= 0)
	{
		kAmat[Eind] = Ke;
		oAmat[Eind] = We;
	}

	if (Nind >= 0)
	{
		kAmat[Nind] = Kn;
		oAmat[Nind] = Wn;
	}


	kbvec[gid] = Ksrc;
	obvec[gid] = Wsrc;

	if (Nind < 0 || Sind < 0)
	{
		kAmat[Sind] = 0.;
		oAmat[Sind] = 0.;
		kAmat[Wind] = 0.;
		oAmat[Wind] = 0.;
		kAmat[Cind] = 0.;
		oAmat[Cind] = 0.;
		kAmat[Eind] = 0.;
		oAmat[Eind] = 0.;
		kAmat[Nind] = 0.;
		oAmat[Nind] = 0.;
		kbvec[gid] = 1.;
		obvec[gid] = 1.;
	}


}




//__kernel void Update_kOmega_Coeffs_Explicit(__global int *__restrict__ IndArr,
//	__global double *__restrict__ iK,
//	__global double *__restrict__ iO,
//	__global double *__restrict__ iNut,
//	__global double *__restrict__ iDk,
//	__global double *__restrict__ iDo,
//	__global double *__restrict__ iF1,
//	__global double *__restrict__ idkdo,
//	__global double *__restrict__ iSxy,
//	__global double *__restrict__ irho,
//	__global double *__restrict__ ivx,
//	__global double *__restrict__ ivy,
//	__global double *__restrict__ kAmat,
//	__global double *__restrict__ oAmat,
//	__global double *__restrict__ kbvec,
//	__global double *__restrict__ obvec,
//	__global double *__restrict__ dXcur, 
//	__global double *__restrict__ WallD,
//#ifndef DEBUG_ARRAYS
//	double DTFD, double TurbInts, double TurbLScale_Inv)
//#else
//	double DTFD, double TurbInts, double TurbLScale_Inv,
//	__global double *__restrict__ DebugArrs)
//#endif
//{
//
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
//
//	int Eind = IndArr[gid + DIST_SIZE], Wind = IndArr[gid + DIST_SIZE * 2];
//	int Nind = IndArr[gid + DIST_SIZE * 3], Sind = IndArr[gid + DIST_SIZE * 4];
//
//
//	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx = dx_e + dx_w;
//	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy = dy_n + dy_s;
//
//	double Xe_coeff = dx_w / (dx_e*dx), Xw_coeff = -dx_e / (dx_w*dx), Xc_coeff = (dx_e - dx_w) / (dx_e*dx_w);
//	double Yn_coeff = dy_s / (dy_n*dy), Ys_coeff = -dy_n / (dy_s*dy), Yc_coeff = (dy_n - dy_s) / (dy_n*dy_s);
//
//	if (dx_w < 1.)
//	{
//		Xe_coeff = 1. / dx_e;
//		Xc_coeff = -Xe_coeff;
//		Xw_coeff = 0.;
//	}
//	else if (dx_e < 1.)
//	{
//		Xc_coeff = 1. / dx_w;
//		Xw_coeff = -Xc_coeff;
//		Xe_coeff = 0.;
//	}
//
//	if (dy_s < 1.)
//	{
//		Yn_coeff = 1. / dy_n;
//		Yc_coeff = -Yn_coeff;
//		Ys_coeff = 0.;
//	}
//	else if (dy_n < 1.)
//	{
//		Yc_coeff = 1. / dy_s;
//		Ys_coeff = -Yc_coeff;
//		Yn_coeff = 0.;
//	}
//	
//	double Xe2_coeff = 2. / (dx_e*dx), Xw2_coeff = 2. / (dx_w*dx), Xc2_coeff = -2. / (dx_e*dx_w);
//	double Yn2_coeff = 2. / (dy_n*dy), Ys2_coeff = 2. / (dy_s*dy), Yc2_coeff = -2. / (dy_n*dy_s);
//
//#ifdef DEBUG_ARRAYS
//	DebugArrs[gid] = Xe_coeff;
//	DebugArrs[gid + DIST_SIZE] = Xw_coeff;
//	DebugArrs[gid + 2*DIST_SIZE] = Xc_coeff;
//	DebugArrs[gid + 3*DIST_SIZE] = Yn_coeff;
//	DebugArrs[gid + 4*DIST_SIZE] = Ys_coeff;
//	DebugArrs[gid + 5*DIST_SIZE] = Yc_coeff;
//	DebugArrs[gid + 6*DIST_SIZE] = Xe2_coeff;
//	DebugArrs[gid + 7*DIST_SIZE] = Xw2_coeff;
//	DebugArrs[gid + 8*DIST_SIZE] = Xc2_coeff;
//	DebugArrs[gid + 9*DIST_SIZE] = Yn2_coeff;
//	DebugArrs[gid + 10*DIST_SIZE] = Ys2_coeff;
//	DebugArrs[gid + 11*DIST_SIZE] = Yc2_coeff;
//#endif
//
//	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
//	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);
//
//	double ux_w = ivx[gidw], uy_w = ivy[gidw]; ///will never be out of bounds because node (0,0) is always solid
//	double ux_e = ivx[gide], uy_e = ivy[gide];
//	double ux_n = ivx[gid + CHANNEL_LENGTH_FULL], uy_n = ivy[gid + CHANNEL_LENGTH_FULL];
//	double ux_s = ivx[gid - CHANNEL_LENGTH_FULL], uy_s = ivy[gid - CHANNEL_LENGTH_FULL];
//	double Ux = ivx[gid], Uy = ivy[gid];
//
//	double dudx = ux_e * Xe_coeff + ux_w * Xw_coeff + Ux * Xc_coeff;
//	double dvdx = uy_e * Xe_coeff + uy_w * Xw_coeff + Uy * Xc_coeff;
//	double dudy = ux_n * Yn_coeff + ux_s * Ys_coeff + Ux * Yc_coeff;
//	double dvdy = uy_n * Yn_coeff + uy_s * Ys_coeff + Uy * Yc_coeff;
//
//#ifdef DEBUG_ARRAYS
//	DebugArrs[gid + 12 * DIST_SIZE] = dudx;
//	DebugArrs[gid + 13 * DIST_SIZE] = dudy;
//	DebugArrs[gid + 14 * DIST_SIZE] = dvdx;
//	DebugArrs[gid + 15 * DIST_SIZE] = dvdy;
//#endif
//
//
//	double diffk_w = iDk[gidw], diffo_w = iDo[gidw]; ///will never be out of bounds because node (0,0) is always solid
//	double diffk_e = iDk[gide], diffo_e = iDo[gide];
//	double diffk_n = iDk[gid + CHANNEL_LENGTH_FULL], diffo_n = iDo[gid + CHANNEL_LENGTH_FULL];
//	double diffk_s = iDk[gid - CHANNEL_LENGTH_FULL], diffo_s = iDo[gid - CHANNEL_LENGTH_FULL];
//	
//	double diffk = iDk[gid], diffo = iDo[gid];
//	double diffk_dx, diffk_dy, diffo_dx, diffo_dy;
//
//	diffk_dx = Xe_coeff * diffk_e + Xw_coeff * diffk_w + Xc_coeff * diffk;
//	diffo_dx = Xe_coeff * diffo_e + Xw_coeff * diffo_w + Xc_coeff * diffo;
//	diffk_dy = Yn_coeff * diffk_n + Ys_coeff * diffk_s + Yc_coeff * diffk;
//	diffo_dy = Yn_coeff * diffo_n + Ys_coeff * diffo_s + Yc_coeff * diffo;
//
//
//#ifdef DEBUG_ARRAYS
//	DebugArrs[gid + 16 * DIST_SIZE] = diffk_dx;
//	DebugArrs[gid + 17 * DIST_SIZE] = diffk_dy;
//	DebugArrs[gid + 18 * DIST_SIZE] = diffo_dx;
//	DebugArrs[gid + 19 * DIST_SIZE] = diffo_dy;
//#endif
//
//	double omegaval = iO[gid], kval = iK[gid], rhoval = irho[gid];
//	double nutval = iNut[gid];
//
//	
//	// Used to calculate S^2 before being converted to tau_ij
//	double tauxx = iSxy[gid], tauxy = iSxy[gid + DIST_SIZE], tauyy = iSxy[gid + 2 * DIST_SIZE];
//	double S2 = tauxx*tauxx + tauyy*tauyy + 2.*tauxy*tauxy;
//	tauxx = 2.*rhoval * (nutval*tauxx - 1. / 3. * kval);// may be the slightly different version in NASA paper
//	tauxy = 2.*rhoval * (nutval*tauxy);
//	tauyy = 2.*rhoval * (nutval*tauyy - 1. / 3. * kval);
//
//	double Pk = min(tauxx*dudx + tauxy*(dudy + dvdx) + tauyy*dvdy, 10.*KO_BETA_STAR*kval*omegaval);
//
//	double F1val = iF1[gid];
//	double beta_val = KO_BETA1*F1val + KO_BETA2*(1. - F1val);
//	double alpha_val = KO_ALPHA1*F1val + KO_ALPHA2*(1. - F1val);
//
//
//#ifdef DEBUG_ARRAYS
//	DebugArrs[gid + 20 * DIST_SIZE] = S2;
//	DebugArrs[gid + 21 * DIST_SIZE] = tauxx;
//	DebugArrs[gid + 22 * DIST_SIZE] = tauxy;
//	DebugArrs[gid + 23 * DIST_SIZE] = tauyy;
//	DebugArrs[gid + 24 * DIST_SIZE] = Pk;
//#endif
//
//
//	// If using implicit, will need to account for source term linearization.
//
//	double Kc = -DTFD*(Ux*Xc_coeff + Uy*Yc_coeff - ((Xc2_coeff + Yc2_coeff)*diffk + Xc_coeff*diffk_dx + Yc_coeff*diffk_dy)) + 1.0;
//	double Ke = -DTFD*(Ux*Xe_coeff - Xe2_coeff*diffk - Xe_coeff*diffk_dx);
//	double Kw = -DTFD*(Ux*Xw_coeff - Xw2_coeff*diffk - Xw_coeff*diffk_dx);
//	double Kn = -DTFD*(Uy*Yn_coeff - Yn2_coeff*diffk - Yn_coeff*diffk_dy);
//	double Ks = -DTFD*(Uy*Ys_coeff - Ys2_coeff*diffk - Ys_coeff*diffk_dy);
//
//	double Wc = -DTFD*(Ux*Xc_coeff + Uy*Yc_coeff - ((Xc2_coeff + Yc2_coeff)*diffo + Xc_coeff*diffo_dx + Yc_coeff*diffo_dy)) + 1.0;
//	double We = -DTFD*(Ux*Xe_coeff - Xe2_coeff*diffo - Xe_coeff*diffo_dx);
//	double Ww = -DTFD*(Ux*Xw_coeff - Xw2_coeff*diffo - Xw_coeff*diffo_dx);
//	double Wn = -DTFD*(Uy*Yn_coeff - Yn2_coeff*diffo - Yn_coeff*diffo_dy);
//	double Ws = -DTFD*(Uy*Ys_coeff - Ys2_coeff*diffo - Ys_coeff*diffo_dy);
//
//	double Ksrc = Pk - KO_BETA_STAR*omegaval*kval;
//	double Wsrc = (alpha_val)*S2 - beta_val*omegaval*omegaval + idkdo[gid];
//
//	if (Sind >= 0)
//	{
//		kAmat[Sind] = Ks;
//		oAmat[Sind] = Ws;
//	}
//	else
//	{
//		double omegawallval = get_omega_wall(WallD[gid]);
//		Ksrc += Ks*kval;
//		Wsrc += (Ws*omegawallval);
//	}
//
//	if (Wind >= 0)
//	{
//		kAmat[Wind] = Kw;
//		oAmat[Wind] = Ww;
//	}
//	//else
//	//{
//	//	double omegawallval = get_omega_wall(WallD[gid])
//	//	Ksrc += (i == 0) ? (Kw * 1.5*(TurbInts*TurbInts*(Ux*Ux + Uy*Uy))) : (Kw*kval);
//	//	Wsrc += (i == 0) ? (Ww * sqrt(kval)*TurbLScale_Inv / KO_C_MU_0_25) : (Ww*omegaval);
//	//}
//
//	kAmat[Cind] = Kc;
//	oAmat[Cind] = Wc;
//
//	if (Eind >= 0)
//	{
//		kAmat[Eind] = Ke;
//		oAmat[Eind] = We;
//	}
//	else
//	{
//		Ksrc += Ke*kval;
//		Wsrc += We*omegaval;
//	}
//
//	if (Nind >= 0)
//	{
//		kAmat[Nind] = Kn;
//		oAmat[Nind] = Wn;
//	}
//	else
//	{
//		double omegawallval = get_omega_wall(WallD[gid]);
//		Wsrc += (Wn*omegawallval);
//		Ksrc += Kn*kval;
//	}
//
//	kbvec[gid] = Ksrc;
//	obvec[gid] = Wsrc;
//
//	if (Nind < 0 || Sind < 0)
//	{
//		kAmat[Sind] = 0.;
//		oAmat[Sind] = 0.;
//		kAmat[Wind] = 0.;
//		oAmat[Wind] = 0.;
//		kAmat[Cind] = 0.;
//		oAmat[Cind] = 0.;
//		kAmat[Eind] = 0.;
//		oAmat[Eind] = 0.;
//		kAmat[Nind] = 0.;
//		oAmat[Nind] = 0.;
//		kbvec[gid] = K_WALL_VALUE;
//		obvec[gid] = OMEGA_WALL_VALUE;
//	}
//
//#ifdef DEBUG_ARRAYS
//	DebugArrs[gid + 25 * DIST_SIZE] = Kc;
//	DebugArrs[gid + 26 * DIST_SIZE] = Ke;
//	DebugArrs[gid + 27 * DIST_SIZE] = Kw;
//	DebugArrs[gid + 28 * DIST_SIZE] = Kn;
//	DebugArrs[gid + 29 * DIST_SIZE] = Ks;
//	DebugArrs[gid + 30 * DIST_SIZE] = Ksrc;
//	DebugArrs[gid + 31 * DIST_SIZE] = Wc;
//	DebugArrs[gid + 32 * DIST_SIZE] = We;
//	DebugArrs[gid + 33 * DIST_SIZE] = Ww;
//	DebugArrs[gid + 34 * DIST_SIZE] = Wn;
//	DebugArrs[gid + 35 * DIST_SIZE] = Ws;
//	DebugArrs[gid + 36 * DIST_SIZE] = Wsrc;
//#endif
//
//
//
//
//
//
//
//
//}





	



