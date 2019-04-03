

#define DBGARR(dbgnum, dbgval)	DebugArrs[gid + dbgnum*DIST_SIZE] = dbgval



__kernel void Update_kOmega_Diffusivities_Alphat(__global int *__restrict__ map,
	__global double *__restrict__ iK,
	__global double *__restrict__ iO,
	__global double *__restrict__ iNut,
	__global double *__restrict__ WallD,
	__global double *__restrict__ iDk,
	__global double *__restrict__ iDo,
	__global double *__restrict__ iF1,
	__global double *__restrict__ idkdxy,

	__global double *__restrict__ irho,
	__global double *__restrict__ iAlphat,
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
	
	double Xe_coeff = dXcur[gid], Xw_coeff = dXcur[gid + DIST_SIZE], Xc_coeff = dXcur[gid + DIST_SIZE * 2];
	double Yn_coeff = dXcur[gid + DIST_SIZE * 3], Ys_coeff = dXcur[gid + DIST_SIZE * 4], Yc_coeff = dXcur[gid + DIST_SIZE * 5];

	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double omegaval = iO[gid], kval = iK[gid], yval = WallD[gid];
	
	double dOdx = Xe_coeff*iO[gide] + Xw_coeff*iO[gidw] + Xc_coeff*omegaval;
	double dKdx = Xe_coeff*iK[gide] + Xw_coeff*iK[gidw] + Xc_coeff*kval;
	
	double dOdy = Yn_coeff*iO[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iO[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*omegaval;
	double dKdy = Yn_coeff*iK[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iK[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*kval;
	
	double dodx_mult = 2.*KO_SIGMA_O2 / omegaval;
	double dody_mult = dodx_mult*dKdy;
	dodx_mult *= dKdx;

	idkdxy[gid] = dodx_mult;
	idkdxy[gid+DIST_SIZE] = dody_mult;

	double CDkw = max(irho[gid] * (dodx_mult*dOdx + dody_mult*dOdy), 1.0e-10);
	
#ifdef DEBUG_ARRAYS
	DBGARR(0, dKdx);
	DBGARR(1, dKdy);
	DBGARR(2, dOdx);
	DBGARR(3, dOdy);
	DBGARR(4, CDkw);
#endif
	
	double F1val = min(max(sqrt(kval) / (KO_BETA_STAR*omegaval*yval), 500. * MU_NUMBER / (yval*yval*omegaval)), 4.*kval*KO_SIGMA_O2 / (CDkw*yval*yval));
	F1val = tanh(pown(F1val, 4));
	iF1[gid] = F1val;
	double nutval = iNut[gid];
	double sigma = KO_SIGMA_K1*F1val + KO_SIGMA_K2*(1. - F1val);
	iDk[gid] = MU_NUMBER + sigma*nutval;

	sigma = KO_SIGMA_O1*F1val + KO_SIGMA_O2*(1. - F1val);
	iDo[gid] = MU_NUMBER + sigma*nutval;

	iAlphat[gid] = MU_NUMBER / PR_NUMBER + nutval / PR_TURB_NUMBER;
}


__kernel void Update_kOmega_Diffusivities(__global int *__restrict__ map,
	__global double *__restrict__ iK,
	__global double *__restrict__ iO,
	__global double *__restrict__ iNut,
	__global double *__restrict__ WallD,
	__global double *__restrict__ iDk,
	__global double *__restrict__ iDo,
	__global double *__restrict__ iF1,
	__global double *__restrict__ idkdxy,

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

	double Xe_coeff = dXcur[gid], Xw_coeff = dXcur[gid + DIST_SIZE], Xc_coeff = dXcur[gid + DIST_SIZE * 2];
	double Yn_coeff = dXcur[gid + DIST_SIZE * 3], Ys_coeff = dXcur[gid + DIST_SIZE * 4], Yc_coeff = dXcur[gid + DIST_SIZE * 5];

	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double omegaval = iO[gid], kval = iK[gid], yval = WallD[gid];

	double dOdx = Xe_coeff*iO[gide] + Xw_coeff*iO[gidw] + Xc_coeff*omegaval;
	double dKdx = Xe_coeff*iK[gide] + Xw_coeff*iK[gidw] + Xc_coeff*kval;

	double dOdy = Yn_coeff*iO[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iO[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*omegaval;
	double dKdy = Yn_coeff*iK[gid + CHANNEL_LENGTH_FULL] + Ys_coeff*iK[gid - CHANNEL_LENGTH_FULL] + Yc_coeff*kval;

	double dodx_mult = 2.*KO_SIGMA_O2 / omegaval;
	double dody_mult = dodx_mult*dKdy;
	dodx_mult *= dKdx;

	idkdxy[gid] = dodx_mult;
	idkdxy[gid + DIST_SIZE] = dody_mult;

	double CDkw = max(irho[gid] * (dodx_mult*dOdx + dody_mult*dOdy), 1.0e-10);

#ifdef DEBUG_ARRAYS
	DBGARR(0, dKdx);
	DBGARR(1, dKdy);
	DBGARR(2, dOdx);
	DBGARR(3, dOdy);
	DBGARR(4, CDkw);
#endif

	double F1val = min(max(sqrt(kval) / (KO_BETA_STAR*omegaval*yval), 500. * MU_NUMBER / (yval*yval*omegaval)), 4.*kval*KO_SIGMA_O2 / (CDkw*yval*yval));
	F1val = tanh(pown(F1val, 4));
	iF1[gid] = F1val;
	double nutval = iNut[gid];
	double sigma = KO_SIGMA_K1*F1val + KO_SIGMA_K2*(1. - F1val);
	iDk[gid] = MU_NUMBER + sigma*nutval;

	sigma = KO_SIGMA_O1*F1val + KO_SIGMA_O2*(1. - F1val);
	iDo[gid] = MU_NUMBER + sigma*nutval;
}



__kernel void Update_kOmega_Coeffs_Implicit(__global int *__restrict__ IndArr,
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


	double Xe_coeff = dXcur[gid], Xw_coeff = dXcur[gid + DIST_SIZE], Xc_coeff = dXcur[gid + DIST_SIZE * 2];
	double Yn_coeff = dXcur[gid + DIST_SIZE * 3], Ys_coeff = dXcur[gid + DIST_SIZE * 4], Yc_coeff = dXcur[gid + DIST_SIZE * 5];
	double Xe2_coeff = dXcur[gid + DIST_SIZE * 6], Xw2_coeff = dXcur[gid + DIST_SIZE * 7], Xc2_coeff = dXcur[gid + DIST_SIZE * 8];
	double Yn2_coeff = dXcur[gid + DIST_SIZE * 9], Ys2_coeff = dXcur[gid + DIST_SIZE * 10], Yc2_coeff = dXcur[gid + DIST_SIZE * 11];

	int gidw = (i > 0) ? gid - 1 : (CHANNEL_LENGTH - 1 + j*CHANNEL_LENGTH_FULL);
	int gide = (i < CHANNEL_LENGTH - 1) ? gid + 1 : (j*CHANNEL_LENGTH_FULL);

	double diffk_w = iDk[gidw], diffo_w = iDo[gidw]; ///will never be out of bounds because node (0,0) is always solid
	double diffk_e = iDk[gide], diffo_e = iDo[gide];
	double diffk_n = iDk[gid + CHANNEL_LENGTH_FULL], diffo_n = iDo[gid + CHANNEL_LENGTH_FULL];
	double diffk_s = iDk[gid - CHANNEL_LENGTH_FULL], diffo_s = iDo[gid - CHANNEL_LENGTH_FULL];

	double diffk = iDk[gid], diffo = iDo[gid];
	
	double diffk_dx = Xe_coeff * diffk_e + Xw_coeff * diffk_w + Xc_coeff * diffk;
	double diffo_dx = Xe_coeff * diffo_e + Xw_coeff * diffo_w + Xc_coeff * diffo;
	double diffk_dy = Yn_coeff * diffk_n + Ys_coeff * diffk_s + Yc_coeff * diffk;
	double diffo_dy = Yn_coeff * diffo_n + Ys_coeff * diffo_s + Yc_coeff * diffo;



#ifdef DEBUG_ARRAYS
	DBGARR(0, diffk_dx);
	DBGARR(1, diffk_dy);
	DBGARR(2, diffo_dx);
	DBGARR(3, diffo_dy);
#endif




	double Ux = ivx[gid], Uy = ivy[gid];

	double Jx_omega = Ux - diffo_dx - idkdo[gid];
	double Jy_omega = Uy - diffo_dy - idkdo[gid+DIST_SIZE];

	double Jx_k = Ux - diffk_dx;
	double Jy_k = Uy - diffk_dy;

	double omegaval = iO[gid], kval = iK[gid]; 
	double F1val = iF1[gid], nutval = iNut[gid];
	double Sxx = iSxy[gid], Sxy = iSxy[gid + DIST_SIZE], Syy = iSxy[gid + 2 * DIST_SIZE];
	double Pk = 2.*nutval*(Sxx*Sxx + Syy*Syy + Sxy*Sxy) - 2. / 3.*kval * (Sxx + Syy);
	Pk = min(Pk, 10.*KO_BETA_STAR*kval*omegaval);
	double gammaval = (F1val*KO_GAMMA1 + (1. - F1val)*KO_GAMMA2) / nutval;


	double Sc = -KO_BETA_STAR * omegaval;

#ifdef DEBUG_ARRAYS
	DBGARR(4, Jx_omega);
	DBGARR(5, Jy_omega);
	DBGARR(6, Jx_k);
	DBGARR(7, Jy_k);
	DBGARR(8, Pk);
	DBGARR(9, Sc);
#endif
	
	double Kc = 1. + DTFD * (Jx_k*Xc_coeff + Jy_k*Yc_coeff - diffk * (Xc2_coeff + Yc2_coeff)) - Sc;
	double Ke = DTFD * (Jx_k*Xe_coeff - diffk * Xe2_coeff);
	double Kw = DTFD * (Jx_k*Xw_coeff - diffk * Xw2_coeff);
	double Kn = DTFD * (Jy_k*Yn_coeff - diffk * Yn2_coeff);
	double Ks = DTFD * (Jy_k*Ys_coeff - diffk * Ys2_coeff);
	double Ksrc = kval + DTFD*Pk;
	
	double Wc = 1. + DTFD * (Jx_omega*Xc_coeff + Jy_omega*Yc_coeff - diffo * (Xc2_coeff + Yc2_coeff)) - Sc;
	double We = DTFD * (Jx_omega*Xe_coeff - diffo * Xe2_coeff);
	double Ww = DTFD * (Jx_omega*Xw_coeff - diffo * Xw2_coeff);
	double Wn = DTFD * (Jy_omega*Yn_coeff - diffo * Yn2_coeff);
	double Ws = DTFD * (Jy_omega*Ys_coeff - diffo * Ys2_coeff);
	double Wsrc = omegaval + DTFD*(gammaval * Pk);

#ifdef DEBUG_ARRAYS
	DBGARR(10, Kc);
	DBGARR(11, Ke);
	DBGARR(12, Kw);
	DBGARR(13, Kn);
	DBGARR(14, Ks);
	DBGARR(15, Ksrc);
	DBGARR(16, Wc);
	DBGARR(17, We);
	DBGARR(18, Ww);
	DBGARR(19, Wn);
	DBGARR(20, Ws);
	DBGARR(21, Wsrc);
#endif



	if (Sind >= 0 && Wind >= 0 && Eind >= 0 && Nind >= 0)
	{
		kAmat[Sind] = Ks;
		oAmat[Sind] = Ws;
		kAmat[Wind] = Kw;
		oAmat[Wind] = Ww;
		kAmat[Cind] = Kc;
		oAmat[Cind] = Wc;
		kAmat[Eind] = Ke;
		oAmat[Eind] = We;
		kAmat[Nind] = Kn;
		oAmat[Nind] = Wn;

		kbvec[gid] = Ksrc;
		obvec[gid] = Wsrc;
	}
	else
	{
		//if (Sind < 0)
		//{
		//	kAmat[Cind] = Kc + Ks;
		//	kAmat[Nind] = Kn;
		//	oAmat[Nind] = 0.;
		//}
		//else
		//{
		//	oAmat[Sind] = 0.;
		//	kAmat[Sind] = Ks;
		//	kAmat[Cind] = Kc + Kn;
		//}

		//
		//kAmat[Wind] = Kw;
		//oAmat[Wind] = 0.;
		//oAmat[Cind] = 1.;
		//kAmat[Eind] = Ke;
		//oAmat[Eind] = 0.;

		//if (Sind < 0)
		//{
		//	kAmat[Nind] = 1e-3;
		//	oAmat[Nind] = 1e-3;
		//}
		//else
		//{
		//	kAmat[Sind] = 1e-3;
		//	oAmat[Sind] = 1e-3;
		//}

		//kAmat[Wind] = 1e-3;
		//kAmat[Eind] = 1e-3;

		kAmat[Cind] = 1.;
		oAmat[Cind] = 1.;
		
		kbvec[gid] = K_WALL_VALUE;
		obvec[gid] = OMEGA_WALL_VALUE;
	}


}

#undef DBGARR