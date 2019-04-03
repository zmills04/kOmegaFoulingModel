__kernel void Update_T_Coeffs(__global double *__restrict__ Temp,
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


	double Ux = ivx[gid], Uy = ivy[gid];
	double dx_e = dXcur[gid], dx_w = dXcur[gid + DIST_SIZE], dx = dx_e + dx_w;
	double dy_n = dXcur[gid + DIST_SIZE * 2], dy_s = dXcur[gid + DIST_SIZE * 3], dy = dy_n + dy_s;

	double xi_0 = 1.0 / dx_e;
	double xi_1 = Ux*DTFD;
	double xi_2 = xi_0*xi_1;
	double xi_3 = 1.0 / dy_n;
	double xi_4 = Uy*DTFD;
	double xi_5 = xi_3*xi_4;
	double xi_6 = 1.0 / dx_w;
	double xi_7 = xi_1*xi_6;
	double xi_8 = 1.0 / dy_s;
	double xi_9 = xi_4*xi_8;
	double xi_10 = 2.0*alpha*DTFD;
	double xi_11 = 1.0 / dx;
	double xi_12 = 2.0*alpha*DTFD*xi_11;
	double xi_13 = 1.0 / dy;
	double xi_14 = 2.0*alpha*DTFD*xi_13;
	//Tc = -xi_0*xi_10*xi_6 - xi_10*xi_3*xi_8 + xi_2 + xi_5 - xi_7 - xi_9 + 1.0;
	//Te = -dx_w*xi_11*xi_2 + xi_0*xi_12;
	//Tw = dx_e*xi_11*xi_7 + xi_12*xi_6;
	//Tn = -dy_s*xi_13*xi_5 + xi_14*xi_3;
	//Ts = dy_n*xi_13*xi_9 + xi_14*xi_8;

	if (Sind >= 0)
		Amat[Sind] = dy_n*xi_13*xi_9 + xi_14*xi_8;

	if (Wind >= 0)
		Amat[Wind] = dx_e*xi_11*xi_7 + xi_12*xi_6;
	else if (i == 0)
		bvec[gid] = TFD_X_IN_VAL * (dx_e*xi_11*xi_7 + xi_12*xi_6);

	Amat[Cind] = -xi_0*xi_10*xi_6 - xi_10*xi_3*xi_8 + xi_2 + xi_5 - xi_7 - xi_9 + 1.0;

	if (Eind >= 0)
		Amat[Eind] = -dx_w*xi_11*xi_2 + xi_0*xi_12;
	else if (i == CHANNEL_LENGTH - 1)
		bvec[gid] = Temp[gid] * (-dx_w*xi_11*xi_2 + xi_0*xi_12);

	if (Nind >= 0)
		Amat[Nind] = -dy_s*xi_13*xi_5 + xi_14*xi_3;

}


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