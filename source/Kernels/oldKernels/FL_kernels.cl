int testInside_Shift(double2 vd, double2 v0, double2 v1)
{
	double2 v10 = v1 - v0, vd0 = vd - v0, vd1 = vd - v1;
			
	if(fabs(length(v10) - length(vd0) - length(vd1)) < CEPS)
		return 1;
	return 0;
}

int bcFindIntersectionLinePlane_shift(double *dist, double2 vL0, bLinks BL, double2 *vN)
{
	double2 vP1 = BL.vP1;
	double2 vP0 = BL.vP0;
	*vN = BL.vNvec;
	double2 vLd = (vP0 - vL0);
	*dist = dot(*vN, vLd);
	if((*dist) <= 0.)
		return 0;
	double2 vCcut = vL0+(*vN) * (*dist);
	return testInside_Shift(vCcut,vP0,vP1);
}

nodeC calc_tr_coeffs(nodeC NodC, int tnum, double dXcx, double dXcy, double dXx, double dXy, double Ka, double Ks)
{
	switch (tnum)
	{
	case -1:
	{
			   break;
	}
	case 1:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( 1., 0., 0., 0.);
			  NodC.CoeffT10 = (double4)((C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0.);
			  NodC.CoeffT01 = (double4)((C1y - 1.) / dXcy + 1., 0., C2y / dXcy, 0.);
			  NodC.CoeffT11 = (double4)((C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., C2x / dXcx, C2y / dXcy, 0.);

			  NodC.CoeffU00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffU10 = (double4)( 1. - 1. / dXcx, 0., 0., 0. );
			  NodC.CoeffU01 = (double4)( 1. - 1. / dXcy, 0., 0., 0. );
			  NodC.CoeffU11 = (double4)( 1. - 1. / dXcy - 1. / dXcx, 0, 0, 0 );
			  break;
	}
	case 2:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( C2x / dXcx, (C1x - 1.) / dXcx + 1., 0, 0 );
			  NodC.CoeffT10 = (double4)( 0., 1., 0., 0. );
			  NodC.CoeffT01 = (double4)( C2x / dXcx, (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., 0., C2y / dXcy );
			  NodC.CoeffT11 = (double4)( 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy );

			  NodC.CoeffU00 = (double4)( 0., 1. - 1. / dXcx, 0., 0. );
			  NodC.CoeffU10 = (double4)( 0., 1., 0., 0. );
			  NodC.CoeffU01 = (double4)( 0., 1. - 1. / dXcx - 1. / dXcy, 0., 0. );
			  NodC.CoeffU11 = (double4)( 0., 1. - 1. / dXcy, 0., 0. );
			  break;
	}
	case 3:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffT10 = (double4)( 0., 1., 0., 0. );
			  NodC.CoeffT01 = (double4)( (C1x - 1.) / dXcx + 1., 0, C2x / dXcx, 0. );
			  NodC.CoeffT11 = (double4)( 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy );

			  NodC.CoeffU00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffU10 = (double4)( 0., 1., 0., 0. );
			  NodC.CoeffU01 = (double4)( 1. - 1. / dXcx, 0., 0., 0. );
			  NodC.CoeffU11 = (double4)( 0., 1. - 1. / dXcy, 0., 0. );
			  break;
	}
	case 4:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( 0., C2y / dXcy, C2x / dXcx, (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1. );
			  NodC.CoeffT10 = (double4)( 0., C2y / dXcy, 0., (C1y - 1.) / dXcy + 1. );
			  NodC.CoeffT01 = (double4)( 0., 0., C2x / dXcx, (C1x - 1.) / dXcx + 1. );
			  NodC.CoeffT11 = (double4)( 0., 0., 0., 1. );

			  NodC.CoeffU00 = (double4)( 0., 0., 0., 1. - 1. / dXcy - 1. / dXcx );
			  NodC.CoeffU10 = (double4)( 0., 0., 0., 1. - 1. / dXcy );
			  NodC.CoeffU01 = (double4)( 0., 0., 0., 1. - 1. / dXcx );
			  NodC.CoeffU11 = (double4)( 0., 0., 0., 1. );
			  break;
	}
	case 5:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffT10 = (double4)( (C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0. );
			  NodC.CoeffT01 = (double4)( (C1y - 1.) / dXcy + 1., 0., C2y / dXcy, 0. );
			  NodC.CoeffT11 = (double4)( 0., 0., 0., 1. );

			  NodC.CoeffU00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffU10 = (double4)( 1. - 1. / dXcx, 0., 0., 0. );
			  NodC.CoeffU01 = (double4)( 1. - 1. / dXcy, 0., 0., 0. );
			  NodC.CoeffU11 = (double4)( 0., 0., 0., 1. );
			  break;
	}
	case 6:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( C2x / dXcx, (C1x - 1.) / dXcx + 1., 0., 0. );
			  NodC.CoeffT10 = (double4)( 0., 1., 0., 0. );
			  NodC.CoeffT01 = (double4)( 0., 0., C2y / dXcy, (C1y - 1.) / dXcy + 1. );
			  NodC.CoeffT11 = (double4)( 0., 0., 0., 1. );

			  NodC.CoeffU00 = (double4)( 0., 1. - 1. / dXcx, 0., 0. );
			  NodC.CoeffU10 = (double4)( 0., 1., 0., 0. );
			  NodC.CoeffU01 = (double4)( 0., 0., 0., 1. - 1. / dXcy );
			  NodC.CoeffU11 = (double4)( 0., 0., 0., 1. );
			  break;
	}
	case 7:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffT10 = (double4)( 0., 1., 0., 0. );
			  NodC.CoeffT01 = (double4)( 1. - (C1x - 1.) / dXcx, 0., -C2x / dXcx, 0. );
			  NodC.CoeffT11 = (double4)( 0., 0., 0., 1. );

			  NodC.CoeffU00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffU10 = (double4)( 0., 1., 0., 0. );
			  NodC.CoeffU01 = (double4)( 1. - 1. / dXcx, 0., 0., 0. );
			  NodC.CoeffU11 = (double4)( 0., 0., 0., 1. );
			  break;
	}
	case 8:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( C2y / dXcy, 0., (C1y - 1.) / dXcy + 1., 0. );
			  NodC.CoeffT10 = (double4)( C2y / dXcy, 0., (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., C2x / dXcx );
			  NodC.CoeffT01 = (double4)( 0., 0., 1., 0. );
			  NodC.CoeffT11 = (double4)( 0., 0., (C1x - 1.) / dXcx + 1., C2x / dXcx );

			  NodC.CoeffU00 = (double4)( 0, 0, 1 - 1 / dXcy, 0 );
			  NodC.CoeffU10 = (double4)( 0, 0, 1. - 1. / dXcy - 1. / dXcx, 0 );
			  NodC.CoeffU01 = (double4)( 0., 0., 1., 0. );
			  NodC.CoeffU11 = (double4)( 0, 0, 1. - 1. / dXcx, 0 );
			  break;
	}
	case 9:
	{
			  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			  NodC.CoeffT00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffT10 = (double4)( (C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0. );
			  NodC.CoeffT01 = (double4)( 0., 0., 1., 0. );
			  NodC.CoeffT11 = (double4)( 0., 0., (C1y - 1.) / dXcy + 1., C2y / dXcy );

			  NodC.CoeffU00 = (double4)( 1., 0., 0., 0. );
			  NodC.CoeffU10 = (double4)( 1. - 1. / dXcx, 0., 0., 0. );
			  NodC.CoeffU01 = (double4)( 0., 0., 1., 0. );
			  NodC.CoeffU11 = (double4)( 0., 0., 1. - 1. / dXcy, 0. );
			  break;
	}
	case 10:
	{
			   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			   double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			   double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			   NodC.CoeffT00 = (double4)( C2x / dXcx, (C1x - 1.) / dXcx + 1., 0, 0 );
			   NodC.CoeffT10 = (double4)( 0., 1., 0., 0. );
			   NodC.CoeffT01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffT11 = (double4)( 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy );

			   NodC.CoeffU00 = (double4)( 0., 1. - 1. / dXcx, 0., 0. );
			   NodC.CoeffU10 = (double4)( 0., 1., 0., 0. );
			   NodC.CoeffU01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffU11 = (double4)( 0., 1. - 1. / dXcy, 0., 0. );
			   break;
	}
	case 11:
	{
			   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));


			   NodC.CoeffT00 = (double4)( 1., 0., 0., 0. );
			   NodC.CoeffT10 = (double4)( 0., 1., 0., 0. );
			   NodC.CoeffT01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffT11 = (double4)( 0., (C1x - 1.) / dXcx + 1., 0, C2x / dXcx );

			   NodC.CoeffU00 = (double4)( 1., 0., 0., 0. );
			   NodC.CoeffU10 = (double4)( 0., 1., 0., 0. );
			   NodC.CoeffU01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffU11 = (double4)( 0., 1. - 1. / dXcx, 0., 0. );
			   break;
	}
	case 12:
	{
			   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			   double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
			   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
			   double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

			   NodC.CoeffT00 = (double4)( C2x / dXcx, 0., (C1x - 1.) / dXcx + 1., 0 );
			   NodC.CoeffT10 = (double4)( 0., C2y / dXcy, 0., (C1y - 1.) / dXcy + 1. );
			   NodC.CoeffT01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffT11 = (double4)( 0., 0., 0., 1. );

			   NodC.CoeffU00 = (double4)( 0., 0., 1. - 1. / dXcx, 0. );
			   NodC.CoeffU10 = (double4)( 0., 0., 0., 1. - 1. / dXcy );
			   NodC.CoeffU01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffU11 = (double4)( 0., 0., 0., 1. );
			   break;
	}
	case 13:
	{
			   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));

			   NodC.CoeffT00 = (double4)( 1., 0., 0., 0 );
			   NodC.CoeffT10 = (double4)( 0., C2x / dXcx, 0., (C1x - 1.) / dXcy + 1. );
			   NodC.CoeffT01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffT11 = (double4)( 0., 0., 0., 1. );

			   NodC.CoeffU00 = (double4)( 1., 0., 0., 0. );
			   NodC.CoeffU10 = (double4)( 0., 0., 0., 1. - 1. / dXcx );
			   NodC.CoeffU01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffU11 = (double4)( 0., 0., 0., 1. );
			   break;
	}
	case 14:
	{
			   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
			   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));

			   NodC.CoeffT00 = (double4)( C2x / dXcx, 0., (C1x - 1.) / dXcx + 1., 0 );
			   NodC.CoeffT10 = (double4)( 0., 1., 0., 0. );
			   NodC.CoeffT01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffT11 = (double4)( 0., 0., 0., 1. );

			   NodC.CoeffU00 = (double4)( 0., 0., 1. - 1. / dXcx, 0. );
			   NodC.CoeffU10 = (double4)( 0., 1., 0., 0. );
			   NodC.CoeffU01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffU11 = (double4)( 0., 0., 0., 1. );
			   break;
	}
	case 15:
	{
			   NodC.CoeffT00 = (double4)( 1., 0., 0., 0. );
			   NodC.CoeffT10 = (double4)( 0., 1., 0., 0. );
			   NodC.CoeffT01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffT11 = (double4)( 0., 0., 0., 1. );

			   NodC.CoeffU00 = (double4)( 1., 0., 0., 0. );
			   NodC.CoeffU10 = (double4)( 0., 1., 0., 0. );
			   NodC.CoeffU01 = (double4)( 0., 0., 1., 0. );
			   NodC.CoeffU11 = (double4)( 0., 0., 0., 1. );
			   break;
	}
	}
	
	return NodC;
}

int bcFindIntersectionFD(int2 *vCcut0, int2 *vCcut1, double2 *vCcut, double *dist, 
	double2 vL0, double2 vLd, double2 vC0, double2 vC1, __global char *M_s, int tb_flag)
{
	double2 vP10 = vC1 - vC0;
	double2 vN = (double2)(-vP10.y,vP10.x);

	double den = dot(vN, vLd);
	if(den == 0.)
		return 0;

	*dist = dot(vN,(vC0 - vL0)) / den;
	*vCcut = vL0 + vLd * (*dist);
	double2 vd0 = *vCcut - vC0, vd1 = *vCcut - vC1;
	if(fabs(length(vP10) - length(vd0) - length(vd1)) >= CEPS)
		return 0;
		
	double distF = (*dist);
	int n_cut_0 = (int)ceil(distF);
	int n_cut = (int)ceil(distF);
	double dist_0 = 1. - (convert_double(n_cut) - distF);
	*dist = 1. - (convert_double(n_cut) - distF);

	*vCcut1 = convert_int2(vL0) + convert_int2(vLd) * n_cut;
	*vCcut0 = *vCcut1 - convert_int2(vLd);

	int2 ii0 = (int2)(MODFAST((*vCcut0).x, FULLSIZEX), MODFAST((*vCcut0).y, DOMAIN_SIZE_Y));

	if(ii0.x != (*vCcut0).x)
		return 0; 

	if(ii0.y != (*vCcut0).y)
		return 0;

	if(M_s[ii0.x * DOMAIN_SIZE_Y + ii0.y] == 1)
	{
		if(*dist < CEPS)
		{
			*dist = 1.;
			n_cut--;
			(*vCcut1) = convert_int2(vL0) + convert_int2(vLd) * n_cut;
			*vCcut0 = *vCcut1 - convert_int2(vLd);

			if(ii0.x != (*vCcut0).x)
				return 0;

			if(ii0.y != (*vCcut0).y)
				return 0;
		}
		else
			return 0;
	}

	return 1;
}


//Alpha: C,E,W,N,S        dX:E,W,N,S
void bcSetInterfaceNode(int2 ii0, int2 ii1, int dir, double dist,
	int xyz, __global char *M, __global char *M_sf, __global char *M_o,
	__global double *dX, __global double *dX_cur,
	__global int *Stor)
{
	if (ii0.x >= FULLSIZEX || ii0.x < 0)
		return;
	if (ii0.y >= DOMAIN_SIZE_Y || ii0.y < 0)
		return;

	if (ii1.x >= FULLSIZEX || ii1.x < 0)
		return;
	if (ii1.y >= DOMAIN_SIZE_Y || ii1.y < 0)
		return;

	if (M[ii0.x*DOMAIN_SIZE_Y + ii0.y] == 0)
		return;
	if (M[ii1.x*DOMAIN_SIZE_Y + ii1.y] == 1)
		return;

	int ii0_1D = ii0.x * DOMAIN_SIZE_Y + ii0.y;
	int ii1_1D = ii1.x * DOMAIN_SIZE_Y + ii1.y;


	int yred0 = Stor[ii0_1D];
	int yred1 = Stor[ii1_1D];

	int ii0red = ii0.x * FULLSIZEY + yred0;
	int ii1red = ii1.x * FULLSIZEY + yred1;
	if (xyz == 0)
	{
		if (CDirX[dir] == -1)
		{
			dX_cur[ii0red * 4 + 1] = dist;
			if (M_o[ii1_1D] == 1 && yred1 > -1)
				dX_cur[ii1red * 4] = 1. - dist;

		}
		else
		{
			dX_cur[ii0red * 4] = dist;
			if (M_o[ii1_1D] == 1 && yred1 > -1)
				dX_cur[ii1red * 4 + 1] = 1. - dist;
		}
	}
	else
	{
		if (CDirY[dir] == -1)
		{
			dX_cur[ii0red * 4 + 3] = dist;
			if (M_o[ii1_1D] == 1 && yred1 > -1)
				dX_cur[ii1red * 4 + 2] = 1. - dist;
		}
		else
		{
			dX_cur[ii0red * 4 + 2] = dist;
			if (M_o[ii1_1D] == 1 && yred1 > -1)
				dX_cur[ii1red * 4 + 3] = 1. - dist;
		}
	}
}


double Cross2(double2 v1, double2 v2)
{//absolute value of cross product
	return fabs(v1.x*v2.y - v1.y*v2.x);
}

int test_inside(double2 P0, double2 P1, double2 P2, double2 L0)
{
	double2 v10 = P1 - P0, v20 = P2 - P0, v21 = P2 - P1, vd0 = L0 - P0, vd1 = L0 - P1;
	double stot = Cross2(v10,v20), s0 = Cross2(v10, vd0), s1 = Cross2(vd0,v20), s2 = Cross2(v21,vd1);

	if(fabs(stot - s0 - s1 - s2) < CEPS)
		return 1;
	return 0;
}

int2 min2(double2 v0, double2 v1)
{
	if(v1.x < v0.x) v0.x = v1.x;
	if(v1.y < v0.y) v0.y = v1.y;
	return convert_int2(v0);
}
				
int2 max2(double2 v0, double2 v1)
{
	if(v1.x > v0.x) v0.x = v1.x;
	if(v1.y > v0.y) v0.y = v1.y;
	
	return convert_int2(ceil(v0));
}

int2 min4(double2 v0, double2 v1, double2 v2, double2 v3)
{
	if(v1.x < v0.x) v0.x = v1.x;
	if(v2.x < v0.x) v0.x = v2.x;
	if(v3.x < v0.x) v0.x = v3.x;
	if(v1.y < v0.y) v0.y = v1.y;
	if(v2.y < v0.y) v0.y = v2.y;
	if(v3.y < v0.y) v0.y = v3.y;

	return convert_int2(v0);
}
				
int2 max4(double2 v0, double2 v1, double2 v2, double2 v3)
{
	if(v1.x > v0.x) v0.x = v1.x;
	if(v2.x > v0.x) v0.x = v2.x;
	if(v3.x > v0.x) v0.x = v3.x;
	if(v1.y > v0.y) v0.y = v1.y;
	if(v2.y > v0.y) v0.y = v2.y;
	if(v3.y > v0.y) v0.y = v3.y;

	return convert_int2(ceil(v0));
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_FL_SHIFT, 1, 1)))
void Shift_walls(__global uint *BLdep,
	__global double *BLdep_tot,
	__global double2 *C,
	__global double2 *C0,
	__global foulI *FI,
	int num_BLs)
{
	int i = get_global_id(0);


	if (i >= num_BLs)
		return;
		
	foulI FItemp = FI[i];
	uint BLstart = FItemp.BL_ind;
	uint8 Depi;
	uint CurBLstart;
	double8 Depf;
	double disp = 0.;
	double4 WeightL = FItemp.WeightsL;
	double4 WeightR = FItemp.WeightsR;
	uint Cind = FItemp.C_ind;
	

	for (int j = 0; j < NUM_PAR_SIZES; j++)
	{
		CurBLstart = BLstart*NUM_PAR_SIZES + j;
		Depi = (uint8)(BLdep[CurBLstart], BLdep[CurBLstart+NUM_PAR_SIZES],
			BLdep[CurBLstart+2*NUM_PAR_SIZES], BLdep[CurBLstart+3*NUM_PAR_SIZES],
			BLdep[CurBLstart+4*NUM_PAR_SIZES], BLdep[CurBLstart+5*NUM_PAR_SIZES],
			BLdep[CurBLstart+6*NUM_PAR_SIZES], BLdep[CurBLstart+7*NUM_PAR_SIZES]);


		Depf = convert_double8(Depi);
		double disttemp = dot(Depf.lo, WeightL) + dot(Depf.hi, WeightR);
		disttemp += BLdep_tot[i*NUM_PAR_SIZES + j];
		BLdep_tot[i*NUM_PAR_SIZES + j] = disttemp;
		disp += (disttemp * Par_multiplier[j]);
	}
	FI[i].disp = disp;
	C[Cind] = C0[Cind] + disp*FItemp.vN;	
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_FL_SHIFT, 1, 1)))
void Smooth_walls1(__global double *BLdep_tot,
__global double *BLdep_tot2,
__global foulI *FI,
int num_BLs,
int2 SSbottom,
int2 SStop)
{
	int i = get_global_id(0);


	if (i >= num_BLs)
		return;

	int2 startstop = (i < num_BLs / 2) ? (SSbottom) : (SStop);

	int kkstart = MAX(i - NEIGHS_PER_SIDE_SMOOTHING, startstop.x);
	int kkstop = MIN(i + 1 + NEIGHS_PER_SIDE_SMOOTHING, startstop.y);
	double num_locs = convert_double(kkstop - kkstart);
	for (int j = 0; j < NUM_PAR_SIZES; j++)
	{
		double Depf = 0.;
		int k = kkstart;
		while(k < kkstop)
		{
			Depf += BLdep_tot[k*NUM_PAR_SIZES + j] * PERCENT_USED_IN_SMOOTHING;
			k++;
		}
		BLdep_tot2[i*NUM_PAR_SIZES + j] = PERCENT_NOT_USED_IN_SMOOTHING*BLdep_tot[i*NUM_PAR_SIZES + j] + Depf / num_locs;
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_FL_SHIFT, 1, 1)))
void Smooth_walls2(__global double *BLdep_tot,
__global double2 *C,
__global double2 *C0,
__global foulI *FI,
int num_BLs)
{
	int i = get_global_id(0);

	if (i >= num_BLs)
		return;

	foulI FItemp = FI[i];
	double disp = 0.;
	uint Cind = FItemp.C_ind;

	for (int j = 0; j < NUM_PAR_SIZES; j++)
	{
		disp += (BLdep_tot[i*NUM_PAR_SIZES + j] * Par_multiplier[j]);
	}
	
	FI[i].disp = disp;
	C[Cind] = C0[Cind] + disp*FItemp.vN;
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_FL_RAMP, 1, 1)))
	void Ramp_ends(__global rampI *RI,
	__global double *Coeff,
	__global double2 *C)
{
	int gid = get_global_id(0);

	if (gid >= WORKGROUPSIZE_FL_RAMP * 4)
		return;
	rampI RItemp = RI[gid];
	uint IOind = RItemp.IOind;
	uint Cind = RItemp.Cind;
	double Ycur = C[IOind].y;

	double Yshift = Ycur - RItemp.Ybegin;
	C[Cind].y = Ycur - RItemp.Coeff * Yshift;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEM, 1, 1)))
	void Update_M(__global char *M,
	__global char *M_o,
	__global int *Stor,
	__global double2 *C,
	__global double2 *C0,
	int2 Dmax,//vlb.nX-1, vlb.nY-1
	uint nCnodes, //vls.nN - 1 nodes,
	uint SkipNode,
	__global double2 *J,
	__global double *Ro) //vls.nN/2 - 1 
{
	int i = get_global_id(0);

	if(i >= nCnodes || i == SkipNode)
		return;

	double2 P0 = C0[i], P1 = C0[i+1], P2 = C[i], P3 = C[i+1];
	if (P0.x == P2.x && P0.y == P2.y && P1.x == P3.x && P1.y == P3.y)
		return;

	int2 Cmax = max4(P0,P1,P2,P3), Cmin = min4(P0,P1,P2,P3);
		
	if(Cmax.y > Dmax.y)
		Cmax.y = Dmax.y;
	if(Cmin.y > Dmax.y)
		return;

	if(Cmax.y < 0)
		return;
	if(Cmin.y < 0)
		Cmin.y = 0;

	if(Cmax.x > Dmax.x)
		Cmax.x = Dmax.x;
	if(Cmin.x > Dmax.x)
		return;

	if(Cmax.x < 0)
		return;
	if(Cmin.x < 0)
		Cmin.x = 0;
		
	for(int ii = Cmin.x; ii <= Cmax.x; ii++)
	{
		for(int jj = Cmin.y; jj <= Cmax.y; jj++)
		{
			int indy = ii*DOMAIN_SIZE_Y+jj;
			int yred = Stor[indy];
						
			if(M_o[indy] == 0)
				continue;
			if (NEQUALS2D(P0, P2))
			{
				if(test_inside(P0,P1,P2,(double2)(ii,jj)))
				{
					M[indy] = 0;
					if (yred > -1)
					{
						J[ii*FULLSIZEY_UT + yred] = 0.;
						Ro[ii*FULLSIZEY_UT + yred] = 1.;
					}
					continue;
				}
			}
			if (NEQUALS2D(P1, P3))
			{
				if(test_inside(P3,P2,P1,(double2)(ii,jj)))
				{
					M[indy] = 0;
					if (yred > -1)
					{
						J[ii*FULLSIZEY_UT + yred] = 0.;
						Ro[ii*FULLSIZEY_UT + yred] = 1.;
					}
					continue;
				}
			}
		}
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
	void Fill_Sf_kernel(__global char *M,
	__global char *M_o,
	__global char *s,
	__global char *sf)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	
	if(i >= FULLSIZEX || j >= DOMAIN_SIZE_Y)
		return;
		
	int ind = i*DOMAIN_SIZE_Y + j;	
	if(M[ind] == 0)
	{
		s[ind] = 1;
		if(M_o[ind] == 1)
		{
			sf[ind] = 1;
		}
	}	
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RED, 1, 1)))
	void Sum_M_arrays1(__global char *M,
	__global char *sf,
	__global double *output,
	__local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;
	
	sdata[tid * 2] = convert_double(M[stride] + M[stride + 1]);
	sdata[tid*2+1] = convert_double( sf[stride] + sf[stride + 1] );

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid*2] += sdata[(tid + s)*2];
			sdata[tid*2+1] += sdata[(tid + s)*2+1];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid*2] = sdata[0];
		output[bid*2+1] = sdata[1];
	}

}

__kernel __attribute__((reqd_work_group_size(1, 1, 1)))
	void Sum_M_arrays2(__global double *input,
	__global double *output,
	int nBlocksMred)//Fluid mass, fouling layer mass
{
	int gid = get_global_id(0);
	
	double temp = input[gid];
	for (unsigned int s = 1; s < nBlocksMred; s++)
	{
		temp += input[s * 2 + gid];
	}

	output[gid] = temp;

}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
	void Reduce_M_array(__global char *M,
	__global char *Mred,
	__global int *Act)
{
	int ix = get_global_id(0);
	int iy = get_global_id(1);
	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
		return;
	
	int red_ind = ix*FULLSIZEY + iy;
	int yact = Act[ix*FULLSIZEY + iy];
	Mred[red_ind] = M[ix * DOMAIN_SIZE_Y + yact]; 
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
	void update_LB_pos(__global char *M,
	__global int *Prop_loc,
	__global int *Act,
	__global int *Stor)
{
	int xloc = get_global_id(0);
	int yred = get_global_id(1);
	
	if (xloc >= FULLSIZEX || yred >= FULLSIZEY)
		return;
		
	int yloc = Act[xloc*FULLSIZEY + yred];
	
	int2 ii = (int2)(xloc, yloc);

	for (int k = 0; k < 8; k++)
	{
		int dirval = k;
		int2 CC = (int2)(CDirX[dirval], CDirY[dirval]);
		int2 ii1 = ii + CC;

#ifdef INLET_OUTLET_BC
		if (ii1.x == -1 || ii1.x == FULLSIZEX)
		{
			ii1 = ii;
			dirval = RevDir[dirval];
		}
#else
		ii1.x = MODFAST(ii1.x, FULLSIZEX);
#endif

		if (M[ii1.x * DOMAIN_SIZE_Y + ii1.y] == 0)
		{
			ii1 = ii;
			dirval = RevDir[dirval];
		}


		ii1.y = Stor[ii1.x * DOMAIN_SIZE_Y + ii1.y];
		int loc_array = xloc*FULLSIZEY + yred;   ///location int Dir, Prop arrays to store values
		int block = ii1.x*FULLSIZEY + ii1.y;
#ifdef SOA_STORAGE
		Prop_loc[loc_array*8+k] = block + dirval*FULLSIZEXY;
#else
		Prop_loc[loc_array * 8 + k] = block*8 + dirval;
#endif
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEFD, 1, 1)))
	void Update_FD_dir_arrays(__global bLinks *BL,
	__global double2 *C,
	__global char *M,
	__global char *M_s,
	__global char *M_sf,
	__global char *M_o,
	__global double *dX,
	__global double *dX_cur,
	__global int *Stor,
	int2 Dmax,
	int NUMBLINKS,
	__global double2 *Cnot)
{
	int i = get_global_id(0);
	if (i >= NUMBLINKS)
		return;

	bLinks BLtemp = BL[i];
	int n0 = BLtemp.P0ind, n1 = BLtemp.P1ind;

	double2 vC0 = C[n0], vC1 = C[n1];
	double2 vC0not = Cnot[n0], vC1not = Cnot[n1];
	BLtemp.vP0 = vC0;
	BLtemp.vP1 = vC1;

	double2 vT = vC1 - vC0;
	double2 vTnot = vC1not - vC0not;
	
	double2 CenterDist = vC0 + vT / 2. - vC0not - vTnot/2.;
	BLtemp.int_type = (length(CenterDist) >= FOUL_SIZE_SWITCH_SOOT2) ? (1) : (0);
	
	BLtemp.blLen = length(vT);
	BLtemp.vTvec = normalize(vT);
	int tb_flag = 0;
	if (i < NUMBLINKS / 2)
	{
		BLtemp.vNvec = (double2)(-BLtemp.vTvec.y, BLtemp.vTvec.x);
	}
	else
	{
		tb_flag = 1;
		BLtemp.vNvec = (double2)(BLtemp.vTvec.y, -BLtemp.vTvec.x);
	}

	double2 centpos = BLtemp.vP0 + BLtemp.vTvec * BLtemp.blLen / 2.;
	int2 blind = convert_int2(centpos);

	if (blind.x >= 0 && blind.x < FULLSIZEX)
		BLtemp.Node_loc = blind.x * FULLSIZEY_TR + blind.y;
	else
		BLtemp.Node_loc = -1;

	double2 vCn = BLtemp.vNvec;

	BL[i] = BLtemp;

	int2 vCmin = min2(vC0, vC1), vCmax = max2(vC0, vC1);

	if (vCmin.x < 0) vCmin.x = 0;
	if (vCmin.x > Dmax.x) return;
	if (vCmax.x < 0) return;
	if (vCmax.x > Dmax.x) vCmax.x = Dmax.x;

	if (vCmin.y < 0) vCmin.y = 0;
	if (vCmin.y > Dmax.y) return;
	if (vCmax.y < 0) return;
	if (vCmax.y > Dmax.y) vCmax.y = Dmax.y;

	int dir = 0;
	double dist;
	double2 Cdir = (double2)(CDirX[dir], CDirY[dir]);
	double Cndir = dot(vCn, Cdir);
	
	int2 vCcut0, vCcut1;
	double2 vCcut;
	int kk, j;

	if (Cndir != 0.)
	{
		if (Cndir > 0.)
		{
			dir = RevDir[dir];
			Cdir = (double2)(CDirX[dir], CDirY[dir]);
		}

		if (CDirX[dir] == 1) kk = vCmin.x;
		if (CDirX[dir] == -1) kk = vCmax.x;

		for (j = vCmin.y; j <= vCmax.y; j++)
		{
			if (bcFindIntersectionFD(&vCcut0, &vCcut1, &vCcut, &dist,
				convert_double2((int2)(kk, j)), Cdir, vC0, vC1, M_s, tb_flag))
			{

				bcSetInterfaceNode(vCcut0, vCcut1, dir, dist, 0, M, M_sf, M_o, dX, dX_cur, Stor);
			}
		}
	}

	dir = 2;
	Cdir = (double2)(CDirX[dir], CDirY[dir]);
	Cndir = dot(vCn, Cdir);
	if (Cndir != 0.)
	{
		if (Cndir > 0.)
		{
			dir = RevDir[dir];
			Cdir = (double2)(CDirX[dir], CDirY[dir]);
		}

		if (CDirY[dir] == 1) j = vCmin.y;
		if (CDirY[dir] == -1) j = vCmax.y;
		for (kk = vCmin.x; kk <= vCmax.x; kk++)
		{
			if (bcFindIntersectionFD(&vCcut0, &vCcut1, &vCcut, &dist,
				convert_double2((int2)(kk, j)), Cdir, vC0, vC1, M_s, tb_flag))
			{
				bcSetInterfaceNode(vCcut0, vCcut1, dir, dist, 1, M, M_sf, M_o, dX, dX_cur, Stor);
			}
		}
	}
}

int test_cross(double2 vL0, double2 vLd, double2 vP0, double2 vP1)
{
	double2 vP10 = vP1 - vP0;
	double2 vNm = (double2)(-1.*vP10.y, vP10.x);
	double den = dot(vNm, vLd);
	if (den == 0.)
		return 0;
	double dist = dot(vNm, (vP0 - vL0)) / den;
	double2 vC = vL0 + vLd * dist;
	double2 vd0 = vC - vP0, vd1 = vC - vP1;
	if (dist >= 0. && dist <= 1.)
	{
		if (fabs(length(vP10) - length(vd0) - length(vd1)) < CEPS)
		{
			return 1;
		}
	}
	return 0;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
	void Update_BD_Wall_Nodes(__global bLinks *BL,
	__global nodeI *Nod,
	__global int *BLind_ind,
	__global int *nwallnodes,
	BLbound BLb,
	uint xstart,
	uint xstop,
	__global int *Active_Inds)
{
	//zero BLind_ind before calling this kernel
	int i = get_global_id(0);
	int wallflagval;

	if (i > BLb.MIN_BL_BOT && i <= BLb.MAX_BL_BOT)
		wallflagval = 1;
	else if (i > BLb.MIN_BL_TOP && i <= BLb.MAX_BL_TOP)
		wallflagval = 2;
	else
		return;

	double2 c0 = BL[i].vP0 + BL[i].vTvec / 2.;

	int2 c0i = convert_int2(c0)- 1;

	for (int i0 = c0i.x; i0 < c0i.x+3; i0++)
	{
		for (int j0 = c0i.y; j0 < c0i.y+3; j0++)
		{
			int nind = i0*FULLSIZEY_TR + j0;
			if (Active_Inds[nind] != 1)
				continue;
			int BLind_temp = atomic_add(&BLind_ind[nind], 1);
			if (BLind_temp == 0)
			{
				atomic_add(&nwallnodes[0], 1);
			}
			if (BLind_temp == MAX_BL_PER_NODE)
				continue;

			Nod[nind].Wall_Flag = wallflagval;
			Nod[nind].BLind[BLind_temp] = i;
		}
	}
}
//
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEWALL, 1, 1)))
//	void find_wall_nodes(__global nodeI *NodI,
//	__global int *Winds,
//	__global int *cur_loc,
//	__global int2 *ActiveNodes, int nActiveNodes)
//{
//	int i = get_global_id(0);
//	if (i >= nActiveNodes)
//		return;
//
//	int2 ii = ActiveNodes[i];
//
//	if (NodI[ii.x*FULLSIZEY_TR+ii.y].Wall_Flag != 0)
//	{
//		int ind = atomic_add(&cur_loc[0], 1);
//		Winds[ind] = i;
//	}
//}

__kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void find_wall_nodes(__global nodeI *NodI,
	__global int *Winds,
	__global int2 *ActiveNodes,
	int nActiveNodes)
{
	int count = 0;
	for (int i = 0; i < nActiveNodes; i++)
	{
		int2 ii = ActiveNodes[i];

		if (NodI[ii.x*FULLSIZEY_TR + ii.y].Wall_Flag != 0)
		{
			Winds[count++] = i;
		}
	}
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATE_GL, 1, 1)))
	void update_GL_wall(__global double2 *C,
	__global float2 *BotWall,
	__global float2 *TopWall,
	uint num_nodes)
{
	int i = get_global_id(0);
	if (i >= num_nodes)
		return;

	int itop = 2 * num_nodes - 1 - i;
	TopWall[i] = convert_float2(C[itop]);
	BotWall[i] = convert_float2(C[i]);
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RERELEASE, 1, 1)))
	void Shift_deposited_particles(__global par *P,
	__global bLinks *BLlist,
	uint offset,
	uint maxel)
{
	int i = get_global_id(0);

	if (i >= maxel)
		return;

	i += offset;

	par Pcur = P[i];

	if (Pcur.Dep_Flag < 0)
		return;

	bLinks BLtemp = BLlist[Pcur.Dep_Flag];
	Pcur.pos = BLtemp.vP0 + BLtemp.vTvec * BLtemp.blLen/2.;

	if (Pcur.pos.x >= X_MAX_VAL)
	{
		Pcur.Dep_Flag = -2;
		Pcur.loc = -2;
	}
	
	P[i] = Pcur;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SHIFT_PAR, 1, 1)))
	void Shift_particles(__global par *P,
	__global nodeI *nI,
	__global bLinks *BLlist,
	double X_release, double Ymin, double Ymax)
{

	int i = get_global_id(0);
	if(i >= TRC_NUM_TRACERS)
		return;

	par Pcur = P[i];
	double2 C2 = Pcur.pos;

	if ((C2.x < X_release) && (C2.y > Ymax || C2.y < Ymin))
	{
		Pcur.Dep_Flag = -2;
		Pcur.loc = -2;
		P[i] = Pcur;
		return;
	}

	nodeI Ncur = nI[Pcur.loc];
	if(Ncur.Wall_Flag == 0)
		return;
	
	double dist, dist_use = 100.;							//||dCc||/||dC|| where dCc is the distance between  
	double2 vN, vN_use;			//Intersection pt, Normal vector, Velocity of surf at vCcut

	int bl, bl_use = -1;									//boundary link index (value stored in each linked list element

	for (int k = 0; k < MAX_BL_PER_NODE; k++)
	{
		bl = Ncur.BLind[k];
		if(bl == -1)
			continue;
			
		if(bcFindIntersectionLinePlane_shift(&dist, C2, BLlist[bl], &vN))
		{
			if(dist < dist_use)
			{
				bl_use = bl;
				dist_use = dist;
				vN_use = vN;
			}
		}
	}


	if(bl_use != -1)
	{
		Pcur.pos += vN_use * 2. * dist_use;
		int2 Posi = convert_int2(Pcur.pos);
		Pcur.loc = Posi.x*FULLSIZEY_TR + Posi.y;
		if (Pcur.pos.x >= X_MAX_VAL)
		{
			Pcur.Dep_Flag = -2;
			Pcur.loc = -2;
		}
		P[i] = Pcur;
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SHIFT_PAR, 1, 1)))
void Clear_OOB_Particles(__global par *P,
double X_release, double Ymin, double Ymax)
{
	int i = get_global_id(0);
	if (i >= TRC_NUM_TRACERS)
		return;

	par Pcur = P[i];
	if (Pcur.pos.x >= X_release)
		return;
	
	if (Pcur.pos.y <= Ymax && Pcur.pos.y >= Ymin)
		return;
	
	Pcur.Dep_Flag = -2;
	Pcur.loc = -2;
	P[i] = Pcur;
}


//Ynum channel height+2
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR, 1, 1)))
	void Update_BD_Nodes(__global char *M,
	__global nodeC *nC,
	__global double *dX_cur,
	__global double *dX,
	__global int *Stor,
	__global int2 *ActiveNodes,
	int nActiveNodes)
{
	int ii = get_global_id(0);

	
	if (ii >= nActiveNodes)
		return;
	
	int2 Nod_xy = ActiveNodes[ii];
	int Nod_ind = Nod_xy.x * FULLSIZEY_TR + Nod_xy.y;
	
	//Nact Ntemp = Neighs[Nod_ind];
	uint ii00 = Nod_xy.x * DOMAIN_SIZE_Y + Nod_xy.y;
	uint ii10 = (Nod_xy.x + 1) * DOMAIN_SIZE_Y + Nod_xy.y;
	uint ii01 = Nod_xy.x * DOMAIN_SIZE_Y + (Nod_xy.y + 1);
	uint ii11 = (Nod_xy.x + 1) * DOMAIN_SIZE_Y + (Nod_xy.y + 1);
	
	uint ii00r = Nod_xy.x * FULLSIZEY + Stor[ii00];
	uint ii10r = (Nod_xy.x + 1) * FULLSIZEY + Stor[ii10];
	uint ii01r = Nod_xy.x * FULLSIZEY + Stor[ii01];
	uint ii11r = (Nod_xy.x + 1) * FULLSIZEY + Stor[ii11];
	
	int tnum = 0;
	if(M[ii00] == 1)
		tnum += 1;
	if(M[ii10] == 1)
		tnum += 2;
	if(M[ii11] == 1)
		tnum += 4;
	if(M[ii01] == 1)
		tnum += 8;
		
	if(tnum == 0)
		tnum = -1;
	
	nodeC NodC = nC[Nod_ind];
		
	double dXcx, dXcy, dXx, dXy;
	
	switch(tnum)
	{
	case -1:
	{
		dXcx = 1.;
		dXcy = 1.;
		dXx = 1.;
		dXy = 1.;
		break;
	}
	case 1:
	{
		dXcx = dX_cur[ii00r*4];
		dXcy = dX_cur[ii00r*4+2];
		dXx = dX[ii00r*4];
		dXy = dX[ii00r*4+2];
		break;
	}
	case 2:
	{
		dXcx = dX_cur[ii10r*4+1];
		dXcy = dX_cur[ii10r*4+2];
		dXx = dX[ii10r*4+1];
		dXy = dX[ii10r*4+2];
		break;
	}
	case 3:
	{
		dXcx = dX_cur[ii00r*4+2];
		dXcy = dX_cur[ii10r*4+2];
		dXx = dX[ii00r*4+2];
		dXy = dX[ii10r*4+2];
		break;
	}
	case 4:
	{
		dXcx = dX_cur[ii11r*4+1];
		dXcy = dX_cur[ii11r*4+3];
		dXx = dX[ii11r*4+1];
		dXy = dX[ii11r*4+3];
		break;
	}
	case 5:
	{
		dXcx = dX_cur[ii00r*4];
		dXcy = dX_cur[ii00r*4+2];
		dXx = dX[ii00r*4];
		dXy = dX[ii00r*4+2];
		break;
	}
	case 6:
	{
		dXcx = dX_cur[ii10r*4+1];
		dXcy = dX_cur[ii11r*4+1];
		dXx = dX[ii10r*4+1];
		dXy = dX[ii11r*4+1];
		break;
	}
	case 7:
	{
		dXcx = dX_cur[ii00r*4+2];
		dXcy = 1.;
		dXx = dX[ii00r*4+2];
		dXy = 1.;
		break;
	}
	case 8:
	{
		dXcx = dX_cur[ii01r*4];
		dXcy = dX_cur[ii01r*4+3];
		dXx = dX[ii01r*4];
		dXy = dX[ii01r*4+3];
		break;
	}
	case 9:
	{
		dXcx = dX_cur[ii00r*4];
		dXcy = dX_cur[ii01r*4];
		dXx = dX[ii00r*4];
		dXy = dX[ii01r*4];
		break;
	}
	case 10:
	{
		dXcx = dX_cur[ii10r*4+1];
		dXcy = dX_cur[ii10r*4+2];
		dXx = dX[ii10r*4+1];
		dXy = dX[ii10r*4+2];
		break;
	}
	case 11:
	{
		dXcx = dX_cur[ii10r*4+2];
		dXcy = 1.;
		dXx = dX[ii10r*4+2];
		dXy = 1.;
		break;
	}
	case 12:
	{
		dXcx = dX_cur[ii01r*4+3];
		dXcy = dX_cur[ii11r*4+3];
		dXx = dX[ii01r*4+3];
		dXy = dX[ii11r*4+3];
		break;
	}
	case 13:
	{
		dXcx = dX_cur[ii11r*4+3];
		dXcy = 1.;
		dXx = dX[ii11r*4+3];
		dXy = 1.;
		break;
	}
	case 14:
	{
		dXcx = dX_cur[ii01r*4+3];
		dXcy = 1.;
		dXx = dX[ii01r*4+3];
		dXy = 1.;
		break;
	}
	case 15:
	{
		dXcx = 1.;
		dXcy = 1.;
		dXx = 1.;
		dXy = 1.;
	}
	}
	
	NodC = calc_tr_coeffs(NodC, tnum, dXcx, dXcy, dXx, dXy, ALPHA_FLUID, ALPHA_FOUL);
	nC[Nod_ind] = NodC;
	
}



__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1))) 
void LB_Update_IBB(__global int2 *ibb_loc,
	__global double2 *ibb_coeff,
	__global int2 *ii0_array,
	__global int2 *iis1_array,
	__global int *dir_array,
	__global int *Stor,
	__global double *D_array,
	int max_el)//7
{
	int i = get_global_id(0);
	if(i >= max_el)
		return;
	
	int2 ii0 = ii0_array[i];
	int yred = Stor[ii0.x*DOMAIN_SIZE_Y + ii0.y];
	double D = D_array[i];
	int Vdir = dir_array[i]-1;
	int Vdirp = RevDir[Vdir];
	int ival = ii0.x * FULLSIZEY + yred;
	int loc1 = ival * 8 + Vdirp;
	ibb_loc[i].x = loc1;

	if (D < 0.5)
	{
		int loc2 = ival * 8 + Vdir;
		ibb_loc[i].y = loc2;
		ibb_coeff[i].x = 2.*D;
		ibb_coeff[i].y = 1. - 2.*D;
	}
	else
	{
		int2 ii1 = iis1_array[i];
		ii1.y = Stor[ii1.x*DOMAIN_SIZE_Y + ii1.y];
		int i2val = ii1.x*FULLSIZEY +ii1.y;
		int loc2 = i2val * 8 + Vdirp;
		ibb_loc[i].y = loc2;
		ibb_coeff[i].x = 1. / (2.*D);
		ibb_coeff[i].y = (2.*D - 1.) / (2.*D);
	}
}