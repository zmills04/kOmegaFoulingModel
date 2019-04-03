
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
void Calc_Ro_J(__global double *Ro_array,
__global double2 *J_array,
__global double *F0,
__global double4 *FA,
__global char *Mred)
{

	int ix = get_global_id(0);
	int iy = get_global_id(1);
	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
		return;

	int jj = ix * FULLSIZEY + iy;
	if (Mred[jj] == 0)
		return;

	int indUT = ix * FULLSIZEY_UT + iy;


	int ii = jj * 2;

	double Fi = F0[jj];
	double4 FAlo = FA[ii];
	double4 FAhi = FA[ii + 1];

	double Ro = Fi + dot(FAhi, (double4)(1.)) + dot(FAlo, (double4)(1.));
	double Jx = dot(CCxlo, FAlo) + dot(CCxhi, FAhi);
	double Jy = dot(CCylo, FAlo) + dot(CCyhi, FAhi);

	Ro_array[indUT] = Ro;
	J_array[indUT] = (double2)(Jx, Jy) / Ro;
}

inline double CalcDeviation(double* FAA, double* Fneq)
{
	double dev = 0.;
	__attribute__((opencl_unroll_hint(9)))
		for (int i = 0; i < 9; i++)
		{
			dev = fmax(dev, Fneq[i] / FAA[i]);
			if (dev > 0.01f)
				return dev;
		}
	return dev;
}

inline double Calc_EntropyIneq(double *FAA, double* Fneq, double alp, double *dent)
{
	double t = FAA[0] + alp*Fneq[0];
	double h = log(t) - log(4./9.);
	double ent = t*h;
	*dent += Fneq[0] * (h + 1.);
	
	t = FAA[1] + alp*Fneq[1];
	h = log(t) - log(1. / 9.);
	ent += t*h;
	*dent += Fneq[1] * (h + 1.);

	t = FAA[2] + alp*Fneq[2];
	h = log(t) - log(1. / 9.);
	ent += t*h;
	*dent += Fneq[2] * (h + 1.);

	t = FAA[3] + alp*Fneq[3];
	h = log(t) - log(1. / 9.);
	ent += t*h;
	*dent += Fneq[3] * (h + 1.);

	t = FAA[4] + alp*Fneq[4];
	h = log(t) - log(1. / 9.);
	ent += t*h;
	*dent += Fneq[4] * (h + 1.);

	t = FAA[5] + alp*Fneq[5];
	h = log(t) - log(1. / 36.);
	ent += t*h;
	*dent += Fneq[5] * (h + 1.);

	t = FAA[6] + alp*Fneq[6];
	h = log(t) - log(1. / 36.);
	ent += t*h;
	*dent += Fneq[6] * (h + 1.);

	t = FAA[7] + alp*Fneq[7];
	h = log(t) - log(1. / 36.);
	ent += t*h;
	*dent += Fneq[7] * (h + 1.);

	t = FAA[8] + alp*Fneq[8];
	h = log(t) - log(1. / 36.);
	ent += t*h;
	*dent += Fneq[8] * (h + 1.);
	
	return ent;
}

inline double Calc_Entropy(double* FAA)
{
	double ent = FAA[0] * (log(FAA[0]) - log(4./9.));
	ent += FAA[1] * (log(FAA[1]) - log(1. / 9.));
	ent += FAA[2] * (log(FAA[2]) - log(1. / 9.));
	ent += FAA[3] * (log(FAA[3]) - log(1. / 9.));
	ent += FAA[4] * (log(FAA[4]) - log(1. / 9.));
	ent += FAA[5] * (log(FAA[5]) - log(1. / 36.));
	ent += FAA[6] * (log(FAA[6]) - log(1. / 36.));
	ent += FAA[7] * (log(FAA[7]) - log(1. / 36.));
	ent += FAA[8] * (log(FAA[8]) - log(1. / 36.));
	return ent;
}

inline double find_max_alpha(double* FAA, double* Fneq)
{
	double max_alpha = 1000.;
	__attribute__((opencl_unroll_hint(9)))
		for (int k = 0; k < 9; k++)
			if (FAA[k] < 0. || Fneq[k] < 0.)
				max_alpha = fmin(max_alpha, -FAA[k] / Fneq[k]);

	return max_alpha;
}

double Calc_Alpha_Root(double *FAA, double* Fneq, double alp)
{
	double ent = Calc_Entropy(FAA);
	int ii = 0;

	double max_alpha = find_max_alpha(FAA, Fneq);

	while (1)
	{
		double delta_ent_der = 0.;
		double ent_ineq = Calc_EntropyIneq(FAA, Fneq, alp, &delta_ent_der);

		if (isnan(ent_ineq) && alp != 1.1)
		{
			alp = 1.1;
			continue;
		}

		double ent_inc = ent_ineq - ent;

		if (fabs(ent_inc) < entropy_tolerance)
			break;

		double new_alpha = alp - ent_inc / delta_ent_der;

		if (new_alpha >= max_alpha)
			new_alpha = (alp + max_alpha) * 0.5;

		if (fabs(new_alpha - alp) < alpha_tolerance)
			break;

		alp = new_alpha;
		ii++;
		if (ii > 1000)
		{
			return 2.;
		}
	}

	if (alp < 1.0 || !isfinite(alp))
	{
		return 2.;
	}
	return alp;
	return 2.;
}

double EstimateAlpha(double *FAA, double *Fneq)
{
	double a1 = 0.;
	double a2 = 0.;
	double a3 = 0.;
	double a4 = 0.;
	for (int k = 0; k < 9; k++)
	{
		double t = Fneq[k];
		double inv = 1. / FAA[k];
		double p = t*t*inv;
		t *= inv;
		a1 += p;
		p *= t;
		a2 += p;
		p *= t;
		a3 += p;
		p *= t;
		a4 += p;
	}
	a1 /= 2.;
	a2 /= -6.;
	a3 /= 12.;
	a4 /= -20.;
	return 2. - 4.*a2 / a1 + 16.*a2*a2 / a1 / a1 - 8.*a3 / a1 + 80.*a2*a3 / a1 / a1 - 80.*a2*a2*a2 / a1 / a1 / a1 - 16.*a4 / a1;
}

double Calc_alpha(double alpha_i, double *FAA, double *Fneq)
{
	double dev = CalcDeviation(FAA, Fneq);
	if (dev < 1.e-6f)
	{
		return 2.0;
	}
	else if (dev < 0.01f)
	{
		return EstimateAlpha(FAA, Fneq);
	}
	else
	{
		return Calc_Alpha_Root(FAA, Fneq, alpha_i);
	}
}

inline void calculate_Feq_variables(double rho, double2 UU, double *Psi, double *Bx, double *By, double *Bxinv, double *Byinv)
{
	double2 op3u = sqrt(1. + 3.*UU*UU);

	*Psi = rho * (2. - op3u.x) * (2. - op3u.y) / 9.;
	*Bx = (2. * UU.x + op3u.x) / (1. - UU.x);
	*Bxinv = 1. / (*Bx);
	*By = (2. * UU.y + op3u.y) / (1. - UU.y);
	*Byinv = 1. / (*By);
}

double calc_scalar_product(double *X, double *Y, double *Fieq)
{
	double ret = 0.;
	__attribute__((opencl_unroll_hint(9)))
		for (int i = 0; i < 9; i++)
		{
			ret += X[i] * Y[i] / Fieq[i];
		}
	return ret;
}

int Max_Hdel(double *Hdel)
{
	__attribute__((opencl_unroll_hint(9)))
		for (int k = 0; k < 9; k++)
		{
			if (fabs(Hdel[k]) > 1e-12)
				return 1;
		}
	return 0;
}

inline double calc_Gamma(double *Sdel, double *Hdel, double *Fieq)
{

	return (Max_Hdel(Hdel) == 1) ? (BETA_INV - (2. - BETA_INV)*calc_scalar_product(Sdel, Hdel, Fieq) / calc_scalar_product(Hdel, Hdel, Fieq)) : (BETA_INV);
}

inline void calc_Hdel(double *Fi, double *Fieq, double *Sdel, double *Hdel)
{
	__attribute__((opencl_unroll_hint(9)))
		for (int i = 0; i < 9; i++)
		{
			Hdel[i] = Fi[i] - Fieq[i] - Sdel[i];
		}
}

inline void calc_Sdist(double *S, double rho, double3 TNP)
{
	S[0] = -TNP.x;
	S[1] = 0.25 * (TNP.x + TNP.y);
	S[2] = 0.25 * (TNP.x + TNP.y);
	S[3] = 0.25 * (TNP.x - TNP.y);
	S[4] = 0.25 * (TNP.x - TNP.y);
	S[5] = 0.25 * TNP.z;
	S[6] = 0.25 * TNP.z;
	S[7] = -0.25 * TNP.z;
	S[8] = -0.25 * TNP.z;
}
//inline void calc_Sdist(double *S, double rho, double3 TNP)
//{
//	S[0] = rho * -TNP.x;
//	S[1] = 0.25 * rho * (TNP.x + TNP.y);
//	S[2] = 0.25 * rho * (TNP.x + TNP.y);
//	S[3] = 0.25 * rho * (TNP.x - TNP.y);
//	S[4] = 0.25 * rho * (TNP.x - TNP.y);
//	S[5] = 0.25 * rho * TNP.z;
//	S[6] = 0.25 * rho * TNP.z;
//	S[7] = -0.25 * rho * TNP.z;
//	S[8] = -0.25 * rho * TNP.z;
//}

inline void calc_Sdel(double *Sdel, double rho, double3 TNP, double3 TNP_eq)
{
	double S[9];
	double Seq[9];
	calc_Sdist(S, rho, TNP);
	calc_Sdist(Seq, rho, TNP_eq);
	__attribute__((opencl_unroll_hint(9)))
	for (int k = 0; k < 9; k++)
	{
		Sdel[k] = S[k] - Seq[k];
	}
}

inline void calculate_moments(double *FA, double *FAA, double *Rho, double2 *UU, double3 *TNP, double3 *TNP_eq)
{
	(*Rho) = FAA[0];
	(*UU) = 0.;
	double M11 = 0., M20 = 0., M02 = 0.;
	__attribute__((opencl_unroll_hint(8)))
		for (int k = 1; k < 9; k++)
		{
			FAA[k] = FA[k - 1];
			(*Rho) += FAA[k];
			(*UU) += CXY_DOUBLE[k - 1] * FAA[k];
			M20 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].x * FAA[k];
			M02 += CXY_DOUBLE[k - 1].y * CXY_DOUBLE[k - 1].y * FAA[k];
			M11 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].y * FAA[k];
		}

	double rhoinv = 1. / (*Rho);
	(*UU) *= rhoinv;
	(*TNP) = (double3)(M20 + M02, M20 - M02, M11)*rhoinv;

	double2 M22eq = (*UU)*(*UU) + 1. / 3.;
	(*TNP_eq) = (double3)(M22eq.x + M22eq.y, M22eq.x - M22eq.y, (*UU).x*(*UU).y);
}



//inline void calculate_moments_2pop(double *FA, double *FAA, double *GA, double *GAA, double *Rho, double2 *UU, double3 *TNP, double3 *TNP_eq,
//	double *RoE2, double2 *qq, double2 *qq_eq, double3 *Req)
//{
//	(*UU) = 0.;
//	(*qq) = 0.;
//	double M11 = 0., M20 = 0., M02 = 0.;
//	__attribute__((opencl_unroll_hint(8)))
//		for (int k = 1; k < 9; k++)
//		{
//			FAA[k] = FA[k - 1];
//			GAA[k] = GA[k - 1];
//			(*Rho) += FAA[k];
//			(*RoE2) += GAA[k];
//			(*UU) += CXY_DOUBLE[k-1] * FAA[k];
//			(*qq) += CXY_DOUBLE[k-1] * GAA[k];
//			M20 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].x * FAA[k];
//			M02 += CXY_DOUBLE[k - 1].y * CXY_DOUBLE[k - 1].y * FAA[k];
//			M11 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].y * FAA[k];
//		}
//
//	double rhoinv = 1. / (*Rho);
//	(*UU) *= rhoinv;
//	(*TNP) = (double3)(M20 + M02, M20 - M02, M11)*rhoinv;
//
//	double2 M22eq = 1. / 3. + (*UU)*(*UU);
//	(*TNP_eq) = (double3)(M22eq.x + M22eq.y, M22eq.x - M22eq.y, (*UU).x*(*UU).y);
//	
//	(*qq_eq) = ((*RoE2) + 2. / 3. * (*Rho)) * (*UU);
//	
//	double3 Pneq = (double3)( M20 - M22eq.x*(*Rho), M02 - M22eq.y*(*Rho), M11 - (*TNP_eq).z * (*Rho) );
//
//	(*qq).x -= 2. * dot((*UU), Pneq.xz);// ((*UU).x*Pneq.x + (*UU).y*Pneq.z);
//	(*qq).y -= 2. * dot((*UU), Pneq.zy);// ((*UU).x*Pneq[2] + (*UU).y*Pneq[1]);
//	
//	(*Req).xy = ((*RoE2) + 2.*(*Rho)/3.) * M22eq + 2.*(*Rho) / 3. * (*UU)*(*UU) - (*RoE2) / 3.;
//	(*Req).z = ((*RoE2) + 4.*(*Rho) / 3.) * (*TNP_eq).z;
//
//
////	(*Req).x = (*RoE2) * M22eq.x + 2.*(*Rho) / 3. * (M22eq.x + (*UU).x*(*UU).x) - (*RoE2) / 3.;
//	//(*Req).y = (*RoE2) * M22eq.y + 2.*(*Rho) / 3. * (M22eq.y + (*UU).y*(*UU).y) - (*RoE2) / 3.;
//	//(*Req).z = ((*RoE2) + 4.*(*Rho) / 3.) * (*TNP_eq).z;
//}

inline void calculate_moments_SRT_2pop(double *FA, double *FAA, double *GA, double *GAA, double *Rho, double2 *UU,
	double *RoE2, double2 *qq, double2 *qq_eq, double3 *Req)
{
	(*UU) = 0.;
	(*qq) = 0.;
	double M11 = 0., M20 = 0., M02 = 0.;
	__attribute__((opencl_unroll_hint(8)))
		for (int k = 1; k < 9; k++)
		{
			FAA[k] = FA[k - 1];
			GAA[k] = GA[k - 1];
			(*Rho) += FAA[k];
			(*RoE2) += GAA[k];
			(*UU) += CXY_DOUBLE[k - 1] * FAA[k];
			(*qq) += CXY_DOUBLE[k - 1] * GAA[k];
			M20 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].x * FAA[k];
			M02 += CXY_DOUBLE[k - 1].y * CXY_DOUBLE[k - 1].y * FAA[k];
			M11 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].y * FAA[k];
		}

	(*UU) /= (*Rho);
	double2 M22eq = (1. / 3. + (*UU)*(*UU));

	(*qq_eq) = ((*RoE2) + 2. * (*Rho) / 3.) * (*UU);

	double3 Pneq = (double3)(M20 - M22eq.x*(*Rho), M02 - M22eq.y*(*Rho), M11 - (*Rho)*(*UU).x*(*UU).y);

	(*qq).x -= 2. * dot((*UU), Pneq.xz);// ((*UU).x*Pneq.x + (*UU).y*Pneq.z);
	(*qq).y -= 2. * dot((*UU), Pneq.zy);// ((*UU).x*Pneq[2] + (*UU).y*Pneq[1]);

	(*Req).xy = ((*RoE2) + 2. * (*Rho) / 3.) * M22eq + 2. * (*Rho) / 3. * (*UU)*(*UU) - (*RoE2) / 3.;
	(*Req).z = ((*RoE2) + 4. * (*Rho) / 3.) * (*UU).x*(*UU).y;
}

inline void calculate_moments_SRT_2pop_dff(double *FA, double *FAA, double *GA, double *GAA, double *Rho, double2 *UU,
	double *RoE2, double2 *qq, double2 *qq_eq, double3 *Req, double dff)
{
	(*UU) = 0.;
	(*qq) = 0.;
	double M11 = 0., M20 = 0., M02 = 0.;
	__attribute__((opencl_unroll_hint(8)))
		for (int k = 1; k < 9; k++)
		{
			double dff1 = (k < 5) ? dff : dff/4.;
			FAA[k] = FA[k - 1] + dff1;
			GAA[k] = GA[k - 1];
			(*Rho) += FAA[k];
			(*RoE2) += GAA[k];
			(*UU) += CXY_DOUBLE[k - 1] * FAA[k];
			(*qq) += CXY_DOUBLE[k - 1] * GAA[k];
			M20 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].x * FAA[k];
			M02 += CXY_DOUBLE[k - 1].y * CXY_DOUBLE[k - 1].y * FAA[k];
			M11 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].y * FAA[k];
		}

	(*UU) /= (*Rho);
	double2 M22eq = (1. / 3. + (*UU)*(*UU));

	(*qq_eq) = ((*RoE2) + 2. * (*Rho) / 3.) * (*UU);

	double3 Pneq = (double3)(M20 - M22eq.x*(*Rho), M02 - M22eq.y*(*Rho), M11 - (*Rho)*(*UU).x*(*UU).y);

	(*qq).x -= 2. * dot((*UU), Pneq.xz);// ((*UU).x*Pneq.x + (*UU).y*Pneq.z);
	(*qq).y -= 2. * dot((*UU), Pneq.zy);// ((*UU).x*Pneq[2] + (*UU).y*Pneq[1]);

	(*Req).xy = ((*RoE2) + 2. * (*Rho) / 3.) * M22eq + 2. * (*Rho) / 3. * (*UU)*(*UU) - (*RoE2) / 3.;
	(*Req).z = ((*RoE2) + 4. * (*Rho) / 3.) * (*UU).x*(*UU).y;
}

inline void calculate_moments_2pop(double *FA, double *FAA, double *GA, double *GAA, double *Rho, double2 *UU, double3 *TNP, double3 *TNP_eq,
	double *RoE2, double2 *qq, double2 *qq_eq, double3 *Req)
{
	(*UU) = 0.;
	(*qq) = 0.;
	double M11 = 0., M20 = 0., M02 = 0.;
	__attribute__((opencl_unroll_hint(8)))
		for (int k = 1; k < 9; k++)
		{
			FAA[k] = FA[k - 1];
			GAA[k] = GA[k - 1];
			(*Rho) += FAA[k];
			(*RoE2) += GAA[k];
			(*UU) += CXY_DOUBLE[k - 1] * FAA[k];
			(*qq) += CXY_DOUBLE[k - 1] * GAA[k];
			M20 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].x * FAA[k];
			M02 += CXY_DOUBLE[k - 1].y * CXY_DOUBLE[k - 1].y * FAA[k];
			M11 += CXY_DOUBLE[k - 1].x * CXY_DOUBLE[k - 1].y * FAA[k];
		}

	double rhoinv = 1. / (*Rho);
	(*TNP) = (double3)(M20 + M02, M20 - M02, M11);

	double2 M22eq = (*Rho) / 3. + (*UU)*(*UU);
	(*TNP_eq) = (double3)(M22eq.x + M22eq.y, M22eq.x - M22eq.y, (*UU).x*(*UU).y);

	(*qq_eq) = ((*RoE2)/(*Rho) + 2. / 3.) * (*UU);

	double3 Pneq = (double3)(M20 - M22eq.x, M02 - M22eq.y, M11 - (*TNP_eq).z);

	(*qq).x -= 2. * dot((*UU), Pneq.xz) / (*Rho);
	(*qq).y -= 2. * dot((*UU), Pneq.zy) / (*Rho);

	(*Req).xy = ((*RoE2)/(*Rho) + 2. / 3.) * M22eq + 2. / 3. * (*UU)*(*UU) - (*RoE2) / 3.;
	(*Req).z = ((*RoE2)/(*Rho) + 4. / 3.) * (*TNP_eq).z;
}

inline void calculate_Feq(double *Feq, double rho, double2 UU)
{
	double2 op3u = sqrt(1. + 3.*UU*UU);
	double Psi = rho * (2. - op3u.x) * (2. - op3u.y) / 9.;
	double Bx = (2. * UU.x + op3u.x) / (1. - UU.x);
	double Bxinv = 1. / Bx;
	double By = (2. * UU.y + op3u.y) / (1. - UU.y);
	double Byinv = 1. / By;

	Feq[0] = 4. * Psi;
	Feq[1] = Psi * Bx;
	Feq[2] = Psi * Bxinv;
	Feq[3] = Psi * By;
	Feq[4] = Psi * Byinv;
	Feq[5] = 0.25 * Psi * Bx * By;
	Feq[6] = 0.25 * Psi * Bxinv * Byinv;
	Feq[7] = 0.25 * Psi * Bx * Byinv;
	Feq[8] = 0.25 * Psi * Bxinv * By;
}

inline void calculate_Fdel(double *Feq, double *FAA, double rho, double2 UU)
{
	double2 op3u = sqrt(1. + 3.*UU*UU);
	double Psi = rho * (2. - op3u.x) * (2. - op3u.y) / 9.;
	double Bx = (2. * UU.x + op3u.x) / (1. - UU.x);
	double Bxinv = 1. / Bx;
	double By = (2. * UU.y + op3u.y) / (1. - UU.y);
	double Byinv = 1. / By;

	Feq[0] = 4. * Psi - FAA[0];
	Feq[1] = Psi * Bx - FAA[1];
	Feq[2] = Psi * Bxinv - FAA[2];
	Feq[3] = Psi * By - FAA[3];
	Feq[4] = Psi * Byinv - FAA[4];
	Feq[5] = 0.25 * Psi * Bx * By - FAA[5];
	Feq[6] = 0.25 * Psi * Bxinv * Byinv - FAA[6];
	Feq[7] = 0.25 * Psi * Bx * Byinv - FAA[7];
	Feq[8] = 0.25 * Psi * Bxinv * By - FAA[8];
}

inline void calculate_Fdel_Fterm(double *Fdel, double *Fterm, double *FAA, double rho, double2 UU)
{
	double2 op3u = sqrt(1. + 3.*UU*UU);
	double dU = UU.x + FTERM_VAL;
	double op3du = sqrt(1. + 3.*dU*dU);

	double Psi = rho * (2. - op3u.x) * (2. - op3u.y) / 9.;
	double Psi_du = rho * (2. - op3du) * (2. - op3u.y) / 9.;

	double Bx = (2. * UU.x + op3u.x) / (1. - UU.x);
	double Bxinv = 1. / Bx;

	double Bx_du = (2. * dU + op3du) / (1. - dU);
	double Bxinv_du = 1. / Bx_du;

	double By = (2. * UU.y + op3u.y) / (1. - UU.y);
	double Byinv = 1. / By;

	double feq = 4. * Psi;
	Fdel[0] = feq - FAA[0];
	Fterm[0] = 4.*Psi_du - feq;


	feq = Psi * Bx;
	Fdel[1] = feq - FAA[1];
	Fterm[1] = Psi_du*Bx_du - feq;

	feq = Psi * Bxinv;
	Fdel[2] = feq - FAA[2];
	Fterm[2] = Psi_du*Bxinv_du - feq;

	feq = Psi * By;
	Fdel[3] = feq - FAA[3];
	Fterm[3] = Psi_du*By - feq;

	feq = Psi * Byinv;
	Fdel[4] = feq - FAA[4];
	Fterm[4] = Psi_du*Byinv - feq;
	
	feq = 0.25 * Psi * Bx * By;
	Fdel[5] = feq - FAA[5];
	Fterm[5] = 0.25*Psi_du*Bx_du*By - feq;

	feq = 0.25 * Psi * Bxinv * Byinv;
	Fdel[6] = feq - FAA[6];
	Fterm[6] = 0.25*Psi_du*Bxinv_du*Byinv - feq;

	feq = 0.25 * Psi * Bx * Byinv;
	Fdel[7] = feq - FAA[7];
	Fterm[7] = 0.25*Psi_du*Bx_du*Byinv - feq;

	feq = 0.25 * Psi * Bxinv * By;
	Fdel[8] = feq - FAA[8];
	Fterm[8] = 0.25*Psi_du*Bxinv_du*By - feq;
}

inline void calc_Geq_Gstar(double *Geq, double *Gstar, double RoE2, double2 qq, double2 qq_eq, double3 Req)
{
	Geq[0] = 4. / 9. * (RoE2 + 4.5 * dot((double2)(-1. / 3., -1. / 3.),Req.xy));
	Geq[1] = 1. / 9. * (RoE2 + qq_eq.x*3. + 4.5 * dot((double2)(2. / 3., -1. / 3.), Req.xy));
	Geq[2] = 1. / 9. * (RoE2 - qq_eq.x*3. + 4.5 * dot((double2)(2. / 3., -1. / 3.), Req.xy));
	Geq[3] = 1. / 9. * (RoE2 + qq_eq.y*3. + 4.5 * dot((double2)(-1. / 3., 2. / 3.), Req.xy));
	Geq[4] = 1. / 9. * (RoE2 - qq_eq.y*3. + 4.5 * dot((double2)(-1. / 3., 2. / 3.), Req.xy));
	Geq[5] = 1. / 36. * (RoE2 + (qq_eq.x + qq_eq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
	Geq[6] = 1. / 36. * (RoE2 - (qq_eq.x + qq_eq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
	Geq[7] = 1. / 36. * (RoE2 + (qq_eq.x - qq_eq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
	Geq[8] = 1. / 36. * (RoE2 + (qq_eq.y - qq_eq.x)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));

	Gstar[0] = 4. / 9. * (RoE2 + 4.5 * dot((double2)(-1. / 3., -1. / 3.), Req.xy));
	Gstar[1] = 1. / 9. * (RoE2 + qq.x*3. + 4.5 * dot((double2)(2. / 3., -1. / 3.), Req.xy));
	Gstar[2] = 1. / 9. * (RoE2 - qq.x*3. + 4.5 * dot((double2)(2. / 3., -1. / 3.), Req.xy));
	Gstar[3] = 1. / 9. * (RoE2 + qq.y*3. + 4.5 * dot((double2)(-1. / 3., 2. / 3.), Req.xy));
	Gstar[4] = 1. / 9. * (RoE2 - qq.y*3. + 4.5 * dot((double2)(-1. / 3., 2. / 3.), Req.xy));
	Gstar[5] = 1. / 36. * (RoE2 + (qq.x + qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
	Gstar[6] = 1. / 36. * (RoE2 - (qq.x + qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
	Gstar[7] = 1. / 36. * (RoE2 + (qq.x - qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
	Gstar[8] = 1. / 36. * (RoE2 + (qq.y - qq.x)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
}

double calculate_Feq_indiv(int ind, double Psi, double Bx, double By, double Bxinv, double Byinv)
{
	switch (ind)
	{
	case 0:
	{
		return 4. * Psi;
	}
	case 1:
	{
		return Psi * Bx;
	}
	case 2:
	{
		return Psi * Bxinv;
	}
	case 3:
	{
		return Psi * By;
	}
	case 4:
	{
		return Psi * Byinv;
	}
	case 5:
	{
		return 0.25 * Psi * Bx * By;
	}
	case 6:
	{
		return 0.25 * Psi * Bxinv * Byinv;
	}
	case 7:
	{
		return 0.25 * Psi * Bx * Byinv;
	}
	case 8:
	{
		return 0.25 * Psi * Bxinv * By;
	}
	}
	return 0.;
}

//Collision and Propagation step (two instances created on host side with 
//FA_buffer = set as FA argument and FB_buffer set as FB argument and vice versa)
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
void LB_collision_ELBM(__global double *Ro_array,
__global double2 *J_array,
__global double *T_array,
__global int *Pos,
__global double *F0,
__global double *FA,
__global double *FB,
__global double *G0,
__global double *GA,
__global double *GB,
__global char *Mred)
{

	int ix = get_global_id(0);
	int iy = get_global_id(1);
	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
		return;

	int jj = ix * FULLSIZEY + iy;
	if (Mred[jj] == 0)
		return;

	int indUT = ix * FULLSIZEY_UT + iy;

	int ii = jj * 8;

	double rho, RoE2;
	double2 UU, qq, qq_eq;
	double3 Req, TNP, TNP_eq;
	double Feq[9], FAA[9], GAA[9], Geq[9], Gstar[9], Sdel[9], Hdel[9];
	FAA[0] = F0[jj];
	GAA[0] = G0[jj];
	rho = FAA[0];
	RoE2 = GAA[0];
	calculate_moments_2pop(&FA[ii], FAA, &GA[ii], GAA, &rho, &UU, &TNP, &TNP_eq, &RoE2, &qq, &qq_eq, &Req);

	Ro_array[indUT] = rho;
	J_array[indUT] = UU/rho;
	T_array[indUT] = (RoE2 / rho - dot(UU,UU)) / 2.;

	calculate_Feq(Feq, rho, UU);
	calc_Sdel(Sdel, rho, TNP, TNP_eq);
	calc_Hdel(FAA, Feq, Sdel, Hdel);
	calc_Geq_Gstar(Geq, Gstar, RoE2, qq, qq_eq, Req);
	double gamma = calc_Gamma(Sdel, Hdel, Feq);

	F0[jj] = FAA[0] - omega_f * (Sdel[0] + 0.5*gamma*Hdel[0]);
	G0[jj] = omega_f*(Geq[0] - GAA[0]) + (omega_f - omega_1f)*(Gstar[0] - Geq[0]) + GAA[0];

	__attribute__((opencl_unroll_hint(8)))
	for (int k = 0; k < 8; k++)
	{
		FB[Pos[ii + k]] = FAA[k + 1] - omega_f * (Sdel[k + 1] + 0.5*gamma*Hdel[k + 1]);
		GB[Pos[ii + k]] = omega_f*(Geq[k + 1] - GAA[k + 1]) + (omega_f - omega_1f)*(Gstar[k + 1] - Geq[k + 1]) + GAA[k + 1];
	}
}

//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
//void LB_collision_SRT(__global double *Ro_array,
//__global double2 *J_array,
//__global double *T_array,
//__global int *Pos,
//__global double *F0,
//__global double *FA,
//__global double *FB,
//__global double *G0,
//__global double *GA,
//__global double *GB,
//__global char *Mred,
//__global double *Alpha_array,
//__constant double *dFF)
//{
//
//	int ix = get_global_id(0);
//	int iy = get_global_id(1);
//	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
//		return;
//
//	int jj = ix * FULLSIZEY + iy;
//	if (Mred[jj] == 0)
//		return;
//
//	int indUT = ix * FULLSIZEY_UT + iy;
//
//	int ii = jj * 8;
//	double rho, RoE2;
//	double2 UU, qq, qq_eq;
//	double3 Req;
//	double Fdel[9], FAA[9], GAA[9], Geq[9], Gstar[9];
//	FAA[0] = F0[jj] + dFF[0]*4.;
//	GAA[0] = G0[jj];
//	rho = FAA[0];
//	RoE2 = GAA[0];
//	calculate_moments_SRT_2pop_dff(&FA[ii], FAA, &GA[ii], GAA, &rho, &UU, &RoE2, &qq, &qq_eq, &Req, dFF[0]);
//
//	Ro_array[indUT] = rho;
//	J_array[indUT] = UU;
//	T_array[indUT] = (RoE2 / rho - dot(UU, UU)) / 2.;
//
//	calculate_Fdel(Fdel, FAA, rho, UU);
//	calc_Geq_Gstar(Geq, Gstar, RoE2, qq, qq_eq, Req);
//	double alpha_temp = Calc_alpha(Alpha_array[jj], FAA, Fdel);
//	Alpha_array[jj] = alpha_temp;
//	alpha_temp *= BETA_VAL;
//	F0[jj] = FAA[0] + alpha_temp * Fdel[0];
//	G0[jj] = omega_f*(Geq[0] - GAA[0]) + (omega_f - omega_1f)*(Gstar[0] - Geq[0]) + GAA[0];
//
//	__attribute__((opencl_unroll_hint(8)))
//		for (int k = 0; k < 8; k++)
//		{
//			FB[Pos[ii + k]] = FAA[k + 1] + alpha_temp * Fdel[k+1];
//			GB[Pos[ii + k]] = omega_f*(Geq[k + 1] - GAA[k + 1]) + (omega_f - omega_1f)*(Gstar[k + 1] - Geq[k + 1]) + GAA[k + 1];
//		}
//}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
void LB_collision_SRT(__global double *Ro_array,
__global double2 *J_array,
__global double *T_array,
__global int *Pos,
__global double *F0,
__global double *FA,
__global double *FB,
__global double *G0,
__global double *GA,
__global double *GB,
__global char *Mred,
__global double *Alpha_array)
{

	int ix = get_global_id(0);
	int iy = get_global_id(1);
	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
		return;

	int jj = ix * FULLSIZEY + iy;
	if (Mred[jj] == 0)
		return;

	int indUT = ix * FULLSIZEY_UT + iy;

	int ii = jj * 8;
	double rho, RoE2;
	double2 UU, qq, qq_eq;
	double3 Req;
	double Fdel[9], Fterm[9], FAA[9], GAA[9], Geq[9], Gstar[9];
	FAA[0] = F0[jj];
	GAA[0] = G0[jj];
	rho = FAA[0];
	RoE2 = GAA[0];
	calculate_moments_SRT_2pop(&FA[ii], FAA, &GA[ii], GAA, &rho, &UU, &RoE2, &qq, &qq_eq, &Req);

	Ro_array[indUT] = rho;
	J_array[indUT] = UU;
	T_array[indUT] = (RoE2 / rho - dot(UU, UU)) / 2.;
	calculate_Fdel_Fterm(Fdel, Fterm, FAA, rho, UU);
	calc_Geq_Gstar(Geq, Gstar, RoE2, qq, qq_eq, Req);
	double alpha_temp = Calc_alpha(Alpha_array[jj], FAA, Fdel);
	Alpha_array[jj] = alpha_temp;
	alpha_temp *= BETA_VAL;
	F0[jj] = FAA[0] + alpha_temp * Fdel[0] + Fterm[0];
	G0[jj] = omega_f*(Geq[0] - GAA[0]) + (omega_f - omega_1f)*(Gstar[0] - Geq[0]) + GAA[0];

	__attribute__((opencl_unroll_hint(8)))
		for (int k = 0; k < 8; k++)
		{
			FB[Pos[ii + k]] = FAA[k + 1] + alpha_temp * Fdel[k + 1] + Fterm[k+1];
			GB[Pos[ii + k]] = omega_f*(Geq[k + 1] - GAA[k + 1]) + (omega_f - omega_1f)*(Gstar[k + 1] - Geq[k + 1]) + GAA[k + 1];
		}
}


//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
//void LB_collision_SRT(__global double *Ro_array,
//__global double2 *J_array,
//__global double *T_array,
//__global int *Pos,
//__global double *F0,
//__global double *FA,
//__global double *FB,
//__global double *G0,
//__global double *GA,
//__global double *GB,
//__global char *Mred,
//__global double *Alpha_array)
//{
//
//	int ix = get_global_id(0);
//	int iy = get_global_id(1);
//	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
//		return;
//
//	int jj = ix * FULLSIZEY + iy;
//	if (Mred[jj] == 0)
//		return;
//
//	int indUT = ix * FULLSIZEY_UT + iy;
//
//	int ii = jj * 8;
//	double rho, RoE2;
//	double2 UU, qq, qq_eq;
//	double3 Req;
//	double Fdel[9], FAA[9], GAA[9], Geq[9], Gstar[9];
//	FAA[0] = F0[jj];
//	GAA[0] = G0[jj];
//	rho = FAA[0];
//	RoE2 = GAA[0];
//	calculate_moments_SRT_2pop(&FA[ii], FAA, &GA[ii], GAA, &rho, &UU, &RoE2, &qq, &qq_eq, &Req);
//
//	Ro_array[indUT] = rho;
//	J_array[indUT] = UU;
//	T_array[indUT] = (RoE2 / rho - dot(UU, UU)) / 2.;
//	calculate_Fdel(Fdel, FAA, rho, UU);
//	calc_Geq_Gstar(Geq, Gstar, RoE2, qq, qq_eq, Req);
//	double alpha_temp = Calc_alpha(Alpha_array[jj], FAA, Fdel);
//	Alpha_array[jj] = alpha_temp;
//	alpha_temp *= BETA_VAL;
//	F0[jj] = FAA[0] + alpha_temp * Fdel[0];
//	G0[jj] = omega_f*(Geq[0] - GAA[0]) + (omega_f - omega_1f)*(Gstar[0] - Geq[0]) + GAA[0];
//
//	__attribute__((opencl_unroll_hint(8)))
//		for (int k = 0; k < 8; k++)
//		{
//			FB[Pos[ii + k]] = FAA[k + 1] + alpha_temp * Fdel[k + 1];
//			GB[Pos[ii + k]] = omega_f*(Geq[k + 1] - GAA[k + 1]) + (omega_f - omega_1f)*(Gstar[k + 1] - Geq[k + 1]) + GAA[k + 1];
//		}
//}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1)))
void LB_IBB_ELBM_No_IBB(__global int2 *ibb_loc,
__global double3 *ibb_coeff,
__global int *node_loc,
__global int *neighs,
__global double2 *dxdy,
__global double *FB,
__global double *GB,
__global double2 *U_array,
__global double *Ro_array,
__global double *T_array,
int max_el)
{
	int i = get_global_id(0);

	if (i >= max_el)
		return;

	int nloc = node_loc[i];
	double2 U_tgt = U_array[nloc];
	double rho_tgt = Ro_array[nloc];

	double T_tgt = 0.;
	int count = 0;
	__attribute__((opencl_unroll_hint(8)))
		for (int k = 0; k < 8; k++)
		{
			if (ibb_loc[i * 8 + k].y != -1)
			{
				T_tgt = (ibb_coeff[i * 8 + k].z*T_array[neighs[i * 8 + k]] + T_array[neighs[i * 8 + RevDir[k]]]) / (1. + ibb_coeff[i * 8 + k].z);
				count++;
			}
		}

	double2 dUdy = (U_array[neighs[i * 8 + 2]] - U_array[neighs[i * 8 + 3]]) / dxdy[i].y;
	double2 dUdx = (U_array[neighs[i * 8]] - U_array[neighs[i * 8 + 1]]) / dxdy[i].x;
	double2 dTdxy = (double2)((T_array[neighs[i * 8 + 2]] - T_array[neighs[i * 8 + 3]]), (T_array[neighs[i * 8]] - T_array[neighs[i * 8 + 1]])) / dxdy[i];

	double3 Pneq = (double3)(dUdx.x, dUdy.y, 0.5*(dUdy.x + dUdx.y)) * 2. * rho_tgt / (-3.*omega_f);

	double2 qq = -2. / omega_1f * rho_tgt * dTdxy / 3. + 2. * (U_tgt.x*Pneq.xz + U_tgt.y*Pneq.zy);
	double3 Rneq = (double3)(2.*dUdx.x, 2.*dUdy.y, (dUdy.x + dUdx.y)) * (2. / 3. + T_tgt);
	double RoE2 = 2.*rho_tgt*T_tgt + rho_tgt*dot(U_tgt, U_tgt);
	double2 qq_eq = (RoE2 + 2. * rho_tgt / 3.) * U_tgt;
	double3 Req = (RoE2 + 4.*rho_tgt / 3.);
	Req.xy = Req.xy * U_tgt*U_tgt + 2.*rho_tgt / 9.;
	Req.z *= U_tgt.x*U_tgt.y;
	Rneq.xy += 2. * U_tgt*dTdxy;
	Rneq.z += dot(U_tgt, dTdxy);
	Rneq *= (-2. / omega_1f / 3. * rho_tgt);


	if (ibb_loc[i * 8].y != -1)
	{
		double gi_eq = (RoE2 + qq_eq.x*3. + 1.5 * dot((double2)(2., -1.), Req.xy));
		double gi_neq = qq.x*3. + 1.5 * dot((double2)(2., -1.), Rneq.xy);
		GB[ibb_loc[i * 8].x] = (gi_eq + gi_neq) / 9.;
	}
	if (ibb_loc[i * 8 + 1].y != -1)
	{
		double gi_eq = (RoE2 - qq_eq.x*3. + 1.5 * dot((double2)(2., -1.), Req.xy));
		double gi_neq = -qq.x*3. + 1.5 * dot((double2)(2., -1.), Rneq.xy);
		GB[ibb_loc[i * 8 + 1].x] = (gi_eq + gi_neq) / 9.;
	}
	if (ibb_loc[i * 8 + 2].y != -1)
	{
		double gi_eq = (RoE2 + qq_eq.y*3. + 1.5 * dot((double2)(-1., 2.), Req.xy));
		double gi_neq = qq.y*3. + 1.5 * dot((double2)(-1., 2.), Rneq.xy);
		GB[ibb_loc[i * 8 + 2].x] = (gi_eq + gi_neq) / 9.;
	}
	if (ibb_loc[i * 8 + 3].y != -1)
	{
		double gi_eq = (RoE2 - qq_eq.y*3. + 1.5 * dot((double2)(-1., 2.), Req.xy));
		double gi_neq = -qq.y*3. + 1.5 * dot((double2)(-1., 2.), Rneq.xy);
		GB[ibb_loc[i * 8 + 3].x] = (gi_eq + gi_neq) / 9.;
	}
	if (ibb_loc[i * 8 + 4].y != -1)
	{
		double gi_eq = (RoE2 + (qq_eq.x + qq_eq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
		double gi_neq = (qq.x + qq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Rneq);
		GB[ibb_loc[i * 8 + 4].x] = (gi_eq + gi_neq) / 36.;
	}
	if (ibb_loc[i * 8 + 5].y != -1)
	{
		double gi_eq = (RoE2 - (qq_eq.x + qq_eq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
		double gi_neq = -(qq.x + qq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Rneq);
		GB[ibb_loc[i * 8 + 5].x] = (gi_eq + gi_neq) / 36.;
	}
	if (ibb_loc[i * 8 + 6].y != -1)
	{
		double gi_eq = (RoE2 + (qq_eq.x - qq_eq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
		double gi_neq = (qq.x - qq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Rneq);
		GB[ibb_loc[i * 8 + 6].x] = (gi_eq + gi_neq) / 36.;
	}
	if (ibb_loc[i * 8 + 7].y != -1)
	{
		double gi_eq = (RoE2 + (qq_eq.y - qq_eq.x) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
		double gi_neq = (qq.y - qq.x) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Rneq);
		GB[ibb_loc[i * 8 + 7].x] = (gi_eq + gi_neq) / 36.;
	}
}

inline void calculate_Geq_variables(double rho, double2 UU, double T, double *RoE2, double2 *qq, double3 *Req)
{
	(*RoE2) = 2.*rho*T + rho*dot(UU,UU);
	(*qq) = ((*RoE2) + 2. * rho / 3.) * UU;
	(*Req).xy = ((*RoE2) + 4. / 3.*rho) * UU*UU + 2. / 9.*rho; 
	(*Req).z = ((*RoE2) + 4./3.*rho) * (UU.x*UU.y);
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_INLET, 1, 1)))
void LB_Outflow_ELBM(__global double2 *U_array,
__global double *Ro_array,
__global double *T_array,
__global double *FB,
__global double *GB)
{
	int i = get_global_id(0);

	if (i >= FULLSIZEY)
		return;

	int UTii = (FULLSIZEX - 1)*(FULLSIZEY_UT)+i;
	int ii = i * 8 + (FULLSIZEX - 1) * 8 * FULLSIZEY;

	double rho = Ro_array[UTii];
	double2 UU = U_array[UTii];
	double TT = T_array[UTii];
	double Psi, Bx, By, Bxinv, Byinv;
	double RoE2;
	double2 qq;
	double3 Req;
	calculate_Feq_variables(rho, UU, &Psi, &Bx, &By, &Bxinv, &Byinv);
	calculate_Geq_variables(rho, UU, TT, &RoE2, &qq, &Req);

	FB[ii + 1] = Psi*Bxinv;
	GB[ii + 1] = 1. / 9. * (RoE2 - qq.x*3. + 3.*Req.x - 1.5*Req.y);

	if (i > 0)
	{
		FB[ii + 7] = 0.25 * Psi * Bxinv * By;
		GB[ii + 7] = 1. / 36. * (RoE2 + 3.*(qq.y - qq.x + Req.x + Req.y) - 9. * Req.z);
	}
	
	if (i < FULLSIZEY - 1)
	{
		FB[ii + 5] = 0.25 * Psi * Bxinv * Byinv;
		GB[ii + 5] = 1. / 36. * (RoE2 - 3.*(qq.x + qq.y - Req.x - Req.y) + 9. * Req.z);
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_INLET, 1, 1)))
void LB_Inlet_ELBM(__global double *Ro_array,
	__global double *F0,
	__global double *FB,
	__global double *G0,
	__global double *GB)
{
	int i = get_global_id(0);

	if (i >= FULLSIZEY)
		return;

	int ii = i * 8;
	double rho = Ro_array[i + FULLSIZEY_UT];
	double2 UU = (double2)(UX_INLET_VAL, 0.);
	double RoE2;
	double2 qq;
	double3 Req;
	calculate_Geq_variables(rho, UU, TFD_X_IN_VAL, &RoE2, &qq, &Req);

	double op3ux = sqrt(1. + 3.*UX_INLET_VAL*UX_INLET_VAL);

	double Psi = rho * (2. - op3ux) / 9.;
	double Bx = (2. * UX_INLET_VAL + op3ux) / (1. - UX_INLET_VAL);
	double Bxinv = 1. / (Bx);
	
	F0[i] = 4. * Psi;
	FB[ii] = Psi * Bx;
	FB[ii+1] = Psi * Bxinv;
	FB[ii+2] = Psi;
	FB[ii+3] = Psi;
	FB[ii+4] = 0.25 * Psi * Bx;
	FB[ii+5] = 0.25 * Psi * Bxinv;
	FB[ii+6] = 0.25 * Psi * Bx;
	FB[ii+7] = 0.25 * Psi * Bxinv;

	G0[i] = 4. / 9. * (RoE2 + 4.5 * dot((double2)(-1. / 3., -1. / 3.), Req.xy));
	GB[ii] = 1. / 9. * (RoE2 + qq.x*3. + 4.5 * dot((double2)(2. / 3., -1. / 3.), Req.xy));
	GB[ii+1] = 1. / 9. * (RoE2 - qq.x*3. + 4.5 * dot((double2)(2. / 3., -1. / 3.), Req.xy));
	GB[ii+2] = 1. / 9. * (RoE2 + qq.y*3. + 4.5 * dot((double2)(-1. / 3., 2. / 3.), Req.xy));
	GB[ii+3] = 1. / 9. * (RoE2 - qq.y*3. + 4.5 * dot((double2)(-1. / 3., 2. / 3.), Req.xy));
	GB[ii+4] = 1. / 36. * (RoE2 + (qq.x + qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
	GB[ii+5] = 1. / 36. * (RoE2 - (qq.x + qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
	GB[ii+6] = 1. / 36. * (RoE2 + (qq.x - qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
	GB[ii+7] = 1. / 36. * (RoE2 + (qq.y - qq.x)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1)))
void LB_IBB_ELBM(__global int2 *ibb_loc,
__global double3 *ibb_coeff,
__global int *node_loc,
__global int *neighs,
__global double2 *dxdy,
__global double *FB,
__global double *GB,
__global double2 *U_array,
__global double *Ro_array,
__global double *T_array,
int max_el)
{
	int i = get_global_id(0);

	if (i >= max_el)
		return;

	int nloc = node_loc[i];
	double2 U_tgt = U_array[nloc];
	double rho_tgt = Ro_array[nloc];
	
	double T_tgt = 0.;
	int count = 0;
	__attribute__((opencl_unroll_hint(8)))
	for (int k = 0; k < 8; k++)
	{
		if (ibb_loc[i * 8 + k].y != -1)
		{
			double Fnew = ibb_coeff[i * 8 + k].x * FB[ibb_loc[i * 8 + k].x] + ibb_coeff[i * 8 + k].y * FB[ibb_loc[i * 8 + k].y];
			FB[ibb_loc[i * 8 + k].x] = Fnew;
			T_tgt = (ibb_coeff[i * 8 + k].z*T_array[neighs[i*8 + k]] + T_array[neighs[i*8 + RevDir[k]]]) / (1. + ibb_coeff[i * 8 + k].z);
			count++;
		}
	}

	double2 dUdy = (U_array[neighs[i*8+2]] - U_array[neighs[i*8+3]])/dxdy[i].y;
	double2 dUdx = (U_array[neighs[i*8]] - U_array[neighs[i*8+1]])/dxdy[i].x;
	double2 dTdxy = (double2)((T_array[neighs[i*8+2]] - T_array[neighs[i*8+3]]),(T_array[neighs[i*8]] - T_array[neighs[i*8+1]]))/dxdy[i];

	double3 Pneq = (double3) (dUdx.x, dUdy.y, 0.5*(dUdy.x + dUdx.y)) * 2. * rho_tgt / (-3.*omega_f);
	
	double2 qq = -2./omega_1f * rho_tgt * dTdxy / 3. + 2. * (U_tgt.x*Pneq.xz + U_tgt.y*Pneq.zy);
	double3 Rneq = (double3)(2.*dUdx.x, 2.*dUdy.y, (dUdy.x + dUdx.y)) * (2./3. + T_tgt);
	double RoE2 =  2.*rho_tgt*T_tgt + rho_tgt*dot(U_tgt,U_tgt);
	double2 qq_eq = (RoE2 + 2. * rho_tgt / 3.) * U_tgt;
	double3 Req = (RoE2 + 4.*rho_tgt/3.);
	Req.xy = Req.xy * U_tgt*U_tgt + 2.*rho_tgt / 9.;
	Req.z *= U_tgt.x*U_tgt.y; 
	Rneq.xy += 2. * U_tgt*dTdxy;
	Rneq.z += dot(U_tgt,dTdxy);
	Rneq *= (-2. / omega_1f / 3. * rho_tgt);


	if (ibb_loc[i * 8].y != -1)
	{
		double gi_eq = (RoE2 + qq_eq.x*3. + 1.5 * dot((double2)(2.,-1.), Req.xy));
		double gi_neq = qq.x*3. + 1.5 * dot((double2)(2.,-1.), Rneq.xy);
		GB[ibb_loc[i * 8].x] = (gi_eq + gi_neq)/9.;
	}
	if (ibb_loc[i * 8 + 1].y != -1)
	{
		double gi_eq = (RoE2 - qq_eq.x*3. + 1.5 * dot((double2)(2., -1.), Req.xy));
		double gi_neq = -qq.x*3. + 1.5 * dot((double2)(2.,-1.), Rneq.xy);
		GB[ibb_loc[i * 8+1].x] = (gi_eq + gi_neq)/9.;
	}
	if (ibb_loc[i * 8 + 2].y != -1)
	{
		double gi_eq = (RoE2 + qq_eq.y*3. + 1.5 * dot((double2)(-1., 2.), Req.xy));
		double gi_neq = qq.y*3. + 1.5 * dot((double2)(-1.,2.), Rneq.xy);
		GB[ibb_loc[i * 8+2].x] = (gi_eq + gi_neq)/9.;
	}
	if (ibb_loc[i * 8 + 3].y != -1)
	{
		double gi_eq = (RoE2 - qq_eq.y*3. + 1.5 * dot((double2)(-1., 2.), Req.xy));
		double gi_neq = -qq.y*3. + 1.5 * dot((double2)(-1.,2.), Rneq.xy);
		GB[ibb_loc[i * 8+3].x] = (gi_eq + gi_neq)/9.;
	}
	if (ibb_loc[i * 8 + 4].y != -1)
	{
		double gi_eq = (RoE2 + (qq_eq.x + qq_eq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
		double gi_neq = (qq.x + qq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Rneq);
		GB[ibb_loc[i * 8+4].x] = (gi_eq + gi_neq)/36.;
	}
	if (ibb_loc[i * 8 + 5].y != -1)
	{
		double gi_eq = (RoE2 - (qq_eq.x + qq_eq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
		double gi_neq = -(qq.x + qq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Rneq);
		GB[ibb_loc[i * 8+5].x] = (gi_eq + gi_neq)/36.;
	}
	if (ibb_loc[i * 8 + 6].y != -1)
	{
		double gi_eq = (RoE2 + (qq_eq.x - qq_eq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
		double gi_neq = (qq.x - qq.y) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Rneq);
		GB[ibb_loc[i * 8+6].x] = (gi_eq + gi_neq)/36.;
	}
	if (ibb_loc[i * 8 + 7].y != -1)
	{
		double gi_eq = (RoE2 + (qq_eq.y - qq_eq.x) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
		double gi_neq = (qq.y - qq.x) * 3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Rneq);
		GB[ibb_loc[i * 8+7].x] = (gi_eq + gi_neq)/36.;
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_INLET, 1, 1)))
void TLB_Inlet_Outlet(__global double2 *U_array,
__global double *Ro_array,
__global double *G0,
__global double *GA,
__global double *GB)
{
	int i = get_global_id(0);

	if (i >= FULLSIZEY)
		return;

	int ii = i * 8;
	double rho = Ro_array[i];
	double2 UU = U_array[i];
	double RoE2;
	double2 qq;
	double3 Req;
	calculate_Geq_variables(rho, UU, TFD_X_IN_VAL, &RoE2, &qq, &Req);

	G0[i] = 4. / 9. * (RoE2 + 4.5 * dot((double2)(-1. / 3., -1. / 3.), Req.xy));
	GB[ii] = 1. / 9. * (RoE2 + qq.x*3. + 4.5 * dot((double2)(2. / 3., -1. / 3.), Req.xy));
	GB[ii + 1] = 1. / 9. * (RoE2 - qq.x*3. + 4.5 * dot((double2)(2. / 3., -1. / 3.), Req.xy));
	GB[ii + 2] = 1. / 9. * (RoE2 + qq.y*3. + 4.5 * dot((double2)(-1. / 3., 2. / 3.), Req.xy));
	GB[ii + 3] = 1. / 9. * (RoE2 - qq.y*3. + 4.5 * dot((double2)(-1. / 3., 2. / 3.), Req.xy));
	GB[ii + 4] = 1. / 36. * (RoE2 + (qq.x + qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
	GB[ii + 5] = 1. / 36. * (RoE2 - (qq.x + qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., 2.), Req));
	GB[ii + 6] = 1. / 36. * (RoE2 + (qq.x - qq.y)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));
	GB[ii + 7] = 1. / 36. * (RoE2 + (qq.y - qq.x)*3. + 4.5 * dot((double3)(2. / 3., 2. / 3., -2.), Req));

	ii = 8 * i + 8 * (FULLSIZEX - 1)*FULLSIZEY;

	GB[ii + 1] = GA[ii + 1];
	if (i > 0)
	{
		GB[ii + 7] = GA[ii + 7];
	}
	if (i < FULLSIZEY - 1)
	{
		GB[ii + 5] = GA[ii + 5];
	}
}

//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_INLET, 1, 1)))
//void LB_Outflow(__global int *loc,
//__global double *FA,
//__global double *FB,
//__global double *GA,
//__global double *GB,
//int max_el)
//{
//	int ii = get_global_id(0);
//	if (ii >= max_el)
//		return;
//	int i = loc[ii];
//	GB[i] = GA[i];
//	FB[i] = FA[i];
//}

//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_INLET, 1, 1)))
//void TLB_Inlet_Outlet(__global int *loc,
//__global double *FA,
//__global double *FB,
//__global double *GA,
//__global double *GB,
//int max_el)
//{
//	int ii = get_global_id(0);
//	if (ii >= max_el)
//		return;
//	int i = loc[ii];
//	GB[i] = GA[i];
//	FB[i] = FA[i];
//}






//Collision and Propagation step (two instances created on host side with 
//FA_buffer = set as FA argument and FB_buffer set as FB argument and vice versa)
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
void LB_collision(__global double *Ro_array,
__global double2 *J_array,
__global int8 *Pos,
__global double *F0,
__global double4 *FA,
__global double *FB,
double Lamda,
double Fbody,
__constant double *dFF,
__global char *Mred)
{

	int ix = get_global_id(0);
	int iy = get_global_id(1);
	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
		return;

	int jj = ix * FULLSIZEY + iy;
	if (Mred[jj] == 0)
		return;

	int indUT = ix * FULLSIZEY_UT + iy;

	int8 Ploc = Pos[jj];
	int ii = jj * 2;


#ifdef OPENCL_VERSION_1_2
	int lid = get_local_id(0) * WORKGROUPSIZEY_LB + get_local_id(1);
#else
	int lid = get_local_linear_id();
#endif

	__local double PPy[WORKGROUPSIZEX_LB*WORKGROUPSIZEY_LB];
	__local double PPx[WORKGROUPSIZEX_LB*WORKGROUPSIZEY_LB];
	__local double PPxy[WORKGROUPSIZEX_LB*WORKGROUPSIZEY_LB];
	__local double4 FAlo[WORKGROUPSIZEX_LB*WORKGROUPSIZEY_LB];

	FAlo[lid] = FA[ii + 1] + dFF[0];
	double Ro = dot(FAlo[lid], (double4)(1.));
	PPxy[lid] = dot(RoUUPF_constxy_hi, FAlo[lid]);
	double Jx = dot(CCxhi, FAlo[lid]);
	double Jy = dot(CCyhi, FAlo[lid]);

	FAlo[lid] = FA[ii] + dFF[0] / 4.;
	PPx[lid] = dot(RoUUPF_constxy_hi.lo, FAlo[lid].lo) + Ro;
	PPy[lid] = dot(RoUUPF_constxy_hi.lo, FAlo[lid].hi) + Ro;
	Ro += F0[jj] + dFF[0] * 4. + dot(FAlo[lid], (double4)(1.));

	Jx += dot(CCxlo.lo, FAlo[lid].lo);
	Jy += dot(CCylo.hi, FAlo[lid].hi);

	Ro_array[indUT] = Ro;
	J_array[indUT].x = Jx / Ro;
	J_array[indUT].y = Jy / Ro;

	double Jx2 = Jx * Jx / Ro, Jy2 = Jy * Jy / Ro;
	double RoUUPFx = (1.5 - Lamda)*Jx2 + Jy2*Lamda / 2. + (PPx[lid] - PPy[lid] / 2. - Ro / 6.) * Lamda + 3.*Fbody * Jx / Ro;
	double RoUUPFy = (1.5 - Lamda)*Jy2 + Jx2*Lamda / 2. + (PPy[lid] - PPx[lid] / 2. - Ro / 6.) * Lamda;
	double RoUUPFz = 1.5 / Ro * (Fbody * Jy + Jx*Jy*(1. - Lamda) + Lamda * PPxy[lid] * Ro);

	Jx = Jx + Fbody;
	Jx = Jx*3.;
	Jy = Jy*3.;
	
	F0[jj] = (Ro - (RoUUPFx + RoUUPFy)) * 4. / 9.;

	FAlo[lid] = (CCxlo * Jx) + (CCylo * Jy) + (RoUUPF_constx_lo * RoUUPFx) + (RoUUPF_constx_lo.zwyx * RoUUPFy) + Ro;
	FAlo[lid] /= 9.;

	FB[Ploc.s0] = (FAlo[lid].s0);
	FB[Ploc.s1] = (FAlo[lid].s1);
	FB[Ploc.s2] = (FAlo[lid].s2);
	FB[Ploc.s3] = (FAlo[lid].s3);

	FAlo[lid] = (CCxhi * Jx) + (CCyhi * Jy) + 2. * (RoUUPFy + RoUUPFx) + (RoUUPF_constxy_hi * RoUUPFz)*6. + Ro;
	FAlo[lid] /= 36.;

	FB[Ploc.s4] = (FAlo[lid].s0);
	FB[Ploc.s5] = (FAlo[lid].s1);
	FB[Ploc.s6] = (FAlo[lid].s2);
	FB[Ploc.s7] = (FAlo[lid].s3);
}

////Collision and Propagation step (two instances created on host side with 
////FA_buffer = set as FA argument and FB_buffer set as FB argument and vice versa)
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1))) 
//void LB_collision(__global double *Ro_array,
//__global double2 *J_array,
//__global int *Pos,
//__global double *F0,
//__global double4 *FA,
//__global double *FB,
//double Lamda,
//double Fbody,
//__global double *dFF,
//__global char *Mred,
//__local double3 *PP, 
//__local double4 *FAlo)
//{
//
//	int ix = get_global_id(0);
//	int iy = get_global_id(1);
//	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
//		return;
//
//	int jj = ix * FULLSIZEY + iy;
//	if (Mred[jj] == 0)
//		return;
//
//	int indUT = ix * FULLSIZEY_UT + iy;
//
//	int lid = get_local_linear_id();
//	
//	int ii = jj * 2;
//
//
//	FAlo[lid] = FA[ii + 1] + dFF[0]/4.;
//	double Ro = dot(FAlo[lid], (double4)(1.));
//	PP[lid].z = dot(RoUUPF_constxy_hi, FAlo[lid]);
//	double Jx = dot(CCxhi, FAlo[lid]);
//	double Jy = dot(CCyhi, FAlo[lid]);
//
//
//	FAlo[lid] = FA[ii] + dFF[0];
//
//	PP[lid].x = dot(RoUUPF_constxy_hi.lo, FAlo[lid].lo) + Ro;
//	PP[lid].y = dot(RoUUPF_constxy_hi.lo, FAlo[lid].hi) + Ro;
//
//	Ro += F0[jj] + dFF[0] * 4. + dot(FAlo[lid], (double4)(1.));
//
//
//	Jx += dot(CCxlo.lo, FAlo[lid].lo);
//	Jy += dot(CCylo.hi, FAlo[lid].hi);
//
//	Ro_array[indUT] = Ro;
//	J_array[indUT].x = Jx / Ro;
//	J_array[indUT].y = Jy / Ro;
//
//	double Jx2 = Jx * Jx / Ro, Jy2 = Jy * Jy / Ro;
//	double newtemp = (PP[lid].x - PP[lid].y / 2. - Ro / 6.) * Lamda;
//	PP[lid].y = (1.5 - Lamda)*Jy2 + Jx2*Lamda / 2. + (PP[lid].y - PP[lid].x / 2. - Ro / 6.) * Lamda;
//	PP[lid].x = (1.5 - Lamda)*Jx2 + Jy2*Lamda / 2. + newtemp + 3.*Fbody * Jx / Ro;
//	PP[lid].z = 1.5 / Ro * (Fbody * Jy + Jx*Jy*(1. - Lamda) + Lamda * PP[lid].z * Ro);
//
//
//	Jx += Fbody;
//	Jx *= 3.;
//	Jy *= 3.;
//
//
//	int Ploc0 = Pos[jj * 8 + 0];
//	int Ploc1 = Pos[jj * 8 + 1];
//	int Ploc2 = Pos[jj * 8 + 2];
//	int Ploc3 = Pos[jj * 8 + 3];
//	int Ploc4 = Pos[jj * 8 + 4];
//	int Ploc5 = Pos[jj * 8 + 5];
//	int Ploc6 = Pos[jj * 8 + 6];
//	int Ploc7 = Pos[jj * 8 + 7];
//
//	F0[jj] = (Ro - (PP[lid].x + PP[lid].y)) * 4. / 9.;
//
//
//	FAlo[lid] = (CCxlo * Jx) + (CCylo * Jy) + (RoUUPF_constx_lo * PP[lid].x) + (RoUUPF_constx_lo.zwyx * PP[lid].y) + Ro;
//
//	FAlo[lid] /= 9.;
//
//	FB[Ploc0] = (FAlo[lid].s0);
//	FB[Ploc1] = (FAlo[lid].s1);
//	FB[Ploc2] = (FAlo[lid].s2);
//	FB[Ploc3] = (FAlo[lid].s3);
//
//	
//	FAlo[lid] = (CCxhi * Jx) + (CCyhi * Jy) + 2. * (PP[lid].x + PP[lid].y) + (RoUUPF_constxy_hi * PP[lid].z)*6. + Ro;
//	FAlo[lid] /= 36.;
//
//	FB[Ploc4] = (FAlo[lid].s0);
//	FB[Ploc5] = (FAlo[lid].s1);
//	FB[Ploc6] = (FAlo[lid].s2);
//	FB[Ploc7] = (FAlo[lid].s3);
//}



//Interpolated bounce back kernel
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1))) 
void LB_IBB(__global int2 *ibb_loc,
	__global double2 *ibb_coeff,
	__global double *FB,
	int max_el)
{
	int i = get_global_id(0);

	if (i >= max_el)
		return;

	double Fold = FB[ibb_loc[i].x];
	double Fnew = ibb_coeff[i].x * Fold + ibb_coeff[i].y * FB[ibb_loc[i].y];
	FB[ibb_loc[i].x] = Fnew;
}

//
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_INLET, 1, 1)))
//void LB_Inlet_Outlet(__global double *F0,
//__global double *FA,
//__global double *FB,
//__global double *Uvals,
//__global double *Ro_array,
//__global double2 *J_array)
//{
//	int i = get_global_id(0);
//
//	if (i >= FULLSIZEY)
//		return;
//
//	int ii = i * 8;
//	int nXm = ((FULLSIZEX - 1)*FULLSIZEY + i) * 8;
//	int nXm1 = ((FULLSIZEX - 2)*FULLSIZEY + i) * 8;
//	int nXm2 = ((FULLSIZEX - 3)*FULLSIZEY + i) * 8;
//
//	double U0 = Uvals[i];
//
//	FB[ii] += U0 / 9.;
//	FB[ii + 4] += (i > 0) ? (U0 / 36.) : (0.);
//	FB[ii + 6] += (i < FULLSIZEY - 1) ? (U0 / 36.) : (0.);
//
//	double2 Uval = J_array[(FULLSIZEX - 1)*FULLSIZEY_UT + i];
//	double Roval = Ro_array[(FULLSIZEX - 1)*FULLSIZEY_UT + i];
//	Uval *= Roval;
//	double Feq = 1. - Uval.x * 3. + 4.5*Uval.x*Uval.x - 1.5 * dot(Uval, Uval);
//	double feq = Feq - 1. + Roval + 6.*Uval.x;
//	double fneq = FA[nXm] - feq / 9.;
//
//	FB[nXm + 1] = Feq / 9. + fneq;
//	FB[nXm + 5] = 2.*FB[nXm1 + 5] - FB[nXm2 + 5];
//	FB[nXm + 7] = 2.*FB[nXm1 + 7] - FB[nXm2 + 7];
//
//}
//

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_INLET, 1, 1)))
void LB_Inlet_Outlet(__global double *F0,
__global double *FA,
__global double *FB,
__global double *Uvals,
__global double *Ro_array,
__global double2 *J_array,
double tau_multiplier)
{
	int i = get_global_id(0);

	if (i >= FULLSIZEY)
		return;

	int ii = i * 8;
	int nXm = ((FULLSIZEX - 1)*FULLSIZEY + i) * 8;

	double U0 = Uvals[i];

	FB[ii] += U0 / 9.;
	FB[ii + 4] += (i > 0) ? (U0 / 36.) : (0.);
	FB[ii + 6] += (i < FULLSIZEY - 1) ? (U0 / 36.) : (0.);

	double2 Uval = J_array[(FULLSIZEX - 1)*FULLSIZEY_UT + i];
	double Roval = Ro_array[(FULLSIZEX - 1)*FULLSIZEY_UT + i];
	double Feq = 1. - Uval.x * 3. + 4.5*Uval.x*Uval.x - 1.5 * dot(Uval, Uval);
	Uval *= Roval;
	double feq = Feq - 1. + Roval + 6.*Uval.x;
	double fneq = FA[nXm] - feq / 9.;

	FB[nXm + 1] = Feq / 9. + tau_multiplier * fneq;

	double Uy = J_array[(FULLSIZEX - 2)*FULLSIZEY_UT + i].y * Roval;
	if (i > 0)
	{
		FB[nXm + 7] = FB[nXm + 6] - (Uval.x - Uy) / 6.;
	}
	if (i < FULLSIZEY - 1)
	{
		FB[nXm + 5] = FB[nXm + 4] - (Uval.x + Uy) / 6.;
	}
}

