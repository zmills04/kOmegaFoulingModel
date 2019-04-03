
//__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
//void LB_collision_SRT_Fluid(__global double *Ro_array,
//__global double2 *J_array,
//__global int8 *Pos,
//__global float *F0,
//__global float *FA,
//__global float *FB,
//__global char *Mred,
//__global float *Alpha_array)
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
//#ifdef OPENCL_VERSION_1_2
//	short lid = get_local_id(0) * WORKGROUPSIZEY_LB + get_local_id(1);
//#else
//	short lid = get_local_linear_id();
//#endif
//	short lid9 = lid * 9;
//	int indUT = ix * FULLSIZEY_UT + iy;
//
//	float2 UU;
//	float rho;
//
//	__local float Fterm[WORKGROUPSIZEX_LB*WORKGROUPSIZEY_LB * 9];
//
//	float FAA[9];
//	FAA[0] = F0[jj];
//	FAA[1] = FA[jj];
//	FAA[2] = FA[jj + FULLSIZEXY];
//	FAA[3] = FA[jj + 2*FULLSIZEXY];
//	FAA[4] = FA[jj + 3*FULLSIZEXY];
//	FAA[5] = FA[jj + 4*FULLSIZEXY];
//	FAA[6] = FA[jj + 5*FULLSIZEXY];
//	FAA[7] = FA[jj + 6*FULLSIZEXY];
//	FAA[8] = FA[jj + 7*FULLSIZEXY];
//
//	UU.x = FAA[1] - FAA[2] + FAA[5] - FAA[6] + FAA[7] - FAA[8];
//	UU.y = FAA[3] - FAA[4] + FAA[5] - FAA[6] - FAA[7] + FAA[8];
//	rho = FAA[0] + FAA[1] + FAA[2] + FAA[3] + FAA[4] + FAA[5] + FAA[6] + FAA[7] + FAA[8];
//
//	UU /= rho;
//
//	Ro_array[indUT] = convert_double(rho);
//	J_array[indUT] = convert_double2(UU);
//
//	float2 op3u = sqrt(1.f + 3.f*UU*UU);
//	float Psi = rho * (2.f - op3u.x)*(2.f - op3u.y) / 9.f;
//	float2 B = (2.f * UU + op3u) / (1.f - UU);
//	float2 Binv = 1.f / B;
//
//	op3u.x = sqrt(1.f + 3.f*(UU.x + FTERM_VAL)*(UU.x + FTERM_VAL));
//	float Psi_du = rho * (2.f - op3u.x)*(2.f - op3u.y) / 9.f;
//	float Bx_du = (2.f * (UU.x + FTERM_VAL) + op3u.x) / (1.f - (UU.x + FTERM_VAL));
//	float Bxinv_du = 1.f / Bx_du;
//	
//	float Fdel[9];
//	
//	float feq = 4.f * Psi;
//	Fdel[0] = feq - FAA[0];
//	Fterm[lid9] = fma(4.f, Psi_du, -feq);
//
//	feq = Psi * B.x;
//	Fdel[1] = feq - FAA[1];
//	Fterm[lid9 + 1] = fma(Psi_du, Bx_du, -feq);
//
//	feq = Psi * Binv.x;
//	Fdel[2] = feq - FAA[2];
//	Fterm[lid9 + 2] = fma(Psi_du, Bxinv_du, -feq);
//
//	feq = Psi * B.y;
//	Fdel[3] = feq - FAA[3];
//	Fterm[lid9 + 3] = fma(Psi_du, B.y, -feq);
//
//	feq = Psi * Binv.y;
//	Fdel[4] = feq - FAA[4];
//	Fterm[lid9 + 4] = fma(Psi_du, Binv.y, -feq);
//
//	Psi_du *= 0.25f;
//	Psi *= 0.25f;
//
//	feq = Psi * B.x * B.y;
//	Fdel[5] = feq - FAA[5];
//	Fterm[lid9 + 5] = fma(Psi_du, Bx_du * B.y, -feq);
//
//	feq = Psi * Binv.x * Binv.y;
//	Fdel[6] = feq - FAA[6];
//	Fterm[lid9 + 6] = fma(Psi_du, Bxinv_du * Binv.y, -feq);
//
//	feq = Psi * B.x * Binv.y;
//	Fdel[7] = feq - FAA[7];
//	Fterm[lid9 + 7] = fma(Psi_du, Bx_du * Binv.y, -feq);
//
//	feq = Psi * Binv.x * B.y;
//	Fdel[8] = feq - FAA[8];
//	Fterm[lid9 + 8] = fma(Psi_du, Bxinv_du * B.y, -feq);
//
//
//	float alpha_temp = Calc_alpha(Alpha_array[jj], FAA, Fdel);
//	Alpha_array[jj] = alpha_temp;
//	alpha_temp *= BETA_VAL_FLOAT;
//
//	int8 ii = Pos[jj];
//
//	F0[jj] += alpha_temp * Fdel[0] + Fterm[lid9];
//	FB[ii.s0] = FAA[1] + alpha_temp * Fdel[1] + Fterm[lid9 + 1];
//	FB[ii.s1] = FAA[2] + alpha_temp * Fdel[2] + Fterm[lid9 + 2];
//	FB[ii.s2] = FAA[3] + alpha_temp * Fdel[3] + Fterm[lid9 + 3];
//	FB[ii.s3] = FAA[4] + alpha_temp * Fdel[4] + Fterm[lid9 + 4];
//	FB[ii.s4] = FAA[5] + alpha_temp * Fdel[5] + Fterm[lid9 + 5];
//	FB[ii.s5] = FAA[6] + alpha_temp * Fdel[6] + Fterm[lid9 + 6];
//	FB[ii.s6] = FAA[7] + alpha_temp * Fdel[7] + Fterm[lid9 + 7];
//	FB[ii.s7] = FAA[8] + alpha_temp * Fdel[8] + Fterm[lid9 + 8];
//}










inline float CalcDeviation(float *FAA, float* Fneq)
{
	float dev = Fneq[0] / FAA[0];
	if (dev > 0.01f)
		return dev;
	dev = fmax(dev, Fneq[1] / FAA[1]);
	if (dev > 0.01f)
		return dev;
	dev = fmax(dev, Fneq[2] / FAA[2]);
	if (dev > 0.01f)
		return dev;
	dev = fmax(dev, Fneq[3] / FAA[3]);
	if (dev > 0.01f)
		return dev;
	dev = fmax(dev, Fneq[4] / FAA[4]);
	if (dev > 0.01f)
		return dev;
	dev = fmax(dev, Fneq[5] / FAA[5]);
	if (dev > 0.01f)
		return dev;
	dev = fmax(dev, Fneq[6] / FAA[6]);
	if (dev > 0.01f)
		return dev;
	dev = fmax(dev, Fneq[7] / FAA[7]);
	if (dev > 0.01f)
		return dev;
	dev = fmax(dev, Fneq[8] / FAA[8]);
	return dev;
}

inline float Calc_EntropyIneq(float *FAA, float* Fneq, float alp, float *dent)
{
	float t = fma(alp, Fneq[0], FAA[0]);
	float h = log(t) - LOG49;
	float ent = t*h;
	*dent += Fneq[0] * (h + 1.f);

	t = fma(alp, Fneq[1], FAA[1]);
	h = log(t) - LOG19;
	ent += t*h;
	*dent += Fneq[1] * (h + 1.f);

	t = fma(alp, Fneq[2], FAA[2]);
	h = log(t) - LOG19;
	ent += t*h;
	*dent += Fneq[2] * (h + 1.f);

	t = fma(alp, Fneq[3], FAA[3]);
	h = log(t) - LOG19;
	ent += t*h;
	*dent += Fneq[3] * (h + 1.f);

	t = fma(alp, Fneq[4], FAA[4]);
	h = log(t) - LOG19;
	ent += t*h;
	*dent += Fneq[4] * (h + 1.f);

	t = fma(alp, Fneq[5], FAA[5]);
	h = log(t) - LOG136;
	ent += t*h;
	*dent += Fneq[5] * (h + 1.f);

	t = fma(alp, Fneq[6], FAA[6]);
	h = log(t) - LOG136;
	ent += t*h;
	*dent += Fneq[6] * (h + 1.f);

	t = fma(alp, Fneq[7], FAA[7]);
	h = log(t) - LOG136;
	ent += t*h;
	*dent += Fneq[7] * (h + 1.f);

	t = fma(alp, Fneq[8], FAA[8]);
	h = log(t) - LOG136;
	ent += t*h;
	*dent += Fneq[8] * (h + 1.f);

	return ent;
}

inline float Calc_Entropy(float *FAA)
{
	float ent = FAA[0] * (log(FAA[0]) - LOG49);
	ent += FAA[1] * (log(FAA[1]) - LOG19);
	ent += FAA[2] * (log(FAA[2]) - LOG19);
	ent += FAA[3] * (log(FAA[3]) - LOG19);
	ent += FAA[4] * (log(FAA[4]) - LOG19);
	ent += FAA[5] * (log(FAA[5]) - LOG136);
	ent += FAA[6] * (log(FAA[6]) - LOG136);
	ent += FAA[7] * (log(FAA[7]) - LOG136);
	ent += FAA[8] * (log(FAA[8]) - LOG136);
	return ent;
}

inline float find_max_alpha(float *FAA, float* Fneq)
{
	float max_alpha = 1000.f;
	if (FAA[0] < 0. || Fneq[0] < 0.)
		max_alpha = fmin(max_alpha, -FAA[0] / Fneq[0]);
	if (FAA[1] < 0. || Fneq[1] < 0.)
		max_alpha = fmin(max_alpha, -FAA[1] / Fneq[1]);
	if (FAA[2] < 0. || Fneq[2] < 0.)
		max_alpha = fmin(max_alpha, -FAA[2] / Fneq[2]);
	if (FAA[3] < 0. || Fneq[3] < 0.)
		max_alpha = fmin(max_alpha, -FAA[3] / Fneq[3]);
	if (FAA[4] < 0. || Fneq[4] < 0.)
		max_alpha = fmin(max_alpha, -FAA[4] / Fneq[4]);
	if (FAA[5] < 0. || Fneq[5] < 0.)
		max_alpha = fmin(max_alpha, -FAA[5] / Fneq[5]);
	if (FAA[6] < 0. || Fneq[6] < 0.)
		max_alpha = fmin(max_alpha, -FAA[6] / Fneq[6]);
	if (FAA[7] < 0. || Fneq[7] < 0.)
		max_alpha = fmin(max_alpha, -FAA[7] / Fneq[7]);
	if (FAA[8] < 0. || Fneq[8] < 0.)
		max_alpha = fmin(max_alpha, -FAA[8] / Fneq[8]);
	return max_alpha;
}

float Calc_Alpha_Root(float *FAA, float* Fneq, float alp)
{
	float ent = Calc_Entropy(FAA);
	int ii = 0;

	float max_alpha = find_max_alpha(FAA, Fneq);

	while (1)
	{
		float delta_ent_der = 0.f;
		float ent_ineq = Calc_EntropyIneq(FAA, Fneq, alp, &delta_ent_der);

		if (isnan(ent_ineq) && alp != 1.1f)
		{
			alp = 1.1f;
			continue;
		}

		float ent_inc = ent_ineq - ent;

		if (fabs(ent_inc) < entropy_tolerance)
			break;

		float new_alpha = alp - ent_inc / delta_ent_der;

		if (new_alpha >= max_alpha)
			new_alpha = (alp + max_alpha) * 0.5f;

		if (fabs(new_alpha - alp) < alpha_tolerance)
			break;

		alp = new_alpha;
		ii++;
		if (ii > 1000)
		{
			return 2.f;
		}
	}

	if (alp < 1.0 || !isfinite(alp))
	{
		return 2.f;
	}
	return alp;
}

float EstimateAlpha(float *FAA, float *Fneq)
{
	float t = Fneq[0];
	float inv = 1.f / FAA[0];
	float p = t*t*inv;
	t *= inv;
	float a1 = p;
	p *= t;
	float a2 = p;
	p *= t;
	float a3 = p;
	p *= t;
	float a4 = p;

	t = Fneq[1];
	inv = 1.f / FAA[1];
	p = t*t*inv;
	t *= inv;
	a1 += p;
	p *= t;
	a2 += p;
	p *= t;
	a3 += p;
	p *= t;
	a4 += p;

	t = Fneq[2];
	inv = 1.f / FAA[2];
	p = t*t*inv;
	t *= inv;
	a1 += p;
	p *= t;
	a2 += p;
	p *= t;
	a3 += p;
	p *= t;
	a4 += p;

	t = Fneq[3];
	inv = 1.f / FAA[3];
	p = t*t*inv;
	t *= inv;
	a1 += p;
	p *= t;
	a2 += p;
	p *= t;
	a3 += p;
	p *= t;
	a4 += p;

	t = Fneq[4];
	inv = 1.f / FAA[4];
	p = t*t*inv;
	t *= inv;
	a1 += p;
	p *= t;
	a2 += p;
	p *= t;
	a3 += p;
	p *= t;
	a4 += p;

	t = Fneq[5];
	inv = 1.f / FAA[5];
	p = t*t*inv;
	t *= inv;
	a1 += p;
	p *= t;
	a2 += p;
	p *= t;
	a3 += p;
	p *= t;
	a4 += p;

	t = Fneq[6];
	inv = 1.f / FAA[6];
	p = t*t*inv;
	t *= inv;
	a1 += p;
	p *= t;
	a2 += p;
	p *= t;
	a3 += p;
	p *= t;
	a4 += p;

	t = Fneq[7];
	inv = 1.f / FAA[7];
	p = t*t*inv;
	t *= inv;
	a1 += p;
	p *= t;
	a2 += p;
	p *= t;
	a3 += p;
	p *= t;
	a4 += p;

	t = Fneq[8];
	inv = 1.f / FAA[8];
	p = t*t*inv;
	t *= inv;
	a1 += p;
	p *= t;
	a2 += p;
	p *= t;
	a3 += p;
	p *= t;
	a4 += p;

	a1 /= 2.f;
	a2 /= -6.f;
	a3 /= 12.f;
	a4 /= -20.f;
	return 2.f - 4.f*a2 / a1 + 16.f*a2*a2 / a1 / a1 - 8.f*a3 / a1 + 80.f*a2*a3 / a1 / a1 - 80.f*a2*a2*a2 / a1 / a1 / a1 - 16.f*a4 / a1;
}



































double computeEntropy(double f[])
{
    double entropy = 0.;
    for (int iPop = 0; iPop < 9; iPop++)
	{
       entropy += f[iPop]*log(f[iPop]/W[iPop]);
    }
   
    return entropy;
}

double computeEntropyGrowth(double f[], double fNeq[], double alpha)
{
    double fAlphaFneq[9];
    for (int iPop = 0; iPop < 9; iPop++)
    {
		fAlphaFneq[iPop] = f[iPop] - alpha*fNeq[iPop];
    }
	return computeEntropy(f) - computeEntropy(fAlphaFneq);
}

double computeEntropyGrowthDerivative(double f[], double fNeq[], double alpha)
{
    double entropyGrowthDerivative = 0.;
    for (int iPop = 0; iPop < 9; iPop++)
    {
		double tmp = f[iPop] - alpha*fNeq[iPop];
		entropyGrowthDerivative += fNeq[iPop]*(log(tmp/W[iPop]));
    }
   
    return entropyGrowthDerivative;
}

double computeF1(double f[], double del[], double alpha)
{
    double entropyGrowthDerivative = 0.;
    for (int iPop = 0; iPop < 9; iPop++)
    {
		double tmp = f[iPop] + alpha*del[iPop];
		entropyGrowthDerivative += del[iPop]*(log(tmp/W[iPop]) + 1.);
    }
   
    return entropyGrowthDerivative;
}

double computeF(double f[], double del[], double alpha)
{
    double entropyGrowthDerivative = 0.;
    for (int k = 0; k < 9; k++)
    {
		double tmp = f[k] + alpha*del[k];
		entropyGrowthDerivative += tmp*log(tmp/W[k]) - f[k]*log(f[k]/W[k]);
    }
   
    return entropyGrowthDerivative;
}

double computeF2(double f[], double del[], double alpha)
{
    double entropyGrowthDerivative = 0.;
    for (int k = 0; k < 9; k++)
    {
		entropyGrowthDerivative += del[k]*del[k]/(alpha*del[k]+f[k]);
	}
   
    return entropyGrowthDerivative;
}

double getAlpha(double alpha, double f[], double fNeq[])
{
    double epsilon = EPSILON;
    double alphaGuess = 0.;
    double var = 100.0;
    double errorMax = epsilon*var;
    double error = 1.0;
    for (int count = 0; count < 10000; count++)
    {
		double entGrowth = computeEntropyGrowth(f,fNeq,alpha);
		double entGrowthDerivative = computeEntropyGrowthDerivative(f,fNeq,alpha);
		if ((error < errorMax) || (fabs(entGrowth) < var*epsilon))
		{
         return alpha;
		}
		alphaGuess = alpha - entGrowth / entGrowthDerivative;
		error = fabs(alpha-alphaGuess);
		alpha = alphaGuess;
    }
    return alpha;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1))) void Entropic_Collision_Kernel(__global double *Ro_array,
	__global double2 *J_array,
	__global int8 *Pos,
	__global double *F0,
	__global double *FA,
	__global double *FB,
	__global double *Alpha,
	double tau,
	double Fbody)
{

	int ix = get_global_id(0);
	int iy = get_global_id(1);
	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
		return;

	int jj = ix * FULLSIZEY + iy;
	int ii = jj * 8;
	int indUT = ix * FULLSIZEY_UT + iy;

	int8 Ploc = Pos[jj];
	double Fprev[9];
	Fprev[0] = F0[jj];
	double Ro = Fprev[0];
	double Jx = 0.;
	double Jy = 0.;
	double alphaij = Alpha[jj];

	for (int k = 0; k < 8; k++)
	{
		Fprev[k+1] = FA[ii + k];
		Ro += Fprev[k+1];
		Jx += CCxt[k] * Fprev[k+1];
		Jy += CCyt[k] * Fprev[k+1];
	}

	Jx /= Ro;
	Jy /= Ro;

	Ro_array[indUT] = Ro;
	J_array[indUT].x = Jx;
	J_array[indUT].y = Jy;

#ifndef USE_SPLIT_FORCE
	Jx += Fbody/Ro;
#endif

	double UC1 = 2. - sqrt(1. + 3.*Jx*Jx);
	double UC2 = (2.*Jx + sqrt(1. + 3.*Jx*Jx))/(1. - Jx);
	double VC1 = 2. - sqrt(1. + 3.*Jy*Jy);
	double VC2 = (2.*Jy + sqrt(1. + 3.*Jy*Jy))/(1. - Jy);
	double Feq[9], Del[9];
	
	Feq[0] = W[0] * Ro * UC1 * VC1;
	Feq[1] = W[1] * Ro * UC1 * UC2 * VC1;
	Feq[2] = W[2] * Ro * UC1 / UC2 * VC1;
	Feq[3] = W[3] * Ro * UC1 * VC1 * VC2;
	Feq[4] = W[4] * Ro * UC1 * VC1 / VC2;
	Feq[5] = W[5] * Ro * UC1 * UC2 * VC1 * VC2;
	Feq[6] = W[6] * Ro * UC1 / UC2 * VC1 / VC2;
	Feq[7] = W[7] * Ro * UC1 * UC2 * VC1 / VC2;
	Feq[8] = W[8] * Ro * UC1 / UC2 * VC1 * VC2;

#ifdef USE_SPLIT_FORCE
	double Jx_f = Jx + Fbody;
	double Feq_f[9];
	UC1 = 2. - sqrt(1. + 3.*Jx_f*Jx_f);
	UC2 = (2.*Jx_f + sqrt(1. + 3.*Jx_f*Jx_f))/(1. - Jx_f);
	Feq_f[0] = W[0] * Ro * UC1 * VC1;
	Feq_f[1] = W[1] * Ro * UC1 * UC2 * VC1;
	Feq_f[2] = W[2] * Ro * UC1 / UC2 * VC1;
	Feq_f[3] = W[3] * Ro * UC1 * VC1 * VC2;
	Feq_f[4] = W[4] * Ro * UC1 * VC1 / VC2;
	Feq_f[5] = W[5] * Ro * UC1 * UC2 * VC1 * VC2;
	Feq_f[6] = W[6] * Ro * UC1 / UC2 * VC1 / VC2;
	Feq_f[7] = W[7] * Ro * UC1 * UC2 * VC1 / VC2;
	Feq_f[8] = W[8] * Ro * UC1 / UC2 * VC1 * VC2;
#endif
	
 	double delmin = 10;
	for(int k = 0; k < 9; k++)
	{
		Del[k] = Feq[k] - Fprev[k];
		if(fabs(Del[k]) < delmin)
			delmin = fabs(Del[k]);
	}
	double F = computeF(Fprev, Del, alphaij);
	double F1 = computeF1(Fprev, Del, alphaij);
	//double F2 = computeF2(Fprev, Del, alphaij);
	double Alpha_new;
	if(delmin < 1.0e-15)
		Alpha_new = 2.;
	else
	{
		Alpha_new = alphaij - F/F1;
		if(computeF(Fprev, Del, Alpha_new) < 0.)
		{
			Alpha_new -= 2.*computeF(Fprev, Del, Alpha_new) / computeF1(Fprev, Del, Alpha_new);
		}
	}

	
	Alpha[jj] = Alpha_new;
#ifdef USE_SPLIT_FORCE
	F0[jj] = Fprev[0] - Alpha_new*tau*(Del[0]) + (Feq_f[0] - Feq[0]);
	FB[Ploc.s0] = Fprev[1] - Alpha_new*tau*Del[1] + (Feq_f[1] - Feq[1]);
	FB[Ploc.s1] = Fprev[2] - Alpha_new*tau*Del[2] + (Feq_f[2] - Feq[2]);
	FB[Ploc.s2] = Fprev[3] - Alpha_new*tau*Del[3] + (Feq_f[3] - Feq[3]);
	FB[Ploc.s3] = Fprev[4] - Alpha_new*tau*Del[4] + (Feq_f[4] - Feq[4]);
	FB[Ploc.s4] = Fprev[5] - Alpha_new*tau*Del[5] + (Feq_f[5] - Feq[5]);
	FB[Ploc.s5] = Fprev[6] - Alpha_new*tau*Del[6] + (Feq_f[6] - Feq[6]);
	FB[Ploc.s6] = Fprev[7] - Alpha_new*tau*Del[7] + (Feq_f[7] - Feq[7]);
	FB[Ploc.s7] = Fprev[8] - Alpha_new*tau*Del[8] + (Feq_f[8] - Feq[8]);
#else
	F0[jj] = Fprev[0] + Alpha_new*tau*(Del[0]);
	FB[Ploc.s0] = Fprev[1] + Alpha_new*tau*Del[1];
	FB[Ploc.s1] = Fprev[2] + Alpha_new*tau*Del[2];
	FB[Ploc.s2] = Fprev[3] + Alpha_new*tau*Del[3];
	FB[Ploc.s3] = Fprev[4] + Alpha_new*tau*Del[4];
	FB[Ploc.s4] = Fprev[5] + Alpha_new*tau*Del[5];
	FB[Ploc.s5] = Fprev[6] + Alpha_new*tau*Del[6];
	FB[Ploc.s6] = Fprev[7] + Alpha_new*tau*Del[7];
	FB[Ploc.s7] = Fprev[8] + Alpha_new*tau*Del[8];
#endif
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1))) void ibb_kernel(__global int2 *ibb_loc, __global double2 *ibb_coeff, __global double *FB, int max_el)
{
	int i = get_global_id(0);

	if (i >= max_el)
		return;

	double Fold = FB[ibb_loc[i].x];
	double Fnew = ibb_coeff[i].x * Fold + ibb_coeff[i].y * FB[ibb_loc[i].y];
	FB[ibb_loc[i].x] = Fnew;
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RED, 1, 1))) void reduce_1(__global double *input, __global double *output, __local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;


	sdata[tid] = input[stride] + input[stride + 1];

	barrier(CLK_LOCAL_MEM_FENCE);
	// do reduction in shared mem
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// write result for this block to global mem
	if (tid == 0)
	{
		output[bid] = sdata[0];
	}

}

__kernel __attribute__((reqd_work_group_size(1,1,1))) void reduce_2(__global double *input, __global double *output)
{
	double temp = input[0];
	for (unsigned int s = 1; s < NUMBLOCKS; s++)
	{
		temp += input[s];
	}
	output[0] = (temp);
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RED, 1, 1))) void reduce_both_1(__global double *inputRo, __global double2 *inputJ, __global double *output, __local double *sdata)
{
	// load shared mem
	unsigned int tid = get_local_id(0);
	unsigned int bid = get_group_id(0);
	unsigned int gid = get_global_id(0);

	unsigned int localSize = get_local_size(0);
	unsigned int stride = gid * 2;

	sdata[tid*2] = inputRo[stride] + inputRo[(stride + 1)];
	sdata[tid*2+1] = inputJ[stride].x + inputJ[(stride + 1)].x;
	
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

__kernel __attribute__((reqd_work_group_size(1, 1, 1))) void reduce_both_2(__global double *input, __global double *output, __global int *savestep)
{
	int gid = get_global_id(0);
	unsigned int step_offset = savestep[0];
	savestep[gid] = step_offset+1;

	double temp = input[gid];
	for (unsigned int s = 1; s < NUMBLOCKS; s++)
	{
		temp += input[s * 2 + gid];
	}

	output[step_offset*2 + gid] = temp;

}

__kernel __attribute__((reqd_work_group_size(1, 1, 1))) void reduce_four_2(__global double *input, __global double *output, __global double *T, __global double2 *J, __global int *savestep)
{
	int gid = get_global_id(0);
	unsigned int step_offset = savestep[0];
	savestep[gid] = step_offset+1;

	if(gid < 2)
	{
		double temp = input[gid];
		for (unsigned int s = 1; s < NUMBLOCKS; s++)
		{
			temp += input[s * 2 + gid];
		}

		output[step_offset*4 + gid] = temp;
	}
	else if(gid == 2)
	{
		double temp1 = 0., temp2 = 0.;
		for (unsigned int s = 0; s < FULLSIZEY; s++)
		{
			temp1 += fabs(J[s].x*T[s]);
			temp2 += fabs(J[s].x);
		}
		output[step_offset*4 + gid] = temp1/temp2;
	}
	else
	{
		double temp1 = 0., temp2 = 0.;
		for (unsigned int s = 0; s < FULLSIZEY; s++)
		{
			int ind = s + (FULLSIZEX - 1) * FULLSIZEY_UT;
			temp1 += fabs(J[ind].x*T[ind]);
			temp2 += fabs(J[ind].x);
		}
		output[step_offset*4 + gid] = temp1/temp2;
	}
}


