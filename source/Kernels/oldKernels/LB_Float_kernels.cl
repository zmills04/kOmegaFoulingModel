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
	h = log(t) + LOG19;
	ent += t*h;
	*dent += Fneq[1] * (h + 1.f);

	t = fma(alp, Fneq[2], FAA[2]);
	h = log(t) + LOG19;
	ent += t*h;
	*dent += Fneq[2] * (h + 1.f);

	t = fma(alp, Fneq[3], FAA[3]);
	h = log(t) + LOG19;
	ent += t*h;
	*dent += Fneq[3] * (h + 1.f);

	t = fma(alp, Fneq[4], FAA[4]);
	h = log(t) + LOG19;
	ent += t*h;
	*dent += Fneq[4] * (h + 1.f);

	t = fma(alp, Fneq[5], FAA[5]);
	h = log(t) + LOG136;
	ent += t*h;
	*dent += Fneq[5] * (h + 1.f);

	t = fma(alp, Fneq[6], FAA[6]);
	h = log(t) + LOG136;
	ent += t*h;
	*dent += Fneq[6] * (h + 1.f);

	t = fma(alp, Fneq[7], FAA[7]);
	h = log(t) + LOG136;
	ent += t*h;
	*dent += Fneq[7] * (h + 1.f);

	t = fma(alp, Fneq[8], FAA[8]);
	h = log(t) + LOG136;
	ent += t*h;
	*dent += Fneq[8] * (h + 1.f);

	return ent;
}

inline float Calc_Entropy(float *FAA)
{
	float ent = FAA[0] * (log(FAA[0]) + LOG49);
	ent += FAA[1] * (log(FAA[1]) + LOG19);
	ent += FAA[2] * (log(FAA[2]) + LOG19);
	ent += FAA[3] * (log(FAA[3]) + LOG19);
	ent += FAA[4] * (log(FAA[4]) + LOG19);
	ent += FAA[5] * (log(FAA[5]) + LOG136);
	ent += FAA[6] * (log(FAA[6]) + LOG136);
	ent += FAA[7] * (log(FAA[7]) + LOG136);
	ent += FAA[8] * (log(FAA[8]) + LOG136);
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

float Calc_alpha(float alpha_i, float *FAA, float *Fneq)
{
	float dev = CalcDeviation(FAA, Fneq);
	if (dev < 1.e-6f)
	{
		return 2.0f;
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

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_FD, WORKGROUPSIZEY_FD, 1)))
void LB_collision_SRT_Heat(__global double *Ro_array,
__global double2 *J_array,
__global double *T_array,
__global int8 *Pos,
__global float *FA,
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

#ifdef OPENCL_VERSION_1_2
	short lid = get_local_id(0) * WORKGROUPSIZEY_LB + get_local_id(1);
#else
	short lid = get_local_linear_id();
#endif
	int indUT = ix * FULLSIZEY_UT + iy;

	double rho = Ro_array[indUT];
	double2 UU = J_array[indUT];

#ifdef SOA_STORAGE
	double GAA[8];
	GAA[0] = GA[jj];
	GAA[1] = GA[jj + FULLSIZEXY];
	GAA[2] = GA[jj + 2 * FULLSIZEXY];
	GAA[3] = GA[jj + 3 * FULLSIZEXY];
	GAA[4] = GA[jj + 4 * FULLSIZEXY];
	GAA[5] = GA[jj + 5 * FULLSIZEXY];
	GAA[6] = GA[jj + 6 * FULLSIZEXY];
	GAA[7] = GA[jj + 7 * FULLSIZEXY];

	double Pneqx, Pneqy, Pneqz;

	Pneqx = UU.x*(rho - (FA[jj + 2 * FULLSIZEXY] + FA[jj + 3 * FULLSIZEXY]) - rho*(1. / 3. + UU.x*UU.x)); //Pneq.y
	Pneqy = UU.y*(rho - (FA[jj] + FA[jj + FULLSIZEXY]) - rho*(1. / 3. + UU.y*UU.y)); //Pneq.x
	Pneqz = (FA[jj + 4 * FULLSIZEXY] + FA[jj + 5 * FULLSIZEXY] - FA[jj + 6 * FULLSIZEXY] - FA[jj + 7 * FULLSIZEXY]) - rho * UU.x*UU.y; //Pneq.z

#else
	double GAA[8];
	GAA[0] = GA[jj * 8];
	GAA[1] = GA[jj * 8 + 1];
	GAA[2] = GA[jj * 8 + 2];
	GAA[3] = GA[jj * 8 + 3];
	GAA[4] = GA[jj * 8 + 4];
	GAA[5] = GA[jj * 8 + 5];
	GAA[6] = GA[jj * 8 + 6];
	GAA[7] = GA[jj * 8 + 7];

	double Pneqx, Pneqy, Pneqz;

	Pneqx = UU.x*(rho - convert_double(FA[jj * 8 + 2] + FA[jj * 8 + 3]) - rho*(1. / 3. + UU.x*UU.x)); //Pneq.y
	Pneqy = UU.y*(rho - convert_double(FA[jj * 8] + FA[jj * 8 + 1]) - rho*(1. / 3. + UU.y*UU.y)); //Pneq.x
	Pneqz = convert_double(FA[jj * 8 + 4] + FA[jj * 8 + 5] - FA[jj * 8 + 6] - FA[jj * 8 + 7]) - rho * UU.x*UU.y; //Pneq.z

#endif

	double RoE2 = G0[jj] + GAA[0] + GAA[1] + GAA[2] + GAA[3] + GAA[4] + GAA[5] + GAA[6] + GAA[7];

	double2 qq, qq_eq = (RoE2 * 3. + 2. * rho) * UU;

	qq.x = 3.*(GAA[0] - GAA[1] + GAA[4] - GAA[5] + GAA[6] - GAA[7]) - 6. * fma(UU.x, Pneqx, UU.y*Pneqz) - qq_eq.x;
	qq.y = 3.*(GAA[2] - GAA[3] + GAA[4] - GAA[5] - GAA[6] + GAA[7]) - 6. * fma(UU.x, Pneqz, UU.y*Pneqy) - qq_eq.y;

	qq *= omega_diff;


	double2 Reqxy = (RoE2 + 4. / 3.*rho) * UU * UU + 2. / 9.*rho;
	double Reqz = 9.*UU.x*UU.y * (RoE2 + 4. / 3.*rho);

	T_array[indUT] = (RoE2 / rho - dot(UU, UU)) / 2. - 1.;

	G0[jj] = G0[jj] * (1. - omega_f) + omega_f * 4. / 9. * (RoE2 - 1.5 * (Reqxy.x + Reqxy.y));

	int8 ii = Pos[jj];

	GB[ii.s0] = GAA[0] * (1. - omega_f) + (qq.x + omega_f * (RoE2 + qq_eq.x + 3.*Reqxy.x - 1.5*Reqxy.y)) / 9.;
	GB[ii.s1] = GAA[1] * (1. - omega_f) - (qq.x - omega_f * (RoE2 - qq_eq.x + 3.*Reqxy.x - 1.5*Reqxy.y)) / 9.;
	GB[ii.s2] = GAA[2] * (1. - omega_f) + (qq.y + omega_f * (RoE2 + qq_eq.y + 3.*Reqxy.y - 1.5*Reqxy.x)) / 9.;
	GB[ii.s3] = GAA[3] * (1. - omega_f) - (qq.y - omega_f * (RoE2 - qq_eq.y + 3.*Reqxy.y - 1.5*Reqxy.x)) / 9.;

	RoE2 += 3.*(Reqxy.y + Reqxy.x);

	GB[ii.s4] = GAA[4] * (1. - omega_f) + ((qq.x + qq.y) + omega_f * (RoE2 + (qq_eq.x + qq_eq.y) + Reqz)) / 36.;
	GB[ii.s5] = GAA[5] * (1. - omega_f) - ((qq.x + qq.y) - omega_f * (RoE2 - (qq_eq.x + qq_eq.y) + Reqz)) / 36.;
	GB[ii.s6] = GAA[6] * (1. - omega_f) + ((qq.x - qq.y) + omega_f * (RoE2 + (qq_eq.x - qq_eq.y) - Reqz)) / 36.;
	GB[ii.s7] = GAA[7] * (1. - omega_f) + ((qq.y - qq.x) + omega_f * (RoE2 + (qq_eq.y - qq_eq.x) - Reqz)) / 36.;
}



__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZEX_LB, WORKGROUPSIZEY_LB, 1)))
void LB_collision_SRT_Fluid(__global double *Ro_array,
__global double2 *J_array,
__global int8 *Pos,
__global float *F0,
__global float *FA,
__global float *FB,
__global char *Mred,
__global float *Alpha_array)
{

	int ix = get_global_id(0);
	int iy = get_global_id(1);
	if (ix >= FULLSIZEX || iy >= FULLSIZEY)
		return;

	int jj = ix * FULLSIZEY + iy;
	if (Mred[jj] == 0)
		return;

#ifdef OPENCL_VERSION_1_2
	short lid = get_local_id(0) * WORKGROUPSIZEY_LB + get_local_id(1);
#else
	short lid = get_local_linear_id();
#endif
	short lid9 = lid * 9;
	int indUT = ix * FULLSIZEY_UT + iy;

	float2 UU;
	float rho;

	__local float Fterm[WORKGROUPSIZEX_LB*WORKGROUPSIZEY_LB * 9];

#ifdef SOA_STORAGE
	float FAA[9];
	FAA[0] = F0[jj];
	FAA[1] = FA[jj];
	FAA[2] = FA[jj + FULLSIZEXY];
	FAA[3] = FA[jj + 2 * FULLSIZEXY];
	FAA[4] = FA[jj + 3 * FULLSIZEXY];
	FAA[5] = FA[jj + 4 * FULLSIZEXY];
	FAA[6] = FA[jj + 5 * FULLSIZEXY];
	FAA[7] = FA[jj + 6 * FULLSIZEXY];
	FAA[8] = FA[jj + 7 * FULLSIZEXY];
#else
	float FAA[9];
	FAA[0] = F0[jj];
	FAA[1] = FA[jj * 8];
	FAA[2] = FA[jj * 8 + 1];
	FAA[3] = FA[jj * 8 + 2];
	FAA[4] = FA[jj * 8 + 3];
	FAA[5] = FA[jj * 8 + 4];
	FAA[6] = FA[jj * 8 + 5];
	FAA[7] = FA[jj * 8 + 6];
	FAA[8] = FA[jj * 8 + 7];
#endif

	UU.x = FAA[1] - FAA[2] + FAA[5] - FAA[6] + FAA[7] - FAA[8];
	UU.y = FAA[3] - FAA[4] + FAA[5] - FAA[6] - FAA[7] + FAA[8];
	rho = FAA[0] + FAA[1] + FAA[2] + FAA[3] + FAA[4] + FAA[5] + FAA[6] + FAA[7] + FAA[8];

	UU /= rho;

	Ro_array[indUT] = convert_double(rho);
	J_array[indUT] = convert_double2(UU);

	float2 op3u = sqrt(1.f + 3.f*UU*UU);
	float Psi = rho * (2.f - op3u.x)*(2.f - op3u.y) / 9.f;
	float2 B = (2.f * UU + op3u) / (1.f - UU);
	float2 Binv = 1.f / B;

	op3u.x = sqrt(1.f + 3.f*(UU.x + FTERM_VAL)*(UU.x + FTERM_VAL));
	float Psi_du = rho * (2.f - op3u.x)*(2.f - op3u.y) / 9.f;
	float Bx_du = (2.f * (UU.x + FTERM_VAL) + op3u.x) / (1.f - (UU.x + FTERM_VAL));
	float Bxinv_du = 1.f / Bx_du;
	
	float Fdel[9];
	
	float feq = 4.f * Psi;
	Fdel[0] = feq - FAA[0];
	Fterm[lid9] = fma(4.f, Psi_du, -feq);

	feq = Psi * B.x;
	Fdel[1] = feq - FAA[1];
	Fterm[lid9 + 1] = fma(Psi_du, Bx_du, -feq);

	feq = Psi * Binv.x;
	Fdel[2] = feq - FAA[2];
	Fterm[lid9 + 2] = fma(Psi_du, Bxinv_du, -feq);

	feq = Psi * B.y;
	Fdel[3] = feq - FAA[3];
	Fterm[lid9 + 3] = fma(Psi_du, B.y, -feq);

	feq = Psi * Binv.y;
	Fdel[4] = feq - FAA[4];
	Fterm[lid9 + 4] = fma(Psi_du, Binv.y, -feq);

	Psi_du *= 0.25f;
	Psi *= 0.25f;

	feq = Psi * B.x * B.y;
	Fdel[5] = feq - FAA[5];
	Fterm[lid9 + 5] = fma(Psi_du, Bx_du * B.y, -feq);

	feq = Psi * Binv.x * Binv.y;
	Fdel[6] = feq - FAA[6];
	Fterm[lid9 + 6] = fma(Psi_du, Bxinv_du * Binv.y, -feq);

	feq = Psi * B.x * Binv.y;
	Fdel[7] = feq - FAA[7];
	Fterm[lid9 + 7] = fma(Psi_du, Bx_du * Binv.y, -feq);

	feq = Psi * Binv.x * B.y;
	Fdel[8] = feq - FAA[8];
	Fterm[lid9 + 8] = fma(Psi_du, Bxinv_du * B.y, -feq);


	float alpha_temp = Calc_alpha(convert_float(Alpha_array[jj]), FAA, Fdel);
	Alpha_array[jj] = alpha_temp;
	alpha_temp *= BETA_VAL_FLOAT;

	int8 ii = Pos[jj];

	F0[jj] += alpha_temp * Fdel[0] + Fterm[lid9];
	FB[ii.s0] = FAA[1] + alpha_temp * Fdel[1] + Fterm[lid9 + 1];
	FB[ii.s1] = FAA[2] + alpha_temp * Fdel[2] + Fterm[lid9 + 2];
	FB[ii.s2] = FAA[3] + alpha_temp * Fdel[3] + Fterm[lid9 + 3];
	FB[ii.s3] = FAA[4] + alpha_temp * Fdel[4] + Fterm[lid9 + 4];
	FB[ii.s4] = FAA[5] + alpha_temp * Fdel[5] + Fterm[lid9 + 5];
	FB[ii.s5] = FAA[6] + alpha_temp * Fdel[6] + Fterm[lid9 + 6];
	FB[ii.s6] = FAA[7] + alpha_temp * Fdel[7] + Fterm[lid9 + 7];
	FB[ii.s7] = FAA[8] + alpha_temp * Fdel[8] + Fterm[lid9 + 8];
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_INLET, 1, 1)))
void TLB_Inlet_Outlet(__global double2 *U_array,
__global double *Ro_array,
__global double *T_array,
__global double *G0,
__global double *GA,
__global double *GB)
{
	int i = get_global_id(0);

	if (i >= FULLSIZEY)
		return;

	double rho = Ro_array[i];
	double2 UU = U_array[i];
	double2 qq;
	double4 Req;

	Req.w = 2.*TFD_X_IN_VAL*rho + rho*dot(UU, UU);
	qq = fma(2., rho, 3.*Req.w) * UU;
	Req.xyz = fma(4. / 3., rho, Req.w);
	Req.xy = Req.xy*UU*UU + 2. / 9.*rho;
	Req.z *= 9. * (UU.x*UU.y);

	G0[i] = 4. / 9. * (Req.w - 1.5 * (Req.x + Req.y));
#ifdef SOA_STORAGE
	GB[i] = 1. / 9. * (Req.w + qq.x + 3.*Req.x - 1.5*Req.y);
	GB[i + FULLSIZEXY] = 1. / 9. * (Req.w - qq.x + 3.*Req.x - 1.5*Req.y);
	GB[i + 2 * FULLSIZEXY] = 1. / 9. * (Req.w + qq.y + 3.*Req.y - 1.5*Req.x);
	GB[i + 3 * FULLSIZEXY] = 1. / 9. * (Req.w - qq.y + 3.*Req.y - 1.5*Req.x);
	Req.w += 3.*(Req.x + Req.y);
	GB[i + 4 * FULLSIZEXY] = 1. / 36. * (Req.w + (qq.x + qq.y) + Req.z);
	GB[i + 5 * FULLSIZEXY] = 1. / 36. * (Req.w - (qq.x + qq.y) + Req.z);
	GB[i + 6 * FULLSIZEXY] = 1. / 36. * (Req.w + (qq.x - qq.y) - Req.z);
	GB[i + 7 * FULLSIZEXY] = 1. / 36. * (Req.w + (qq.y - qq.x) - Req.z);
#else
	GB[i * 8] = 1. / 9. * (Req.w + qq.x + 3.*Req.x - 1.5*Req.y);
	GB[i * 8 + 1] = 1. / 9. * (Req.w - qq.x + 3.*Req.x - 1.5*Req.y);
	GB[i * 8 + 2] = 1. / 9. * (Req.w + qq.y + 3.*Req.y - 1.5*Req.x);
	GB[i * 8 + 3] = 1. / 9. * (Req.w - qq.y + 3.*Req.y - 1.5*Req.x);
	Req.w += 3.*(Req.x + Req.y);
	GB[i * 8 + 4] = 1. / 36. * (Req.w + (qq.x + qq.y) + Req.z);
	GB[i * 8 + 5] = 1. / 36. * (Req.w - (qq.x + qq.y) + Req.z);
	GB[i * 8 + 6] = 1. / 36. * (Req.w + (qq.x - qq.y) - Req.z);
	GB[i * 8 + 7] = 1. / 36. * (Req.w + (qq.y - qq.x) - Req.z);
#endif
	rho = Ro_array[i + (FULLSIZEX - 1)*FULLSIZEY_UT];
	UU = U_array[i + (FULLSIZEX - 1)*FULLSIZEY_UT];
	double Tval = 1. + MAX(0., T_array[i + (FULLSIZEX - 1)*FULLSIZEY_UT]);
	Req.w = 2. * Tval * rho + rho*dot(UU, UU);

	i += (FULLSIZEX - 1)*FULLSIZEY;
	qq = fma(2., rho, 3.*Req.w) * UU;
	Req.xyz = fma(4. / 3., rho, Req.w);
	Req.xy = Req.xy*UU*UU + 2. / 9.*rho;
	Req.z *= 9. * (UU.x*UU.y);

#ifdef SOA_STORAGE
	GB[i + FULLSIZEXY] = 1. / 9. * (Req.w - qq.x + 3.*Req.x - 1.5*Req.y);
	Req.w += 3.*(Req.x + Req.y);
	if (i > 0)
	{
		GB[i + 7 * FULLSIZEXY] = 1. / 36. * (Req.w + (qq.y - qq.x) - Req.z);
	}
	if (i < FULLSIZEY - 1)
	{
		GB[i + 5 * FULLSIZEXY] = 1. / 36. * (Req.w - (qq.x + qq.y) + Req.z);
	}
#else
	GB[i * 8 + 1] = 1. / 9. * (Req.w - qq.x + 3.*Req.x - 1.5*Req.y);
	Req.w += 3.*(Req.x + Req.y);
	if (i > 0)
	{
		GB[i * 8 + 7] = 1. / 36. * (Req.w + (qq.y - qq.x) - Req.z);
	}
	if (i < FULLSIZEY - 1)
	{
		GB[i * 8 + 5] = 1. / 36. * (Req.w - (qq.x + qq.y) + Req.z);
	}
#endif
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1)))
void LB_IBB_ELBM_Fluid(__global short *ibb_flag,
__global int *ibb_loc,
__global double *ibb_coeff,
__global int *neighs,
__global double2 *dxdy,
__global float *F0,
__global float *FB,
__global double2 *U_array,
int max_el)
{
	int i = get_global_id(0);

	if (i >= max_el)
		return;

	int nloc = ibb_loc[i];
	double2 UUdouble = 0.;
	float rho = F0[nloc];
	char flag = ibb_flag[i];
	float count = 0.;
	__attribute__((opencl_unroll_hint(8)))
		for (int k = 0; k < 8; k++)
		{
			if (flag & IBB_Flag[k])
			{
				UUdouble += (ibb_coeff[i * 8 + k])*(U_array[neighs[i * 8 + RevDir[k]]] - U_array[neighs[i * 8 + k]]) + U_array[neighs[i * 8 + k]];
				count += 1.f;
			}
#ifdef SOA_STORAGE
			rho += FB[nloc + FULLSIZEXY * k];
#else
			rho += FB[nloc * 8 + k];
#endif
		}

	float2 UU = convert_float2(UUdouble) / count;
	float2 dUdy = convert_float2((U_array[neighs[i * 8 + 2]] - U_array[neighs[i * 8 + 3]]) / dxdy[i].y);
	float2 dUdx = convert_float2((U_array[neighs[i * 8]] - U_array[neighs[i * 8 + 1]]) / dxdy[i].x);

	float3 Pneq = (float3)(2.f*dUdx.x, 2.f*dUdy.y, 9.f*(dUdy.x + dUdx.y)) * rho / (-27.f*omega_f_float);

	float2 op3u = sqrt(1.f + 3.f*UU*UU);
	float Psi = rho * (2.f - op3u.x)*(2.f - op3u.y) / 9.f;
	float2 B = (2.f * UU + op3u) / (1.f - UU);
	float2 Binv = 1.f / B;
#ifdef SOA_STORAGE		
	if (flag & 0x1)
	{
		FB[nloc] = Psi*B.x + (3.f*Pneq.x - 1.5f*Pneq.y);
	}
	if (flag & 0x2)
	{
		FB[nloc + FULLSIZEXY] = Psi*Binv.x + (3.f*Pneq.x - 1.5f*Pneq.y);
	}
	if (flag & 0x4)
	{
		FB[nloc + 2 * FULLSIZEXY] = Psi*B.y + (3.f*Pneq.y - 1.5f*Pneq.x);
	}
	if (flag & 0x8)
	{
		FB[nloc + 3 * FULLSIZEXY] = Psi*Binv.y + (3.f*Pneq.y - 1.5f*Pneq.x);
	}
	if (flag & 0x10)
	{
		FB[nloc + 4 * FULLSIZEXY] = (Psi*B.x*B.y + (3.f*(Pneq.x + Pneq.y) + Pneq.z)) / 4.f;
	}
	if (flag & 0x20)
	{
		FB[nloc + 5 * FULLSIZEXY] = (Psi*Binv.x*Binv.y + (3.f*(Pneq.x + Pneq.y) + Pneq.z)) / 4.f;
	}
	if (flag & 0x40)
	{
		FB[nloc + 6 * FULLSIZEXY] = (Psi*B.x*Binv.y + (3.f*(Pneq.x + Pneq.y) - Pneq.z)) / 4.f;
	}
	if (flag & 0x80)
	{
		FB[nloc + 7 * FULLSIZEXY] = (Psi*Binv.x*B.y + (3.f*(Pneq.x + Pneq.y) - Pneq.z)) / 4.f;
	}
#else
	nloc *= 8;
	if (flag & 0x1)
	{
		FB[nloc] = Psi*B.x + (3.f*Pneq.x - 1.5f*Pneq.y);
	}
	if (flag & 0x2)
	{
		FB[nloc + 1] = Psi*Binv.x + (3.f*Pneq.x - 1.5f*Pneq.y);
	}
	if (flag & 0x4)
	{
		FB[nloc + 2] = Psi*B.y + (3.f*Pneq.y - 1.5f*Pneq.x);
	}
	if (flag & 0x8)
	{
		FB[nloc + 3] = Psi*Binv.y + (3.f*Pneq.y - 1.5f*Pneq.x);
	}
	if (flag & 0x10)
	{
		FB[nloc + 4] = (Psi*B.x*B.y + (3.f*(Pneq.x + Pneq.y) + Pneq.z)) / 4.f;
	}
	if (flag & 0x20)
	{
		FB[nloc + 5] = (Psi*Binv.x*Binv.y + (3.f*(Pneq.x + Pneq.y) + Pneq.z)) / 4.f;
	}
	if (flag & 0x40)
	{
		FB[nloc + 6] = (Psi*B.x*Binv.y + (3.f*(Pneq.x + Pneq.y) - Pneq.z)) / 4.f;
	}
	if (flag & 0x80)
	{
		FB[nloc + 7] = (Psi*Binv.x*B.y + (3.f*(Pneq.x + Pneq.y) - Pneq.z)) / 4.f;
	}
#endif
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1)))
void LB_IBB_ELBM_Heat(__global short *ibb_flag,
__global int *ibb_loc,
__global double *ibb_coeff,
__global int *macro_loc,
__global int *neighs,
__global double2 *dxdy,
__global double *GB,
__global double2 *U_array,
__global double *Ro_array,
__global double *T_array,
int max_el)
{
	int i = get_global_id(0);

	if (i >= max_el)
		return;

	int nloc = macro_loc[i];
	double2 U_tgt = U_array[nloc];
	double rho_tgt = Ro_array[nloc];

	double T_tgt = 0.;
	double count = 0.;
	char flag = ibb_flag[i];
	__attribute__((opencl_unroll_hint(8)))
		for (int k = 0; k < 8; k++)
		{
			if (flag & IBB_Flag[k])
			{
				T_tgt += 1. + (ibb_coeff[i * 8 + k]) * (T_array[neighs[i * 8 + RevDir[k]]] - T_array[neighs[i * 8 + k]]) + T_array[neighs[i * 8 + k]];
				count += 1.;
			}
		}

	T_tgt /= count;
	double2 dUdy = (U_array[neighs[i * 8 + 2]] - U_array[neighs[i * 8 + 3]]) / dxdy[i].y;
	double2 dUdx = (U_array[neighs[i * 8]] - U_array[neighs[i * 8 + 1]]) / dxdy[i].x;
	double2 dTdxy = (double2)((T_array[neighs[i * 8]] - T_array[neighs[i * 8 + 1]]), (T_array[neighs[i * 8 + 2]] - T_array[neighs[i * 8 + 3]])) / dxdy[i];

	double3 Pneq = (double3)(2.*dUdx.x, 2.*dUdy.y, (dUdy.x + dUdx.y)) * rho_tgt / (-3.*omega_f);

	double2 qq = -2. / omega_1f * rho_tgt / 3. * dTdxy + 2.*(double2)(U_tgt.x*Pneq.x + U_tgt.y*Pneq.z, U_tgt.x*Pneq.z + U_tgt.y*Pneq.y);

	double3 Rneq;
	Rneq.xy = 2.*U_tgt*dTdxy + 2.*(double2)(dUdx.x, dUdy.y) * (2. / 3. + T_tgt);
	Rneq.z = U_tgt.x*dTdxy.y + U_tgt.y*dTdxy.x + (dUdy.x + dUdx.y)*(2. / 3. + T_tgt);
	Rneq *= -2. / omega_1f / 3.*rho_tgt;

	double RoE2 = 2.*rho_tgt*T_tgt + rho_tgt*dot(U_tgt, U_tgt);
	qq += (RoE2 + 2.*rho_tgt / 3.) * U_tgt;
	qq *= 3.;
	Rneq += RoE2*(double3)(U_tgt.x*U_tgt.x, U_tgt.y*U_tgt.y, U_tgt.x*U_tgt.y) + 2.*rho_tgt / 3.*(double3)(1. / 3. + 2.*U_tgt.x*U_tgt.x, 1. / 3. + 2.*U_tgt.y*U_tgt.y, 2.*U_tgt.x*U_tgt.y);

#ifdef SOA_STORAGE
	int iloc = ibb_loc[i];
	if (flag & 0x1)
	{
		GB[iloc] = (RoE2 + qq.x + 3.*Rneq.x - 1.5*Rneq.y) / 9.;
	}
	if (flag & 0x2)
	{
		GB[iloc + FULLSIZEXY] = (RoE2 - qq.x + 3.*Rneq.x - 1.5*Rneq.y) / 9.;
	}
	if (flag & 0x4)
	{
		GB[iloc + 2 * FULLSIZEXY] = (RoE2 + qq.y + 3.*Rneq.y - 1.5*Rneq.x) / 9.;
	}
	if (flag & 0x8)
	{
		GB[iloc + 3 * FULLSIZEXY] = (RoE2 - qq.y + 3.*Rneq.y - 1.5*Rneq.x) / 9.;
	}

	RoE2 += 3.*(Rneq.x + Rneq.y);

	if (flag & 0x10)
	{
		GB[iloc + 4 * FULLSIZEXY] = (RoE2 + qq.x + qq.y + 9.*Rneq.z) / 36.;
	}
	if (flag & 0x20)
	{
		GB[iloc + 5 * FULLSIZEXY] = (RoE2 - qq.x - qq.y + 9.*Rneq.z) / 36.;
	}
	if (flag & 0x40)
	{
		GB[iloc + 6 * FULLSIZEXY] = (RoE2 + qq.x - qq.y - 9.*Rneq.z) / 36.;
	}
	if (flag & 0x80)
	{
		GB[iloc + 7 * FULLSIZEXY] = (RoE2 - qq.x + qq.y - 9.*Rneq.z) / 36.;
	}
#else
	int iloc = ibb_loc[i] * 8;
	if (flag & 0x1)
	{
		GB[iloc] = (RoE2 + qq.x + 3.*Rneq.x - 1.5*Rneq.y) / 9.;
	}
	if (flag & 0x2)
	{
		GB[iloc + 1] = (RoE2 - qq.x + 3.*Rneq.x - 1.5*Rneq.y) / 9.;
	}
	if (flag & 0x4)
	{
		GB[iloc + 2] = (RoE2 + qq.y + 3.*Rneq.y - 1.5*Rneq.x) / 9.;
	}
	if (flag & 0x8)
	{
		GB[iloc + 3] = (RoE2 - qq.y + 3.*Rneq.y - 1.5*Rneq.x) / 9.;
	}

	RoE2 += 3.*(Rneq.x + Rneq.y);

	if (flag & 0x10)
	{
		GB[iloc + 4] = (RoE2 + qq.x + qq.y + 9.*Rneq.z) / 36.;
	}
	if (flag & 0x20)
	{
		GB[iloc + 5] = (RoE2 - qq.x - qq.y + 9.*Rneq.z) / 36.;
	}
	if (flag & 0x40)
	{
		GB[iloc + 6] = (RoE2 + qq.x - qq.y - 9.*Rneq.z) / 36.;
	}
	if (flag & 0x80)
	{
		GB[iloc + 7] = (RoE2 - qq.x + qq.y - 9.*Rneq.z) / 36.;
	}
#endif
}

//Interpolated bounce back kernel
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_IBB, 1, 1)))
void LB_IBB(__global int2 *ibb_loc,
__global double2 *ibb_coeff,
__global float *FB,
int max_el)
{
	int i = get_global_id(0);

	if (i >= max_el)
		return;
	float2 coeff = convert_float2(ibb_coeff[i]);
	FB[ibb_loc[i].x] = coeff.x * FB[ibb_loc[i].x] + coeff.y * FB[ibb_loc[i].y];
}