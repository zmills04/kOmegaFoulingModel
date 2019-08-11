// Updates variables associated with tracer methods.
//
//
//
///////////////////////////////////////////////////////////////////////////
void calc_tr_coeffs(int linInd, int tnum, 
	double dXcx, double dXcy, double dXx, double dXy,
	__global double4* __restrict__ CoeffT00,
	__global double4* __restrict__ CoeffT10,
	__global double4* __restrict__ CoeffT01,
	__global double4* __restrict__ CoeffT11,
	__global double4* __restrict__ CoeffU00,
	__global double4* __restrict__ CoeffU10,
	__global double4* __restrict__ CoeffU01,
	__global double4* __restrict__ CoeffU11)
{
	switch (tnum)
	{
	case -1:
	{
		break;
	}
	case 1:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { 1., 0., 0., 0. } };
		CoeffT10[linInd] = { { (C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0. } };
		CoeffT01[linInd] = { { (C1y - 1.) / dXcy + 1., 0., C2y / dXcy, 0. } };
		CoeffT11[linInd] = { { (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., C2x / dXcx, C2y / dXcy, 0. } };


		CoeffU00[linInd] = { { 1., 0., 0., 0. } };
		CoeffU10[linInd] = { { 1. - 1. / dXcx, 0., 0., 0. } };
		CoeffU01[linInd] = { { 1. - 1. / dXcy, 0., 0., 0. } };
		CoeffU11[linInd] = { { 1. - 1. / dXcy - 1. / dXcx, 0, 0, 0 } };

		break;
	}
	case 2:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { C2x / dXcx, (C1x - 1.) / dXcx + 1., 0, 0 } };
		CoeffT10[linInd] = { { 0., 1., 0., 0. } };
		CoeffT01[linInd] = { { C2x / dXcx, (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., 0., C2y / dXcy } };
		CoeffT11[linInd] = { { 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy } };

		CoeffU00[linInd] = { { 0., 1. - 1. / dXcx, 0., 0. } };
		CoeffU10[linInd] = { { 0., 1., 0., 0. } };
		CoeffU01[linInd] = { { 0., 1. - 1. / dXcx - 1. / dXcy, 0., 0. } };
		CoeffU11[linInd] = { { 0., 1. - 1. / dXcy, 0., 0. } };
		break;
	}
	case 3:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { 1., 0., 0., 0. } };
		CoeffT10[linInd] = { { 0., 1., 0., 0. } };
		CoeffT01[linInd] = { { (C1x - 1.) / dXcx + 1., 0, C2x / dXcx, 0. } };
		CoeffT11[linInd] = { { 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy } };

		CoeffU00[linInd] = { { 1., 0., 0., 0. } };
		CoeffU10[linInd] = { { 0., 1., 0., 0. } };
		CoeffU01[linInd] = { { 1. - 1. / dXcx, 0., 0., 0. } };
		CoeffU11[linInd] = { { 0., 1. - 1. / dXcy, 0., 0. } };
		break;
	}
	case 4:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { C2y / dXcy, 0., (C1y - 1.) / dXcy + 1., 0. } };
		CoeffT10[linInd] = { { C2y / dXcy, 0., (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., C2x / dXcx } };
		CoeffT01[linInd] = { { 0., 0., 1., 0. } };
		CoeffT11[linInd] = { { 0., 0., (C1x - 1.) / dXcx + 1., C2x / dXcx } };

		CoeffU00[linInd] = { { 0, 0, 1 - 1 / dXcy, 0 } };
		CoeffU10[linInd] = { { 0, 0, 1. - 1. / dXcy - 1. / dXcx, 0 } };
		CoeffU01[linInd] = { { 0., 0., 1., 0. } };
		CoeffU11[linInd] = { { 0, 0, 1. - 1. / dXcx, 0 } };
		break;
	}
	case 5:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { 1., 0., 0., 0. } };
		CoeffT10[linInd] = { { (C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0. } };
		CoeffT01[linInd] = { { 0., 0., 1., 0. } };
		CoeffT11[linInd] = { { 0., 0., (C1y - 1.) / dXcy + 1., C2y / dXcy } };

		CoeffU00[linInd] = { { 1., 0., 0., 0. } };
		CoeffU10[linInd] = { { 1. - 1. / dXcx, 0., 0., 0. } };
		CoeffU01[linInd] = { { 0., 0., 1., 0. } };
		CoeffU11[linInd] = { { 0., 0., 1. - 1. / dXcy, 0. } };
		break;
	}
	case 6:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { C2x / dXcx, (C1x - 1.) / dXcx + 1., 0, 0 } };
		CoeffT10[linInd] = { { 0., 1., 0., 0. } };
		CoeffT01[linInd] = { { 0., 0., 1., 0. } };
		CoeffT11[linInd] = { { 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy } };

		CoeffU00[linInd] = { { 0., 1. - 1. / dXcx, 0., 0. } };
		CoeffU10[linInd] = { { 0., 1., 0., 0. } };
		CoeffU01[linInd] = { { 0., 0., 1., 0. } };
		CoeffU11[linInd] = { { 0., 1. - 1. / dXcy, 0., 0. } };
		break;
	}
	case 7:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));


		CoeffT00[linInd] = { { 1., 0., 0., 0. } };
		CoeffT10[linInd] = { { 0., 1., 0., 0. } };
		CoeffT01[linInd] = { { 0., 0., 1., 0. } };
		CoeffT11[linInd] = { { 0., (C1x - 1.) / dXcx + 1., 0, C2x / dXcx } };

		CoeffU00[linInd] = { { 1., 0., 0., 0. } };
		CoeffU10[linInd] = { { 0., 1., 0., 0. } };
		CoeffU01[linInd] = { { 0., 0., 1., 0. } };
		CoeffU11[linInd] = { { 0., 1. - 1. / dXcx, 0., 0. } };
		break;
	}
	case 8:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { 0., C2y / dXcy, C2x / dXcx, (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1. } };
		CoeffT10[linInd] = { { 0., C2y / dXcy, 0., (C1y - 1.) / dXcy + 1. } };
		CoeffT01[linInd] = { { 0., 0., C2x / dXcx, (C1x - 1.) / dXcx + 1. } };
		CoeffT11[linInd] = { { 0., 0., 0., 1. } };

		CoeffU00[linInd] = { { 0., 0., 0., 1. - 1. / dXcy - 1. / dXcx } };
		CoeffU10[linInd] = { { 0., 0., 0., 1. - 1. / dXcy } };
		CoeffU01[linInd] = { { 0., 0., 0., 1. - 1. / dXcx } };
		CoeffU11[linInd] = { { 0., 0., 0., 1. } };
		break;
	}
	case 9:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { 1., 0., 0., 0. } };
		CoeffT10[linInd] = { { (C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0. } };
		CoeffT01[linInd] = { { (C1y - 1.) / dXcy + 1., 0., C2y / dXcy, 0. } };
		CoeffT11[linInd] = { { 0., 0., 0., 1. } };

		CoeffU00[linInd] = { { 1., 0., 0., 0. } };
		CoeffU10[linInd] = { { 1. - 1. / dXcx, 0., 0., 0. } };
		CoeffU01[linInd] = { { 1. - 1. / dXcy, 0., 0., 0. } };
		CoeffU11[linInd] = { { 0., 0., 0., 1. } };
		break;
	}
	case 10:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { C2x / dXcx, (C1x - 1.) / dXcx + 1., 0., 0. } };
		CoeffT10[linInd] = { { 0., 1., 0., 0. } };
		CoeffT01[linInd] = { { 0., 0., C2y / dXcy, (C1y - 1.) / dXcy + 1. } };
		CoeffT11[linInd] = { { 0., 0., 0., 1. } };

		CoeffU00[linInd] = { { 0., 1. - 1. / dXcx, 0., 0. } };
		CoeffU10[linInd] = { { 0., 1., 0., 0. } };
		CoeffU01[linInd] = { { 0., 0., 0., 1. - 1. / dXcy } };
		CoeffU11[linInd] = { { 0., 0., 0., 1. } };
		break;
	}
	case 11:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { 1., 0., 0., 0. } };
		CoeffT10[linInd] = { { 0., 1., 0., 0. } };
		CoeffT01[linInd] = { { 1. - (C1x - 1.) / dXcx, 0., -C2x / dXcx, 0. } };
		CoeffT11[linInd] = { { 0., 0., 0., 1. } };

		CoeffU00[linInd] = { { 1., 0., 0., 0. } };
		CoeffU10[linInd] = { { 0., 1., 0., 0. } };
		CoeffU01[linInd] = { { 1. - 1. / dXcx, 0., 0., 0. } };
		CoeffU11[linInd] = { { 0., 0., 0., 1. } };
		break;
	}
	case 12:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C1y = -(K_AIR_LB * (dXcy - dXy)) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2y = (K_SOOT_LB * dXcy) / (K_SOOT_LB * dXcy - K_AIR_LB * (dXcy - dXy));

		CoeffT00[linInd] = { { C2x / dXcx, 0., (C1x - 1.) / dXcx + 1., 0 } };
		CoeffT10[linInd] = { { 0., C2y / dXcy, 0., (C1y - 1.) / dXcy + 1. } };
		CoeffT01[linInd] = { { 0., 0., 1., 0. } };
		CoeffT11[linInd] = { { 0., 0., 0., 1. } };

		CoeffU00[linInd] = { { 0., 0., 1. - 1. / dXcx, 0. } };
		CoeffU10[linInd] = { { 0., 0., 0., 1. - 1. / dXcy } };
		CoeffU01[linInd] = { { 0., 0., 1., 0. } };
		CoeffU11[linInd] = { { 0., 0., 0., 1. } };
		break;
	}
	case 13:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));

		CoeffT00[linInd] = { { 1., 0., 0., 0 } };
		CoeffT10[linInd] = { { 0., C2x / dXcx, 0., (C1x - 1.) / dXcy + 1. } };
		CoeffT01[linInd] = { { 0., 0., 1., 0. } };
		CoeffT11[linInd] = { { 0., 0., 0., 1. } };

		CoeffU00[linInd] = { { 1., 0., 0., 0. } };
		CoeffU10[linInd] = { { 0., 0., 0., 1. - 1. / dXcx } };
		CoeffU01[linInd] = { { 0., 0., 1., 0. } };
		CoeffU11[linInd] = { { 0., 0., 0., 1. } };
		break;
	}
	case 14:
	{
		double C1x = -(K_AIR_LB * (dXcx - dXx)) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));
		double C2x = (K_SOOT_LB * dXcx) / (K_SOOT_LB * dXcx - K_AIR_LB * (dXcx - dXx));

		CoeffT00[linInd] = { { C2x / dXcx, 0., (C1x - 1.) / dXcx + 1., 0 } };
		CoeffT10[linInd] = { { 0., 1., 0., 0. } };
		CoeffT01[linInd] = { { 0., 0., 1., 0. } };
		CoeffT11[linInd] = { { 0., 0., 0., 1. } };

		CoeffU00[linInd] = { { 0., 0., 1. - 1. / dXcx, 0. } };
		CoeffU10[linInd] = { { 0., 1., 0., 0. } };
		CoeffU01[linInd] = { { 0., 0., 1., 0. } };
		CoeffU11[linInd] = { { 0., 0., 0., 1. } };
		break;
	}
	case 15:
	{
		CoeffT00[linInd] = { { 1., 0., 0., 0. } };
		CoeffT10[linInd] = { { 0., 1., 0., 0. } };
		CoeffT01[linInd] = { { 0., 0., 1., 0. } };
		CoeffT11[linInd] = { { 0., 0., 0., 1. } };

		CoeffU00[linInd] = { { 1., 0., 0., 0. } };
		CoeffU10[linInd] = { { 0., 1., 0., 0. } };
		CoeffU01[linInd] = { { 0., 0., 1., 0. } };
		CoeffU11[linInd] = { { 0., 0., 0., 1. } };
		break;
	}
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR, 1, 1)))
void updateTRCoeffs(__global NTYPE_TYPE* nType,
	__global double4* __restrict__ CoeffT00,
	__global double4* __restrict__ CoeffT10,
	__global double4* __restrict__ CoeffT01,
	__global double4* __restrict__ CoeffT11,
	__global double4* __restrict__ CoeffU00,
	__global double4* __restrict__ CoeffU10,
	__global double4* __restrict__ CoeffU01,
	__global double4* __restrict__ CoeffU11,
	__global short* __restrict__ niWallFlag,
	__global double* __restrict__ dX_cur,
	__global double* __restrict__ dX,
	__global int* __restrict__ ActiveNodes,
	int nActiveNodes)
{
	int ii = get_global_id(0);

	if (ii >= nActiveNodes)
		return;

	// get index of nodes aways from wall
	int iLin = activeNodes[i0];

	int i, j;
	decodeGlobalIdx(iLin, &i, &j);

	i += TR_X_IND_START;

	//Nact Ntemp = Neighs[Nod_ind];
	uint ii00 = i + j * CHANNEL_LENGTH_FULL;
	uint ii10 = ii00 + 1;
	uint ii01 = ii00 + CHANNEL_LENGTH_FULL;
	uint ii11 = ii10 + CHANNEL_LENGTH_FULL;

	short wfTemp = WF_EMPTY;
	int tnum = 0;

	NTYPE_TYPE ntype = nType[ii00];
	if (ntype & M_SOLID_NODE)
		wfTemp |= WF_00_SOLID;
	if (ntype & M_FLUID_NODE)
		tnum += 1;
	
	ntype = nType[ii10];
	if (ntype & M_SOLID_NODE)
		wfTemp |= WF_10_SOLID;
	if (ntype & M_FLUID_NODE)
		tnum += 2;

	ntype = nType[ii01];
	if (ntype & M_SOLID_NODE)
		wfTemp |= WF_01_SOLID;
	if (ntype & M_FLUID_NODE)
		tnum += 4;

	ntype = nType[ii11];
	if (ntype & M_SOLID_NODE)
		wfTemp |= WF_11_SOLID;
	if (ntype & M_FLUID_NODE)
		tnum += 8;

	// if all are solid, set node as solid
	if (wfTemp == WF_TEST_ALL_SOLID)
	{
		wfTemp |= WF_SOLID;
	}

	// if any are not solid, set as fluid node
	if (wfTemp ^= WF_TEST_ALL_SOLID)
	{
		wfTemp |= WF_FLUID;
	}

	niWallFlag[iLin] = wfTemp;

	if (tnum == 0)
		tnum = -1;

	double dXcx, dXcy, dXx, dXy;

	switch (tnum)
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
		// only 00 is a fluid node, rest are solid
		// use east and north directions from ii00 node
		dXcx = dX_cur[ii00 * 4];
		dXcy = dX_cur[ii00 * 4 + 2];
		dXx = dX[ii00 * 4];
		dXy = dX[ii00 * 4 + 2];
		break;
	}
	case 2:
	{
		// only 10 is a fluid node, rest are solid
		// use west and north directions from ii10 node
		dXcx = dX_cur[ii10 * 4 + 1];
		dXcy = dX_cur[ii10 * 4 + 2];
		dXx = dX[ii10 * 4 + 1];
		dXy = dX[ii10 * 4 + 2];
		break;
	}
	case 3:
	{
		// 00 and 10 are fluid nodes, rest are solid
		// use north from both ii00 and ii10 nodes
		dXcx = dX_cur[ii00 * 4 + 2];
		dXcy = dX_cur[ii10 * 4 + 2];
		dXx = dX[ii00 * 4 + 2];
		dXy = dX[ii10 * 4 + 2];
		break;
	}
	case 4:
	{
		// only 01 is fluid node
		// use east and south from ii01 node
		dXcx = dX_cur[ii01 * 4];
		dXcy = dX_cur[ii01 * 4 + 3];
		dXx = dX[ii01 * 4];
		dXy = dX[ii01 * 4 + 3];
		break;
	}
	case 5:
	{
		// 00 and 01 are fluid nodes, rest are solid
		// use east from both ii00 and ii01 nodes
		dXcx = dX_cur[ii00 * 4];
		dXcy = dX_cur[ii01 * 4];
		dXx = dX[ii00 * 4];
		dXy = dX[ii01 * 4];
		break;
	}
	case 6:
	{
		// 10 and 01 are fluid nodes, rest are solid
		// this should rarely if ever occur, will just treat it the same
		// as tnum == 2
		dXcx = dX_cur[ii10 * 4 + 1];
		dXcy = dX_cur[ii10 * 4 + 2];
		dXx = dX[ii10 * 4 + 1];
		dXy = dX[ii10 * 4 + 2];
		break;
	}
	case 7:
	{
		// 00, 01 and 10 are fluid nodes, 11 is solid
		// using north from ii10 for only dx value.
		dXcx = dX_cur[ii10 * 4 + 2];
		dXcy = 1.;
		dXx = dX[ii10 * 4 + 2];
		dXy = 1.;
		break;
	}
	case 8:
	{
		// 11 is fluid node, rest are solid
		// using west and south from ii11
		dXcx = dX_cur[ii11 * 4 + 1];
		dXcy = dX_cur[ii11 * 4 + 3];
		dXx = dX[ii11 * 4 + 1];
		dXy = dX[ii11 * 4 + 3];
		break;
	}
	case 9:
	{
		// 00 and 11 are fluid nodes, rest are solid
		// this should rarely if ever occur, will just treat it the same
		// as tnum == 1
		dXcx = dX_cur[ii00 * 4];
		dXcy = dX_cur[ii00 * 4 + 2];
		dXx = dX[ii00 * 4];
		dXy = dX[ii00 * 4 + 2];
		break;
	}
	case 10:
	{
		// 10 and 11 are fluid nodes, rest are solid
		// using west from both ii10 and ii11 nodes
		dXcx = dX_cur[ii10 * 4 + 1];
		dXcy = dX_cur[ii11 * 4 + 1];
		dXx = dX[ii10 * 4 + 1];
		dXy = dX[ii11 * 4 + 1];
		break;
	}
	case 11:
	{
		// 00 and 10 and 11 are fluid nodes, 01 is solid
		// using north from ii00 for only dx value
		dXcx = dX_cur[ii00 * 4 + 2];
		dXcy = 1.;
		dXx = dX[ii00 * 4 + 2];
		dXy = 1.;
		break;
	}
	case 12:
	{
		// 01 and 11 are fluid nodes, 00 and 10 are solid
		// using south from both 01 and 11 nodes
		dXcx = dX_cur[ii01 * 4 + 3];
		dXcy = dX_cur[ii11 * 4 + 3];
		dXx = dX[ii01 * 4 + 3];
		dXy = dX[ii11 * 4 + 3];
		break;
	}
	case 13:
	{
		// 00 and 01 and 11 are fluid nodes, 10 is solid
		// using south from ii11 for only dx value
		dXcx = dX_cur[ii11 * 4 + 3];
		dXcy = 1.;
		dXx = dX[ii11 * 4 + 3];
		dXy = 1.;
		break;
	}
	case 14:
	{
		// 00 and 01 and 11 are fluid nodes, 10 is solid
		// using south from ii01 node for only dx value
		dXcx = dX_cur[ii01 * 4 + 3];
		dXcy = 1.;
		dXx = dX[ii01 * 4 + 3];
		dXy = 1.;
		break;
	}
	case 15:
	{
		// all fluid nodes, so all dx values are 1
		dXcx = 1.;
		dXcy = 1.;
		dXx = 1.;
		dXy = 1.;
	}
	}

	calc_tr_coeffs(iLin, tnum, dXcx, dXcy, dXx, dXy,
		CoeffT00, CoeffT10, CoeffT01, CoeffT11,
		CoeffU00, CoeffU10, CoeffU01, CoeffU11);
}

// TODO: create dist array to track distances, and fill with 1000. before
//		calling this function every time
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_WALL, 1, 1)))
void updateTRWallNodes(__global ushort2* __restrict__ blP01,
	__global double2* __restrict__ Cvals,
	__global double* __restrict__ distVals,
	__global short* __restrict__ niWallFlag,
	__global int* __restrict__ niBLInd)
{
	//zero BLind_ind before calling this kernel
	int i = get_global_id(0);
	short wallflagval;

	if (i > MIN_BL_BOT && i <= MAX_BL_BOT)
		wallflagval = WF_BOT_WALL;
	else if (i > MIN_BL_TOP && i <= MAX_BL_TOP)
		wallflagval = WF_TOP_WALL;
	else
		return;

	ushort2 p01ind = blP01[i];

	double2 c1Val = Cvals[p01ind.x], c2Val = Cvals[p01ind.y];
	double2 cCenter = 0.5*c1Val + 0.5*c2Val;

	double2 vT = normalize(c2Val - c1Val);


	c1Val -= 0.05 * vT;
	c2Val += 0.05 * vT;

	int2 Cmin = min2(c0, c1);
	int2 Cmax = max2(c0, c1);

	Cmin.x -= TR_X_IND_START;
	Cmax.x -= TR_X_IND_START;
	
	Cmin.x = max(0, Cmin.x);
	Cmax.x = min(TR_X_IND_STOP - 1, Cmax.x);

	Cmin.y = max(0, Cmin.y);
	Cmax.y = min(FULLSIZEY_TR - 1, Cmax.y);


	if (Cmax.x < 0 || Cmin.x >= TR_X_IND_STOP ||
		Cmax.y < 0 || Cmin.y >= FULLSIZEY_TR)
		return;

	for (int i0 = Cmin.x; i0 < Cmax.x; i0++)
	{
		for (int j0 = Cmin.y; j0 < Cmax.y; j0++)
		{
			int linInd = i0 + FULLSIZEX_TR_PADDED * j0;
			double centerVal = convert_double2((int2)(i0 + TR_X_IND_START, j0));
			centerVal += 0.5 - cCenter;

			double distTemp = length(centerVal);
			if (distTemp > 2.)
				continue;

			if (AtomicMin(distVals[linInd], distTemp))
			{
				niBLInd[linInd] = i;
				niWallFlag[linInd] |= wallflagval;
			}
		}
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_UPDATEWALL, 1, 1)))
	void findTRWallNodes(__global short* __restrict__ niWallFlag,
	__global int* __restrict__ Winds,
	__global int* __restrict__ cur_loc,
	__global int* __restrict__ ActiveNodes, 
	int nActiveNodes, int windsCurLength)
{
	int i = get_global_id(0);
	if (i >= nActiveNodes)
		return;

	int ii = ActiveNodes[i];

	if (niWallFlag[ii] & WF_WALL)
	{
		int ind = atomic_add(&cur_loc[0], 1);
		ind = min(windsCurLength - 1, ind);
		Winds[ind] = i;
	}
}



int testInside_Shift(double2 vd, double2 v0, double2 v1)
{
	double2 v10 = v1 - v0, vd0 = vd - v0, vd1 = vd - v1;

	if (fabs(length(v10) - length(vd0) - length(vd1)) < CEPS)
		return 1;
	return 0;
}

int bcFindIntersectionLinePlane_shift(double* dist, double2 vL0, double2 vP0, double2 vP1, double2 * vN)
{
	*vN = normalize(double2(vP0.y - vP1.y, vP1.x - vP0.x));
	double2 vLd = (vP0 - vL0);
	*dist = dot(*vN, vLd);
	if ((*dist) <= 0.)
		return 0;
	double2 vCcut = vL0 + (*vN) * (*dist);
	return testInside_Shift(vCcut, vP0, vP1);
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_RERELEASE, 1, 1)))
void Shift_deposited_particles(__global int* __restrict__ pDepFlag,
	__global int* __restrict__ pLoc,
	__global double2 * __restrict__ pPos,
	__global ushort2 * __restrict__ blP01,
	__global double2 * __restrict__ Cvals,
	uint offset,
	uint maxel)
{
	int i = get_global_id(0);

	if (i >= maxel)
		return;

	i += offset;

	int depBL = pDepFlag[i];

	if (depBL < 0)
		return;

	ushort2 p01ind = blP01[depBL];
	double2 pPosNew = 0.5 * (Cvals[p01ind.x] + Cvals[p01ind.y]);

	if (pPosNew.x >= X_MAX_VAL)
	{
		pDepFlag[i] = -2;
		pLoc[i] = -2;
	}

	pPos[i] = pPosNew;
}


__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SHIFT_PAR, 1, 1)))
void Shift_particles(__global int* __restrict__ pDepFlag,
	__global int* __restrict__ pLoc,
	__global double2 * __restrict__ pPos,
	__global ushort2 * __restrict__ blP01,
	__global double2 * __restrict__ Cvals,
	__global short* __restrict__ niWallFlag,
	__global int* __restrict__ niBLInd)
{

	int i = get_global_id(0);
	if (i >= TRC_NUM_TRACERS)
		return;

	double2 C2 = pPos[i];

	if ((C2.x < X_release) && (C2.y > FL_YMAX || C2.y < FL_YMIN))
	{
		pDepFlag[i] = -2;
		pLoc[i] = -2;
		return;
	}

	trLoc = pLoc[i];


	short wallFlag = niWallFlag[trLoc];
	if (!(wallFlag & WF_WALL))
		return;

	double dist, dist_use = 100.;	//||dCc||/||dC|| where dCc is the distance between  
	double2 vN, vN_use;				//Intersection pt, Normal vector, Velocity of surf at vCcut

	int bl, bl_use = -1;			//boundary link index (value stored in each linked list element
	int blInd = niBLInd[trLoc];

	int minBl = (wallFlag & WF_BOT_WALL) ? MIN_BL_BOT : MIN_BL_TOP;
	int maxBl = (wallFlag & WF_BOT_WALL) ? MAX_BL_BOT : MAX_BL_TOP;

	int indStart = max(minBL, blInd - 5);
	int indStop = min(maxBl, blInd + 6);


	for (int k = indStart; k < indStop; k++)
	{
		ushort2 p01ind = blP01[k];
		double2 vC0 = Cvals[p01ind.x], vC1 = Cvals[p01ind.y];
		if (bcFindIntersectionLinePlane_shift(&dist, C2, vC0, vC1, &vN))
		{
			if (dist < dist_use)
			{
				bl_use = bl;
				dist_use = dist;
				vN_use = vN;
			}
		}
	}


	if (bl_use != -1)
	{
		C2 += vN_use * 2. * dist_use;
		int2 Posi = convert_int2(C2);
		trLoc = Posi.x + FULLSIZEX_TR_PADDED * Posi.y;
		if (C2.x >= X_MAX_VAL)
		{
			pDepFlag[i] = -2;
			trLoc = -2;
		}
		pLoc[i] = trLoc;
		pPos[i] = C2;
	}
}







//
//__kernel __attribute__((reqd_work_group_size(1, 1, 1)))
//void find_wall_nodes(__global nodeI* NodI,
//	__global int* Winds,
//	__global int2* ActiveNodes,
//	int nActiveNodes)
//{
//	int count = 0;
//	for (int i = 0; i < nActiveNodes; i++)
//	{
//		int2 ii = ActiveNodes[i];
//
//		if (NodI[ii.x * FULLSIZEY_TR + ii.y].Wall_Flag != 0)
//		{
//			Winds[count++] = i;
//		}
//	}
//}
