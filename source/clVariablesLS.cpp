// clVariables.cpp: implementation of the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clVariablesFL.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"


void clVariablesLS::allocateArrays()
{
	ssArrIndMap.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	nType.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	BL.zeros(nBL, 2);
	C0.zeros(nN);
	C.zeros(nN);
	dXArr.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 4, 4);
	dXArr0.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 4, 4);

	if (vlb.kOmegaClass.kOmegaSolverFlag || vfd.thermalSolverFlag)
		dXCoeffs.zeros(p.XsizeFull, p.nY, 12);


	ssArrInds.zeros(p.nX);
	lsMap.zeros(p.nX, p.XsizeFull, 2, 2);

	Masses.allocate(2);

	ERROR_CHECKING(BL.load("load" SLASH "lsbl") == false,
		"Could not load lsbl.txt", ERROR_INITIALIZING_VLS);

	ERROR_CHECKING(C0.load("load" SLASH "lsc0") == false,
		"No lsc load file (initial locations provided in load/lsc0)",
		ERROR_INITIALIZING_VLS);
}

void clVariablesLS::allocateBuffers()
{
	//M.allocate_buffer_w_copy();
	nType.allocate_buffer_w_copy();
	BL.allocate_buffer_w_copy();
	C0.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	C.allocate_buffer_w_copy();
	dXArr.allocate_buffer_w_copy();
	dXArr.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	Masses.allocate_buffer_w_copy(2);
	ibbArr.allocate_buffer(CL_MEM_READ_ONLY);
	ibbArr.copy_dynamic_to_buffer();
	ssArrInds.allocate_buffer_w_copy();
	ssArrIndMap.allocate_buffer_w_copy();
	lsMap.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	if (vlb.kOmegaClass.kOmegaSolverFlag || vfd.thermalSolverFlag)
		dXCoeffs.allocate_buffer_w_copy();

}

void clVariablesLS::bcFindDirection(int dir, int bl, cl_double2 vC0,
	cl_double2 vC1, cl_double2 vCn, cl_int2 vCmin, cl_int2 vCmax, int findtype)
{
	double dist;
	cl_double2 Cdir = vlb.Cxy_double[dir];
	double Cndir = DOT_PROD(vCn, Cdir);
	if (Cndir == 0.)
		return;
	if (Cndir > 0.)
	{
		dir = vlb.rev[dir];
		Cdir = vlb.Cxy_double[dir];
	}
	cl_int2 vCcut0, vCcut1;
	cl_double2 vCcut;
	if (vlb.Cx[dir] != 0)
	{
		int i = 0;
		if (vlb.Cx[dir] == 1) i = vCmin.x;
		if (vlb.Cx[dir] == -1) i = vCmax.x;

		for (int j = vCmin.y; j <= vCmax.y; j++)
		{
			if (bcFindIntersection(vCcut0, vCcut1, vCcut, dist, { { (double)i, (double)j } }, Cdir, vC0, vC1, vCn) == true)
			{
				if (findtype == find_dx0)
					bcSetdXArr0(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
				else if (findtype == find_dx)
					bcSetdXArr(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
			}
		}
	}
	if (vlb.Cy[dir] != 0)
	{
		int j = 0;
		if (vlb.Cy[dir] == 1) j = vCmin.y;
		if (vlb.Cy[dir] == -1) j = vCmax.y;

		for (int i = vCmin.x; i <= vCmax.x; i++)
		{
			if (bcFindIntersection(vCcut0, vCcut1, vCcut, dist, { { (double)i, (double)j } }, Cdir, vC0, vC1, vCn) == true)
			{
				if (findtype == find_dx0)
					bcSetdXArr0(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
				else if (findtype == find_dx)
					bcSetdXArr(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
			}
		}
	}
}


bool clVariablesLS::bcFindIntersection(cl_int2 & vCcut0, cl_int2 & vCcut1,
	cl_double2 & vCcut, double& dist, cl_double2 vL0, cl_double2 vLd,
	cl_double2 vC0, cl_double2 vC1, cl_double2 vCn)
{
	double distF;
	if (bcFindIntersectionLinePlane(vCcut, distF, vL0, vLd, vC0, vC1) == false)
		return false;

	if (testInside(vCcut, vC0, vC1) == false)
		return false;

	int n_cut_0 = (int)ceil(distF);
	int n_cut = (int)ceil(distF);
	double dist_0 = 1. - (n_cut - (double)distF);
	dist = 1. - (n_cut - (double)distF);

	vCcut1 = { { (int)(vL0.x + vLd.x * n_cut), (int)(vL0.y + vLd.y * n_cut) } };
	vCcut0 = { { vCcut1.x - (int)vLd.x, vCcut1.y - (int)vLd.y } };

	cl_int2 ii0 = { { MOD(vCcut0.x, p.nX), MOD(vCcut0.y, p.nY) } };

	if (ii0.x != vCcut0.x)
	{
		return false;
	}

	if (ii0.y != vCcut0.y)
		return false;

	if (vls.nType(ii0.x, ii0.y) & M_SOLID_NODE)
	{
		if (dist < p.eps)
		{
			dist = 1.;
			n_cut--;
			vCcut1 = { { (int)(vL0.x + vLd.x * n_cut),
				(int)(vL0.y + vLd.y * n_cut) } };
			vCcut0 = { { vCcut1.x - (int)vLd.x,
				vCcut1.y - (int)vLd.y } };

			if (ii0.x != vCcut0.x)
				return false;

			if (ii0.y != vCcut0.y)
				return false;
		}
		else
			return false;
	}

	return true;
}

bool clVariablesLS::bcFindIntersectionLinePlane(cl_double2& vC, double& dist,
	cl_double2 vL0, cl_double2 vLd, cl_double2 vP0, cl_double2 vP1)
{
	cl_double2 vP10 = Subtract2(vP1, vP0);
	cl_double2 vN = { { -vP10.y, vP10.x } };
	double den = DOT_PROD(vN, vLd);
	if (den == 0.)
		return false;
	cl_double2 vPL = Subtract2(vP0, vL0);
	dist = DOT_PROD(vN, vPL) / den;
	vC = { { vL0.x + vLd.x * dist, vL0.y + vLd.y * dist } };
	return true;
}

void clVariablesLS::bcFindNodes(int bl, int findtype)
{
	int n0 = BL(bl, 0), n1 = BL(bl, 1);
	if (n0 < 0)
		return;
	cl_double2 vC0 = { { C[n0].x, C[n0].y } };
	cl_double2 vC1 = { { C[n1].x, C[n1].y } };

	if (findtype == find_dx0)
	{
		vC0 = { { C0[n0].x, C0[n0].y } };
		vC1 = { { C0[n1].x, C0[n1].y } };
	}

	cl_double2 vCn = { { vC0.y - vC1.y, vC1.x - vC0.x } };
	cl_int2 vCmin = min2(vC0, vC1), vCmax = max2(vC0, vC1);

	if (vCmin.x < 0) vCmin.x = 0;
	if (vCmin.x > p.nX - 1) return;
	if (vCmax.x < 0) return;
	if (vCmax.x > p.nX - 1) vCmax.x = p.nX - 1;

	if (vCmin.y < 0) vCmin.y = 0;
	if (vCmin.y > p.nY - 1) return;
	if (vCmax.y < 0) return;
	if (vCmax.y > p.nY - 1) vCmax.y = p.nY - 1;


	for (int dir = 1; dir < 9; dir += 2)
		bcFindDirection(dir, bl, vC0, vC1, vCn, vCmin, vCmax, findtype);
}

void clVariablesLS::bcSetdXArr(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
	int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vC2)
{
	if (ii0.x >= p.nX || ii0.x < 0)
		return;
	if (ii0.y >= p.nY || ii0.y < 0)
		return;

	ii1.x = MOD(ii1.x, p.nX);
	if (ii1.y >= p.nY || ii1.y < 0)
		return;

	if (nType(ii0.x, ii0.y) & M_SOLID_NODE)
		return;
	if ((nType(ii1.x, ii1.y) & M_SOLID_NODE) == 0)
		return;

	dXArr(ii0.x, ii0.y, dir - 1) = MIN(dist, dXArr(ii0.x, ii0.y, dir - 1));
}

void clVariablesLS::bcSetdXArr0(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
	int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vC2)
{
	if (ii0.x >= p.nX || ii0.x < 0)
		return;
	if (ii0.y >= p.nY || ii0.y < 0)
		return;

	ii1.x = MOD(ii1.x, p.nX);
	if (ii1.y >= p.nY || ii1.y < 0)
		return;

	if (vls.nType(ii0.x, ii0.y) & M0_SOLID_NODE)
		return;
	if ((vls.nType(ii1.x, ii1.y) & M0_SOLID_NODE) == 0)
		return;
	
	dXArr0(ii0.x, ii0.y, dir - 1) = MIN(dist, dXArr0(ii0.x, ii0.y, dir - 1));
}


void clVariablesLS::calcDxTerms()
{
	/// Doing it this way will allow for more complicated corrections to be made
	/// without worrying about effects on speed.

	

	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (nType(i, j) & M_SOLID_NODE)
				continue;


			double dx_e = dXArr(i, j, 0), dx_w = dXArr(i, j, 1), dx = dx_e + dx_w;
			double dy_n = dXArr(i, j, 2), dy_s = dXArr(i, j, 3), dy = dy_n + dy_s;

			double Xe_coeff = dx_w / (dx_e * dx), Xw_coeff = -dx_e / (dx_w * dx), Xc_coeff = (dx_e - dx_w) / (dx_e * dx_w);
			double Yn_coeff = dy_s / (dy_n * dy), Ys_coeff = -dy_n / (dy_s * dy), Yc_coeff = (dy_n - dy_s) / (dy_n * dy_s);

			double Xe2_coeff = 2. / (dx_e * dx), Xw2_coeff = 2. / (dx_w * dx), Xc2_coeff = -2. / (dx_e * dx_w);
			double Yn2_coeff = 2. / (dy_n * dy), Ys2_coeff = 2. / (dy_s * dy), Yc2_coeff = -2. / (dy_n * dy_s);


			dXCoeffs(i, j, dxind_e) = Xe_coeff;
			dXCoeffs(i, j, dxind_w) = Xw_coeff;
			dXCoeffs(i, j, dxind_c) = Xc_coeff;
			dXCoeffs(i, j, dyind_n) = Yn_coeff;
			dXCoeffs(i, j, dyind_s) = Ys_coeff;
			dXCoeffs(i, j, dyind_c) = Yc_coeff;

			dXCoeffs(i, j, dx2ind_e) = Xe2_coeff;
			dXCoeffs(i, j, dx2ind_w) = Xw2_coeff;
			dXCoeffs(i, j, dx2ind_c) = Xc2_coeff;
			dXCoeffs(i, j, dy2ind_n) = Yn2_coeff;
			dXCoeffs(i, j, dy2ind_s) = Ys2_coeff;
			dXCoeffs(i, j, dy2ind_c) = Yc2_coeff;
		}
	}
}


void clVariablesLS::createKernels()
{

}

#ifndef IN_KERNEL_IBB
cl_bool2 clVariablesLS::fillBoundaryArray()
{
	bFlag.fastfill(0);
	cl_bool2 resizeFlag = { { CL_FALSE, CL_FALSE} };
	ssArrIndMap.fastfill(-1);
	for (int i = 0; i < p.nX; i++)
	{
		ssArrInds(i) = ssArr.curSize();
		for (int j = 0; j < p.nY; j++)
		{
			if (nType(i, j) & M_SOLID_NODE)
				continue;

			for (int k = 0; k < 8; k++)
			{
				if (dXArr(i, j, k) != 1. && dXArr(i, j, k) > 0.)
				{
					resizeFlag.x |= ibbArr.append(i + j * p.XsizeFull +
						p.distSize * (k + 1));

					if (k < 4 && (i >= minSSArrNodeX) &&
						(i <= maxSSArrNodeX) && ((bFlag(i, j) & SHEAR_NODE_EXISTS) == 0))
					{
						ssArrIndMap(i, j) = ssArr.curSize();
						resizeFlag.y |= ssArr.append(i + j * p.XsizeFull);
						bFlag(i, j) |= SHEAR_NODE_EXISTS;
					}
				}
			}
		}
	}

	return resizeFlag;
}
#else
bool clVariablesLS::fillBoundaryArray()
{
	bFlag.fastfill(0);
	bool resizeFlag = false;
	ssArrIndMap.fastfill(-1);
	for (int i = 0; i < p.nX; i++)
	{
		ssArrInds(i) = ssArr.curSize();
		for (int j = 0; j < p.nY; j++)
		{
			if (nType(i, j) & M_SOLID_NODE)
				continue;

			for (int k = 0; k < 8; k++)
			{
				if (dXArr(i, j, k) != 1. && dXArr(i, j, k) > 0.)
				{
					if (k < 4 && (i >= minSSArrNodeX) &&
						(i <= maxSSArrNodeX) && ((bFlag(i, j) & SHEAR_NODE_EXISTS) == 0))
					{
						ssArrIndMap(i, j) = ssArr.curSize();
						resizeFlag |= ssArr.append(i + j * p.XsizeFull);
						bFlag(i, j) |= SHEAR_NODE_EXISTS;
					}
				}
			}
		}
	}

	return resizeFlag;
}
#endif

//TODO: figure out what arrays can be added here
void clVariablesLS::freeHostArrays()
{

}

void clVariablesLS::ini()
{
	if(!restartRunFlag)
	{
		for (int i = 0; i < nN; i++)
			C(i) = C0(i);
	}
	// save C0 (only needs to be done once)
	C0.savetxt();
		
	iniFillMap();
	iniFillMap0();
	iniCountSolid();
	inidXArrays();
	fillBoundaryArray();
	iniNodeBoundaryInfo();
	iniLSMap();

	if (vlb.kOmegaClass.kOmegaSolverFlag || vfd.thermalSolverFlag)
		calcDxTerms();


	allocateBuffers();

	setSourceDefines();
	std::function<void(void)> createKerPtr = std::bind(&clVariablesLS::createKernels, this);
	std::function<void(void)> setArgsPtr = std::bind(&clVariablesLS::setKernelArgs, this);

	sourceGenerator::SourceInstance()->addIniFunction(createKerPtr, setArgsPtr);

	
	if (saveMacroStart & !restartRunFlag)
		save2file();
}


void clVariablesLS::iniCountSolid()
{
	int nS = 0;
	int nSf = 0;


	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (nType(i, j) & M_SOLID_NODE)
			{
				nS++;
				if (nType(i, j) & M0_FLUID_NODE)
				{
					nSf++;
				}
			}
		}
	}

	int nFluid = p.Channel_Height * p.nX - nSf;
	Masses(0) = (double)nFluid;
	Masses(1) = (double)nSf;
}



void clVariablesLS::iniFillMap()
{
	// Set all M0_FLUID_NODES to M_FLUID_NODE
	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if(nType(i, j) & M0_FLUID_NODE)
				nType(i, j) |= M_FLUID_NODE;
			else if (restartRunFlag == false)
			{// if not a restart, we can set all nodes
			// same as those found in inifillMap0 and
			// dont need to complete rest of function
				nType(i, j) |= M_SOLID_NODE;
			}

		}
	}

	// Return if restarting run
	if (restartRunFlag == false)
	{
		return;
	}

	for (int n = 0; n < nBL; n++)
	{
		int n0 = BL(n, 0), n1 = BL(n, 1);
		cl_double2 c0 = { { vls.C(n0).x, vls.C(n0).y } }, c1 = { { vls.C[n1].x, vls.C[n1].y } };
		cl_double2 cn = { { c0.y - c1.y, c1.x - c0.x } };
		cn = VNORM(cn);
		cl_int2 Cmax = max2(c0, c1), Cmin = min2(c0, c1);
		double maxx = c0.x;
		if (c1.x > maxx) maxx = c1.x;
		double minx = c0.x;
		if (c1.x < minx) minx = c1.x;

		if (Cmax.x < 0 || Cmin.x >= p.nX)
			continue;

		if (Cmin.x < 0)
			Cmin.x = 0;

		if (Cmax.x >= p.nX)
			Cmax.x = p.nX - 1;

		if (Cmin.y < 0)
			Cmin.y = 0;

		if (Cmax.y >= p.nY)
			Cmax.y = p.nY - 1;




		for (int i = Cmin.x; i < Cmax.x + 1; i++)
		{
			for (int j = Cmin.y; j < Cmax.y + 1; j++)
			{
				cl_double2 L0 = { { (double)i, (double)j } };
				if (i < minx || i > maxx)
					continue;

				cl_double2 vL0 = Subtract2(L0, c0);
				//if dot product is less than 0, this node is a solid
				if (DOT_PROD(vL0, cn) < 0)
				{	// If solid now, we need to unset M_FLUID_NODE bit,
					// then set M_SOLID_NODE bit
					nType(i, j) &= !M_FLUID_NODE;
					nType(i, j) |= M_SOLID_NODE;
				}
			}
		}
	}

	for (int i = 0; i < p.nX; i++)
	{
		int flag = 0;
		int ind = 0;
		
		// starting at y = 0, find first M_SOLID_NODE in column i
		while (flag == 0)
		{
			// break when reaching first M_SOLID_NODE
			if (nType(i, ind) & M_SOLID_NODE)
				flag = 1;
			else
			{	// Set all nodes to M_SOLID_NODE until reaching 
				// first solid node set along wall. again need
				// to make sure to remove fluid flag
				nType(i, ind) &= !M_FLUID_NODE;
				nType(i, ind) |= M_SOLID_NODE;
			}

			ind++;
		}

		// Skip ahead 5 nodes in y direction
		ind += 5;
		flag = 0;
		// continue testing nodes until next M_SOLID_NODE is 
		// found, which will be located on outside of upper wall
		while (flag == 0 || ind >= p.nY)
		{
			// Once next M_SOLID_NODE is found, break
			if (nType(i, ind) & M_SOLID_NODE)
				flag = 1;

			ind++;
		}

		// for remaining nodes above this M_SOLID_NODE,
		// set to M_SOLID_NODE, while removing M_FLUID_NODE
		// bit
		for (int j = ind; j < p.nY; j++)
		{
			nType(i, j) &= !M_FLUID_NODE;
			nType(i, j) = M_SOLID_NODE;
		}
	}
}

void clVariablesLS::iniFillMap0()
{
	nType.fill(M0_FLUID_NODE);


	for (int n = 0; n < nBL; n++)
	{
		int n0 = BL(n, 0), n1 = BL(n, 1);
		cl_double2 c0 = { { vls.C0(n0).x, vls.C0(n0).y } }, c1 = { { vls.C0[n1].x, vls.C0[n1].y } };
		cl_double2 cn = { { c0.y - c1.y, c1.x - c0.x } };
		cn = VNORM(cn);
		cl_int2 Cmax = max2(c0, c1), Cmin = min2(c0, c1);
		double maxx = c0.x;
		if (c1.x > maxx) maxx = c1.x;
		double minx = c0.x;
		if (c1.x < minx) minx = c1.x;

		if (Cmax.x < 0 || Cmin.x >= p.nX)
			continue;

		if (Cmin.x < 0)
			Cmin.x = 0;

		if (Cmax.x >= p.nX)
			Cmax.x = p.nX - 1;

		if (Cmin.y < 0)
			Cmin.y = 0;

		if (Cmax.y >= p.nY)
			Cmax.y = p.nY - 1;


		for (int i = Cmin.x; i < Cmax.x + 1; i++)
		{
			for (int j = Cmin.y; j < Cmax.y + 1; j++)
			{
				cl_double2 L0 = { { (double)i, (double)j } };
				if (i < minx || i > maxx)
					continue;

				cl_double2 vL0 = Subtract2(L0, c0);
				if (DOT_PROD(vL0, cn) < 0)
					vls.nType(i, j) = M0_SOLID_NODE;
			}
		}
	}

	for (int i = 0; i < p.nX; i++)
	{
		int flag = 0;
		int ind = 0;
		while (flag == 0)
		{
			if (nType(i, ind) == M0_SOLID_NODE)
				flag = 1;
			else
				nType(i, ind) = M0_SOLID_NODE;

			ind++;
		}

		ind += 5;
		flag = 0;
		while (flag == 0 || ind >= p.nY)
		{
			if (nType(i, ind) == M0_SOLID_NODE)
				flag = 1;

			ind++;
		}

		for (int j = ind; j < p.nY; j++)
			nType(i, j) = M0_SOLID_NODE;
	}
}


void clVariablesLS::iniLSMap()
{
	for (int ibot = 0; ibot < nBL/2; ibot++)
	{
		// if part of BL falls inside domain,
		// calculate the i node that falls in between two
		// ends of BL and set the BL index to the corresponding
		// location in lsMap.
		double cx1 = C0(BL(ibot, 0)).x, cx2 = C0(BL(ibot, 1)).x;
		if (!(cx2 < 0. || cx1 >= (double)p.nX - 1.))
		{
			int iloc = (int)ceil(cx1);
			lsMap(iloc, 0) = (cl_ushort)ibot;
		}

		int itop = ibot + nBL / 2;

		// Do same calculation but for top wall
		cx2 = C0(BL(itop, 0)).x, cx1 = C0(BL(itop, 1)).x;
		if (!(cx2 < 0. || cx1 >= (double)p.nX - 1.))
		{
			int iloc = (int)ceil(cx1);
			lsMap(iloc, 1) = (cl_ushort)itop;
		}
	}
}

void clVariablesLS::iniNodeBoundaryInfo()
{
	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (nType(i, j) & M_SOLID_NODE)
				continue;

			for (int m = 1; m < 9; m++)
			{
				int iis = MOD(i + vlb.Cx[m], p.nX);
				int jjs = j + vlb.Cy[m];

				if (jjs < 0 || jjs >= p.nY)
				{// want to avoid indexing out of bounds
					nType(i, j) |= vlb.boundsArr[m];
#ifdef IN_KERNEL_IBB
					// at edge (which should never happen), dy <= 0.5
					nType(i,j) |= vlb.boundsArrT1[m];
#endif
					continue;
				}

				if (vls.nType(iis, jjs) & M_SOLID_NODE)
				{
					nType(i, j) |= vlb.boundsArr[m];
#ifdef IN_KERNEL_IBB
					if (vls.dXArr(i, j, m - 1) <= 0.5)
						nType(i,j) |= vlb.boundsArrT1[m];
					else if (vls.dXArr(i, j, m - 1) > 0.5)
						nType(i,j) |= vlb.boundsArrT2[m];
					else
						printf("dir at (%d,%d) points to solid node, but dx = 1\n", i, j);
#else
					if (vls.dXArr(i, j, m - 1) == 1.0)
						printf("dir at (%d,%d) points to solid node, but dx = 1\n", i, j);
#endif
				}
			}
		}
	}
}


void clVariablesLS::loadParams()
{
	nN = p.getParameter<int>("nN");
	nBL = nN - 2;

	restartRunFlag = p.getParameter<bool>("Restart Run", false);
	lsSpacing = p.getParameter<double>("LS Spacing", LS_SPACING);
	saveMacroStart = p.getParameter<bool>("Save Macro Start", true);
	restartRunFlag &= testRestartRun();
}

cl_int2 clVariablesLS::min2(cl_double2 v0, cl_double2 v1)
{
	double x = v0.x, y = v0.y;
	if (v1.x < x) x = v1.x;
	if (v1.y < y) y = v1.y;
	return{ { (int)floor(x), (int)floor(y) } };
}

cl_int2 clVariablesLS::max2(cl_double2 v0, cl_double2 v1)
{
	double x = v0.x, y = v0.y;
	if (v1.x > x) x = v1.x;
	if (v1.y > y) y = v1.y;
	return{ { (int)ceil(x), (int)ceil(y) } };
}

void clVariablesLS::renameSaveFiles()
{
	C.RenameTxtFile();
	nType.RenameTxtFile();

}

void clVariablesLS::saveRestartFiles()
{
	C.save_bin_from_device("lsc");
}

void clVariablesLS::save2file()
{
	C.save_txt_from_device();
	nType.save_txt_from_device();
	//M.save_txt_from_device();
}

void clVariablesLS::saveDebug()
{
	ibbArr.savetxt();
	dXArr.savetxt();
	dXArr0.savetxt();
}

void clVariablesLS::saveParams()
{
	p.setParameter("nN", nN);
	p.setParameter("LS Spacing", lsSpacing);
	if (p.Time > 0)
		p.setParameter<bool>("Restart Run", true);
	else
		p.setParameter<bool>("Restart Run", false);
	p.getParameter<bool>("Save Macro Start", saveMacroStart);
}

void clVariablesLS::setKernelArgs()
{

}


void clVariablesLS::setSourceDefines()
{
	int ind_ = nN / 2 - 10;
	while (true)
	{
		if ((C(ind_).y - C(ind_ - 1).y) > p.Pipe_radius)
			break;
		ind_++;
	}

	minSSArrNodeX = MAX(vtr.xReleasePos - 10, 0);
	maxSSArrNodeX = vtr.xStopPos + 10;

	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "MIN_SS_ARR_NODE_X", minSSArrNodeX);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "MAX_SS_ARR_NODE_X", maxSSArrNodeX);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "LSC_START_TOP", ind_);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "LSC_NN", nN);
#ifdef IN_KERNEL_IBB
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NTYPE_TYPE", "int");
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "IN_KERNEL_IBB");
#else
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NTYPE_TYPE", "short");
#endif

}

bool clVariablesLS::testInside(cl_double2 vd, cl_double2 v0, cl_double2 v1)
{
	cl_double2 v10 = Subtract2(v1, v0);
	cl_double2 vd0 = Subtract2(vd, v0);
	cl_double2 vd1 = Subtract2(vd, v1);

	if (fabs(GETLEN(v10) - GETLEN(vd0) - GETLEN(vd1)) < p.eps)
		return true;
	return false;
}

bool clVariablesLS::testPeriodicIndex(cl_int2 i0, cl_int2 i1)
{
	if (i0.y != i1.y)
		return false;
	if (i0.x != i1.x)
		return false;

	return true;
}


bool clVariablesLS::testRestartRun()
{
	allocateArrays();
	bool ret = C.load("load" SLASH "lsc");
	return ret;
}

// TODO: Replace with kernel based on legacy code
//		implementation
void clVariablesLS::updatedXArr()
{
	dXArr.fill(0.);

	for (int j = 0; j < p.nY; j++)
	{
		for (int i = 0; i < p.nX; i++)
		{
			if (nType(i, j) & M_FLUID_NODE)
			{
				for (int k = 0; k < 8; k++)
				{
					dXArr(i, j, k) = 1.;
				}
			}
		}
	}

	for (int i = 0; i < nBL; i++)
	{
		bcFindNodes(i, find_dx);
	}
	dXArr.copy_to_buffer();
}
// TODO: make sure that all kernel arguments do not need to be
//		reset when reallocating buffer size
// TODO: find out which kernels need to be updated with size of
//		ibbArr.
// TODO: if dXArr is updated on device, can try to update IBBArr
//		as well. If not, make sure dXArr update kernel is called
//		before calling this function, and copy dXArr to host.
void clVariablesLS::updateBoundaryArrays()
{
	dXArr.read_from_buffer();

	ssArr.zeros();

#ifndef IN_KERNEL_IBB
	ibbArr.zeros();
	cl_bool2 reallocFlag = fillBoundaryArray();
	ssArrIndMap.copy_to_buffer();
	if (reallocFlag.x)
	{
		ibbArr.reallocate_device_dynamic();
	}
	if (reallocFlag.y)
	{
		ssArr.reallocate_device_dynamic();
	}

	ibbArr.copy_dynamic_to_buffer();
	vlb.updateIBBArrays(reallocFlag.x);
	ssArr.copy_dynamic_to_buffer();
	ssArrInds.copy_to_buffer();
	vtr.wallShear.updateShearArrays(reallocFlag.y);
#else
	bool reallocFlag = fillBoundaryArray();
	ssArrIndMap.copy_to_buffer();
	if (reallocFlag)
	{
		ssArr.reallocate_device_dynamic();
	}
	ssArr.copy_dynamic_to_buffer();
	ssArrInds.copy_to_buffer();
	vtr.wallShear.updateShearArrays(reallocFlag);
#endif
}



const int clVariablesLS::LF[9] = {	LB_NO_BOUNDARY_LINK,
									LB_BOUNDARY_LINK_1,
									LB_BOUNDARY_LINK_2,
									LB_BOUNDARY_LINK_3,
									LB_BOUNDARY_LINK_4,
									LB_BOUNDARY_LINK_5,
									LB_BOUNDARY_LINK_6,
									LB_BOUNDARY_LINK_7,
									LB_BOUNDARY_LINK_8 };

