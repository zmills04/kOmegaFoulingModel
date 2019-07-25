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



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////                                           //////////////
/////////////           ALLOCATION FUNCTIONS            //////////////
/////////////                                           //////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


void clVariablesLS::allocateArrays()
{
	ssArrIndMap.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	nType.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	nTypePrev.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	BL.zeros(nBL, 2);
	C0.zeros(nN);
	C.zeros(nN);
	dXArr.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 4, 4);
	dXArr0.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 4, 4);

	ssArrInds.zeros(p.nX);
	lsMap.zeros(p.nX, p.XsizeFull, 2, 2);

	Masses.allocate(2);

#ifndef IN_KERNEL_IBB
	ibbArrCurIndex.zeros(1);
#endif

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
	ssArrInds.allocate_buffer_w_copy();
	ssArrIndMap.allocate_buffer_w_copy();
	lsMap.allocate_buffer_w_copy(CL_MEM_READ_ONLY);

#ifndef IN_KERNEL_IBB
	ibbArr.allocate_buffer();
	ibbArr.copy_dynamic_to_buffer();
	ibbDistArr.allocate_buffer_w_copy();
	ibbArrCurIndex.allocate_buffer();
#endif
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////                                           //////////////
/////////////               MISC FUNCTIONS              //////////////
/////////////                                           //////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



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

bool clVariablesLS::fillShearArray()
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

			for (int k = 0; k < 4; k++)
			{
				int ii = MODFAST(i + vlb.Cx[k + 1], p.nX);
				int jj = j + vlb.Cy[k + 1];

				if (nType(ii, jj) & M_SOLID_NODE)
				{
					if ((i >= minSSArrNodeX) &&
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

//TODO: figure out what arrays can be added here
void clVariablesLS::freeHostArrays()
{

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



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////                                           //////////////
/////////////         KERNEL CREATION FUNCTIONS         //////////////
/////////////                                           //////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


void clVariablesLS::createKernels()
{
	updateNType.create_kernel(GetSourceProgram, LBQUEUE_REF, "updateNType");
	updateNType.set_size(nN - 1, WORKGROUPSIZE_UPDATEM);

	updateBoundArr.create_kernel(GetSourceProgram, LBQUEUE_REF, "updateBoundaryNodes");
	updateBoundArr.set_size(nBL, WORKGROUPSIZE_UPDATEM);

#ifndef IN_KERNEL_IBB
	updateIBBOnly.create_kernel(GetSourceProgram, LBQUEUE_REF, "updateIBBOnly");
	updateIBBOnly.set_size(nN - 1, WORKGROUPSIZE_UPDATEM);
#endif
}

void clVariablesLS::setKernelArgs()
{
	updateNType.set_argument(0, nType.get_buf_add());
	updateNType.set_argument(1, C.get_buf_add());
	updateNType.set_argument(2, C0.get_buf_add());
	updateNType.set_argument(3, vlb.Ux_array.get_buf_add());
	updateNType.set_argument(4, vlb.Uy_array.get_buf_add());
	updateNType.set_argument(5, vlb.Ro_array.get_buf_add());
	updateNType.set_argument(6, vlb.kOmegaClass.Kappa.get_add_A());
	updateNType.set_argument(7, vlb.kOmegaClass.Kappa.get_add_b());
	updateNType.set_argument(8, vlb.kOmegaClass.Omega.get_add_A());
	updateNType.set_argument(9, vlb.kOmegaClass.Omega.get_add_b());
	updateNType.set_argument(10, vlb.kOmegaClass.Omega.get_add_IndArr());


	int ind = 0;
	BLinks::arrName BLArrList[] = { BLinks::P01Arr, BLinks::vNArr,
		BLinks::lenArr, BLinks::nodLocArr, BLinks::typeArr };
	vtr.BL.setBuffers(updateBoundArr, ind, BLArrList, 5);
	updateBoundArr.set_argument(5, C.get_buf_add());
	updateBoundArr.set_argument(6, C0.get_buf_add());
	updateBoundArr.set_argument(7, nType.get_buf_add());
	updateBoundArr.set_argument(8, dXArr.get_buf_add());
#ifndef IN_KERNEL_IBB
	updateBoundArr.set_argument(9, ibbArr.get_buf_add());
	updateBoundArr.set_argument(10, ibbDistArr.get_buf_add());
	updateBoundArr.set_argument(11, ibbArrCurIndex.get_buf_add());
	updateBoundArr.setOptionInd(12);
	int ibbarrsize_ = ibbArr.getBufferFullSize();
	updateBoundArr.set_argument(12, &ibbarrsize_);

	updateIBBOnly.set_argument(0, C.get_buf_add());
	updateIBBOnly.set_argument(1, nType.get_buf_add());
	updateIBBOnly.set_argument(2, ibbArr.get_buf_add());
	updateIBBOnly.set_argument(3, ibbDistArr.get_buf_add());
	updateIBBOnly.set_argument(4, ibbArrCurIndex.get_buf_add());
	updateIBBOnly.setOptionInd(5);
	updateIBBOnly.set_argument(5, &ibbarrsize_);
#endif


	// Calls kernel, so must go here after kernel has been initialized
	iniIBBArrays();

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

	// These are used by updateM kernel for indexing
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "LSC_UPDATE_M_NN", nN - 1);
	int nodeskip = (nN - 1) / 2 - 1;
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "LSC_UPDATE_M_SKIP", nodeskip);



#ifdef IN_KERNEL_IBB
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NTYPE_TYPE", "int");
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "IN_KERNEL_IBB");
#else
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NTYPE_TYPE", "short");
#endif
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////                                           //////////////
/////////////             PARAMETER FUNCTIONS           //////////////
/////////////                                           //////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void clVariablesLS::loadParams()
{
	nN = p.getParameter<int>("nN");
	nBL = nN - 2;

	restartRunFlag = p.getParameter<bool>("Restart Run", false);
	lsSpacing = p.getParameter<double>("LS Spacing", LS_SPACING);
	saveMacroStart = p.getParameter<bool>("Save Macro Start", true);
	restartRunFlag &= testRestartRun();
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

bool clVariablesLS::testRestartRun()
{
	allocateArrays();
	bool ret = C.load("load" SLASH "lsc");
	return ret;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////                                           //////////////
/////////////              OUTPUT FUNCTIONS             //////////////
/////////////                                           //////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

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


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////                                           //////////////
/////////////         INITIALIZATION FUNCTIONS          //////////////
/////////////                                           //////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void clVariablesLS::ini()
{
	if(!restartRunFlag)
	{
		for (int i = 0; i < nN; i++)
			C(i) = C0(i);
	}
	// save C0 (only needs to be done once)
	C0.savetxt();

	iniFillMap0();
	iniFillMap();

	// Copy contents of nType to nTypePrev before calling iniNodeBoundaryInfo
	memcpy(nTypePrev.get_array(), nType.get_array(), nType.getFullSize() * sizeof(double));

	iniCountSolid();
	inidXArrays();
	iniShearArray();
	iniNodeBoundaryInfo();
	iniLSMap();

	allocateBuffers();

	// This uses a kernel to initialize arrays, so must be called from
	// setKernelArgs function.
	//iniIBBArrays();

	sumFluid.ini(nType, restartRunFlag, "redUx");

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


void clVariablesLS::inidXArrays()
{
	dXArr0.fill(0.);

	for (int j = 0; j < p.nY; j++)
	{
		for (int i = 0; i < p.nX; i++)
		{
			if (nType(i, j) & M_FLUID_NODE)
			{
				for (int k = 0; k < 4; k++)
				{
					dXArr0(i, j, k) = 1.;
				}
			}
		}
	}

	for (int i = 0; i < nBL; i++)
	{
		bcFindNodes(i, find_dx0);
	}

	// if new run, dXArr = dXArr0 so memory can be copied
	if (!restartRunFlag)
	{
		memcpy(dXArr.get_array(), dXArr0.get_array(), dXArr.getFullSize() * sizeof(double));
		return;
	}

	//otherwise, use updatedXArr()
	updatedXArr();
}



void clVariablesLS::iniFillMap()
{
	// Set all M0_FLUID_NODES to M_FLUID_NODE
	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			// Fills all M0_FLUID_NODES as 
			// M_FLUID_NODES as well
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
				{	// If solid, make sure that M_FLUID_NODE bit is set to
					// zero (will be 1 if M0_FLUID_NODE is 1) then set
					// M_SOLID_NODE bit to 1
					RESET_NODE_TO_M_SOLID(nType(i, j));
				}
			}
		}
	}

	// Now that boundary nodes in solid are set to M_SOLID_NODE, set remaining 
	// nodes inside the solid to M_SOLID_NODE while unsetting M_FLUID_NODE bit to 0

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
				RESET_NODE_TO_M_SOLID(nType(i, ind));
			}

			ind++;
		}

		// Skip ahead 5 nodes in y direction
		ind += 5;
		flag = 0;
		// continue testing nodes until next M_SOLID_NODE is 
		// found, which will be located on outside of upper wall
		while (flag == 0 && ind < p.nY)
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
			RESET_NODE_TO_M_SOLID(nType(i, j));
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
				{
					RESET_NODE_TO_M0_SOLID(vls.nType(i, j));
				}
			}
		}
	}

	for (int i = 0; i < p.nX; i++)
	{
		int flag = 0;
		int ind = 0;
		while (flag == 0)
		{
			if (nType(i, ind) & M0_SOLID_NODE)
				flag = 1;
			else
			{
				RESET_NODE_TO_M0_SOLID(vls.nType(i, ind));
			}
			ind++;
		}

		ind += 5;
		flag = 0;
		while (flag == 0 && ind < p.nY)
		{
			if (nType(i, ind) & M0_SOLID_NODE)
				flag = 1;

			ind++;
		}

		for (int j = ind; j < p.nY; j++)
		{
			RESET_NODE_TO_M0_SOLID(vls.nType(i, j));
		}
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
	// fill boundary info for unfouled domain first
	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (nType(i, j) & M0_SOLID_NODE)
				continue;

			for (int m = 1; m < 9; m++)
			{
				int iis = MOD(i + vlb.Cx[m], p.nX);
				int jjs = j + vlb.Cy[m];

				ERROR_CHECKING((jjs < 0 || jjs >= p.nY), "Fluid node located "\
					"at either y = 0 or y = nY-1. These two rows must only "\
					"contain solid nodes", ERROR_INITIALIZING_VLS);

				if (vls.nType(iis, jjs) & M0_SOLID_NODE)
				{
					nType(i, j) |= BOUNDARY_NODE0;
					if (m < 5)
					{
						nType(iis, jjs) |= SOLID_BOUNDARY_NODE0;
						nType(i, j) |= NESW_BOUNDARY_NODE0;
					}
				}
			}
		}
	}
	   	  
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

				ERROR_CHECKING((jjs < 0 || jjs >= p.nY), "Fluid node located "\
					"at either y = 0 or y = nY-1. These two rows must only "\
					"contain solid nodes", ERROR_INITIALIZING_VLS);

				if (vls.nType(iis, jjs) & M_SOLID_NODE)
				{
					nType(i, j) |= vlb.boundsArr[m];
					if (m < 5)
					{
						nType(iis, jjs) |= SOLID_BOUNDARY_NODE;
					}
#ifdef IN_KERNEL_IBB
					if (vls.dXArr(i, j, m - 1) <= 0.5)
						nType(i,j) |= vlb.boundsArrT1[m];
					else if (vls.dXArr(i, j, m - 1) > 0.5 && vls.dXArr(i, j, m - 1) < 1.0)
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


void clVariablesLS::iniShearArray()
{
	ssArr.reset();
	fillShearArray();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////                                           //////////////
/////////////              UPDATE FUNCTIONS             //////////////
/////////////                                           //////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void clVariablesLS::updateBoundaryArrays()
{
	updateBoundArr.call_kernel();
	ibbArrCurIndex.read_from_buffer(LBQUEUE_REF);

	bool resizeFlag = false;
	if (ibbArrCurIndex(0) >= ibbArr.getBufferFullSize())
	{
		ibbArrCurIndex.FillBuffer(0, LBQUEUE_REF);
		resizeFlag = true;
		ibbArr.resize_device_dynamic();
		ibbDistArr.resize_device_dynamic();

		updateIBBOnly.set_argument(2, ibbArr.get_buf_add());
		updateIBBOnly.set_argument(3, ibbDistArr.get_buf_add());
		int ibbarrsize_ = ibbArr.getBufferFullSize();
		updateIBBOnly.set_argument(5, &ibbarrsize_);

		updateIBBOnly.call_kernel();

		ibbArrCurIndex.read_from_buffer(LBQUEUE_REF);
	}


	int ibbnumel_ = ibbArrCurIndex(0);

	vlb.ibbKernel.set_argument(3, &ibbnumel_);

	// dont need to update these kernels until after updateIBBOnly
	// has finished updating ibb arrays since they will not be
	// called again until next update step.
	if (resizeFlag)
	{
		updateBoundArr.set_argument(9, ibbArr.get_buf_add());
		updateBoundArr.set_argument(10, ibbDistArr.get_buf_add());
		int ibbarrsize_ = ibbArr.getBufferFullSize();
		updateBoundArr.set_argument(12, &ibbarrsize_);


		vlb.ibbKernel.set_argument(0, vls.ibbArr.get_buf_add());
		vlb.ibbKernel.set_argument(1, vls.ibbDistArr.get_buf_add());
	}
}

void clVariablesLS::iniIBBArrays()
{
	ibbArrCurIndex.FillBuffer(0, LBQUEUE_REF);

	updateIBBOnly.call_kernel();
	
	ibbArrCurIndex.read_from_buffer(LBQUEUE_REF);

	bool resizeFlag = false;
	if (ibbArrCurIndex(0) >= ibbArr.getBufferFullSize())
	{
		ibbArrCurIndex.FillBuffer(0, LBQUEUE_REF);
		resizeFlag = true;
		ibbArr.resize_device_dynamic();
		ibbDistArr.resize_device_dynamic();

		updateIBBOnly.set_argument(2, ibbArr.get_buf_add());
		updateIBBOnly.set_argument(3, ibbDistArr.get_buf_add());
		int ibbarrsize_ = ibbArr.getBufferFullSize();
		updateIBBOnly.set_argument(5, &ibbarrsize_);

		updateIBBOnly.call_kernel();

		ibbArrCurIndex.read_from_buffer(LBQUEUE_REF);
	}

	int ibbnumel_ = ibbArrCurIndex(0);

	vlb.ibbKernel.set_argument(3, &ibbnumel_);

	// dont need to update these kernels until after updateIBBOnly
	// has finished updating ibb arrays since they will not be
	// called again until next update step.
	if (resizeFlag)
	{
		updateBoundArr.set_argument(9, ibbArr.get_buf_add());
		updateBoundArr.set_argument(10, ibbDistArr.get_buf_add());
		int ibbarrsize_ = ibbArr.getBufferFullSize();
		updateBoundArr.set_argument(12, &ibbarrsize_);


		vlb.ibbKernel.set_argument(0, vls.ibbArr.get_buf_add());
		vlb.ibbKernel.set_argument(1, vls.ibbDistArr.get_buf_add());
	}
}

void clVariablesLS::updateLS()
{
	// Prepares dXArr and nType for updating.
	nType.enqueue_copy_to_buffer(nTypePrev.get_buffer(), -1, LBQUEUE_REF);
	dXArr.enqueue_copy_to_buffer(dXArr0.get_buffer(), -1, LBQUEUE_REF);
	ibbArrCurIndex.FillBuffer(0, LBQUEUE_REF);

	// updates nType
	updateNType.call_kernel();

	// Saves updated nType to nTypePrev to use in next update step
	nTypePrev.enqueue_copy_to_buffer(nType.get_buffer(), -1, LBQUEUE_REF);

	updateBoundaryArrays();


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
}

void clVariablesLS::updateShearArrays()
{
	nType.read_from_buffer();
	// TODO: check if it is necessary to set this to zeros, or if
	// resetting the current index is sufficient
	ssArr.reset();

	bool reallocFlag = fillShearArray();
	ssArrIndMap.copy_to_buffer();
	if (reallocFlag)
	{
		ssArr.reallocate_device_dynamic();
	}
	ssArr.copy_dynamic_to_buffer();
	ssArrInds.copy_to_buffer();
	vtr.wallShear.updateShearArrays(reallocFlag);
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

