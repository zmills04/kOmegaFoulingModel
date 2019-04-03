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
	sf.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	s.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	M_o.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	M.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	BL.zeros(nBL, 2);
	
	C0.zeros(nN);
	C.zeros(nN);
}


void clVariablesLS::ini()
{
	ERROR_CHECKING(BL.load("load" SLASH "lsbl") == FALSE,
		"Could not load lsbl.txt", ERROR_INITIALIZING_VLS);
		
	ERROR_CHECKING(C0.load("load" SLASH "lsc0") == FALSE,
		"No lsc load file (initial locations provided in load/lsc0",
		ERROR_INITIALIZING_VLS);

	if (C.load("load" SLASH "lsc") == FALSE)
	{
		restartRunFlag = false;
		for (int i = 0; i < nN; i++)
			C(i) = C0(i);
	}



	C0.savetxt("lsc0");
	Masses.allocate(2);

	
	inifillMap();
	inifillMap0();
	inifillS();

	ini_dX_arrays();

	M_o.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	M.allocate_buffer_size(CL_MEM_READ_WRITE, p.FullSize);
	M.FillBuffer(0, p.FullSize);	// must be power of two for reduction step
	M.copy_to_buffer();
	sf.allocate_buffer_size(CL_MEM_READ_WRITE, p.FullSize);
	sf.FillBuffer(0, p.FullSize);
	sf.copy_to_buffer();
	s.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
	C.allocate_buffer_w_copy();
	C0.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	Masses.allocate_buffer_w_copy();
	dXArr.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	dXArrCur.allocate_buffer_w_copy(CL_MEM_READ_WRITE);

	if (saveMacroStart)
		save2file();
}



void clVariablesLS::inifillMap()
{
	M.fill(1);


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
				if (DOT_PROD(vL0, cn) < 0)
					vls.M(i, j) = LB_SOLID;
			}
		}
	}

	for (int i = 0; i < p.nX; i++)
	{
		int flag = 0;
		int ind = 0;
		while (flag == 0)
		{
			if (M(i, ind) == LB_SOLID)
				flag = 1;
			else
				M(i, ind) = LB_SOLID;

			ind++;
		}

		ind += 5;
		flag = 0;
		while (flag == 0 || ind >= p.nY)
		{
			if (M(i, ind) == LB_SOLID)
				flag = 1;

			ind++;
		}

		for (int j = ind; j < p.nY; j++)
			M(i, j) = LB_SOLID;

	}
}

void clVariablesLS::inifillMap0()
{
	M_o.fill(1);


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
					vls.M_o(i, j) = LB_SOLID;
			}
		}
	}

	for (int i = 0; i < p.nX; i++)
	{
		int flag = 0;
		int ind = 0;
		while (flag == 0)
		{
			if (M_o(i, ind) == LB_SOLID)
				flag = 1;
			else
				M_o(i, ind) = LB_SOLID;

			ind++;
		}

		ind += 5;
		flag = 0;
		while (flag == 0 || ind >= p.nY)
		{
			if (M_o(i, ind) == LB_SOLID)
				flag = 1;

			ind++;
		}

		for (int j = ind; j < p.nY; j++)
			M_o(i, j) = LB_SOLID;

	}
}

void clVariablesLS::inifillS()
{
	int nS = 0;
	int nSf = 0;


	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (M(i, j) == LB_SOLID)
			{
				s(i, j) = TRUE;
				nS++;
				if (M_o(i, j) != LB_SOLID)
				{
					sf(i, j) = TRUE;
					nSf++;
				}
			}
		}
	}

	int nL = p.Channel_Height*p.nX - nSf;
	Masses(0) = (double)nL;
	Masses(1) = (double)nSf;
}

void clVariablesLS::iniGL()
{
//	LSt_vbo.allocate(vls.nN);
//	LSt_vbo.AllocateColor(vls.nN / 2 * 3);
//	LSt_vbo.set_color_array({ { 0.f, 0.f, 0.f } });
//
//	LSb_vbo.allocate(vls.nN);
//	LSb_vbo.AllocateColor(vls.nN / 2 * 3);
//	LSb_vbo.set_color_array({ { 0.f, 0.f, 0.f } });
//
//	LSt0_vbo.allocate(vls.nN);
//	LSb0_vbo.allocate(vls.nN);
//
//	for (int i = 0; i < vls.nN / 2; i++)
//	{
//		LSb_vbo[i * 2] = vls.C[i].x;
//		LSb_vbo[i * 2 + 1] = vls.C[i].y;
//		LSb0_vbo[i * 2] = vls.C0[i].x;
//		LSb0_vbo[i * 2 + 1] = vls.C0[i].y;
//
//		LSt_vbo[i * 2] = vls.C[vls.nN - 1 - i].x;
//		LSt_vbo[i * 2 + 1] = vls.C[vls.nN - 1 - i].y;
//		LSt0_vbo[i * 2] = vls.C0[vls.nN - 1 - i].x;
//		LSt0_vbo[i * 2 + 1] = vls.C0[vls.nN - 1 - i].y;
//	}
//
//	vtr.BL.read_from_buffer(p.IOqueue);
//
//	for (int i = 0; i < vls.nBL / 2; i++)
//	{
//		vtr.BL[i].Color_ind = i;
//		vtr.BL[vls.nBL / 2 + i].Color_ind = -i;
//	}
//
//	vtr.BL.copy_to_buffer(p.IOqueue);
//
//	LSb_vbo.CreateVBO_position(LINE_VBO, p.context, CL_MEM_READ_WRITE);
//	LSb_vbo.CreateVBO_color(p.context, CL_MEM_READ_WRITE);
//
//	LSt_vbo.CreateVBO_position(LINE_VBO, p.context, CL_MEM_READ_WRITE);
//	LSt_vbo.CreateVBO_color(p.context, CL_MEM_READ_WRITE);
//
//	LSb0_vbo.CreateVBO_position(LINE_VBO);
//	LSb0_vbo.set_color_vector({ { 0.f, 0.f, 0.f } });
//
//	LSt0_vbo.CreateVBO_position(LINE_VBO);
//	LSt0_vbo.set_color_vector({ { 0.f, 0.f, 0.f } });
//
//	int gsize_GL = (int)ceil((double)vls.nN / 2. / WORKGROUPSIZE_UPDATE_GL) * WORKGROUPSIZE_UPDATE_GL;
//	vfl.update_GL_kernel.create_kernel(p.program, &p.LBqueue, "update_GL_wall");
//	vfl.update_GL_kernel.set_size_1D(gsize_GL, WORKGROUPSIZE_UPDATE_GL);
//
//	int ind = 0;
//	int num_nodes = vls.nN / 2;
//	vfl.update_GL_kernel.set_argument(ind++, sizeof(cl_mem), vls.C.get_buf_add());
//	vfl.update_GL_kernel.set_argument(ind++, sizeof(cl_mem), vls.LSb_vbo.get_buf_add());
//	vfl.update_GL_kernel.set_argument(ind++, sizeof(cl_mem), vls.LSt_vbo.get_buf_add());
//	vfl.update_GL_kernel.set_argument(ind++, sizeof(int), (void *)&num_nodes);
//
//#ifdef OPENGL_GRIDLINES
//	int numlinesH = p.nX + 1;
//	int numlinesV = p.nY + 1;
//	LinesV.allocate(numlinesV * 4);
//	LinesH.allocate(numlinesH * 4);
//	for (int i = 0; i < numlinesV; i++)
//	{
//		LinesV(i * 4) = 0.f + (float)i;
//		LinesV(i * 4 + 1) = 0.f;
//		LinesV(i * 4 + 2) = 0.f + (float)i;
//		LinesV(i * 4 + 3) = (float)p.nY;
//	}
//	for (int i = 0; i < numlinesH; i++)
//	{
//		LinesH(i * 4) = 0.;
//		LinesH(i * 4 + 1) = (float)i;
//		LinesH(i * 4 + 2) = (float)p.nX;
//		LinesH(i * 4 + 3) = (float)i;
//	}
//
//	LinesH.set_color_vector({ { 0.f, 0.f, 0.f } });
//	LinesV.set_color_vector({ { 0.f, 0.f, 0.f } });
//	LinesH.CreateVBO_position(LINE_VBO);
//	LinesV.CreateVBO_position(LINE_VBO);
//	LinesH.FreeHost();
//	LinesV.FreeHost();
//#endif
//
//	LSt_vbo.FreeHostColor();
//	LSt_vbo.FreeHost();
//	LSb_vbo.FreeHostColor();
//	LSb_vbo.FreeHost();
//	LSt0_vbo.FreeHost();
//	LSb0_vbo.FreeHost();
}


void clVariablesLS::freeHostMem()
{
	// No need to worry about host memory
}

void clVariablesLS::loadParams()
{
	restartRunFlag = p.getParameter<bool>("Restart Run", false);
	nN = p.getParameter<int>("nN");
	lsSpacing = p.getParameter<double>("LS Spacing", LS_SPACING);
	saveMacroStart = p.getParameter<bool>("Save Macro Start", true);
	nBL = nN - 2;
	allocateArrays();
}

void clVariablesLS::releaseObjects()
{
	//classes take care of this
}

void clVariablesLS::saveRestartFiles()
{
	C.save_bin_from_device("lsc");
}

void clVariablesLS::save2file()
{
	C.save_txt_from_device("lsc");
	M.save_txt_from_device("lbm");

//#ifdef _DEBUG
//	dXArr.savetxt_as_multi2D("dxArr", TRUE);
//#endif
}

void clVariablesLS::saveDebug()
{
	s.save_from_device("lbms");
	M_o.save_from_device("lbm0");
	sf.save_txt_from_device("lbmsf");
	//ii0_array.savetxt("ii0_array");
	//dir_array.savetxt("dir_array");
	//IBB_flag_fluid.savetxt("Ibb_flag");
	//ibb_coeff_fluid.savetxt("ibb_coeff");
	//Node_loc_fluid.savetxt("ibb_loc");
	//Neighs_elbm_fluid.savetxt("ibb_neighs");
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

/////////////Min Max Functions/////////////////

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

//////IBB functions//////////////////

//void clVariablesLS::ini_IBB_arrays()
//{
//	Bflag.zeros(p.nX, p.nY);
//
//	cur_el = 0;
//	num_el = nN;
//	length_ibb = nN;
//	fullsize_ibb_arrays = nN;
//	fullsize_Bnodes = nN / 2;
//	ii0_array.zeros(nN);
//	iis1_array.zeros(nN);
//	dir_array.zeros(nN);
//	D_array.zeros(nN);
//
//	ii0_array.allocate_buffer_w_copy(p.context, p.IOqueue, CL_MEM_READ_ONLY);
//	iis1_array.allocate_buffer_w_copy(p.context, p.IOqueue, CL_MEM_READ_ONLY);
//	dir_array.allocate_buffer_w_copy(p.context, p.IOqueue, CL_MEM_READ_ONLY);
//	D_array.allocate_buffer_w_copy(p.context, p.IOqueue, CL_MEM_READ_ONLY);
//
//	Bnodes.zeros(nN / 2);
//	Bnodes.allocate_buffer_w_copy(p.context, p.IOqueue, CL_MEM_READ_ONLY);
//	lengthBnodes = nN / 2;
//
//	update_IBB_arrays();
//
//	fullsize_ibb_arrays = length_ibb * 2;
//	fullsize_Bnodes = lengthBnodes * 2;
//
//
//	ii0_array.reallocate(fullsize_ibb_arrays, p.context, p.IOqueue);
//	iis1_array.reallocate(fullsize_ibb_arrays, p.context, p.IOqueue);
//	dir_array.reallocate(fullsize_ibb_arrays, p.context, p.IOqueue);
//	D_array.reallocate(fullsize_ibb_arrays, p.context, p.IOqueue);
//	Bnodes.reallocate(fullsize_Bnodes, p.context, p.IOqueue);
//
//	ii0_array.copy_to_buffer(p.IOqueue);
//	iis1_array.copy_to_buffer(p.IOqueue);
//	dir_array.copy_to_buffer(p.IOqueue);
//	D_array.copy_to_buffer(p.IOqueue);
//	Bnodes.copy_to_buffer(p.IOqueue);
//}

//void clVariablesLS::update_IBB_arrays()
//{
//	cur_el = 0;
//	curElBnodes = 0;
//	Bflag.fill(0);
//
//	for (int i = 0; i < nBL; i++)
//		bcFindNodes(i);
//
//	num_el = cur_el;
//	nBnodes = curElBnodes;
//
//	length_ibb = (int)ceil((double)num_el / WORKGROUPSIZE_IBB)*WORKGROUPSIZE_IBB;
//	lengthBnodes = (int)ceil((double)nBnodes / WORKGROUPSIZE_TR_SHEAR)*WORKGROUPSIZE_TR_SHEAR;
//}

void clVariablesLS::ini_dX_arrays()
{
	dXArr.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 8, 8);
	dXArrCur.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 8, 8);

	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (M_o(i, j) == LB_LIQUID1)
			{
				for (int k = 0; k < 8; k++)
				{
					dXArr(i, j, k) = 1.;
					dXArrCur(i, j, k) = 1.;
				}
			}

		}
	}

	for (int i = 0; i < nBL; i++)
	{
		if (i == 510 || i == 0 || i == 1)
			int a = 3;
		bcFindNodes(i, find_dx);
		bcFindNodes(i, find_dx_cur);
	}
}


void clVariablesLS::update_dXCur_array()
{
	for (int i = 0; i < nBL; i++)
	{
		bcFindNodes(i, find_dx);
		bcFindNodes(i, find_dx_cur);
	}
}


void clVariablesLS::bcFindNodes(int bl, int findtype)
{
	int n0 = BL(bl, 0), n1 = BL(bl, 1);
	if (n0 < 0)
		return;
	cl_double2 vC0 = { { C[n0].x, C[n0].y } };
	cl_double2 vC1 = { { C[n1].x, C[n1].y } };;
	if (findtype == find_dx)
	{
		vC0 = { { C0[n0].x, C0[n0].y } };
		vC1 = { { C0[n1].x, C0[n1].y } };
	}

	cl_double2 vCn = { { vC0.y - vC1.y, vC1.x - vC0.x } };
	cl_int2 vCmin = min2(vC0, vC1), vCmax = max2(vC0, vC1);

	if(vCmin.x < 0) vCmin.x = 0;
	if(vCmin.x > p.nX-1) return;
	if(vCmax.x < 0) return;
	if(vCmax.x > p.nX-1) vCmax.x = p.nX-1;

	if(vCmin.y < 0) vCmin.y = 0;
	if(vCmin.y > p.nY-1) return;
	if(vCmax.y < 0) return;
	if(vCmax.y > p.nY-1) vCmax.y = p.nY-1;


	for(int dir = 1; dir < 9; dir+=2)
		bcFindDirection(dir, bl, vC0, vC1, vCn, vCmin, vCmax, findtype);
}

void clVariablesLS::bcFindDirection(int dir, int bl, cl_double2 vC0, 
	cl_double2 vC1, cl_double2 vCn, cl_int2 vCmin, cl_int2 vCmax, int findtype)
{
	double dist;
	cl_double2 Cdir = vlb.Cxy_double[dir];
	double Cndir = DOT_PROD(vCn,Cdir);
	if(Cndir == 0.)
		return;
	if(Cndir > 0.)
	{
		dir = vlb.rev[dir];
		Cdir = vlb.Cxy_double[dir];
	}
	cl_int2 vCcut0, vCcut1;
	cl_double2 vCcut;
	if(vlb.Cx[dir] != 0)
	{
		int i = 0;
		if(vlb.Cx[dir] == 1) i = vCmin.x;
		if(vlb.Cx[dir] == -1) i = vCmax.x;

		for(int j = vCmin.y; j <= vCmax.y; j++)
		{
			if (bcFindIntersection(vCcut0, vCcut1, vCcut, dist, { { (double)i, (double)j } }, Cdir, vC0, vC1, vCn) == TRUE)
			{
				if (findtype == find_ibb)
					bcSetBoundaryNode(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
				if (findtype == find_dx)
					bcSetdXArr(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
				if (findtype == find_dx_cur)
					bcSetdXCurArr(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
			}
		}
	}
	if(vlb.Cy[dir] != 0)
	{
		int j = 0;
		if(vlb.Cy[dir] == 1) j = vCmin.y;
		if(vlb.Cy[dir] == -1) j = vCmax.y;

		for(int i = vCmin.x; i <= vCmax.x; i++)
		{
			if (bcFindIntersection(vCcut0, vCcut1, vCcut, dist, { { (double)i, (double)j } }, Cdir, vC0, vC1, vCn) == TRUE)
			{
				if (findtype == find_ibb)
					bcSetBoundaryNode(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
				if (findtype == find_dx)
					bcSetdXArr(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
				if (findtype == find_dx_cur)
					bcSetdXCurArr(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn);
			}
		}
	}
}


BOOL clVariablesLS::bcFindIntersection(cl_int2 &vCcut0, cl_int2 &vCcut1,
	cl_double2 &vCcut, double &dist, cl_double2 vL0, cl_double2 vLd,
	cl_double2 vC0, cl_double2 vC1, cl_double2 vCn)
{
	double distF;
	if(bcFindIntersectionLinePlane(vCcut, distF, vL0, vLd, vC0, vC1) == FALSE)
		return FALSE;

	if(testInside(vCcut, vC0, vC1) == FALSE)
		return FALSE;

	int n_cut_0 = (int)ceil(distF);
	int n_cut = (int)ceil(distF);
	double dist_0 = 1. - (n_cut - (double)distF);
	dist = 1. - (n_cut - (double)distF);

	vCcut1 = { { (int)(vL0.x + vLd.x * n_cut), (int)(vL0.y + vLd.y * n_cut) } };
	vCcut0 = { { vCcut1.x - (int)vLd.x, vCcut1.y - (int)vLd.y } };

	cl_int2 ii0 = { { MOD(vCcut0.x, p.nX), MOD(vCcut0.y, p.nY) } };

	if(ii0.x != vCcut0.x)
	{
		return FALSE; 
	}

	if(ii0.y != vCcut0.y)
		return FALSE;

	if(vls.M(ii0.x,ii0.y) == LB_SOLID)
	{
		if(dist < p.eps)
		{
			dist = 1.;
			n_cut--;
			vCcut1 = { { (int)(vL0.x + vLd.x * n_cut),
				(int)(vL0.y + vLd.y * n_cut) } };
			vCcut0 = { { vCcut1.x - (int)vLd.x,
				vCcut1.y - (int)vLd.y } };

			if(ii0.x != vCcut0.x)
				return FALSE;

			if(ii0.y != vCcut0.y)
				return FALSE;
		}
		else
			return FALSE;
	}

	return TRUE;
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

	if (vls.M(ii0.x, ii0.y) == LB_SOLID)
		return;
	if (vls.M(ii1.x, ii1.y) != LB_SOLID)
		return;

	dXArr(ii0.x, ii0.y, dir - 1) = dist;
}

void clVariablesLS::bcSetdXCurArr(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
	int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vC2)
{
	if (ii0.x >= p.nX || ii0.x < 0)
		return;
	if (ii0.y >= p.nY || ii0.y < 0)
		return;

	ii1.x = MOD(ii1.x, p.nX);
	if (ii1.y >= p.nY || ii1.y < 0)
		return;

	if (vls.M_o(ii0.x, ii0.y) == LB_SOLID)
		return;
	if (vls.M_o(ii1.x, ii1.y) != LB_SOLID)
		return;

	dXArrCur(ii0.x, ii0.y, dir - 1) = dist;
}





void clVariablesLS::bcSetBoundaryNode(cl_int2 ii0, cl_int2 ii1, int dir,
	double dist, int bl, cl_double2 vCc, cl_double2 vC0,
	cl_double2 vC1, cl_double2 vCn)
{
	if (cur_el == fullsize_ibb_arrays)
	{
		fullsize_ibb_arrays *= 1.3;
		ii0_array.reallocate(fullsize_ibb_arrays);
		iis1_array.reallocate(fullsize_ibb_arrays);
		dir_array.reallocate(fullsize_ibb_arrays);
		D_array.reallocate(fullsize_ibb_arrays);

	}

	if (curElBnodes == fullsize_Bnodes)
	{
		fullsize_Bnodes *= 1.3;
		Bnodes.reallocate(fullsize_Bnodes);
	}

	ii0 = MOD2(ii0, p.nn);

	int BMflag = vls.Bflag(ii0.x, ii0.y);


	if (!(BMflag & LFNearest) && dir < 5)
	{
		if (ii0.x >= imin_Bnode && ii0.x <= imax_Bnode)
		{
			if (bl < nBL / 2)
			{
				bot_ind_end = curElBnodes;
			}
			Bnodes(curElBnodes++) = { { ii0.x, ii0.y, cur_el } };
		}
	}


	if (BMflag & vlb.LF[dir])
		return;
	vls.Bflag(ii0.x, ii0.y) = BMflag | vlb.LF[dir];
	cl_int2 CC = vlb.Cxy[dir];
	char map = vls.M(ii0.x, ii0.y);
	int dirp = vlb.rev[dir];
	cl_int2 iis1t = Subtract2(ii0, CC);
	cl_int2 iis1 = MOD2(iis1t, p.nn);

	if (map != vls.M(iis1.x, iis1.y) || testperiodicindex(iis1t, iis1) == FALSE)
	{
		dist = 0.5;
	}

#ifdef INLET_OUTLET_BC
	if (iis1t.x != iis1.x)
		return;
#endif

	ii0_array(cur_el) = ii0;
	dir_array(cur_el) = dir;
	D_array(cur_el) = dist;
	iis1_array(cur_el) = iis1;
	cur_el++;
}

BOOL clVariablesLS::testperiodicindex(cl_int2 i0, cl_int2 i1)
{
	if(i0.y != i1.y)
		return FALSE;
	if (i0.x != i1.x)
		return FALSE;



	return TRUE;
}

BOOL clVariablesLS::bcFindIntersectionLinePlane(cl_double2 &vC, double &dist,
	cl_double2 vL0, cl_double2 vLd, cl_double2 vP0, cl_double2 vP1)
{
	cl_double2 vP10 = Subtract2(vP1, vP0);
	cl_double2 vN = { { -vP10.y, vP10.x } };
	double den = DOT_PROD(vN, vLd);
	if(den == 0.)
		return FALSE;
	cl_double2 vPL = Subtract2(vP0, vL0);
	dist = DOT_PROD(vN, vPL) / den;
	vC = { { vL0.x + vLd.x * dist, vL0.y + vLd.y * dist } };
	return TRUE;
}


BOOL clVariablesLS::testInside(cl_double2 vd, cl_double2 v0, cl_double2 v1)
{
	cl_double2 v10 = Subtract2(v1, v0);
	cl_double2 vd0 = Subtract2(vd, v0);
	cl_double2 vd1 = Subtract2(vd, v1);

	if(fabs(GETLEN(v10) - GETLEN(vd0) - GETLEN(vd1)) < p.eps)
		return TRUE;
	return FALSE;
}
