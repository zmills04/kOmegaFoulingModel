// clVariablesFD.cpp: implementation of the clVariables class.
//
// (c) Zachary Grant Mills, 2016 
//////////////////////////////////////////////////////////////////////

#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"


void clVariablesFD::allocateArrays()
{
	Alphat.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	Temp_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY);

	//dX_cur.allocate_buffer_w_copy(p.context, p.IOqueue);
	//dX.allocate_buffer_w_copy(p.context, p.IOqueue, CL_MEM_READ_ONLY);
	//Neighs.allocate_buffer_w_copy(p.context, p.IOqueue);
}

double clVariablesFD::calcAverageT()
{
	return sumTemp.reduceSingle();
}


void clVariablesFD::createKernels()
{
	if (vlb.kOmegaSolverFlag)
	{
		//TempUpdateCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF, "Update_T_Coeffs_Implicit_Turbulent");
		TempUpdateCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF, "Update_T_Coeffs_Implicit_Turbulent");
		TempUpdateCoeffs.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY, WORKGROUPSIZEY_LB);
	}
	else
	{
		TempUpdateCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF, "Update_T_Coeffs_Implicit");
		TempUpdateCoeffs.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY, WORKGROUPSIZEY_LB);
	}

//	SetSSCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF, "Set_SS_Temp_Coeffs");
//	SetSSCoeffs.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY, WORKGROUPSIZEY_LB);
}


void clVariablesFD::freeHostMem()
{
}


void clVariablesFD::ini()
{
	allocateArrays();
	setSourceDefines();

	if (thermalSolverFlag == false)
		return;

	sourceGenerator::SourceInstance()->addFile2Kernel("tfdKernels.cl");

	TempInds.ini(p.nX, p.XsizeFull, p.nY, p.nY, &vls.M_o);

	if (tempLoadedFlag == false)
	{
		Temp_array.fill(ROE0);
	}


	Temp.CreateSolver(&Temp_array, &TempInds, FDQUEUE_REF,
		tempMaxIters, tempMaxRelTol, tempMaxAbsTol);

	sumTemp.ini(*Temp.getMacroArray(), "redTemp");
	
	std::function<void(void)> createKerPtr = std::bind(&clVariablesFD::createKernels, this);
	std::function<void(void)> setArgsPtr = std::bind(&clVariablesFD::setKernelArgs, this);

	sourceGenerator::SourceInstance()->addIniFunction(createKerPtr, setArgsPtr);



//	if (!Restart_Run())
//	{
//		Restart_Run_Flag = 0;
//
//		ini_bounds();
//		for (int i = 0; i < p.nX; i++)
//		{
//			for (int j = 0; j < p.nY; j++)
//			{
//				for (int k = 0; k < 4; k++)
//				{
//					dX_cur_full(i, j, k) = dX_full(i, j, k);
//				}
//			}
//		}
//		Create_Coeffs();
//	}
//
//	Find_Neighbors();
//
//
//#ifdef CREATE_BIN_FILES
//	dX.savebin("fddx");
//#endif

	if (saveMacroStart)
		save2file();


}

void clVariablesFD::iniAlpha()
{
}

void clVariablesFD::loadParams()
{
	thermalSolverFlag = p.getParameter<bool>("Thermal Solver", USE_THERMAL_SOLVER);
	calcNuFlag = p.getParameter<bool>("Calculate Nu", CALC_NUSSELT);
	saveMacroStart = p.getParameter<bool>("Save Macros On Start", THERMAL_SAVE_MACROS_ON_START);
	calcSSTempFlag = p.getParameter<bool>("Solve Steady Temp", SOLVE_SS_TEMP);
	
	testRestartRun();
	
	double PrNum = p.getParameter<double>("Pr Number", PR_NUMBER);
	PrTurbNum = p.getParameter<double>("Pr Number", PR_TURB_NUMBER);
	Alpha_fluid = vlb.MuVal / PrNum;
	kSoot = p.getParameter<double>("K Soot", THERMAL_CONDUCTIVITY_FOUL);
	kAir = p.getParameter<double>("K air", THERMAL_CONDUCTIVITY_AIR);
	rhoSoot = p.getParameter<double>("Rho Soot", DENSITY_SOOT);
	cpSoot = p.getParameter<double>("Cp Soot", HEAT_CAPACITY_SOOT);
	Alpha_foul = (kSoot / rhoSoot / cpSoot)*(p.DELTA_L*p.DELTA_L / p.DELTA_T);
	
	T_Actual_Max = p.getParameter<double>("Tinlet Actual", LBT_ACTUAL_TEMP_MAX);
	T_Actual_Min = p.getParameter<double>("Twalls Actual", LBT_ACTUAL_TEMP_MIN);
	T_Actual_Diff = T_Actual_Max - T_Actual_Min;
	ROE0 = p.getParameter<double>("Tinitial", TFD_INI_TEMP);
	ROE_INX = p.getParameter<double>("Tinlet Sim", TFD_X_IN);

	
	double DENSE_SOOT = rhoSoot * (vlb.RhoVal / vlb.rhoAir);
	double HCAP_SOOT = cpSoot * p.DELTA_F*p.DELTA_L / p.DELTA_M;

	Alpha_den_soot = DENSE_SOOT*HCAP_SOOT;
	Alpha_num_air = (kAir * p.DELTA_F / p.DELTA_T);
	Alpha_num_soot = (kSoot * p.DELTA_F / p.DELTA_T);
	Alpha_den_air = Alpha_num_air / Alpha_fluid;
	Alpha_den_soot = DENSE_SOOT*HCAP_SOOT;

	tempMaxRelTol = p.getParameter<double>("Temp Max Rel Tol", THERMAL_MAX_REL_TOL);
	tempMaxAbsTol = p.getParameter<double>("Temp Max Abs Tol", THERMAL_MAX_ABS_TOL);
	tempMaxIters = p.getParameter<int>("Temp Max Iterations", THERMAL_MAX_ITERS);
}

//Releases buffers
void clVariablesFD::releaseObjects()
{

}

BOOL clVariablesFD::restartRun()
{
	//dX.zeros(p.FullSize, 4);
	//dX_cur.zeros(p.FullSize, 4);
	//dX_full.zeros(p.nX, p.nY, 4);
	//dX_cur_full.zeros(p.nX, p.nY, 4);
	//LoadTempBin = 0;
	//Restart_Run_Flag = 0;

	//LoadTempBin = 1;

	//if (dX.load("load" SLASH "fddx") == FALSE)
	//	return FALSE;

	//if (dX_cur.load("load" SLASH "fddx_cur") == FALSE)
	//	return FALSE;

	//if (c.flContinue == 0)
	//	return FALSE;

	//Restart_Run_Flag = 1;

	//for (int i = 0; i < p.nX; i++)
	//{
	//	for (int j = 0; j < p.Channel_Height; j++)
	//	{
	//		int j0 = vlb.Act(i, j);
	//		int indc = i*p.Channel_Height + j;
	//		for (int k = 0; k < 4; k++)
	//		{
	//			dX_full(i, j0, k) = dX(indc, k);
	//			dX_cur_full(i, j0, k) = dX_cur(indc, k);
	//		}
	//	}
	//}

	return TRUE;
}




void clVariablesFD::reset_Nu_out()
{
	Save_loc_Nu = 0;
}


void clVariablesFD::save2file()
{
	//Temp.saveAxbCSR_from_device();
	Temp.savetxt_from_device();
}

void clVariablesFD::saveRestartFiles()
{
	Temp.saveCheckPoint();
//	dX_cur.save_bin_from_device("fddx_cur", p.IOqueue);
}

void clVariablesFD::saveDebug()
{
}

void clVariablesFD::setKernelArgs()
{
	cl_int ind = 0;
	TempUpdateCoeffs.set_argument(ind++, Temp.get_add_Macro());
	TempUpdateCoeffs.set_argument(ind++, vlb.Ux_array.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, vlb.Uy_array.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, Temp.get_add_A());
	TempUpdateCoeffs.set_argument(ind++, Temp.get_add_b());
	//TempUpdateCoeffs.set_argument(ind++, vls.dXArr.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, vlb.dXCoeffs.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, Temp.get_add_IndArr());
	if (vlb.kOmegaSolverFlag)
		TempUpdateCoeffs.set_argument(ind++, Alphat.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, &Alpha_fluid);
	TempUpdateCoeffs.set_argument(ind++, &p.dTtfd);

	//ind = 0;
	//vlb.Update_output_kernel[1].set_argument(ind++, vlb.reduceBoth.get_buf_add());
	//vlb.Update_output_kernel[1].set_argument(ind++, vlb.RoJout.get_buf_add());
	//vlb.Update_output_kernel[1].set_argument(ind++, vlb.Tlb_array.get_buf_add());
	//vlb.Update_output_kernel[1].set_argument(ind++, vlb.J_array.get_buf_add());


}


void clVariablesFD::setSourceDefines()
{
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "TIN", ROE_INX);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "TMIN", T_Actual_Min);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "TDIFF", T_Actual_Diff);

	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "FULLSIZEXM1", p.nX - 1);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_MULTIPLIER", 4.*p.Pipe_radius);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "LOCAL_SOLVE_SIZE_FD", 
		WORKGROUPSIZEX_FD*WORKGROUPSIZEY_FD);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "PR_TURB_NUMBER", PrTurbNum);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "PR_NUMBER", PrNum);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "TFD_X_IN_VAL", TFD_X_IN);
}


// Not implemented yet
void clVariablesFD::solveSSThermal()
{
	if (vfd.calcSSTempFlag == false)
		return;
	vfd.calcSSTempFlag = true;
}


void clVariablesFD::solveTemp()
{
	TempUpdateCoeffs.call_kernel();
	clFinish(FDQUEUE);

	Temp.saveCSR_row_col_val("temp", TRUE);
	Temp.bVec.save_txt_from_device("temp_bvec");
	Temp.xVec->save_txt_from_device_full("temp_xvec");




	Temp.solve();
}

void clVariablesFD::saveParams()
{
	p.setParameter("Restart Run", restartRunFlag);
	p.setParameter("Thermal Solver", thermalSolverFlag);
	p.setParameter("Calculate Nu", calcNuFlag);
	p.setParameter("Save Macros On Start", saveMacroStart);
	p.setParameter("Solve Steady Temp", calcSSTempFlag);

	p.setParameter("Pr Number", PrNum);
	p.setParameter("Pr Number", PrTurbNum);
	p.setParameter("K Soot", kSoot);
	p.setParameter("K air", kAir);
	p.setParameter("Rho Soot", rhoSoot);
	p.setParameter("Cp Soot", cpSoot);

	p.setParameter("Tinlet Actual", T_Actual_Max);
	p.setParameter("Twalls Actual", T_Actual_Min);

	p.setParameter("Tinitial", ROE0);
	p.setParameter("Tinlet Sim", ROE_INX);

	p.setParameter("Temp Max Rel Tol", tempMaxRelTol);
	p.setParameter("Temp Max Abs Tol", tempMaxAbsTol);
	p.setParameter("Temp Max Iterations", tempMaxIters);
}


void clVariablesFD::testRestartRun()
{
	allocateArrays();
	restartRunFlag = p.getParameter<bool>("Restart Run", false);

	if (Temp_array.load("load" SLASH "temp"))
	{
		tempLoadedFlag = true;
	}
	else
	{
		tempLoadedFlag = false;
		restartRunFlag = false;
	}
}

void clVariablesFD::updateAlphaArray()
{
	//for (int i = 0; i < vls.nBL; i++)
	//	bcFindNodes(i);

	//for (int i = 1; i < p.nY - 1; i++)
	//{
	//	if (vls.M(0, i) != LB_SOLID)
	//	{
	//		dX_cur_full(0, i, 1) = 0.5;
	//	}
	//}
}


void clVariablesFD::UpdateNu()
{
	//#ifdef USE_ORIG_NU
	//	vfd.Reinitialize_Nu_kernel();
	//#else
	//	Nu_kernel.set_argument(6, (void*)&Save_loc_Nu);
	//	Nu_kernel.call_kernel();
	//	clFlush(p.FDqueue);
	//	Save_loc_Nu++;
	//#endif
}



//void clVariablesFD::ini_bounds()
//{
//	//dX_full.fill(1.);
//
//	//for (int i = 0; i < vls.nBL; i++)
//	//	bcFindNodes_C0(i);
//
//	//for (int i = 1; i < p.nY - 1; i++)
//	//{
//	//	if (vls.M(0, i) != LB_SOLID)
//	//	{
//	//		dX_full(0, i, 1) = 0.5;
//	//	}
//	//}
//
//	//ini_alpha();
//	//UpdateAlphaArray();
//}

//void clVariablesFD::Find_Neighbors()
//{
//	//Neighs.allocate(p.nX * p.Channel_Height);
//
//	//for (int i = 0; i < p.nX; i++)
//	//{
//	//	for (int j = 0; j < p.Channel_Height; j++)
//	//	{
//	//		int yact = vlb.Act(i, j);
//	//		cl_int2 n_ind = { { i, yact + 1 } };
//	//		cl_int2 s_ind = { { i, yact - 1 } };
//	//		cl_int2 e_ind = { { i + 1, yact } };
//	//		cl_int2 w_ind = { { i - 1, yact } };
//
//
//	//		int n_red = vlb.Stor(i, yact + 1);
//	//		if (n_red != -1)
//	//		{
//	//			n_red += i * (p.Channel_Height + 1);
//	//		}
//
//	//		int s_red = vlb.Stor(i, yact - 1);
//	//		if (s_red != -1)
//	//		{
//	//			s_red += i * (p.Channel_Height + 1);
//	//		}
//
//	//		int e_red;
//	//		if (i == p.nX - 1)
//	//		{
//	//			e_red = (i - 2)*(p.Channel_Height + 1) + vlb.Stor(i - 2, yact);
//	//		}
//	//		else
//	//		{
//	//			e_red = vlb.Stor(i + 1, yact);
//	//			if (e_red != -1)
//	//			{
//	//				e_red += (i + 1) * (p.Channel_Height + 1);
//	//			}
//	//		}
//
//	//		int w_red;
//	//		if (i == 0)
//	//		{
//	//			w_red = -1;
//	//		}
//	//		else
//	//		{
//	//			w_red = vlb.Stor(i - 1, yact);
//	//			if (w_red != -1)
//	//			{
//	//				w_red += (i - 1) * (p.Channel_Height + 1);
//	//			}
//	//		}
//
//	//		Neighs(i*p.Channel_Height + j) = { { e_red, w_red, n_red, s_red } };
//	//	}
//	//}
//}
//
//void clVariablesFD::bcFindDirection(int dir, int bl, cl_double2 vC0,
//	cl_double2 vC1, cl_double2 vCn, cl_int2 vCmin, cl_int2 vCmax, int fflag)
//{
//	double dist;
//	cl_double2 Cdir = vlb.Cxy_double[dir];
//	double Cndir = DOT_PROD(vCn, Cdir);
//	if (Cndir == 0.)
//		return;
//	if (Cndir > 0.)
//	{
//		dir = vlb.rev[dir];
//		Cdir = vlb.Cxy_double[dir];
//	}
//	cl_int2 vCcut0, vCcut1;
//	cl_double2 vCcut;
//
//	if (vlb.Cx[dir] != 0)
//	{
//		int i = vCmin.x;
//		if (vlb.Cx[dir] == -1) i = vCmax.x;
//
//		for (int j = vCmin.y; j <= vCmax.y; j++)
//		{
//			if (vls.bcFindIntersection(vCcut0, vCcut1, vCcut, dist, { { (double)i, (double)j } },
//				Cdir, vC0, vC1, vCn) == TRUE)
//			{
//				if (fflag == 0)
//					bcSetBoundaryNode(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn, 0);
//				else
//					bcSetInterfaceNode(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn, 0);
//			}
//		}
//	}
//
//	if (vlb.Cy[dir] != 0)
//	{
//		int j = vCmin.y;
//		if (vlb.Cy[dir] == -1) j = vCmax.y;
//
//		for (int i = vCmin.x; i <= vCmax.x; i++)
//		{
//			if (vls.bcFindIntersection(vCcut0, vCcut1, vCcut, dist, { { (double)i, (double)j } },
//				Cdir, vC0, vC1, vCn) == TRUE)
//			{
//				if (fflag == 0)
//					bcSetBoundaryNode(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn, 1);
//				else
//					bcSetInterfaceNode(vCcut0, vCcut1, dir, dist, bl, vCcut, vC0, vC1, vCn, 1);
//			}
//		}
//	}
//}
//
//BOOL clVariablesFD::bcFindIntersectionNormal(cl_double2 &vC, double &dist, cl_double2 vL0,
//	cl_double2 vP0, cl_double2 vP1, cl_double2 &vN)
//{
//	cl_double2 vP10 = Subtract2(vP1, vP0);
//	cl_double2 vLd = Subtract2(vL0, vP0);
//	cl_double2 vNm = { { -vP10.y, vP10.x } };
//	vN = VNORM(vNm);
//	dist = DOT_PROD(vLd, vN);
//	cl_double2 vNdist = { { vN.x * dist, vN.y*dist } };
//
//	vC = Subtract2(vL0, vNdist);
//
//	cl_double2 vd0 = Subtract2(vC, vP0), vd1 = Subtract2(vC, vP1);
//
//	if (fabs(GETLEN(vP10) - GETLEN(vd0) - GETLEN(vd1)) < c.eps)
//		return TRUE;
//	return FALSE;
//}
//
//void clVariablesFD::bcFindNodes(int bl)
//{
//	int n0 = vls.BL(bl, 0), n1 = vls.BL(bl, 1);
//	if (n0 < 0)
//		return;
//	cl_double2 vC0 = { { vls.C[n0].x, vls.C[n0].y } }, vC1 = { { vls.C[n1].x, vls.C[n1].y } };
//	cl_double2 vCn = { { vC0.y - vC1.y, vC1.x - vC0.x } };
//	cl_int2 vCmin = vls.min2(vC0, vC1), vCmax = vls.max2(vC0, vC1);
//
//	if (vCmin.x < 0) vCmin.x = 0;
//	if (vCmin.x > p.nX - 1) return;
//	if (vCmax.x < 0) return;
//	if (vCmax.x > p.nX - 1) vCmax.x = p.nX - 1;
//
//	if (vCmin.y < 0) vCmin.y = 0;
//	if (vCmin.y > p.nY - 1) return;
//	if (vCmax.y < 0) return;
//	if (vCmax.y > p.nY - 1) vCmax.y = p.nY - 1;
//
//	for (int dir = 1; dir < vlb.nL; dir += 2)
//		bcFindDirection(dir, bl, vC0, vC1, vCn, vCmin, vCmax, 1);
//}
//
//void clVariablesFD::bcFindNodes_C0(int bl)
//{
//	int n0 = vls.BL(bl, 0), n1 = vls.BL(bl, 1);
//	if (n0 < 0)
//		return;
//	cl_double2 vC0 = { { vls.C0[n0].x, vls.C0[n0].y } }, vC1 = { { vls.C0[n1].x, vls.C0[n1].y } };
//	cl_double2 vCn = { { vC0.y - vC1.y, vC1.x - vC0.x } };
//	cl_int2 vCmin = vls.min2(vC0, vC1), vCmax = vls.max2(vC0, vC1);
//
//	if (vCmin.x < 0) vCmin.x = 0;
//	if (vCmin.x > p.nX - 1) return;
//	if (vCmax.x < 0) return;
//	if (vCmax.x > p.nX - 1) vCmax.x = p.nX - 1;
//
//	if (vCmin.y < 0) vCmin.y = 0;
//	if (vCmin.y > p.nY - 1) return;
//	if (vCmax.y < 0) return;
//	if (vCmax.y > p.nY - 1) vCmax.y = p.nY - 1;
//
//	for (int dir = 1; dir < vlb.nL; dir += 2)
//		bcFindDirection(dir, bl, vC0, vC1, vCn, vCmin, vCmax, 0);
//}
//
////Alpha: C,E,W,N,S        dX:E,W,N,S
//void clVariablesFD::bcSetBoundaryNode(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
//	int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vC2, int xyz)
//{
//	if (ii0.x >= p.nX || ii0.x < 0)
//		return;
//	if (ii0.y >= p.nY || ii0.y < 0)
//		return;
//
//	if (ii1.x >= p.nX || ii1.x < 0)
//		return;
//	if (ii1.y >= p.nY || ii1.y < 0)
//		return;
//
//	if (vls.M_o(ii0.x, ii0.y) == LB_SOLID)
//		return;
//	if (vls.M_o(ii1.x, ii1.y) != LB_SOLID)
//		return;
//
//	dX_full(ii0.x, ii0.y, dir - 1) = dist;
//
//
//	//if (xyz == 0)
//	//{
//	//	if (vlb.Cx[dir] == -1)
//	//	{
//	//		dX_full(ii0.x, ii0.y, 1) = dist;
//	//	}
//	//	else
//	//	{
//	//		dX_full(ii0.x, ii0.y, 0) = dist;
//	//	}
//	//}
//	//else
//	//{
//	//	if (vlb.Cy[dir] == -1)
//	//	{
//	//		dX_full(ii0.x, ii0.y, 3) = dist;
//	//	}
//	//	else
//	//	{
//	//		dX_full(ii0.x, ii0.y, 2) = dist;
//	//	}
//	//}
//}
//
//void clVariablesFD::bcSetInterfaceNode(cl_int2 ii0, cl_int2 ii1, int dir, double dist,
//	int bl, cl_double2 vCc, cl_double2 vC0, cl_double2 vC1, cl_double2 vC2, int xyz)
//{
//	if (ii0.x >= p.nX || ii0.x < 0)
//		return;
//	if (ii0.y >= p.nY || ii0.y < 0)
//		return;
//
//	if (ii1.x >= p.nX || ii1.x < 0)
//		return;
//	if (ii1.y >= p.nY || ii1.y < 0)
//		return;
//
//	if (vls.M(ii0.x, ii0.y) == LB_SOLID)
//		return;
//	if (vls.sf(ii1.x, ii1.y) == TRUE)
//		return;
//
//	dX_cur_full(ii0.x, ii0.y, dir - 1) = dist;
//
//
//	//if (xyz == 0)
//	//{
//	//	if (vlb.Cx[dir] == -1)
//	//	{
//	//		dX_cur_full(ii0.x, ii0.y, 1) = dist;
//	//	}
//	//	else
//	//	{
//	//		dX_cur_full(ii0.x, ii0.y, 0) = dist;
//	//	}
//	//}
//	//else
//	//{
//	//	if (vlb.Cy[dir] == -1)
//	//	{
//	//		dX_cur_full(ii0.x, ii0.y, 3) = dist;
//	//	}
//	//	else
//	//	{
//	//		dX_cur_full(ii0.x, ii0.y, 2) = dist;
//	//	}
//	//}
//}

//void clVariablesFD::Create_Coeffs()
//{
//	//for (int i = 0; i < p.nX; i++)
//	//{
//	//	int jj = 0;
//	//	for (int j = 1; j < p.nY - 1; j++)
//	//	{
//	//		char map = vls.M_o(i, j);
//	//		if (map != LB_SOLID)
//	//		{
//	//			int yred = vlb.Stor(i, j);
//	//			int indc = i * p.Channel_Height + yred;
//	//			for (int k = 0; k < 4; k++)
//	//			{
//	//				dX(indc, k) = dX_full(i, j, k);
//	//				dX_cur(indc, k) = dX_cur_full(i, j, k);
//	//			}
//	//			jj++;
//	//		}
//	//	}
//	//}
//}
//
//void clVariablesFD::Update_Transient_Coeffs()
//{
//	//for (int j = 0; j < p.Channel_Height; j++)
//	//{
//	//	int yred = j;
//	//
//	//	int indc = (p.nX - 1)*p.Channel_Height + yred;
//	//	C(indc, 0) = -2.;
//	//	C(indc, 2) = 1.;
//	//	A(indc, 0) = 1.;
//	//	A(indc, 2) = 0.;
//	//	B(indc, 0) = 1.;
//	//	B(indc, 2) = -1.;
//	//}
//
//	//A.copy_to_buffer(p.IOqueue, CL_TRUE);
//	//B.copy_to_buffer(p.IOqueue, CL_TRUE);
//	//C.copy_to_buffer(p.IOqueue, CL_TRUE);
//
//	//A.FreeHost();
//	//B.FreeHost();
//	//C.FreeHost();
//	//FD_update_base_kernel.call_kernel();
//	//clFinish(p.FDqueue);
//}


//#ifndef USE_ORIG_NU 
//
//void clVariablesFD::initialize_Nu_kernel()
//{
//	int ind = 0;
//	Save_loc_Nu = 0;
//	num_nu_nodes = vtr.numbl_bounds;
//
//	int Fullsize_nu = WORKGROUPSIZE_NU;
//	while (Fullsize_nu < num_nu_nodes)
//		Fullsize_nu += WORKGROUPSIZE_NU;
//
//	Nu.zeros(OUTPUT_MAX_LINES_NU, num_nu_nodes * 3);
//	Nu.allocate_buffer_w_copy(p.context, p.IOqueue);
//
//	Nu_kernel.create_kernel(GetSourceProgram, FDQUEUE_REF, "FD_calc_Nu");
//	Nu_kernel.set_size_1D(Fullsize_nu, WORKGROUPSIZE_NU);
//
//	Nu_kernel.set_argument(ind++, vtr.BL.get_buf_add());
//	Nu_kernel.set_argument(ind++, vtr.BLindicies.get_buf_add());
//	Nu_kernel.set_argument(ind++, Nu.get_buf_add());
//	Nu_kernel.set_argument(ind++, vtr.NodV.get_buf_add());
//	Nu_kernel.set_argument(ind++, vlb.J_array.get_buf_add());
//	Nu_kernel.set_argument(ind++, Temp.get_buf_add());
//	Nu_kernel.set_argument(6, (void*)&Save_loc_Nu);
//	Nu_kernel.set_argument(7, (void*)&num_nu_nodes);
//}
//
//#else
//
//
//void clVariablesFD::Reinitialize_Nu_kernel()
//{
//	vls.C.read_from_buffer(p.IOqueue);
//	vls.M.read_from_buffer(p.IOqueue);
//
//	LBdist_temp.fill(100.);
//
//	int NodeNum = X_STOP_POS - X_MIN_VAL + 1;
//
//	LBnodey.fill({ { -1, -1 } });
//
//	int shiftind = 0;
//
//	for (int i = 0; i < vls.nBL; i++)
//	{
//		int bl_index = i;
//		if (i >= vls.nBL / 2)
//		{
//			shiftind = NodeNum;
//			bl_index = (vls.nBL - 1 - i) + vls.nBL / 2;
//		}
//
//
//
//		cl_double2 C0t = vls.C[vls.BL(i, 0)], C1t = vls.C[vls.BL(i, 1)];
//		cl_double2 C0 = { { C0t.x, C0t.y } }, C1 = { { C1t.x, C1t.y } };
//		cl_int2 C0i = vls.min2(C0, C1), C1i = vls.max2(C0, C1);
//		C0i = { { C0i.x - 1, C0i.y - 1 } };
//		C1i = { { C1i.x + 1, C1i.y + 1 } };
//		if (C0i.x < X_MIN_VAL) C0i.x = X_MIN_VAL;
//		if (C0i.y < 0) C0i.y = 0;
//		if (C1i.x > X_STOP_POS) C1i.x = X_STOP_POS;
//		if (C1i.y > p.nY - 1) C1i.y = p.nY - 1;
//		if (C0i.x > X_STOP_POS)
//			continue;
//		if (C1i.x < X_MIN_VAL)
//			continue;
//
//		double dist;
//		cl_double2 vCcut, vN;
//
//		for (int ii = C0i.x; ii <= C1i.x; ii++)
//		{
//			for (int jj = C0i.y; jj <= C1i.y; jj++)
//			{
//				char map = vls.M(ii, jj);
//				if (map == LB_SOLID)
//					continue;
//
//				cl_double2 L0 = { { (double)ii, (double)jj } };
//				if (bcFindIntersectionNormal(vCcut, dist, L0, C0, C1, vN) == TRUE)
//				{
//					if (dist < LBdist_temp(ii + shiftind))
//					{
//						LBdist_temp(ii + shiftind) = dist;
//
//						int act_ind = ii - X_MIN_VAL;
//
//						if (shiftind == 0)
//						{
//							LBdist_bot(act_ind) = { { (double)ii - vCcut.x, (double)jj - vCcut.y } };
//							LBnodey(act_ind).x = jj;
//							BLinds(act_ind).x = bl_index;
//						}
//						else
//						{
//							LBdist_top(act_ind) = { { (double)ii - vCcut.x, (double)jj - vCcut.y } };
//							LBnodey(act_ind).y = jj;
//							BLinds(act_ind).y = bl_index;
//
//						}
//					}
//				}
//			}
//		}
//	}
//
//	LBnodey.copy_to_buffer(p.FDqueue);
//	LBdist_bot.copy_to_buffer(p.FDqueue);
//	LBdist_top.copy_to_buffer(p.FDqueue);
//	BLinds.copy_to_buffer(p.FDqueue);
//	clFinish(p.FDqueue);
//	Nu_kernel.set_argument(6, (void*)&Save_loc_Nu);
//	Nu_kernel.call_kernel();
//	clFlush(p.FDqueue);
//	Save_loc_Nu++;
//}
//
//
void clVariablesFD::initialize_Nu_kernel()
{
	//Array1Dv2i LBnode_temp;
	//Array1Dv2d vCcut_array_temp;
	//Array1Di Blinds_temp;
	//LBnode_temp.allocate(2 * p.nX);
	//LBdist_temp.allocate(2 * p.nX);
	//LBdist_temp.fill(100.);
	//vCcut_array_temp.allocate(2 * p.nX);
	//Blinds_temp.allocate(2 * p.nX);

	//int shiftind = 0;

	//for (int i = 0; i < vls.nBL; i++)
	//{
	//	int bl_index = i;
	//	if (i >= vls.nBL / 2)
	//	{
	//		shiftind = p.nX;
	//		bl_index = (vls.nBL - 1 - i) + vls.nBL / 2;
	//	}



	//	cl_double2 C0t = vls.C[vls.BL(i, 0)], C1t = vls.C[vls.BL(i, 1)];
	//	cl_double2 C0 = { { C0t.x, C0t.y } }, C1 = { { C1t.x, C1t.y } };
	//	cl_int2 C0i = vls.min2(C0, C1), C1i = vls.max2(C0, C1);
	//	C0i = { { C0i.x - 1, C0i.y - 1 } };
	//	C1i = { { C1i.x + 1, C1i.y + 1 } };
	//	if (C0i.x < 0) C0i.x = 0;
	//	if (C0i.y < 0) C0i.y = 0;
	//	if (C1i.x > p.nX - 1) C1i.x = p.nX - 1;
	//	if (C1i.y > p.nY - 1) C1i.y = p.nY - 1;
	//	if (C0i.x > p.nX - 1)
	//		continue;
	//	if (C1i.x < 0)
	//		continue;

	//	double dist;
	//	cl_double2 vCcut, vN;

	//	for (int ii = C0i.x; ii <= C1i.x; ii++)
	//	{
	//		for (int jj = C0i.y; jj <= C1i.y; jj++)
	//		{
	//			char map = vls.M_o(ii, jj);
	//			if (map == LB_SOLID)
	//				continue;

	//			cl_double2 L0 = { { (double)ii, (double)jj } };
	//			if (bcFindIntersectionNormal(vCcut, dist, L0, C0, C1, vN) == TRUE)
	//			{
	//				if (dist < LBdist_temp(ii + shiftind))
	//				{
	//					Blinds_temp(ii + shiftind) = bl_index;
	//					LBdist_temp(ii + shiftind) = dist;
	//					LBnode_temp(ii + shiftind) = { { ii, jj } };
	//					vCcut_array_temp(ii + shiftind) = vCcut;
	//				}
	//			}
	//		}
	//	}
	//}


	//int NodeNum = X_STOP_POS - X_MIN_VAL + 1;

	//LBnodey.allocate(NodeNum);
	//LBnodey.fill({ { -1, -1 } });
	//LBdist_bot.allocate(NodeNum);
	//LBdist_top.allocate(NodeNum);
	//BLinds.allocate(NodeNum);


	//for (int i = X_MIN_VAL; i <= X_STOP_POS; i++)
	//{
	//	int act_ind = i - X_MIN_VAL;

	//	int ypos1 = LBnode_temp(i).y;
	//	int ypos2 = LBnode_temp(i + p.nX).y;

	//	cl_double2 Node_Pos1 = { { (double)LBnode_temp(i).x, (double)ypos1 } };
	//	cl_double2 Node_Pos2 = { { (double)LBnode_temp(i).x, (double)ypos2 } };

	//	LBdist_bot(act_ind) = Subtract2(Node_Pos1, vCcut_array_temp(i));
	//	LBdist_top(act_ind) = Subtract2(Node_Pos2, vCcut_array_temp(i + p.nX));

	//	int jj1 = vlb.Stor(LBnode_temp(i));
	//	int jj2 = vlb.Stor(LBnode_temp(i + p.nX));

	//	LBnodey(act_ind).x = (LBdist_temp(i) < 100) ? (LBnode_temp(i).y) : (-1);
	//	LBnodey(act_ind).y = (LBdist_temp(i + p.nX) < 100) ? (LBnode_temp(i + p.nX).y) : (-1);

	//	BLinds(act_ind) = { { Blinds_temp(i), Blinds_temp(i + p.nX) } };
	//}

	//Nu.zeros(OUTPUT_MAX_LINES_NU, NodeNum * 5);
	//Nu.allocate_buffer_w_copy(p.context, p.IOqueue);

	//LBnodey.allocate_buffer_w_copy(p.context, p.IOqueue);
	//LBdist_bot.allocate_buffer_w_copy(p.context, p.IOqueue);
	//LBdist_top.allocate_buffer_w_copy(p.context, p.IOqueue);
	//BLinds.allocate_buffer_w_copy(p.context, p.IOqueue);

	//Save_loc_Nu = 0;
	//num_nu_nodes = NodeNum;

	//int Fullsize_nu = WORKGROUPSIZE_NU;
	//while (Fullsize_nu < num_nu_nodes)
	//	Fullsize_nu += WORKGROUPSIZE_NU;

	//Nu_kernel.create_kernel(GetSourceProgram, FDQUEUE_REF, "FD_calc_Nu_Orig");
	//Nu_kernel.set_size_1D(Fullsize_nu, WORKGROUPSIZE_NU);
	//int NodeX_offset = X_MIN_VAL;
	//int ind = 0;
	//Nu_kernel.set_argument(ind++, vlb.Tlb_array.get_buf_add());
	//Nu_kernel.set_argument(ind++, Nu.get_buf_add());
	//Nu_kernel.set_argument(ind++, vlb.J_array.get_buf_add());
	//Nu_kernel.set_argument(ind++, LBdist_bot.get_buf_add());
	//Nu_kernel.set_argument(ind++, LBdist_top.get_buf_add());
	//Nu_kernel.set_argument(ind++, vtr.BL.get_buf_add());
	//Nu_kernel.set_argument(ind++, (void*)&Save_loc_Nu);
	//Nu_kernel.set_argument(ind++, (void*)&num_nu_nodes);
	//Nu_kernel.set_argument(ind++, BLinds.get_buf_add());
	//Nu_kernel.set_argument(ind++, (void*)&NodeX_offset);
	//Nu_kernel.set_argument(ind++, LBnodey.get_buf_add());
	//Nu_kernel.set_argument(ind++, vfd.dX_cur.get_buf_add());
	//Nu_kernel.set_argument(ind++, vfd.dX.get_buf_add());
	//Nu_kernel.set_argument(ind++, vlb.Stor.get_buf_add());


	//LBdist_temp.zeros(2 * NodeNum);
}
//
//#endif


