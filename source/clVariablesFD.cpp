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
	tempArray.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	alphaArray.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 5, 5);
	
	// Nu intialized in separate function (need to calculate
	// number of nodes before initializing)
	//dX_cur.allocate_buffer_w_copy(p.context, p.IOqueue);
	//dX.allocate_buffer_w_copy(p.context, p.IOqueue, CL_MEM_READ_ONLY);
	//Neighs.allocate_buffer_w_copy(p.context, p.IOqueue);
}

void clVariablesFD::allocateBuffers()
{
	Alphat.allocate_buffer_w_copy();
	alphaArray.allocate_buffer_w_copy();
	//tempArray buffer allocated by BiCGStab class
}


double clVariablesFD::calcAverageT()
{
	return sumTemp.reduceSingle();
}


void clVariablesFD::createKernels()
{
	SetSSCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF, "steady_State_T_Coeffs");
	SetSSCoeffs.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY, WORKGROUPSIZEY_LB);

	if (vlb.kOmegaClass.kOmegaSolverFlag)
	{
		TempUpdateCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF, "Update_T_Coeffs_Implicit_Turbulent");
		TempUpdateCoeffs.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY, WORKGROUPSIZEY_LB);
	}
	else
	{
		TempUpdateCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF,
			"Update_T_Coeffs_Implicit");
		TempUpdateCoeffs.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY,
			WORKGROUPSIZEY_LB);
	}

	// TODO: Correct kernel names, and sizes
	if (calcNuFlag)
	{
		int fullSizeNu = getGlobalSizeMacro(numNuNodes, WORKGROUPSIZE_NU);

		calcNu.create_kernel(GetSourceProgram, FDQUEUE_REF, "FD_calc_Nu");
		calcNu.set_size(fullSizeNu, WORKGROUPSIZE_NU);
		updateNuVars.create_kernel(GetSourceProgram, FDQUEUE_REF, "updateNuVars");
		//updateNuVars.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY, WORKGROUPSIZEY_LB);

	}

}


void clVariablesFD::freeHostArrays()
{
	Alphat.FreeHost();
	alphaArray.FreeHost();
}


void clVariablesFD::ini()
{
	// Set source defines even if temp solver not being used, because
	// opencl source may not compile without variables defined
	setSourceDefines();

	if (thermalSolverFlag == false)
		return;

	// Add source file to string containing source to be compiled
	sourceGenerator::SourceInstance()->addFile2Kernel("tfdKernels.cl");
	
	TempInds.ini(p.nX, p.XsizeFull, p.nY, p.nY, &vls.M_o);

	// allocate device buffers
	allocateBuffers();

	// Create BiCGStab Solver for temperature
	Temp.CreateSolver(&tempArray, &TempInds, FDQUEUE_REF,
		tempMaxIters, tempMaxRelTol, tempMaxAbsTol);

	// Initialize reduction class for temp
	sumTemp.ini(*Temp.getMacroArray(), "redTemp");
	
	// Create function pointers to functions to be called after compiling
	// opencl source code
	std::function<void(void)> createKerPtr = 
		std::bind(&clVariablesFD::createKernels, this);
	std::function<void(void)> setArgsPtr = 
		std::bind(&clVariablesFD::setKernelArgs, this);

	// Pass function pointers to sourceGenerator
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



	if (saveMacroStart)
		save2file();
}

void clVariablesFD::iniAlpha()
{

}




#ifndef USE_ORIG_NU 

void clVariablesFD::iniNuKernel()
{
	int ind = 0;
	Save_loc_Nu = 0;
	numNuNodes = vtr.;

	Nu.zeros(OUTPUT_MAX_LINES_NU, numNuNodes * 3);
	Nu.allocate_buffer_w_copy();
}
#else
// Refer to old version of code for this implementation

#endif

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

bool clVariablesFD::testRestartRun()
{
	allocateArrays();
	restartRunFlag = p.getParameter("Restart Run", false);

	// Temp can still be read from file and if it is not
	// a restarted run
	if (tempArray.load("load" SLASH "temp"))
	{
		tempLoadedFlag = true;
	}
	else
	{
		restartRunFlag = false;
		tempLoadedFlag = false;
		tempArray.fill(ROE0);
	}

	return restartRunFlag;
}

void clVariablesFD::save2file()
{
	Temp.savetxt_from_device();
}

void clVariablesFD::saveRestartFiles()
{
	Temp.saveCheckPoint();
}

void clVariablesFD::saveDebug()
{
	alphaArray.save_txt_from_device_as_multi2D();
	Alphat.save_txt_from_device();
}

// Saves time output data (i.e. avg velocity, shear, Nu, etc)
void clVariablesFD::saveTimeData()
{
	// Append Nu to file
	Save_loc_Nu = 0;
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
	if (vlb.kOmegaClass.kOmegaSolverFlag)
		TempUpdateCoeffs.set_argument(ind++, Alphat.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, &Alpha_fluid);
	TempUpdateCoeffs.set_argument(ind++, &p.dTtfd);


	ind = 0;
	SetSSCoeffs.set_argument(ind++, Temp.get_add_Macro());
	SetSSCoeffs.set_argument(ind++, vlb.Ux_array.get_buf_add());
	SetSSCoeffs.set_argument(ind++, vlb.Uy_array.get_buf_add());
	SetSSCoeffs.set_argument(ind++, Temp.get_add_A());
	SetSSCoeffs.set_argument(ind++, Temp.get_add_b());
	SetSSCoeffs.set_argument(ind++, vlb.dXCoeffs.get_buf_add());
	SetSSCoeffs.set_argument(ind++, Temp.get_add_IndArr());
	SetSSCoeffs.set_argument(ind++, Alphat.get_buf_add());
	SetSSCoeffs.set_argument(ind++, &Alpha_fluid);
	int komegaFlag_ = (vlb.kOmegaClass.kOmegaSolverFlag) ? 1 : 0;
	SetSSCoeffs.set_argument(ind++, &komegaFlag_);
	
	ind = 0;
	calcNu.set_argument(ind++, vtr.BL.get_buf_add());
	calcNu.set_argument(ind++, vtr.BLindicies.get_buf_add());
	calcNu.set_argument(ind++, Nu.get_buf_add());
	calcNu.set_argument(ind++, vtr.NodV.get_buf_add());
	calcNu.set_argument(ind++, vlb.J_array.get_buf_add());
	calcNu.set_argument(ind++, Temp.get_buf_add());
	calcNu.set_argument(6, (void*)& Save_loc_Nu);
	calcNu.set_argument(7, (void*)& num_nu_nodes);



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


void clVariablesFD::solveSSThermal()
{
	if (vfd.calcSSTempFlag == false)
		return;

	SetSSCoeffs.call_kernel();
	
	Temp.setMaxIters(1000);
	Temp.setAbsTol(1.e-12);
	Temp.setRelTol(1.e-5);
	Temp.solve();

	Temp.setMaxIters(tempMaxIters);
	Temp.setAbsTol(tempMaxAbsTol);
	Temp.setRelTol(tempMaxRelTol);
}


void clVariablesFD::Solve()
{
	TempUpdateCoeffs.call_kernel();
	clFinish(FDQUEUE);
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
//				Cdir, vC0, vC1, vCn) == true)
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
//				Cdir, vC0, vC1, vCn) == true)
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
//bool clVariablesFD::bcFindIntersectionNormal(cl_double2 &vC, double &dist, cl_double2 vL0,
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
//		return true;
//	return false;
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
//	if (vls.sf(ii1.x, ii1.y) == true)
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





#ifdef USE_ORIG_NU 
void clVariablesFD::reIniNuKernel()
{
	vls.C.read_from_buffer(p.IOqueue);
	vls.M.read_from_buffer(p.IOqueue);

	LBdist_temp.fill(100.);

	int NodeNum = X_STOP_POS - X_MIN_VAL + 1;

	LBnodey.fill({ { -1, -1 } });

	int shiftind = 0;

	for (int i = 0; i < vls.nBL; i++)
	{
		int bl_index = i;
		if (i >= vls.nBL / 2)
		{
			shiftind = NodeNum;
			bl_index = (vls.nBL - 1 - i) + vls.nBL / 2;
		}



		cl_double2 C0t = vls.C[vls.BL(i, 0)], C1t = vls.C[vls.BL(i, 1)];
		cl_double2 C0 = { { C0t.x, C0t.y } }, C1 = { { C1t.x, C1t.y } };
		cl_int2 C0i = vls.min2(C0, C1), C1i = vls.max2(C0, C1);
		C0i = { { C0i.x - 1, C0i.y - 1 } };
		C1i = { { C1i.x + 1, C1i.y + 1 } };
		if (C0i.x < X_MIN_VAL) C0i.x = X_MIN_VAL;
		if (C0i.y < 0) C0i.y = 0;
		if (C1i.x > X_STOP_POS) C1i.x = X_STOP_POS;
		if (C1i.y > p.nY - 1) C1i.y = p.nY - 1;
		if (C0i.x > X_STOP_POS)
			continue;
		if (C1i.x < X_MIN_VAL)
			continue;

		double dist;
		cl_double2 vCcut, vN;

		for (int ii = C0i.x; ii <= C1i.x; ii++)
		{
			for (int jj = C0i.y; jj <= C1i.y; jj++)
			{
				char map = vls.M(ii, jj);
				if (map == LB_SOLID)
					continue;

				cl_double2 L0 = { { (double)ii, (double)jj } };
				if (bcFindIntersectionNormal(vCcut, dist, L0, C0, C1, vN) == true)
				{
					if (dist < LBdist_temp(ii + shiftind))
					{
						LBdist_temp(ii + shiftind) = dist;

						int act_ind = ii - X_MIN_VAL;

						if (shiftind == 0)
						{
							LBdist_bot(act_ind) = { { (double)ii - vCcut.x, (double)jj - vCcut.y } };
							LBnodey(act_ind).x = jj;
							BLinds(act_ind).x = bl_index;
						}
						else
						{
							LBdist_top(act_ind) = { { (double)ii - vCcut.x, (double)jj - vCcut.y } };
							LBnodey(act_ind).y = jj;
							BLinds(act_ind).y = bl_index;

						}
					}
				}
			}
		}
	}

	LBnodey.copy_to_buffer(p.FDqueue);
	LBdist_bot.copy_to_buffer(p.FDqueue);
	LBdist_top.copy_to_buffer(p.FDqueue);
	BLinds.copy_to_buffer(p.FDqueue);
	clFinish(p.FDqueue);
	Nu_kernel.set_argument(6, (void*)&Save_loc_Nu);
	Nu_kernel.call_kernel();
	clFlush(p.FDqueue);
	Save_loc_Nu++;
}
#endif

void clVariablesFD::updateTimeData()
{
	calcNu.setOptionCallKernel(&Save_loc_Nu);
	clFlush(*calcNu.getDefaultQueue());
	Save_loc_Nu++;
	if (Save_loc_Nu == OUTPUT_MAX_LINES_NU)
		saveTimeData();
}



