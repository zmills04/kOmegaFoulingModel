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
	tempArray.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	alphaDer.allocate(p.nX, p.XsizeFull, p.nY, p.nY);
	alphaDer.fill({ {0.,0.} });
	if (chtCorrectionFlag)
	{
		rhoCpDer.allocate(p.nX, p.XsizeFull, p.nY, p.nY);
		rhoCpDer.fill({ {0.,0.} });
	}
	// Nu intialized in separate function (need to calculate
	// number of nodes before initializing)
	//dX_cur.allocate_buffer_w_copy(p.context, p.IOqueue);
	//dX.allocate_buffer_w_copy(p.context, p.IOqueue, CL_MEM_READ_ONLY);
	//Neighs.allocate_buffer_w_copy(p.context, p.IOqueue);
}

void clVariablesFD::allocateBuffers()
{
	if (chtCorrectionFlag)
		rhoCpDer.allocate_buffer_w_copy();
	alphaDer.allocate_buffer_w_copy();
	//tempArray buffer allocated by BiCGStab class
}

double clVariablesFD::calcAlphaDirection(const int i1, const int j1, const int i2, const int j2, const int dirInd)
{
	double k_c, ro_c, cp_c, k_n, ro_n, cp_n;

	if ((vls.nType(i1, j1) & M_FLUID_NODE))
	{
		k_c = kAirLB;
		cp_c = cpAirLB;
		ro_c = rhoAirLB;
	}
	else
	{
		k_c = kSootLB;
		cp_c = cpSootLB;
		ro_c = rhoSootLB;
	}

	if (vls.nType(i2, j2) & M_FLUID_NODE)
	{
		k_n = kAirLB;
		cp_n = cpAirLB;
		ro_n = rhoAirLB;
	}
	else
	{
		k_n = kSootLB;
		cp_n = cpSootLB;
		ro_n = rhoSootLB;
	}

	double dXcur = vls.dXArr(i1, j1, dirInd), dX0 = vls.dXArr0(i1, j1, dirInd);
	if (dXcur - dX0 < p.eps)
		return k_c / ro_c / cp_c;

	double del = dXcur / dX0;

	double coeffA = k_n * ro_c * cp_c;
	double coeffB = k_n * ro_c * cp_n + k_n * ro_n * cp_c + k_c * ro_c * cp_c;
	double coeffC = k_n * ro_n * cp_n + k_c * ro_c * cp_n + k_c * ro_n * cp_c;
	double coeffD = k_c * ro_n * cp_n;

	return k_c * k_n / (coeffA * del * del * del + coeffB * del * del * (1.0 - del) +
		coeffC * del * (1.0 - del) * (1.0 - del) + coeffD * (1.0 - del) * (1.0 - del) * (1.0 - del));
}

cl_double2 clVariablesFD::calcAlphaRhoCpDirection(const int i1, const int j1, const int i2, const int j2, const int dirInd)
{
	double k_c, ro_c, cp_c, k_n, ro_n, cp_n;

	if ((vls.nType(i1, j1) & M_FLUID_NODE))
	{
		k_c = kAirLB;
		cp_c = cpAirLB;
		ro_c = rhoAirLB;
	}
	else
	{
		k_c = kSootLB;
		cp_c = cpSootLB;
		ro_c = rhoSootLB;
	}

	if (vls.nType(i2, j2) & M_FLUID_NODE)
	{
		k_n = kAirLB;
		cp_n = cpAirLB;
		ro_n = rhoAirLB;
	}
	else
	{
		k_n = kSootLB;
		cp_n = cpSootLB;
		ro_n = rhoSootLB;
	}


	double dXcur = vls.dXArr(i1, j1, dirInd), dX0 = vls.dXArr0(i1, j1, dirInd);
	if (dXcur - dX0 < p.eps)
		return { {k_c / ro_c / cp_c, ro_c * cp_c} };


	double del = dXcur / dX0;
	
	double coeffA = k_n * ro_c * cp_c;
	double coeffB = k_n * ro_c * cp_n + k_n * ro_n * cp_c + k_c * ro_c * cp_c;
	double coeffC = k_n * ro_n * cp_n + k_c * ro_c * cp_n + k_c * ro_n * cp_c;
	double coeffD = k_c * ro_n * cp_n;

	cl_double2 retVal;

	// calculate alpha based on harmonic mean interpolation
	retVal.x = k_c * k_n / (coeffA * del * del * del + coeffB * del * del * (1.0 - del) +
		coeffC * del * (1.0 - del) * (1.0 - del) + coeffD * (1.0 - del) * (1.0 - del) * (1.0 - del));

	// calculate alpha based on linear interp of rho*cp
	retVal.y = ro_c * cp_c * (1. - del) + ro_n * cp_n * del;
	return retVal;
}


double clVariablesFD::calcAverageT()
{
	return sumTemp.reduceSingle();
}


void clVariablesFD::createKernels()
{
	SetSSCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF, "Steady_State_T_Coeffs");
	SetSSCoeffs.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY, WORKGROUPSIZEY_LB);

	TempUpdateCoeffs.create_kernel(GetSourceProgram, FDQUEUE_REF,
		"Update_T_Coeffs_Implicit");
	TempUpdateCoeffs.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY,
			WORKGROUPSIZEY_LB);

	updateDAlphaKernel.create_kernel(GetSourceProgram, FDQUEUE_REF,
		"updateDerivativeArrays");
	updateDAlphaKernel.set_size(p.XsizeFull, WORKGROUPSIZEX_LB, p.nY,
		WORKGROUPSIZEY_LB);

	// Make sure that this can be located here
	calcNuKernel.create_kernel(GetSourceProgram, FDQUEUE_REF, "FD_calc_Nu");
	calcNuKernel.set_size(nuNumNodes, WORKGROUPSIZE_NU);

	int numBLPerSide = vtr.Bounds.MAX_BL_BOT - vtr.Bounds.MIN_BL_BOT;
	updateNuKernel[0].create_kernel(GetSourceProgram, FDQUEUE_REF, "updateNuCoeff1");
	updateNuKernel[0].set_size(numBLPerSide, WORKGROUPSIZE_NU);

	updateNuKernel[1].create_kernel(GetSourceProgram, FDQUEUE_REF, "updateNuCoeff2");
	updateNuKernel[1].set_size(2 * nuNumNodes, WORKGROUPSIZE_NU);

	calcNuMeanKernel.create_kernel(GetSourceProgram, FDQUEUE_REF, "getNuMean");
	wgSizeNuMean = getGlobalSizeMacro(p.nY/2, 32);
	calcNuMeanKernel.set_size(nuMeanNumNodes, 1, wgSizeNuMean, wgSizeNuMean);
}

bool clVariablesFD::bcFindIntersectionNormal(cl_double2& vC, double& dist, cl_double2 vL0,
	cl_double2 vP0, cl_double2 vP1, cl_double2& vN)
{
	cl_double2 vP10 = Subtract2(vP1, vP0);
	cl_double2 vLd = Subtract2(vL0, vP0);
	cl_double2 vNm = { { -vP10.y, vP10.x } };
	vN = VNORM(vNm);
	dist = DOT_PROD(vLd, vN);
	cl_double2 vNdist = { { vN.x * dist, vN.y * dist } };

	vC = Subtract2(vL0, vNdist);

	cl_double2 vd0 = Subtract2(vC, vP0), vd1 = Subtract2(vC, vP1);

	if (fabs(GETLEN(vP10) - GETLEN(vd0) - GETLEN(vd1)) < p.eps)
		return true;
	return false;
}


void clVariablesFD::freeHostArrays()
{
	rhoCpDer.FreeHost();
	alphaDer.FreeHost();
}


void clVariablesFD::ini()
{
	// Set source defines even if temp solver not being used, because
	// opencl source may not compile without variables defined
	setSourceDefines();

	if (!tempLoadedFlag)
	{
		tempArray.fillByNodeType(ROE0, vls.nType, M0_FLUID_NODE);
	}

	if (thermalSolverFlag == false)
		return;

	// Add source file to string containing source to be compiled
	sourceGenerator::SourceInstance()->addFile2Kernel("tfdKernels.cl");
	
	// initialize CSR storage arrays
	TempInds.ini(p.nX, p.XsizeFull, p.nY, p.nY, &vls.nType);

	// initialize alpha and rhoCp derivative arrays
	iniDerivativeArrays();

	// Create BiCGStab Solver for temperature
	Temp.CreateSolver(&tempArray, &TempInds, FDQUEUE_REF,
		tempMaxIters, tempMaxRelTol, tempMaxAbsTol);

	// Initialize reduction class for temp
	sumTemp.ini(*Temp.getMacroArray(), restartRunFlag, "redTemp");
	
	// Create function pointers to functions to be called after compiling
	// opencl source code
	std::function<void(void)> createKerPtr = 
		std::bind(&clVariablesFD::createKernels, this);
	std::function<void(void)> setArgsPtr = 
		std::bind(&clVariablesFD::setKernelArgs, this);

	// Pass function pointers to sourceGenerator
	sourceGenerator::SourceInstance()->addIniFunction(createKerPtr, setArgsPtr);

	// allocate device buffers
	allocateBuffers();

	if (saveMacroStart && !restartRunFlag)
		save2file();

	LOGMESSAGE("vfd initialized");
}

void clVariablesFD::iniDerivativeArrays()
{
	for (int i = 0; i < p.nX; i++)
	{
		int ie = MODFAST(i + 1, p.nX);
		int iw = MODFAST(i - 1, p.nX);

		// no fluid nodes should be found on y = 0 or nY-1
		for (int j = 1; j < p.nY - 1; j++)
		{
			if (vls.nType(i, j) & M0_SOLID_NODE)
				continue;

			int jn = j + 1;
			int js = j - 1;

			double alpha_e, alpha_w, alpha_n, alpha_s;

			// Based on distances between actual nodes, so using dXArr0
			// instead of dXArr
			double dxe = vls.dXArr0(i, j, 0), dxw = vls.dXArr0(i, j, 1);
			double dyn = vls.dXArr0(i, j, 2), dys = vls.dXArr0(i, j, 3);

			// Need to get ro*Cp values between node and its neigh along with 
			// Nu when using turbulence model
			if (chtCorrectionFlag)
			{

				cl_double2 alphaRoCp_e = calcAlphaRhoCpDirection(i, j, ie, j, 0);
				cl_double2 alphaRoCp_w = calcAlphaRhoCpDirection(i, j, iw, j, 1);
				cl_double2 alphaRoCp_n = calcAlphaRhoCpDirection(i, j, i, jn, 2);
				cl_double2 alphaRoCp_s = calcAlphaRhoCpDirection(i, j, i, js, 3);

				double roCp_c = 0.25 * (alphaRoCp_e.y + alphaRoCp_w.y + alphaRoCp_n.y + alphaRoCp_s.y);

				double dRoCpdX = 2. * dxw * alphaRoCp_e.y / (dxe * (dxe + dxw)) - 2. * dxe * alphaRoCp_w.y / (dxw * (dxe + dxw))
					+ 2. * (dxe - dxw) * roCp_c / (dxe * dxw);

				double dRoCpdY = 2. * dys * alphaRoCp_n.y / (dyn * (dyn + dys)) - 2. * dyn * alphaRoCp_s.y / (dys * (dyn + dys))
					+ 2. * (dyn - dys) * roCp_c / (dyn * dys);

				dRoCpdX /= roCp_c;
				dRoCpdY /= roCp_c;

				rhoCpDer(i, j) = { {dRoCpdX, dRoCpdY} };

				alpha_e = alphaRoCp_e.x;
				alpha_w = alphaRoCp_w.x;
				alpha_n = alphaRoCp_n.x;
				alpha_s = alphaRoCp_s.x;

			}
			else
			{
				alpha_e = calcAlphaDirection(i, j, ie, j, 0);
				alpha_w = calcAlphaDirection(i, j, iw, j, 1);
				alpha_n = calcAlphaDirection(i, j, i, jn, 2);
				alpha_s = calcAlphaDirection(i, j, i, js, 3);
			}

			double alpha_c = 0.25 * (alpha_e + alpha_w + alpha_n + alpha_s);

			double dAlphadX = 2. * dxw * alpha_e / (dxe * (dxe + dxw)) - 2. * dxe * alpha_w / (dxw * (dxe + dxw))
				+ 2. * (dxe - dxw) * alpha_c / (dxe * dxw);

			double dAlphadY = 2. * dys * alpha_n / (dyn * (dyn + dys)) - 2. * dyn * alpha_s / (dys * (dyn + dys))
				+ 2. * (dyn - dys) * alpha_c / (dyn * dys);

			alphaDer(i, j) = { {dAlphadX, dAlphadY} };
		}
	}
}


// Utilizes variables from vtr, so must be called after 
// initializing vtr, and will throw error if simulation 
// does not have flag set to use vtr.
void clVariablesFD::iniNuCoeffs()
{
	// Need to make sure info about wavy periods is provided, and throw error if not
	ERROR_CHECKING((p.xWavyStart == -1 || p.wavyPeriodLen == -1 || p.numWavyPeriods == -1),
		"Must provide period length and starting location of wavy section "\
		"when calculating Nusselt number.", ERROR_INITIALIZING_VFD);

	NuMean.iniFile(restartRunFlag);
	nuMeanNumNodes = 2*(p.numWavyPeriods + 1);
	nuMeanNumNodesFull = getGlobalSizeMacro(nuMeanNumNodes, WORKGROUPSIZE_NU);
	NuMean.zeros(nuMeanNumNodes, nuMeanNumNodesFull, maxOutputLinesNu, maxOutputLinesNu);
	NuMean.createTimeArray();

	nuNumNodes = vtr.trDomainSize.x;
	
	nuDist.zeros(nuNumNodes * 2);
	nuDist.fill(100.);
	
	nuDistVec.allocate(nuNumNodes * 2);

	nuYInds.allocate(nuNumNodes * 2);
	nuYInds.fill(-1);

	numNuNodesFull = getGlobalSizeMacro(nuNumNodes*5, WORKGROUPSIZE_NU);
	Nu.iniFile(restartRunFlag); 
	Nu.zeros(nuNumNodes*5, numNuNodesFull, maxOutputLinesNu, maxOutputLinesNu);
	Nu.createTimeArray();

	int numBLPerSide = vtr.Bounds.MAX_BL_BOT - vtr.Bounds.MIN_BL_BOT;

	int offsetInd = vtr.trDomainXBounds.x;
	for (int i = vtr.Bounds.MIN_BL_BOT; i < vtr.Bounds.MAX_BL_BOT; i++)
	{
		nuFindNearestNode(i, offsetInd);
	}

	offsetInd = vtr.trDomainXBounds.x - nuNumNodes;
	for (int i = vtr.Bounds.MIN_BL_TOP; i < vtr.Bounds.MAX_BL_TOP; i++)
	{
		nuFindNearestNode(i, offsetInd);
	}

	for (int i = 0; i < 2*nuNumNodes; i++)
	{
		int jj = nuYInds(i);
		int nodeType = 0x0;
		if (jj != -1)
		{
			// Note: nuDistVec corresponds to direction from BL to node,
			// so nuDistVec.x > 0 indicates that BL is to the left of node
			// and nuDistVec.y > 0 indicates that BL is below node.

			// nodeType defines which values of dXcur should be used when
			// calculating the derivative.
			if (nuDistVec(i).x < 0.)
			{
				nodeType |= 0x1;
			}
			if (nuDistVec(i).y < 0.)
			{
				nodeType |= 0x2;
			}
			
		}
		nuYInds(i) = (nuYInds(i) << 2) | nodeType;
	}

	
	Nu.allocate_buffer_w_copy();
	NuMean.allocate_buffer_w_copy();
	nuYInds.allocate_buffer_w_copy();
	nuDistVec.allocate_buffer_w_copy();
	nuDist.allocate_buffer_w_copy();
	

	// Needs to be set here, because setSourceDefines is called before
	// these variables are initialized
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_NUM_NODES", nuNumNodes);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_NUM_NODES_FULLSIZE", numNuNodesFull);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_BL_PER_SIDE_NU", numBLPerSide);
}


void clVariablesFD::loadParams()
{

	thermalSolverFlag = p.getParameter<bool>("Thermal Solver", USE_THERMAL_SOLVER);
	calcNuFlag = p.getParameter<bool>("Calculate Nu", CALC_NUSSELT);
	saveMacroStart = p.getParameter<bool>("Save Macros On Start", THERMAL_SAVE_MACROS_ON_START);
	calcSSTempFlag = p.getParameter<bool>("Solve Steady Temp", SOLVE_SS_TEMP);
	chtCorrectionFlag = p.getParameter<bool>("Use CHT Correction", CHT_SOURCE_CORRECTION);
	
	testRestartRun();
	
	PrNum = p.getParameter<double>("Pr Number", PR_NUMBER);
	PrTurbNum = p.getParameter<double>("Pr Number", PR_TURB_NUMBER);
	Alpha_fluid = vlb.MuVal / PrNum;
	kSoot = p.getParameter<double>("K Soot", THERMAL_CONDUCTIVITY_FOUL);
	kAir = p.getParameter<double>("K air", THERMAL_CONDUCTIVITY_AIR);
	rhoSoot = p.getParameter<double>("Rho Soot", DENSITY_SOOT);
	cpSoot = p.getParameter<double>("Cp Soot", HEAT_CAPACITY_SOOT);
	
	
	Alpha_foul = (kSoot / rhoSoot / cpSoot)*(p.DELTA_L*p.DELTA_L / p.DELTA_T);

	kSootLB = kSoot * p.DELTA_F / p.DELTA_T;
	kAirLB = kAir * p.DELTA_F / p.DELTA_T;
	rhoSootLB = rhoSoot * p.DELTA_M / p.DELTA_L / p.DELTA_L / p.DELTA_L;
	rhoAirLB = vlb.RhoVal;
	cpSootLB = cpSoot * p.DELTA_L * p.DELTA_L / p.DELTA_T / p.DELTA_T;
	cpAirLB = kAirLB / rhoAirLB / Alpha_fluid;

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

	nuCutoffRadius = p.getParameter<double>("Nu Cutoff Radius", NU_CUTOFF_RADIUS);
	maxOutputLinesNu = p.getParameter<int>("Max Nu Output Lines", OUTPUT_MAX_LINES_NU);
}

void clVariablesFD::nuFindNearestNode(const int i, const int shiftInd)
{
	cl_ushort2 P01i = vtr.BL.P01ind(i);
	cl_double2 C0 = vls.C[P01i.x], C1 = vls.C[P01i.y];
	cl_int2 C0i = vls.min2(C0, C1), C1i = vls.max2(C0, C1);
	C0i = { { C0i.x - 1, C0i.y - 1 } };
	C1i = { { C1i.x + 1, C1i.y + 1 } };
	C0i.x = max(C0i.x, vtr.trDomainXBounds.x);
	C0i.y = max(C0i.y, 0);
	C1i.x = min(C1i.x, vtr.trDomainXBounds.y - 1);
	C1i.y = min(C1i.y, p.nY - 1);
	if (C0i.x >= vtr.trDomainXBounds.y)
		return;
	if (C1i.x < vtr.trDomainXBounds.x)
		return;

	double dist;
	cl_double2 vCcut, vN;

	for (int ii = C0i.x; ii <= C1i.x; ii++)
	{
		for (int jj = C0i.y; jj <= C1i.y; jj++)
		{
			if (vls.nType(ii, jj) & M_SOLID_NODE)
				continue;

			cl_double2 L0 = { { (double)ii, (double)jj } };
			if (bcFindIntersectionNormal(vCcut, dist, L0, C0, C1, vN) == TRUE)
			{
				if (dist < nuDist(ii - shiftInd))
				{
					nuDist(ii - shiftInd) = dist;
					nuDistVec(ii - shiftInd) =
					{ { (double)ii - vCcut.x, (double)jj - vCcut.y } };
					nuYInds(ii - shiftInd) = jj;
				}
			}
		}
	}
}

void clVariablesFD::renameSaveFiles()
{
	return;
	Temp.xVec->RenameTxtFile();
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
	alphaDer.save_txt_from_device();
	rhoCpDer.save_txt_from_device();
	Temp.saveAxb_w_indicies_from_device();
}

// Saves time output data (i.e. avg velocity, shear, Nu, etc)
void clVariablesFD::saveTimeData()
{
	if (calcNuFlag)
	{
		if (Nu.appendData_from_device(TRQUEUE_REF))
			Nu.FillBuffer(-1000., TRQUEUE_REF);
		NuMean.appendData_from_device(TRQUEUE_REF);
	}
}

void clVariablesFD::setKernelArgs()
{
	cl_int ind = 0;

	TempUpdateCoeffs.set_argument(ind++, Temp.get_add_IndArr());
	TempUpdateCoeffs.set_argument(ind++, vls.nType.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, vls.dXArr.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, Temp.get_add_A());
	TempUpdateCoeffs.set_argument(ind++, Temp.get_add_b());
	TempUpdateCoeffs.set_argument(ind++, Temp.get_add_Macro());
	TempUpdateCoeffs.set_argument(ind++, alphaDer.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, vlb.Ux_array.get_buf_add());
	TempUpdateCoeffs.set_argument(ind++, vlb.Uy_array.get_buf_add());
	if (vlb.kOmegaClass.kOmegaSolverFlag)
	{
		TempUpdateCoeffs.set_argument(ind++, vlb.kOmegaClass.Nut_array.get_buf_add());
	}
	if(chtCorrectionFlag)
	{
		TempUpdateCoeffs.set_argument(ind++, rhoCpDer.get_buf_add());
	}
	TempUpdateCoeffs.set_argument(ind++, &p.dTtfd);

	
	ind = 0;
	SetSSCoeffs.set_argument(ind++, Temp.get_add_IndArr());
	SetSSCoeffs.set_argument(ind++, vls.nType.get_buf_add());
	SetSSCoeffs.set_argument(ind++, vls.dXArr.get_buf_add());
	SetSSCoeffs.set_argument(ind++, Temp.get_add_A());
	SetSSCoeffs.set_argument(ind++, Temp.get_add_b());
	SetSSCoeffs.set_argument(ind++, Temp.get_add_Macro());
	SetSSCoeffs.set_argument(ind++, alphaDer.get_buf_add());
	if (vlb.kOmegaClass.kOmegaSolverFlag)
	{
		SetSSCoeffs.set_argument(ind++, vlb.kOmegaClass.Nut_array.get_buf_add());
	}
	if (chtCorrectionFlag)
	{
		SetSSCoeffs.set_argument(ind++, rhoCpDer.get_buf_add());
	}
	SetSSCoeffs.set_argument(ind++, vlb.Ux_array.get_buf_add());
	SetSSCoeffs.set_argument(ind++, vlb.Uy_array.get_buf_add());

	ind = 0;

	updateDAlphaKernel.set_argument(ind++, vls.nType.get_buf_add());
	updateDAlphaKernel.set_argument(ind++, vls.dXArr.get_buf_add());
	updateDAlphaKernel.set_argument(ind++, vls.dXArr0.get_buf_add());
	if(chtCorrectionFlag)
		updateDAlphaKernel.set_argument(ind++, rhoCpDer.get_buf_add());
	updateDAlphaKernel.set_argument(ind++, alphaDer.get_buf_add());



	
	ind = 0;
	calcNuKernel.set_argument(ind++, Temp.get_add_Macro());
	calcNuKernel.set_argument(ind++, Nu.get_buf_add());
	calcNuKernel.set_argument(ind++, vlb.Ux_array.get_buf_add());
	calcNuKernel.set_argument(ind++, vlb.Uy_array.get_buf_add());
	calcNuKernel.set_argument(ind++, nuYInds.get_buf_add());
	calcNuKernel.set_argument(ind++, nuDistVec.get_buf_add());
	calcNuKernel.set_argument(ind++, vls.dXArr.get_buf_add());
	calcNuKernel.set_argument(ind++, vls.dXArr0.get_buf_add());
	calcNuKernel.setOptionInd(ind);
	int zer = 0;
	calcNuKernel.setOption(&zer);

	ind = 0;
	updateNuKernel[0].set_argument(ind++, vls.C.get_buf_add());
	updateNuKernel[0].set_argument(ind++, vtr.BL.P01ind.get_buf_add());
	updateNuKernel[0].set_argument(ind++, vls.nType.get_buf_add());
	updateNuKernel[0].set_argument(ind++, nuYInds.get_buf_add());
	updateNuKernel[0].set_argument(ind++, nuDistVec.get_buf_add());
	updateNuKernel[0].set_argument(ind++, nuDist.get_buf_add());

	ind = 0;
	updateNuKernel[1].set_argument(ind++, nuYInds.get_buf_add());
	updateNuKernel[1].set_argument(ind++, nuDistVec.get_buf_add());


	ind = 0;
	calcNuMeanKernel.set_argument(ind++, vlb.Ux_array.get_buf_add());
	calcNuMeanKernel.set_argument(ind++, vlb.Uy_array.get_buf_add());
	calcNuMeanKernel.set_argument(ind++, Temp.get_add_Macro());
	calcNuMeanKernel.set_argument(ind++, NuMean.get_buf_add());
	calcNuMeanKernel.setOptionInd(ind);
	calcNuMeanKernel.set_argument(ind++, &zer);
#ifdef OPENCL_VERSION_1_2
	calcNuMeanKernel.set_local_memory(ind++, wgSizeNuMean * sizeof(double));
	calcNuMeanKernel.set_local_memory(ind++, wgSizeNuMean * sizeof(double));
#endif

}


void clVariablesFD::setSourceDefines()
{
	if (vfd.thermalSolverFlag)
	{
		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "USING_THERMAL_SOLVER");
	}
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "TIN", ROE_INX);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "TMIN", T_Actual_Min);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "TDIFF", T_Actual_Diff);

	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "FULLSIZEXM1", p.nX - 1);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_MULTIPLIER", 4.*p.Pipe_radius);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "LOCAL_SOLVE_SIZE_FD", 
		WORKGROUPSIZEX_FD*WORKGROUPSIZEY_FD);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "PR_TURB_NUMBER", PrTurbNum);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "PR_TURB_NUMBER_INV", 1./PrTurbNum);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "LB_ALPHA_FLUID", Alpha_fluid);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "LB_ALPHA_FOUL", Alpha_foul);

	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "K_SOOT_LB", kSootLB);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "K_AIR_LB", kAirLB);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "RHO_SOOT_LB", rhoSootLB);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "RHO_AIR_LB", rhoAirLB);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "CP_SOOT_LB", cpSootLB);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "CP_AIR_LB", cpAirLB);


	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "PR_NUMBER", PrNum);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "TFD_X_IN_VAL", ROE_INX);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_CUTOFF_RADIUS", nuCutoffRadius);
	if (chtCorrectionFlag)
	{
		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "USING_CHT_SOURCE_CORRECTION");
	}

	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WAVY_SECTION_START", p.xWavyStart);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WAVY_PERIOD_LENGTH", p.wavyPeriodLen);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NUM_NU_LOCS", p.numWavyPeriods+1);

	// These are set when iniNuCoeffs is called later, but if not using Nu this function will not be called,
	// so dummy values need to be set to ensure that opencl source compiles
	if (!calcNuFlag)
	{
		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_NUM_NODES", 0);
		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_NUM_NODES_FULLSIZE", 0);
		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NU_BL_PER_SIDE_NU", 0);
	}


}


void clVariablesFD::solveSSThermal()
{
	if (calcSSTempFlag == false || !thermalSolverFlag)
		return;

	SetSSCoeffs.call_kernel();

	Temp.setMaxIters(100);
	Temp.setAbsTol(1.e-12);
	Temp.setRelTol(1.e-5);
	Temp.solve();

	Temp.savetxt_from_device();
	Temp.solve();
	Temp.savetxt_from_device();
	Temp.solve();
	Temp.savetxt_from_device();
	Temp.solve();
	Temp.savetxt_from_device();


	Temp.setMaxIters(tempMaxIters);
	Temp.setAbsTol(tempMaxAbsTol);
	Temp.setRelTol(tempMaxRelTol);
}




void clVariablesFD::saveParams()
{
	p.setParameter("ParamsName", thermalParamNum);

	if (!thermalSolverFlag)
		return;

	if (p.Time > 0)
		p.setParameter("Restart Run", true);
	else
		p.setParameter("Restart Run", false);

	p.setParameter("Thermal Solver", thermalSolverFlag);
	p.setParameter("Calculate Nu", calcNuFlag);
	p.setParameter("Save Macros On Start", saveMacroStart);
	p.setParameter("Solve Steady Temp", calcSSTempFlag);
	p.setParameter("Use CHT Correction", chtCorrectionFlag);

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

	p.setParameter("Nu Cutoff Radius", nuCutoffRadius);
	p.setParameter("Max Nu Output Lines", maxOutputLinesNu);
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
	else if (vlb.loadVelTxtFlag)
	{
		restartRunFlag = false;
		if (tempArray.load("restart_files\\tempSolver"))
		{
			tempLoadedFlag = true;
		}
		else
		{
			tempLoadedFlag = false;
			tempArray.fill(ROE0);
		}
	}
	else
	{
		restartRunFlag = false;
		tempLoadedFlag = false;
		tempArray.fill(ROE0);
	}

	return restartRunFlag;
}

void clVariablesFD::update()
{
	updateDerivativeArrays();
}

void clVariablesFD::updateDerivativeArrays()
{
	updateDAlphaKernel.call_kernel();
}

void clVariablesFD::updateNuCoeffs()
{
	nuYInds.FillBuffer(-1, FDQUEUE_REF);
	nuDist.FillBuffer(100., FDQUEUE_REF);


	updateNuKernel[0].call_kernel();
	updateNuKernel[1].call_kernel();
}

void clVariablesFD::updateTimeData()
{
	if (calcNuFlag)
	{
		updateNuCoeffs();

		calcNuKernel.setOptionCallKernel(Nu.getCurIndAdd());

		

		//if (Nu.setTimeAndIncrement(p.Time, FDQUEUE_REF))
		//{
		//	// need to fill buffer with set value in order to know
		//	// which nodes do not have correct values at any given time.
		//	Nu.FillBuffer(-1000.);
		//}
		// SHouldnt need to fill this with -1000, because it will set empty nodes
		// in kernel
		Nu.setTimeAndIncrement(p.Time, FDQUEUE_REF);



		calcNuMeanKernel.setOptionCallKernel(NuMean.getCurIndAdd());
		NuMean.setTimeAndIncrement(p.Time, FDQUEUE_REF);
	}
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

//
//
//
//
//#ifdef USE_ORIG_NU 
//void clVariablesFD::reIniNuKernel()
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
//				if (bcFindIntersectionNormal(vCcut, dist, L0, C0, C1, vN) == true)
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
//#endif

//void clVariablesFD::updateTimeData()
//{
//	calcNu.setOptionCallKernel(&Save_loc_Nu);
//	clFlush(*calcNu.getDefaultQueue());
//	Save_loc_Nu++;
//	if (Save_loc_Nu == OUTPUT_MAX_LINES_NU)
//		saveTimeData();
//}
//


