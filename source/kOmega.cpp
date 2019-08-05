#include "kOmega.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"


// Currently, code uses blending function to calculate omega and U at nodes
// positioned beside the boundary while k is solved for like in other fluid nodes.
// k is set to 0 at boundary nodes inside the solid region (represent the value at the
// solid wall surface. This is correct according to theory (and suggested as a bc for k
// in https://turbmodels.larc.nasa.gov/sst.html). Literature is unclear about value of
// at wall surface, but this is not needed as omega is calculated with blending function
// at all nodes positioned beside wall (these nodes are the only ones that would need 
// value of omega at wall if solving with system of eqns.
// Turbulent viscosity is set to zero at wall (as expected from theory) and calculated from
// k and omega (and other flow parameters as well) at all fluid nodes during collision step.



void kOmega::allocateArrays()
{
	if (!kOmegaSolverFlag)
		return
	Kappa_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	Omega_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY);

	/////// Strain rate tensor elements (stress tensor can be calculated from this with bousinesque approximation)
	Sxy_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 3, 3);
	//////// Turbulent Viscosity (LB calculates this, so no need to initialize)
	Nut_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	////////// Diffusivities of Turb parameters
	Diff_Omega.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	Diff_K.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	//////// Coefficients calculated in first update step
	Fval_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY); //F1
	dKdO_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 2, 2); //2*(1-f1)*sigma_w2/omega*dk/dxi
}

void kOmega::allocateBuffers()
{
	Sxy_array.allocate_buffer_w_copy();
	Nut_array.allocate_buffer_w_copy();
	Diff_Omega.allocate_buffer_w_copy();
	Diff_K.allocate_buffer_w_copy();
	Fval_array.allocate_buffer_w_copy();
	dKdO_array.allocate_buffer_w_copy();
	WallD.allocate_buffer_w_copy();
}

void kOmega::createKernels()
{
	if (!kOmegaSolverFlag)
		return;

	int GlobalX = p.XsizeFull;
	int GlobalY = p.nY;


	kOmegaUpdateDiffCoeffs.create_kernel(GetSourceProgram, LBQUEUE_REF,
		"Update_kOmega_Diffusivities");
	kOmegaUpdateDiffCoeffs.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);

	kOmegaUpdateCoeffs.create_kernel(GetSourceProgram, LBQUEUE_REF,
		"Update_kOmega_Coeffs_Implicit");
	kOmegaUpdateCoeffs.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);

	updateWallDKernel.create_kernel(GetSourceProgram, LBQUEUE_REF,
		"updateWallD");
	updateWallDKernel.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);

}


// Iterates through BL indicies provided by first two arguments and 
// returns the minimum distance between nLoc and the wall.
double kOmega::findMinDist(cl_int2 botSearchInds, cl_int2 topSearchInds, cl_double2& nLoc)
{
	double dmin = 1000.;

	for (int i = botSearchInds.x; i < botSearchInds.x; i++)
	{
		cl_double2 P0 = vls.C(vls.BL(i, 0)), P1 = vls.C(vls.BL(i, 1));
		cl_double2 vT = Subtract2(P1, P0);
		double vTmag = GETLEN(vT);
		double dtemp = fabs(vT.y * nLoc.x - vT.x * nLoc.y + P1.x * P0.y - P1.y * P0.x) / vTmag;


		vT = Divide2(vT, vTmag);
		cl_double2 vN = { { -vT.y, vT.x} };

		cl_double2 Pcut = Multiply2(vN, dtemp);
		Pcut = Subtract2(nLoc, Pcut);
		if (Pcut.x < P0.x)
		{
			cl_double2 P0Loc = Subtract2(nLoc, P0);
			dtemp = GETLEN(P0Loc);
		}
		else if (Pcut.x > P1.x)
		{
			cl_double2 P1Loc = Subtract2(nLoc, P1);
			dtemp = GETLEN(P1Loc);
		}
		dmin = MIN(dmin, dtemp);
	}

	for (int i = topSearchInds.x; i < topSearchInds.x; i++)
	{
		cl_double2 P0 = vls.C(vls.BL(i, 0)), P1 = vls.C(vls.BL(i, 1));
		cl_double2 vT = Subtract2(P1, P0);
		double vTmag = GETLEN(vT);
		double dtemp = fabs(vT.y * nLoc.x - vT.x * nLoc.y + P1.x * P0.y - P1.y * P0.x) / vTmag;


		vT = Divide2(vT, vTmag);
		cl_double2 vN = { { -vT.y, vT.x} };

		cl_double2 Pcut = Multiply2(vN, dtemp);
		Pcut = Subtract2(nLoc, Pcut);
		if (Pcut.x > P0.x)
		{
			cl_double2 P0Loc = Subtract2(nLoc, P0);
			dtemp = GETLEN(P0Loc);
		}
		else if (Pcut.x < P1.x)
		{
			cl_double2 P1Loc = Subtract2(nLoc, P1);
			dtemp = GETLEN(P1Loc);
		}
		dmin = MIN(dmin, dtemp);
	}

	return dmin;
}



void kOmega::freeHostArrays()
{
	Sxy_array.FreeHost();
	Nut_array.FreeHost();
	Diff_Omega.FreeHost();
	Diff_K.FreeHost();
	Fval_array.FreeHost();
	dKdO_array.FreeHost();
}

// TODO: implement wall function values for non-parallel channels
//			this will most likely require a calculation of wall function values
//			which can be done in the kernel

// TODO: Figure out why using the same CSR_Inds instance for both k and omega 
//		solvers is causing error, and see if possible to fix it. 
//		there is no need to have two identical instances instead of a single one.
void kOmega::ini()
{
	if (!kOmegaSolverFlag)
		return;

	sourceGenerator::SourceInstance()->addFile2Kernel("kOmegaKernels.cl");
	


	// initialize WallD array
	iniWallD();

	// initialize k and omega arrays, solvers and reduction classes
	iniKOmegaArrays();

	// fill with visc since these are f(mu) just in case
	Diff_Omega.fillByNodeType(vlb.MuVal, vls.nType, M_FLUID_NODE);
	Diff_K.fillByNodeType(vlb.MuVal, vls.nType, M_FLUID_NODE);

	allocateBuffers();
}

void kOmega::iniKOmegaArrays()
{
	// set initial values of k, omega and nut
	setInitialValues();

	// Initialize solver classes for k and omega
	KappaInds.ini(p.nX, p.XsizeFull, p.nY, p.nY, &vls.nType);
	Kappa.CreateSolver(&Kappa_array, &KappaInds, LBQUEUE_REF,
		kappaMaxIters, kappaMaxRelTol, kappaMaxAbsTol);

	OmegaInds.ini(p.nX, p.XsizeFull, p.nY, p.nY, &vls.nType);
	Omega.CreateSolver(&Omega_array, &OmegaInds, LBQUEUE_REF,
		omegaMaxIters, omegaMaxRelTol, omegaMaxAbsTol);

	// initialize reduce classes for k and omega
	sumOmega.ini(*Omega.getMacroArray(), vlb.restartRunFlag, "redOmega");
	sumKappa.ini(*Kappa.getMacroArray(), vlb.restartRunFlag, "redKappa");
	minOmega.ini(*Omega.getMacroArray(), vlb.restartRunFlag, "minOmega");
	minKappa.ini(*Kappa.getMacroArray(), vlb.restartRunFlag, "minKappa");
}

void kOmega::iniTimeData()
{
	// no need to initialize anything here, since
	// reduce kernels are used elsewhere and therefore
	// will be initialized in the ini function
}


// Fills WallD with distance to nearest wall
void kOmega::iniWallD()
{
	WallD.fill(0.);

	for (int i = 0; i < p.nX; i++)
	{
		cl_int2 botSearchInds = { { MAX((vls.lsMap(i,0) - wallDSearchRar), 0),
			MIN((vls.lsMap(i,0) + wallDSearchRar), vls.nBL / 2)} };
		cl_int2 topSearchInds = { { MAX((vls.lsMap(i,1) - wallDSearchRar), vls.nBL / 2),
			MIN((vls.lsMap(i,1) + wallDSearchRar), vls.nBL) } };

		for (int j = 0; j < p.nY; j++)
		{
			if (vls.nType(i, j) & M_SOLID_NODE)
				continue;
			cl_double2 nLoc = { {(double)i, (double)j} };
			WallD(i,j) = findMinDist(botSearchInds, topSearchInds, nLoc);
		}
	}

#ifdef DEBUG_TURBARR
	WallD.savetxt("WallD");
#endif
}


void kOmega::loadParams()
{
	iniTurbVelocity = p.getParameter<bool>("Ini Turb Velocity", INI_TURB_PROFILE);
	perturbVelField = p.getParameter<bool>("Perturb Vel Field", PERTURB_VEL_FIELD);

	perturbDUPlus = p.getParameter<double>("DUplus Perturb Multiplier", DUPLUS_MULT);
	perturbEpsilon = p.getParameter<double>("Epsilon Perturb Multiplier", EPSILON_MULT);
	perturbBeta = p.getParameter<double>("Beta Perturb Multiplier", BETA_MULT);
	perturbAlpha = p.getParameter<double>("Alpha Perturb Multiplier", ALPHA_MULT);
	perturbSigma = p.getParameter<double>("Sigma Perturb Multiplier", SIGMA_MULT);

	roughnessFactor = p.getParameter<double>("Wall Roughness Factor", ROUGHNESS_FACTOR);

	kappaMaxRelTol = p.getParameter<double>("Kappa Max Rel Tol", KOMEGA_MAX_REL_TOL);
	kappaMaxAbsTol = p.getParameter<double>("Kappa Max Abs Tol", KOMEGA_MAX_ABS_TOL);
	kappaMaxIters = p.getParameter<int>("Kappa Max Iterations", KOMEGA_MAX_ITERS);

	omegaMaxRelTol = p.getParameter<double>("Omega Max Rel Tol", kappaMaxRelTol);
	omegaMaxAbsTol = p.getParameter<double>("Omega Max Abs Tol", kappaMaxAbsTol);
	omegaMaxIters = p.getParameter<int>("Omega Max Iterations", kappaMaxIters);

	wallDSearchRar = p.getParameter<int>("Wall D Search Rad", WALLD_SEARCH_RADIUS);

	kIniVal = p.getParameter<double>("K Initial Val", -1.);
	omegaIniVal = p.getParameter<double>("K Initial Val", -1.);
}


void kOmega::renameSaveFiles()
{
	Kappa.xVec->RenameTxtFile();
	Omega.xVec->RenameTxtFile();
}


void kOmega::save2file()
{
	Kappa.savetxt_from_device();
	Omega.savetxt_from_device();
}

void kOmega::saveDebug(int saveFl)
{
#ifdef DEBUG_TURBARR
	Array2Dd outtemp;
	outtemp.zeros(p.nX, p.nY);

	if (saveFl == koDbgSave || saveFl == DbgSave || saveFl == koDbgSave1)
	{
		std::vector<std::string> koDbgNames1{ "dKdx", "dKdy", "dOdx", "dOdy", "CDkw" };

		koDbgArr1.read_from_buffer();

		for (int k = 0; k < koDbgNames1.size(); k++)
		{
			std::string outname = "koDbg1_";
			outname.append(koDbgNames1[k]);
			for (int i = 0; i < p.nX; i++)
			{
				for (int j = 0; j < p.nY; j++)
				{
					outtemp(i, j) = koDbgArr1(i, j, k);
				}
			}
			outtemp.savetxt(outname.c_str());
		}
	}

	if (saveFl == koDbgSave || saveFl == DbgSave || saveFl == koDbgSave2)
	{
		std::vector<std::string> koDbgNames2{ "diffk_dx", "diffk_dy", "diffo_dx", "diffo_dy",
			"Jx_omega", "Jy_omega", "Jx_k", "Jy_k", "Pk", "Sc", "Kc", "Ke", "Kw", "Kn", "Ks",
			"Ksrc", "Wc", "We", "Ww", "Wn", "Ws", "Wsrc" };

		koDbgArr2.read_from_buffer();

		for (int k = 0; k < koDbgNames2.size(); k++)
		{
			std::string outname = "koDbg2_";
			outname.append(koDbgNames2[k]);
			for (int i = 0; i < p.nX; i++)
			{
				for (int j = 0; j < p.nY; j++)
				{
					outtemp(i, j) = koDbgArr2(i, j, k);
				}
			}
			outtemp.savetxt(outname.c_str());
		}
	}
#endif
	Fval_array.save_txt_from_device("Fvals");
	Diff_Omega.save_txt_from_device("Diff_Omega");
	Diff_K.save_txt_from_device("Diff_K");
}


void kOmega::saveRestartFiles()
{
	Kappa.saveCheckPoint();
	Omega.saveCheckPoint();
}



void kOmega::saveParams()
{
	p.setParameter("Use Turb Model", kOmegaSolverFlag);

	p.setParameter("Kappa Max Rel Tol", kappaMaxRelTol);
	p.setParameter("Kappa Max Abs Tol", kappaMaxAbsTol);
	p.setParameter("Kappa Max Iterations", kappaMaxIters);
	p.setParameter("Omega Max Rel Tol", omegaMaxRelTol);
	p.setParameter("Omega Max Abs Tol", omegaMaxAbsTol);
	p.setParameter("Omega Max Iterations", omegaMaxIters);

	p.setParameter("DUplus Perturb Multiplier", perturbDUPlus);
	p.setParameter("Epsilon Perturb Multiplier", perturbEpsilon);
	p.setParameter("Beta Perturb Multiplier", perturbBeta);
	p.setParameter("Alpha Perturb Multiplier", perturbAlpha);
	p.setParameter("Sigma Perturb Multiplier", perturbSigma);
	p.setParameter("Wall Roughness Factor", roughnessFactor);

	p.setParameter("Wall D Search Rad", wallDSearchRar);

	p.setParameter("K Initial Val", kIniVal);
	p.setParameter("Omega Initial Val", omegaIniVal);

	// These do not need to be save since they relate to initialization
	//p.setParameter("Ini Turb Velocity", false);
	//p.setParameter("Perturb Vel Field", false);

}


void kOmega::setInitialValues()
{
	// Initial values are calculated from turbulent intensity and length scales assuming that
	// velocity is f(Re, geomPenaltyFactor).

	// Calculate actual mean velocity from Ux and rho (in case it is a restart
	// after getting initial values only). If it is zero, Umean is calculated
	// from set Re and modified by geomPenaltyFactor.
	double umean_ini = vlb.Ux_array.reduce() / vlb.Ro_array.reduce();
	if (umean_ini == 0.)
	{
		umean_ini = vlb.Re * vlb.MuVal / p.Pipe_radius / 4. * vlb.geomPenaltyFactor;
	}


	double Reini = p.Pipe_radius * 4. * umean_ini / vlb.MuVal;
	TurbIntensity = 0.16 * pow(Reini, -1. / 8.);
	TurbLScale_inv = 1. / (0.036 * p.Pipe_radius);

	if (kIniVal <= 0.)
		kIniVal = 1.5 * (TurbIntensity * umean_ini) * (TurbIntensity * umean_ini);

	if (omegaIniVal <= 0.)
		omegaIniVal = sqrt(kIniVal) / KO_C_MU_0_25 * TurbLScale_inv;

	ERROR_CHECKING((kIniVal <= 0.), "Initial value for initial tKE <= 0.\n", \
		ERROR_INITIALIZING_KOMEGA);

	ERROR_CHECKING((omegaIniVal <= 0.), "Initial value for initial omega <= 0.\n", \
		ERROR_INITIALIZING_KOMEGA);

	double turbViscIniVal = kIniVal / omegaIniVal;

	ERROR_CHECKING((turbViscIniVal > MAX_TURB_VISC_VALUE), "Initial Turbulent viscosity greater "\
		"than set maximum value\n", ERROR_INITIALIZING_KOMEGA);

	// Fill with initial values if not loaded from bin
	Omega_array.fillByNodeType(omegaIniVal, vls.nType, M_FLUID_NODE);
	Kappa_array.fillByNodeType(kIniVal, vls.nType, M_FLUID_NODE);
	Nut_array.fillByNodeType(turbViscIniVal, vls.nType, M_FLUID_NODE);
}


void kOmega::setKernelArgs()
{

	cl_int ind = 0;
	int zer = 0;

	kOmegaUpdateDiffCoeffs.set_argument(ind++, vls.nType.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Kappa.get_add_Macro());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Omega.get_add_Macro());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Nut_array.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, WallD.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Diff_K.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Diff_Omega.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Fval_array.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, dKdO_array.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, vlb.Ro_array.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, vls.dXArr.get_buf_add());
#ifdef DEBUG_TURBARR
	kOmegaUpdateDiffCoeffs.set_argument(ind++, koDbgArr1.get_buf_add());
#endif


	ind = 0;
	kOmegaUpdateCoeffs.set_argument(ind++, Omega.get_add_IndArr());
	kOmegaUpdateCoeffs.set_argument(ind++, Kappa.get_add_Macro());
	kOmegaUpdateCoeffs.set_argument(ind++, Omega.get_add_Macro());
	kOmegaUpdateCoeffs.set_argument(ind++, Nut_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, Diff_K.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, Diff_Omega.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, Fval_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, dKdO_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, Sxy_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, vlb.Ro_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, vlb.Ux_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, vlb.Uy_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, Kappa.get_add_A());
	kOmegaUpdateCoeffs.set_argument(ind++, Omega.get_add_A());
	kOmegaUpdateCoeffs.set_argument(ind++, Kappa.get_add_b());
	kOmegaUpdateCoeffs.set_argument(ind++, Omega.get_add_b());
	kOmegaUpdateCoeffs.set_argument(ind++, vls.dXArr.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, WallD.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, vls.lsMap.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, vtr.wallShear.Tau.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, vls.nType.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, &p.dTlb);
#ifdef DEBUG_TURBARR
	kOmegaUpdateCoeffs.set_argument(ind++, koDbgArr2.get_buf_add());
#endif


	ind = 0;
	updateWallDKernel.set_argument(ind++, WallD.get_buf_add());
	updateWallDKernel.set_argument(ind++, vls.lsMap.get_buf_add());
	updateWallDKernel.set_argument(ind++, vls.C.get_buf_add());
	updateWallDKernel.set_argument(ind++, vls.nType.get_buf_add());

}


void kOmega::saveTimeData()
{

}

#define SETSOURCEDEFINE SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr()

void kOmega::setSourceDefines()
{
#ifdef DEBUG_TURBARR
	SETSOURCEDEFINE, "DEBUG_ARRAYS");
#endif
#ifdef DEBUG_TURBARR
	SETSOURCEDEFINE, "DISABLE_TURBULENT_VISC");
#endif

	SETSOURCEDEFINE, "RE_TURBULENT", ReTurbVal);
	SETSOURCEDEFINE, "UTAU_VAL", UtauVal);
	SETSOURCEDEFINE, "YPLUS_WALL_NODE", yPlusWallVal);
	SETSOURCEDEFINE, "K_WALL_VALUE", kappaWallVal);
	SETSOURCEDEFINE, "NUT_WALL_VALUE", nutWallVal);
	SETSOURCEDEFINE, "OMEGA_WALL_VALUE", omegaWallVal);
	SETSOURCEDEFINE, "WALLD_SEARCH_RADIUS", wallDSearchRar);
	if(kOmegaSolverFlag)
		SETSOURCEDEFINE, "USING_KOMEGA_SOLVER");

}

#undef SETSOURCEDEFINE


void kOmega::Solve()
{
	kOmegaUpdateDiffCoeffs.call_kernel();
	kOmegaUpdateCoeffs.call_kernel();
	clFinish(LBQUEUE);
	Kappa.solve();
	Omega.solve();
}




bool kOmega::testRestartRun()
{
	// if not solving, return true since we dont want to show as
	// not a restart when kappa and omega dont load
	if (!kOmegaSolverFlag)
		return true;

	allocateArrays();

	bool koBool = Kappa_array.load("load" SLASH "lbkappa") && Omega_array.load("load" SLASH "lbomega");
	if (koBool)
	{
		kOmegaLoadedFlag = true;
	}
	else
	{
		kOmegaLoadedFlag = false;
	}

	return koBool;
}


void kOmega::updateTimeData()
{
}





///////////////////////////////////////////////////////////////////////////////////////
//////////                                                                  ///////////
//////////   Calculation of wall values used for straight channel tests.    ///////////
//////////                                                                  ///////////
///////////////////////////////////////////////////////////////////////////////////////

//// calculation of yplus_lam
//// specific to flow in straight channel
//// based on openfoam's implementation
//// see http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2016/FangqingLiu/openfoamFinal.pdf
//// and openfoam documentation
//
//double ypl = 11.;
//for (int i = 0; i < 11; i++)
//{
//	ypl = log(MAX(roughnessFactor * ypl, 1)) / 0.41;
//}
//
//ReTurbVal = sqrt(vlb.Re * 3.);
//yPlusWallVal = ReTurbVal / p.Pipe_radius / 2.;
//UtauVal = ReTurbVal * vlb.MuVal / p.Pipe_radius;
//if (yPlusWallVal <= ypl)
//{
//	double Cf = (1. / ((yPlusWallVal + 11.) * (yPlusWallVal + 11.)) +
//		2. * yPlusWallVal / 11. / 11. / 11. - 1. / 11. / 11.);
//	kappaWallVal = UtauVal * UtauVal * (2400. / 1.9 / 1.9 * Cf);
//}
//else
//{
//	// Nothing implemented for yPlus < yPlus_laminar
//	exit(ERROR_INITIALIZING_KOMEGA);
//}
//
//nutWallVal = vlb.MuVal * (0.41 * yPlusWallVal / log(roughnessFactor * yPlusWallVal) - 1.);
//double omega_vis = 6. * vlb.MuVal / (0.075 * yPlusWallVal * yPlusWallVal);
//double omega_log = UtauVal / (0.3 * yPlusWallVal * 0.41);
//omegaWallVal = sqrt(omega_vis * omega_vis + omega_log * omega_log);
//omegaWallVal = MAX(omegaWallVal, 0.);
//double umean_ini = vlb.Ux_array.reduce() / vlb.Ro_array.reduce();
//if (umean_ini == 0.)
//{
//	umean_ini = vlb.Re * vlb.MuVal / p.Pipe_radius / 4.;
//}
//
//double Reini = p.Pipe_radius * vlb.Re / vlb.MuVal * 4.;
//TurbIntensity = 0.16 * pow(Reini, -1. / 8.);
//TurbLScale_inv = 1. / (0.036 * p.Pipe_radius);
//
//double kIniVal = 1.5 * (TurbIntensity * umean_ini) * (TurbIntensity * umean_ini);
//
//double omegaIniVal = sqrt(kIniVal) / KO_C_MU_0_25 * TurbLScale_inv;
//double turbViscIniVal = kIniVal / omegaIniVal;
//
//int i = 0;
//while (vls.nType(0, i) & M_SOLID_NODE)
//{
//	i++;
//}
//int wallindlow = i;
//while (vls.nType(0, i) & M_FLUID_NODE)
//{
//	i++;
//}
//
//std::vector<int> wallinds = { wallindlow, i - 1 };
//Kappa.setInitialValueRows(kappaWallVal, wallinds); // set to inlet val
//Omega.setInitialValueRows(omegaWallVal, wallinds); // set to inlet val
