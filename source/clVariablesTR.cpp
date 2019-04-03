// clVariables.cpp: implementation of the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////


//TODO:	save necessary variables, finish organizing functions and variables,
//		setup all load from bin files to allow a restart


#include "StdAfx.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"
#include "clProblem.h"
#include <random>

#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(),
void clVariablesTR::setSourceDefines()
{
	setSrcDefinePrefix "START_THERMO_VEL", startThermoVel);
	setSrcDefinePrefix "MAX_BL_PER_NODE", maxBLPerNode);
	setSrcDefinePrefix "TRC_NUM_TRACERS", nN);
	setSrcDefinePrefix "NUM_PAR_SIZES", Nd);
	setSrcDefinePrefix "FULLSIZEX_TR", TrDomainSize.x);
	setSrcDefinePrefix "FULLSIZEY_TR", TrDomainSize.y);
	setSrcDefinePrefix "MU_NUMBER", vlb.MuVal);
	setSrcDefinePrefix "DTTR", p.dTtr);
	setSrcDefinePrefix "DTTR_WALL", p.dTtr_wall);
	setSrcDefinePrefix "FOUL_SIZE_SWITCH_SOOT2", foulSizeSwitchSoot);
	setSrcDefinePrefix "KAIR", vfd.kAir*p.DELTA_F/p.DELTA_T);
	setSrcDefinePrefix "KSOOT", vfd.kSoot*p.DELTA_F / p.DELTA_T);
	setSrcDefinePrefix "ALPHA_DEN_AIR", vfd.Alpha_den_air);
	setSrcDefinePrefix "ALPHA_DEN_SOOT", vfd.Alpha_den_soot);
	setSrcDefinePrefix "ALPHA_TEST_AIR", ((vfd.kAir * p.DELTA_F / p.DELTA_T) / vfd.Alpha_den_air));
	setSrcDefinePrefix "ALPHA_TEST_SOOT", ((vfd.kSoot * p.DELTA_F / p.DELTA_T) / vfd.Alpha_den_soot));
	setSrcDefinePrefix "PAR_TIMER_START", MIN(timeBeforeReRelease, 20000));
	setSrcDefinePrefix "DEP_TIMER_START", MIN(timeBeforeStick, 20000));
	setSrcDefinePrefix "X_START_VAL_INT", xReleasePos);
	setSrcDefinePrefix "X_MAX_VAL_INT", xStopPos);
	setSrcDefinePrefix "X_START_VAL", (double)xReleasePos);
	setSrcDefinePrefix "X_MIN_VAL", (double)xStopPos - 0.5);
	setSrcDefinePrefix "X_MAX_VAL", xStopPos);
	setSrcDefinePrefix "Y_MIN_VAL", 0.5);
	setSrcDefinePrefix "Y_MAX_VAL", p.nY - 0.5);
	setSrcDefinePrefix "XVALS_HEIGHT", p.Channel_Height+2);
	setSrcDefinePrefix "NUM_PAR_GL_DIV", NUM_PAR_GL_DIV);
	setSrcDefinePrefix "ALPHA_FLUID", vfd.Alpha_fluid);
	setSrcDefinePrefix "ALPHA_FOUL", vfd.Alpha_foul);
	if (MAX_NUM_BL_ROLL != -1)
	{
		setSrcDefinePrefix "MAX_NUM_BL_ROLL", maxNumBLRoll);
	}
	int count = 1;
	int tempvar = 2;
	while (tempvar < WORKGROUPSIZE_SORT)
	{
		count++;
		tempvar *= 2;
	}
		
	setSrcDefinePrefix "SORT_NUM_MERGES", count);
	
	// The defines are added even when trSolver is not being used just in case
	// a kernel references one of these defines (due to methods occasionally
	// overlapping in functions). This variable is not included in files outside
	// of TR_kernels, so it will be fine to not include it in definition file.
	if (trSolverFlag)
	{
		char parMultVal[80];
		sprintf(parMultVal, "%20.16g", V_par[0] / vls.lsSpacing);
		std::string parmult = "const double Par_multiplier[" + std::to_string(Nd) + "] = { " + parMultVal;

		for (int i = 1; i < vtr.Nd; i++)
		{
			sprintf(parMultVal, ", %20.16g, ", V_par[i] / vls.lsSpacing);
			parmult.append(parMultVal);
		}

		parmult.append("};");
		SOURCEINSTANCE->addString2Defines(parmult);
	}
}
#undef setSrcDefinePrefix


// TODO: Recalculate vtr.TrDomainSize in vtr.loadParams function 
void clVariablesTR::ini()
{
	setSourceDefines();
	if (!trSolverFlag)
		return;
	Ploc_inds = { { -1, -1, -1, 0 } };
	rmax = (double)RAND_MAX + 1.;
	X_release = X_RELEASE_POS;



	P.allocate(nN);
	NodV.allocate(p.nX - 1, p.nY - 1);//Stores Velocities and Temps used for calculating particle vels
	NodI.allocate(p.nX - 1, p.nY - 1);//Contains boundary links near lattice site as well as flag for wall 
	NodC.allocate(p.nX - 1, p.nY - 1);//Coefficients used to calculate NodV each time step
	BL.allocate(vls.nBL);//Boundary link info
	Active_flag.zeros(p.nX - 1, p.nY - 1);
	BL_dep.zeros(vls.nBL, Nd);
	RandList.allocate(nN);//rng array
	Ploc.allocate(vtr.TrDomainSize.x*vtr.TrDomainSize.y + 2);
	BLind_ind.zeros(p.nX - 1, p.nY - 1); 

	Array1Di TR_Load_Param;
	TR_Load_Param.zeros(6);
	if (TR_Load_Param.load("load" SLASH "TR_Params") == TRUE)
	{
		nActiveNodes = TR_Load_Param(0);
		Sort_timer = TR_Load_Param(1);
		Bounds.MIN_BL_BOT = TR_Load_Param(2);
		Bounds.MAX_BL_BOT = TR_Load_Param(3);
		Bounds.MIN_BL_TOP = TR_Load_Param(4);
		Bounds.MAX_BL_TOP = TR_Load_Param(5);

		BL_dep.load("load" SLASH "BLdep_temp");
		P.load("load" SLASH "trp");
		Ploc.load("load" SLASH "Ploc");
		ini_trp();
		restartRunFlag = 1;
	}
	else
	{
		
		Sort_timer = 0;
		ini_trp();
		ini_particles();
		Ploc(0) = Ploc_inds.lo;
		Ploc(1) = Ploc_inds.hi;
		restartRunFlag = 0;
	}

	if (RandList.load("load" SLASH "RandList") == FALSE)
		ini_rand();

	ini_node();
	ini_blinks();
	Split_Node_arrays();
		
	//createKernels();

	iniShear();

	Find_node_neighs();
	
	allocate_buffers();
	
	ini_shear_coeffs();
	
	ini_sort();

	setSourceDefines();
	std::function<void(void)> createKerPtr = std::bind(&clVariablesTR::createKernels, this);
	std::function<void(void)> setArgsPtr = std::bind(&clVariablesTR::setKernelArgs, this);

	sourceGenerator::SourceInstance()->addIniFunction(createKerPtr, setArgsPtr);
	//setKernelArgs();

	if (saveMacroStart)
		save2file();

#ifdef SAVE_BIN_IN_LOAD
	save_restart_files_ini();
	saveRestartFiles();
#endif
	Save_loc_SS = 0;

	Winds.save_txt_from_device("Winds");
	Node_neigh.save_txt_from_device("Nneigh");
}

void clVariablesTR::UpdateSS()
{
	Save_shear_kernel.set_argument(4, &Save_loc_SS);
	Save_shear_kernel.call_kernel();
	clFlush(TRQUEUE);
	Save_loc_SS++;
}
void clVariablesTR::reset_SS_out()
{
	Save_loc_SS = 0;
}

void clVariablesTR::Update_IO_Dists()
{

	Save_Dists_kernel.set_argument(6, &Save_IO_Loc);
	Save_Dists_kernel.call_kernel();
	clFinish(TRQUEUE);
	Save_IO_Loc++;
}

void clVariablesTR::Reset_IO_Dists()
{
	Save_IO_Loc = 0;
	IO_dists_save.FillBuffer(0);
}

void clVariablesTR::testRestartRun()
{
	//allocateArrays();
	restartRunFlag = p.getParameter("Restart Run", false);
	if (restartRunFlag == false)
		return;
	if (BL_dep.load("load" SLASH "BLdep_temp") == FALSE)
		restartRunFlag = false;
	if (P.load("load" SLASH "trp") == FALSE)
		restartRunFlag = false;
	if (Ploc.load("load" SLASH "Ploc") == FALSE)
		restartRunFlag = false;
}


void clVariablesTR::loadParams()
{
	trSolverFlag = p.getParameter("Use Par Solver",
		USE_PARTICLE_SOLVER);
	saveMacroStart = p.getParameter("Save Macros On Start", TR_SAVE_MACROS_ON_START);
	testRestartRun();

	nN = p.getParameter("Num Particles", TRC_NUM_TRACERS);
	if (trSolverFlag)
		Nd = p.getParameter<int>("Num Par Sizes");
	else
		Nd = 0;
	numStepsBtwSort = p.getParameter("Num Steps Btw Sort", NUM_STEPS_BTW_SORT);
	sootNumConc = p.getParameter("Soot Num Concentration", SOOT_NUMBER_CONCENTRATION);
	
	parThermalCond = p.getParameter("Par Themal Conductivity", THERMAL_CONDUCTIVITY_PARTICLE);
	
	surfEnergySurf = p.getParameter("Surface Energy Wall", SURF_ENERGY_SURF);
	poissonSurf = p.getParameter("Poisson Number Wall", POISSON_SURF);
	yModSurf = p.getParameter("Youngs Modulus Wall", Y_MOD_SURF);

	surfEnergySoot = p.getParameter("Surface Energy Soot", SURF_ENERGY_SOOT);
	poissonSoot = p.getParameter("Poisson Number Soot", POISSON_SOOT);
	yModSoot = p.getParameter("Youngs Modulus Soot", Y_MOD_SOOT);

	mfpAir = p.getParameter("MFP Air", MEAN_FREE_PATH_AIR);

	hamakerConst = p.getParameter("Hamaker Constant", HAMAKER_CONST);
	depPorosity = p.getParameter("Deposit Porosity", DEP_POROSITY);
	wallCorrection = p.getParameter("Wall Correction", WALL_CORRECTION);
	liftCoeff = p.getParameter("Lift Coefficient", LIFT_COEFFICIENT);

	indRadiusSearch = p.getParameter("Index Radius Search", INDEX_RADIUS_SEARCH);
	massFluxInlet = p.getParameter("Mass Flux Inlet", MASS_FLUX_INLET);
	maxBLPerNode = p.getParameter("Max BL Per Node", MAX_BL_PER_NODE);
	parReleaseTime = p.getParameter("Particle Release Time", PARTICLE_RELEASE_TIME);
	stopDistX = p.getParameter("X Stop Distance", STOP_DIST_X);
	timeBeforeReRelease = p.getParameter("Time Before Re-release", TIME_BEFORE_RERELEASE) / p.trSteps_wall;
	timeBeforeStick = p.getParameter("Time Before Stick", TIME_BEFORE_STICK) / p.trSteps_wall;
	maxNumBLRoll = p.getParameter("Max Num BL Roll", MAX_NUM_BL_ROLL);
	xReleasePos = p.getParameter("X Release Pos", X_RELEASE_POS);
	xStopPos = p.getParameter("X Stop Pos", X_STOP_POS);
	reduceDepStop1 = p.getParameter("Reduce Dep Stop1", REDUCE_DEP_STOP1);
	reduceDepStop2 = p.getParameter("Reduce Dep Stop2", REDUCE_DEP_STOP2);
	startThermoVel = p.getParameter("Start Thermo Pos", START_THERMO_VEL);
	amtReduceDep = p.getParameter("Amount Reduce Dep", AMT_REDUCE_DEP);
	xMinVal = p.getParameter("X Min Val", X_MIN_VAL);
	cutoffRadius = p.getParameter("Cutoff Radius", CUTOFF_RADIUS);
	foulSizeSwitchSoot = p.getParameter("Foul Size Switch Soot", FOUL_SIZE_SWITCH_SOOT2);
	numEachPar = p.getParameter("Num Each Par", NUM_EACH_PAR);
	parVolMultiplier = numEachPar / (1. - depPorosity);
	
	if (trSolverFlag == false)
		return;

	D_p_real.zeros(Nd);
	D_dists.zeros(Nd);
	p.yamlIn["Dp Dists"] >> D_dists;
	p.yamlIn["Par Diameters"] >> D_p_real;

	double Num_per_m2 = sootNumConc / p.DELTA_L / p.DELTA_L;
	double Num_per_m = Num_per_m2 * 2. * p.Pipe_radius;
	double Umean = vlb.Re * (vlb.tau - 0.5) / 3. / p.Pipe_radius;
	double Num_per_s = Num_per_m;

	Par_conc_inlet = (cl_uint)Num_per_s * (cl_uint)numStepsBtwSort; //number of particles released after each sort step 
	//must be multiplied by Umean to get true value

	WofA = 2.*sqrt(surfEnergySurf * surfEnergySoot);
	Estar = 1. / ((1. - poissonSurf*poissonSurf) / yModSurf + (1. - poissonSoot*poissonSoot) / yModSoot);
	Estar_s = yModSoot;
	WofA_s = surfEnergySoot;




	double ConstA = 1.2, ConstB = 0.4, ConstC = 1.1;
	double ConstCs = 1.14, ConstCm = 1.17, ConstCt = 2.18;

	Kth_pars.zeros(Nd);
	D_p.zeros(Nd);
	Q_A_prime.zeros(Nd, 2);
	Q_A.zeros(Nd, 2);
	M_p.zeros(Nd);
	F_po.zeros(Nd);
	R_d.zeros(Nd, 2);
	Tau_crit.zeros(Nd, 2);
	V_par.zeros(Nd);
	Tau_crit_max = 0.;

	for (int i = 0; i < Nd; i++)
	{
		double Knud = mfpAir / D_p_real(i);

		double CCF = 1. + Knud * (ConstA + ConstB * exp(-ConstC / Knud));
		double k_ratio = vfd.kAir / parThermalCond;
		double Kth1 = 2. * ConstCs * CCF / (1. + 3.*ConstCm * Knud);
		double Kth2 = (k_ratio + ConstCt * Knud) / (1. + 2. * k_ratio + 2.*ConstCt*Knud);

		Kth_pars(i) = Kth1*Kth2;

		double a3 = 9.*PI_NUMBER*WofA*D_p_real(i)*D_p_real(i) / Estar / 4.;
		double Rdi = pow(a3, (1. / 3.));
		Rdi *= p.DELTA_L;
		R_d(i, 0) = Rdi;

		double a3_s = 9.*PI_NUMBER*WofA_s*D_p_real(i)*D_p_real(i) / Estar_s / 4.;
		double Rdi_s = pow(a3_s, (1. / 3.));
		Rdi_s *= p.DELTA_L;
		R_d(i, 1) = Rdi_s;

		double QAi = 2.*WofA*PI_NUMBER*R_d(i, 0)*R_d(i, 0);
		QAi *= (p.DELTA_F / p.DELTA_L);
		Q_A(i, 0) = QAi;

		double QAi_s = 2.*WofA_s*PI_NUMBER*R_d(i, 1)*R_d(i, 1);
		QAi *= (p.DELTA_F / p.DELTA_L);
		Q_A(i, 1) = QAi_s;

		double Q_A_pt = pow(WofA, 5.)*pow(D_p_real(i), 4.) / Estar / Estar / 16;
		double QA_pi = 7.09 * pow(Q_A_pt, (1. / 3.));
		QA_pi *= p.DELTA_P;
		Q_A_prime(i, 0) = QA_pi;

		double Q_A_pt_s = pow(WofA_s, 5.)*pow(D_p_real(i), 4.) / Estar_s / Estar_s / 16;
		double QA_pi_s = 7.09 * pow(Q_A_pt_s, (1. / 3.));
		QA_pi_s *= p.DELTA_P;
		Q_A_prime(i, 1) = QA_pi_s;

		double Mpi = PI_NUMBER*D_p_real(i)*D_p_real(i)*vfd.rhoSoot / 4.;
		Mpi *= p.DELTA_M;
		M_p(i) = Mpi;

		double F_poi = 625.*hamakerConst / 3. / D_p_real(i);
		F_poi *= p.DELTA_F;
		F_po(i) = F_poi;

		double Dpi = D_p_real(i)*p.DELTA_L;
		D_p(i) = Dpi;

		double Vpi = D_p(i) * D_p(i) * PI_NUMBER / 4. / (1. - depPorosity);
		V_par(i) = Vpi;

		double Tci = 4.*F_po(i)*R_d(i, 0) / (3. * PI_NUMBER * pow(D_p(i), 3.) * wallCorrection);

		Tau_crit(i, 0) = Tci;
		if (Tau_crit(i, 0) > Tau_crit_max)
			Tau_crit_max = Tau_crit(i, 0);

		double Tci_s = 4.*F_po(i)*R_d(i, 1) / (3. * PI_NUMBER * pow(D_p(i), 3.) * wallCorrection);


		Tau_crit(i, 1) = Tci_s;
		if (Tau_crit(i, 1) > Tau_crit_max)
			Tau_crit_max = Tau_crit(i, 1);
	}

	double Distsum = 0.;
	Mean_index = 0;
	parP.allocate(Nd);
	double dist_max = 0;
	int max_dist_num = 0;
	for (int i = 0; i < vtr.Nd; i++)
	{
		parP(i).Dp = D_p(i);
		parP(i).Q_A.x = Q_A(i, 0);
		parP(i).Q_A.y = Q_A(i, 1);
		parP(i).Q_A_prime.x = Q_A_prime(i, 0);
		parP(i).Q_A_prime.y = Q_A_prime(i, 1);
		parP(i).tau_crit.x = Tau_crit(i, 0);
		parP(i).tau_crit.y = Tau_crit(i, 1);
		parP(i).Mp = M_p(i);
		parP(i).Kth = Kth_pars(i);
		parP(i).D_dist = Distsum;
		if (dist_max < D_dists(i))
			max_dist_num = i;
		Distsum += D_dists(i);
		parP(i).L_coeff = liftCoeff * parP[i].Dp * parP[i].Dp * parP[i].Dp / vlb.MuVal / 8.;
		parP(i).D_coeff = 6 * PI_NUMBER * parP[i].Dp * parP[i].Dp * wallCorrection / 4.;
	}

	Tcrit_color.zeros(2);
	Tcrit_color(0) = parP(max_dist_num).tau_crit.x;
	Tcrit_color(1) = parP(max_dist_num).tau_crit.y;

	save_particle_arrays();	
}


void clVariablesTR::iniShear()
{
	// Sizes of arrays and kernels are function of number of BL's
	// and therefore constant (will not need reallocation at update step)
	numbl_bounds = (Bounds.MAX_BL_BOT - Bounds.MIN_BL_BOT) + (Bounds.MAX_BL_TOP - Bounds.MIN_BL_TOP) + 2;
	int shearsize = (int)ceil((double)numbl_bounds / WORKGROUPSIZE_TR_SHEAR) * WORKGROUPSIZE_TR_SHEAR;
	Shear_kernels[2].set_size(shearsize, WORKGROUPSIZE_TR_SHEAR);
	
	
	Update_SS_kernel[1].set_size(shearsize, WORKGROUPSIZE_TR_SHEAR);

	Save_shear_kernel.set_size(shearsize, WORKGROUPSIZE_TR_SHEAR);

	Sind.allocate(numbl_bounds);

	BLindicies.allocate(numbl_bounds);
	
	Weights.allocate(numbl_bounds);
	
	int curind = 0;

	for (int i = Bounds.MIN_BL_BOT; i <= Bounds.MAX_BL_BOT; i++)
	{
		BLindicies[curind++] = i;
	}

	for (int i = Bounds.MIN_BL_TOP; i <= Bounds.MAX_BL_TOP; i++)
	{
		BLindicies[curind++] = i;
	}
	////////////////////////////////////////////////////////////////////////////////////////

	// Sizes of arrays and kernels are not constant
	// and may need reallocation at update step
	Shear_array_len = vls.fullsize_Bnodes;
	
	Shear_inds.allocate(Shear_array_len);
	Shear_coeffs.allocate(Shear_array_len * 2);
	Tau.allocate(Shear_array_len);
	if (Tau.load("load" SLASH "Tau") == FALSE)
		Tau.fill(0);

	int shearnodesize = (int)ceil((double)vls.nBnodes / WORKGROUPSIZE_TR_SHEAR) * WORKGROUPSIZE_TR_SHEAR;
	Shear_kernels[0].set_size(shearnodesize, WORKGROUPSIZE_TR_SHEAR);
	Shear_kernels[1].set_size(shearnodesize, WORKGROUPSIZE_TR_SHEAR);

	Update_SS_kernel[0].set_size(shearnodesize, WORKGROUPSIZE_TR_SHEAR);
}

void clVariablesTR::createKernels()
{
	TR_Node_kernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_nodes_temp");
	TR_Node_kernel[0].set_size(nActiveNodes, WORKGROUPSIZE_TR);

	TR_Node_kernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_nodes_vel");
	TR_Node_kernel[1].set_size(nActiveNodes, WORKGROUPSIZE_TR);


	// The actual global size of any kernels below that are set with dummy 2*LocalSize
	// will be set before launching, (calling set_size w/ dummy arg to ensure dim is set)
	TR_Wall_Node_kernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_nodes_along_wall_temp");
	TR_Wall_Node_kernel[0].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	TR_Wall_Node_kernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_nodes_along_wall_vel");
	TR_Wall_Node_kernel[1].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	TR_ReRelease_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_release_par");
	TR_ReRelease_kernel.set_size(2*WORKGROUPSIZE_RERELEASE, WORKGROUPSIZE_RERELEASE);
	//Global size is set before enqueing kernel (calling set size w/ dummy arg to ensure dim is set)

	TR_Shear_Removal_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_shear_removal");
	TR_Shear_Removal_kernel.set_size(2 * WORKGROUPSIZE_RERELEASE, WORKGROUPSIZE_RERELEASE);
	//Global size is set before enqueing kernel

	TR_Update_Par_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_no_wall");
	TR_Update_Par_kernel.set_size(nActiveNodes, WORKGROUPSIZE_TR);

	TR_Wall_Par_kernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_along_wall");
	TR_Wall_Par_kernel[0].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);
	//Global size is set after num_wall_nodes is calculated

	TR_Wall_Par_kernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_reflect_particles");
	TR_Wall_Par_kernel[1].set_size(2*WORKGROUPSIZE_TR_WALL_REFLECT,WORKGROUPSIZE_TR_WALL_REFLECT);
	
	TR_Wall_Par_kernel[2].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_contact_wall");
	TR_Wall_Par_kernel[2].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	TR_Wall_Par_kernel[3].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_contact_wall2");
	TR_Wall_Par_kernel[3].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	TR_Wall_Par_kernel[4].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_deposit_particles_on_wall");
	TR_Wall_Par_kernel[4].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	Sort_kernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "Sort_merge_local");
	Sort_kernel[0].set_size(nN, WORKGROUPSIZE_SORT);

	Sort_kernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "Sort_merge_global");
	Sort_kernel[1].set_size(nN, WORKGROUPSIZE_SORT);

	Sort_kernel[2].create_kernel(GetSourceProgram, TRQUEUE_REF, "Sort_update_loc");
	Sort_kernel[2].set_size(nN, WORKGROUPSIZE_PLOC);

	Shear_kernels[0].create_kernel(GetSourceProgram, LBQUEUE_REF, "TR_shear_1");
	//Sizes set after Shear Variables are created

	Shear_kernels[1].create_kernel(GetSourceProgram, LBQUEUE_REF, "TR_shear_1");
	//Sizes set after Shear Variables are created

	Shear_kernels[2].create_kernel(GetSourceProgram, LBQUEUE_REF, "TR_shear_2");
	//Sizes set after Shear Variables are created

	Save_shear_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_shear_save_out");
	//Sizes set after Shear Variables are created

	Get_Umax_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "find_umax");

	Save_Dists_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_save_dists");

	Update_SS_kernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_Find_Shear_Coeffs1");
	Update_SS_kernel[0].set_size(2*WORKGROUPSIZE_TR_SHEAR, WORKGROUPSIZE_TR_SHEAR);

	Update_SS_kernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_Find_Shear_Coeffs2");
	Update_SS_kernel[1].set_size(2 * WORKGROUPSIZE_TR_SHEAR, WORKGROUPSIZE_TR_SHEAR);

	Clump_Particle_kernels[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "Rerelease_Clumps");
	Clump_Particle_kernels[0].set_size(2 * WORKGROUPSIZE_TR, WORKGROUPSIZE_TR);

	Clump_Particle_kernels[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "Test_Bounds_Clumped_Particles");
	Clump_Particle_kernels[1].set_size(2 * WORKGROUPSIZE_TR, WORKGROUPSIZE_TR);
}

void clVariablesTR::save_restart_files_ini()
{
	Active_indicies.save_bin_from_device("Active_indicies");
	Active_flag.save_bin_from_device("Active_flag");
	TR_indicies.save_bin_from_device("TR_indicies");
}

void clVariablesTR::saveRestartFiles()
{
	P.save_bin_from_device("trP");
	//RandList.save_bin_from_device("RandList");            /////Only when debugging
	Ploc.save_bin_from_device("Ploc");

	Array1Di TR_Load_Param;
	TR_Load_Param.zeros(6);
	TR_Load_Param(0) = nActiveNodes;
	TR_Load_Param(1) = Sort_timer;
	TR_Load_Param(2) = Bounds.MIN_BL_BOT;
	TR_Load_Param(3) = Bounds.MAX_BL_BOT;
	TR_Load_Param(4) = Bounds.MIN_BL_TOP;
	TR_Load_Param(5) = Bounds.MAX_BL_TOP;

	TR_Load_Param.save("TR_Params");

}

void clVariablesTR::save2file()
{

}

void clVariablesTR::saveDebug()
{
	NodI.save_txt_from_device("NodI");
	NodC.save_txt_from_device("NodC");
	NodV.save_txt_from_device("NodV");
	BL.save_txt_from_device("BL");
	BLind_ind.save_txt_from_device("BLind_ind");
	Sind.save_txt_from_device("Sind");
	BLindicies.save_txt_from_device("BLindicies");
	Weights.save_txt_from_device("Weights");
	Shear_inds.save_txt_from_device("Shear_inds");
	Shear_coeffs.save_txt_from_device("Shear_coeffs");
	Ploc.save_txt_from_device("Ploc");
	Node_neigh.save_txt_from_device("Node_neigh");
	Winds.save_txt_from_device("Winds");
	TR_indicies.save_txt_from_device("TR_indicies");
	Active_indicies.save_txt_from_device("Active_indicies");
	Active_flag.save_txt_from_device("Active_flag");
}

void clVariablesTR::ini_blinks()
{
	for (int i = 0; i < vls.nBL / 2; i++)
	{
		BL[i].P0ind = vls.BL(i, 0);
		BL[i].P1ind = vls.BL(i, 1);

		BL[i].vP0.x = vls.C[vls.BL(i, 0)].x;
		BL[i].vP0.y = vls.C[vls.BL(i, 0)].y;
		BL[i].vP1.x = vls.C[vls.BL(i, 1)].x;
		BL[i].vP1.y = vls.C[vls.BL(i, 1)].y;

		cl_double2 tanvec = Subtract2(BL[i].vP1, BL[i].vP0);
		cl_double2 normvec = { { -tanvec.y, tanvec.x } };
		BL[i].blLen = GETLEN(tanvec);
		BL[i].vNvec = Divide2(normvec, BL[i].blLen);
		BL[i].vTvec = Divide2(tanvec, BL[i].blLen);
		BL[i].Tau = 0.;
		BL[i].dir = 1;

		cl_double2 vP0not = { { vls.C0[vls.BL(i, 0)].x, vls.C0[vls.BL(i, 0)].y } };
		cl_double2 vP1not = { { vls.C0[vls.BL(i, 1)].x, vls.C0[vls.BL(i, 1)].y } };
		cl_double2 tanvecnot = Subtract2(vP1not, vP0not);
		tanvecnot = Divide2(tanvecnot, 2.);
		cl_double2 centernot = Subtract2(vP1not, tanvecnot);
		
		cl_double2 centpos = { { BL[i].vP0.x + tanvec.x * BL[i].blLen / 2,
			BL[i].vP0.y + tanvec.y * BL[i].blLen / 2 } };

		cl_double2 centerdist = Subtract2(centpos, centernot);

		if (GETLEN(centerdist) <= FOUL_SIZE_SWITCH_SOOT2)
		{
			BL[i].int_type = 0;
		}
		else
		{
			BL[i].int_type = 1;
		}
		int xind = (int)floor(centpos.x);
		int yind = (int)floor(centpos.y);
		
		if (xind >= 0 && xind < p.nX-1)
			BL[i].Node_loc = xind*(p.nY-1)+yind;
		else
			BL[i].Node_loc = -1;
	}

	for (int ii = vls.nBL - 1; ii >= vls.nBL / 2; ii--)
	{
		int i = (vls.nBL - 1 - ii) + vls.nBL / 2;
		BL[i].P0ind = vls.BL(ii, 1);
		BL[i].P1ind = vls.BL(ii, 0);

		BL[i].vP0.x = vls.C[vls.BL(ii, 1)].x;
		BL[i].vP0.y = vls.C[vls.BL(ii, 1)].y;
		BL[i].vP1.x = vls.C[vls.BL(ii, 0)].x;
		BL[i].vP1.y = vls.C[vls.BL(ii, 0)].y;

		cl_double2 tanvec = Subtract2(BL[i].vP1, BL[i].vP0);
		cl_double2 normvec = { { tanvec.y, -tanvec.x } };
		BL[i].blLen = GETLEN(tanvec);
		BL[i].vNvec = Divide2(normvec, BL[i].blLen);
		BL[i].vTvec = Divide2(tanvec, BL[i].blLen);
		BL[i].Tau = 0.;
		BL[i].dir = 1;
		cl_double2 centpos = { { BL[i].vP0.x + tanvec.x * BL[i].blLen / 2, BL[i].vP0.y + tanvec.y * BL[i].blLen / 2 } };
		int xind = (int)floor(centpos.x);
		int yind = (int)floor(centpos.y);
		if (xind >= 0 && xind < p.nX)
			BL[i].Node_loc = xind*(p.nY - 1) + yind;
		else
			BL[i].Node_loc = -1;
	}

	int xstart = (int)floor(X_MIN_VAL);
	double cur_x = (double)xstart + 0.5;
	int xstop = (int)ceil(X_STOP_POS);
	double end_x = (double)xstop + 0.5;
	for (int i = 0; i < vls.nBL / 2; i++)
	{
		cl_double2 P0t = vls.C[vls.BL(i, 0)], P1t = vls.C[vls.BL(i, 1)];
		if (P1t.x >(double)xstart)
		{
			Bounds.MIN_BL_BOT = i;
			Bounds.MIN_BL_TOP = vls.nBL / 2 + i;
			break;
		}
	}

	for (int i = Bounds.MIN_BL_BOT+1; i < vls.nBL / 2; i++)
	{
		cl_double2 P1t = vls.C[vls.BL(i, 1)];
		if (P1t.x >(double)xstop + 1.)
		{
			Bounds.MAX_BL_BOT = i;
			Bounds.MAX_BL_TOP = vls.nBL / 2 + i;
			break;
		}
	}

	for (int i = 0; i < vls.nBL / 2; i++)
	{
		int n0 = BL[i].P0ind, n1 = BL[i].P1ind;
		cl_double2 c0t = vls.C[n0], c1t = vls.C[n1];
		cl_double2 c0 = { { c0t.x - BL[i].vTvec.x*0.05, c0t.y - BL[i].vTvec.y*0.05 } };
		cl_double2 c1 = { { c1t.x + BL[i].vTvec.x*0.05, c1t.y + BL[i].vTvec.y*0.05 } };

		cl_int2 Cmin = vls.min2(c0, c1);
		cl_int2 Cmax = vls.max2(c0, c1);
		
		if (Cmin.x < 0)
			Cmin.x = 0;
		if (Cmax.x >= p.nX)
			Cmax.x = p.nX-1;

		if (Cmin.y < 0)
			Cmin.y = 0;

		if (Cmax.y >= p.nY)
			Cmax.y = p.nY - 1;


		if ((Cmax.x - Cmin.x) == 1 && (Cmax.y - Cmin.y) == 1)
		{
			int nind = Cmin.x * (p.nY - 1) + Cmin.y; 
			if (nind == -1)
			{
				printf("Error with BL number %d", i);
				exit(0);
			}
			if (BLind_ind[nind] == MAX_BL_PER_NODE)
			{
				printf("Too many BL's in node %d\n", i);
				exit(0);

			}
			Nod[nind].Wall_Flag = 1;
			Nod[nind].BLind[BLind_ind[nind]] = i;
			BLind_ind[nind]++;
		}
		else
		{
			for (int i0 = Cmin.x; i0 < Cmax.x; i0++)
			{
				for (int j0 = Cmin.y; j0 < Cmax.y; j0++)
				{
					test_node(i0, j0, c0, c1, i, 1);
				}
			}
		}
	}

	for (int i = vls.nBL/2; i < vls.nBL; i++)
	{
		int n0 = BL[i].P0ind, n1 = BL[i].P1ind;
		cl_double2 c0t = vls.C[n0], c1t = vls.C[n1];
		
		cl_double2 c0 = { { c0t.x - BL[i].vTvec.x*0.05, c0t.y - BL[i].vTvec.y*0.05 } };
		cl_double2 c1 = { { c1t.x + BL[i].vTvec.x*0.05, c1t.y + BL[i].vTvec.y*0.05 } };

		cl_int2 Cmin = vls.min2(c0, c1);
		cl_int2 Cmax = vls.max2(c0, c1);

		if (Cmin.x < 0)
			Cmin.x = 0;
		if (Cmax.x >= p.nX)
			Cmax.x = p.nX-1;

		if (Cmin.y < 0)
			Cmin.y = 0;

		if (Cmax.y >= p.nY)
			Cmax.y = p.nY - 1;


		if ((Cmax.x - Cmin.x) == 1 && (Cmax.y - Cmin.y) == 1)
		{
			int nind = (Cmin.x*(p.nY-1) + Cmin.y);
			if (nind == -1)
			{
				printf("Error with BL number %d", i);
				exit(0);
			}
			if (BLind_ind[nind] == MAX_BL_PER_NODE)
			{
				printf("Too many BL's in node %d\n", i);
				exit(0);

			}
			Nod[nind].Wall_Flag = 2;
			Nod[nind].BLind[BLind_ind[nind]] = i;
			BLind_ind[nind]++;
		}
		else
		{
			for (int i0 = Cmin.x; i0 < Cmax.x; i0++)
			{
				for (int j0 = Cmin.y; j0 < Cmax.y; j0++)
				{
					test_node(i0, j0, c0, c1, i, 2);
				}
			}
		}
	}
}

void clVariablesTR::call_update_wall_particles()
{
	//for (int i = 0; i < 9; i++)
	//{
	//	vtr.TR_Wall_Par_kernel[0].set_argument(13, sizeof(int), (void*)&i);
	//	vtr.TR_Wall_Par_kernel[0].call_kernel();
	//}

	//Reflect_inds.read_from_buffer(TRQUEUE, CL_TRUE);

	//if (Reflect_inds(0) > 0)
	//{
	//	TR_Wall_Par_kernel[1].set_argument(4, sizeof(cl_int), &Reflect_inds(0));
	//	TR_Wall_Par_kernel[1].calc_and_set_global_call_kernel(Reflect_inds(0));
	//	Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 0);
	//}

	//if (Reflect_inds(1) == 0)
	//	return;

	//while (TRUE)
	//{
	//	TR_Wall_Par_kernel[2].set_argument(8, sizeof(cl_int), &Reflect_inds(1));
	//	Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 1);
	//	TR_Wall_Par_kernel[2].calc_and_set_global_call_kernel(Reflect_inds(1));
	//	Reflect_inds.read_from_buffer(TRQUEUE, CL_TRUE);
	//	if (Reflect_inds(1) == 0)
	//		break;
	//	
	//	TR_Wall_Par_kernel[3].set_argument(6, sizeof(cl_int), &Reflect_inds(1));
	//	Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 1);
	//	TR_Wall_Par_kernel[3].calc_and_set_global_call_kernel(Reflect_inds(1));;
	//	Reflect_inds.read_from_buffer(TRQUEUE, CL_TRUE);
	//	if (Reflect_inds(1) == 0)
	//		break;
	//}

	//if (Reflect_inds(0) > 0)
	//{
	//	vtr.TR_Wall_Par_kernel[4].set_argument(3, sizeof(int), (void*)&Reflect_inds(0));
	//	vtr.TR_Wall_Par_kernel[4].calc_and_set_global_call_kernel(Reflect_inds(0));

	//}
	//Reflect_inds.FillBuffer(TRQUEUE, 0);
}

void clVariablesTR::call_update_wall_particles(cl_event *fill_evt)
{
	//int i = 0;
	//vtr.TR_Wall_Par_kernel[0].set_argument(13, sizeof(int), (void*)&i);
	//vtr.TR_Wall_Par_kernel[0].call_kernel(1, fill_evt);
	//for (i = 1; i < 9; i++)
	//{
	//	vtr.TR_Wall_Par_kernel[0].set_argument(13, sizeof(int), (void*)&i);
	//	vtr.TR_Wall_Par_kernel[0].call_kernel();
	//}

	//Reflect_inds.read_from_buffer(TRQUEUE, CL_TRUE);

	//if (Reflect_inds(0) > 0)
	//{
	//	TR_Wall_Par_kernel[1].set_argument(4, sizeof(cl_int), &Reflect_inds(0));
	//	TR_Wall_Par_kernel[1].calc_and_set_global_call_kernel(Reflect_inds(0));
	//	Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 0);
	//}

	//if (Reflect_inds(1) == 0)
	//	return;

	//while (TRUE)
	//{
	//	TR_Wall_Par_kernel[2].set_argument(8, sizeof(cl_int), &Reflect_inds(1));
	//	Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 1);
	//	TR_Wall_Par_kernel[2].calc_and_set_global_call_kernel(Reflect_inds(1));
	//	Reflect_inds.read_from_buffer(TRQUEUE, CL_TRUE);
	//	if (Reflect_inds(1) == 0)
	//		break;

	//	TR_Wall_Par_kernel[3].set_argument(6, sizeof(cl_int), &Reflect_inds(1));
	//	Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 1);
	//	TR_Wall_Par_kernel[3].calc_and_set_global_call_kernel(Reflect_inds(1));;
	//	Reflect_inds.read_from_buffer(TRQUEUE, CL_TRUE);
	//	if (Reflect_inds(1) == 0)
	//		break;
	//}

	//if (Reflect_inds(0) > 0)
	//{
	//	/*P.save_txt_from_device_full("trc", TRQUEUE);*/
	//	vtr.TR_Wall_Par_kernel[4].set_argument(3, sizeof(int), (void*)&Reflect_inds(0));
	//	vtr.TR_Wall_Par_kernel[4].calc_and_set_global_call_kernel(Reflect_inds(0));
	//}
	//Reflect_inds.FillBuffer(TRQUEUE, 0);
}

void clVariablesTR::setKernelArgs()
{
//	cl_int ind = 0;
//	TR_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Tlb_array.get_buf_add());
//	TR_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), NodV.get_buf_add());
//	TR_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), NodI.get_buf_add());
//	TR_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), NodC.get_buf_add());
//	TR_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), TR_indicies.get_buf_add());
//	TR_Node_kernel[0].set_argument(ind++, sizeof(cl_int), (void*)&nActiveNodes);
//
//	ind = 0;
//	TR_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), vlb.J_array.get_buf_add());
//	TR_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), NodV.get_buf_add());
//	TR_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), NodI.get_buf_add());
//	TR_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), NodC.get_buf_add());
//	TR_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), TR_indicies.get_buf_add());
//	TR_Node_kernel[1].set_argument(ind++, sizeof(cl_int), (void*)&nActiveNodes);
//
//	ind = 0;
//	TR_Wall_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Tlb_array.get_buf_add());
//	TR_Wall_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), NodC.get_buf_add());
//	TR_Wall_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), NodV.get_buf_add());
//	TR_Wall_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), TR_indicies.get_buf_add());
//	TR_Wall_Node_kernel[0].set_argument(ind++, sizeof(cl_mem), Winds.get_buf_add());
//	TR_Wall_Node_kernel[0].set_argument(ind++, sizeof(cl_int), (void*)&Num_wall_nodes);
//
//	ind = 0;
//	TR_Wall_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), vlb.J_array.get_buf_add());
//	TR_Wall_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), NodC.get_buf_add());
//	TR_Wall_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), NodV.get_buf_add());
//	TR_Wall_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), TR_indicies.get_buf_add());
//	TR_Wall_Node_kernel[1].set_argument(ind++, sizeof(cl_mem), Winds.get_buf_add());
//	TR_Wall_Node_kernel[1].set_argument(ind++, sizeof(cl_int), (void*)&Num_wall_nodes);
//
//	ind = 0; //re-release
//	TR_ReRelease_kernel.set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	TR_ReRelease_kernel.set_argument(ind++, sizeof(cl_mem), RandList.get_buf_add());
//	TR_ReRelease_kernel.set_argument(ind++, sizeof(Trparam), (void*)&trP);
//	TR_ReRelease_kernel.set_argument(ind++, sizeof(cl_mem), vlb.J_array.get_buf_add());
//	TR_ReRelease_kernel.set_argument(6, sizeof(cl_mem), Umax_val.get_buf_add());
//	TR_ReRelease_kernel.set_argument(7, sizeof(cl_mem), parP.get_buf_add());
//
//
//	ind = 0;//removal
//	TR_Shear_Removal_kernel.set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
//	TR_Shear_Removal_kernel.set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	TR_Shear_Removal_kernel.set_argument(ind++, sizeof(cl_mem), parP.get_buf_add());
//	TR_Shear_Removal_kernel.set_argument(ind++, sizeof(cl_mem), BL_dep.get_buf_add());
//
//	ind = 0; //Update Particles
//	TR_Update_Par_kernel.set_argument(ind++, sizeof(cl_mem), NodV.get_buf_add());
//	TR_Update_Par_kernel.set_argument(ind++, sizeof(cl_mem), Ploc.get_buf_add());
//	TR_Update_Par_kernel.set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	TR_Update_Par_kernel.set_argument(ind++, sizeof(cl_mem), parP.get_buf_add());
//	TR_Update_Par_kernel.set_argument(ind++, sizeof(cl_mem), Node_neigh.get_buf_add());
//	TR_Update_Par_kernel.set_argument(ind++, sizeof(cl_mem), Update_flag.get_buf_add());
//	TR_Update_Par_kernel.set_argument(ind++, sizeof(cl_mem), NodI.get_buf_add());
//	TR_Update_Par_kernel.set_argument(ind++, sizeof(cl_uint), (void*)&nActiveNodes);
//
//
//	ind = 0; //Update Particles
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), NodV.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), Ploc.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), PV.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), parP.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), Node_neigh.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), Update_flag.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), Winds.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_uint), (void*)&Num_wall_nodes);
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), NodI.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), Reflect_inds.get_buf_add());
//	TR_Wall_Par_kernel[0].set_argument(ind++, sizeof(cl_mem), Reflect_info.get_buf_add());
//	//index 13 for neigh_ind
//
//
//	ind = 0; //update Walls
//	TR_Wall_Par_kernel[1].set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
//	TR_Wall_Par_kernel[1].set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	TR_Wall_Par_kernel[1].set_argument(ind++, sizeof(cl_mem), PV.get_buf_add());
//	TR_Wall_Par_kernel[1].set_argument(ind++, sizeof(cl_mem), Reflect_info.get_buf_add());
//	//index 4 for number of pars
//
//	ind = 0; //update Walls
//	TR_Wall_Par_kernel[2].set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
//	TR_Wall_Par_kernel[2].set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	TR_Wall_Par_kernel[2].set_argument(ind++, sizeof(cl_mem), PV.get_buf_add());
//	TR_Wall_Par_kernel[2].set_argument(ind++, sizeof(cl_mem), parP.get_buf_add());
//	TR_Wall_Par_kernel[2].set_argument(ind++, sizeof(cl_mem), RandList.get_buf_add());
//	TR_Wall_Par_kernel[2].set_argument(ind++, sizeof(cl_mem), Reflect_info.get_buf_add());
//	TR_Wall_Par_kernel[2].set_argument(ind++, sizeof(cl_mem), Reflect_inds.get_buf_add());
//	TR_Wall_Par_kernel[2].set_argument(ind++, sizeof(cl_int), &nN);
//	//index 8 for num particles
//
//	ind = 0; //update Walls
//	TR_Wall_Par_kernel[3].set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
//	TR_Wall_Par_kernel[3].set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	TR_Wall_Par_kernel[3].set_argument(ind++, sizeof(cl_mem), PV.get_buf_add());
//	TR_Wall_Par_kernel[3].set_argument(ind++, sizeof(cl_mem), Reflect_info.get_buf_add());
//	TR_Wall_Par_kernel[3].set_argument(ind++, sizeof(cl_mem), Reflect_inds.get_buf_add());
//	TR_Wall_Par_kernel[3].set_argument(ind++, sizeof(cl_int), &nN);
//	//index 6 for num particles
//
//	ind = 0; //update Walls
//	TR_Wall_Par_kernel[4].set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
//	TR_Wall_Par_kernel[4].set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	TR_Wall_Par_kernel[4].set_argument(ind++, sizeof(cl_mem), Reflect_info.get_buf_add());
//	//index 3 for num particles
//
//	ind = 0;
//	Shear_kernels[0].set_argument(ind++, sizeof(cl_mem), vlb.Ro_array.get_buf_add());
//	Shear_kernels[0].set_argument(ind++, sizeof(cl_mem), vlb.J_array.get_buf_add());
//	Shear_kernels[0].set_argument(ind++, sizeof(cl_mem), Shear_coeffs.get_buf_add());
//	Shear_kernels[0].set_argument(ind++, sizeof(cl_mem), vlb.FA.get_buf_add());
//	Shear_kernels[0].set_argument(ind++, sizeof(cl_mem), Tau.get_buf_add());
//	Shear_kernels[0].set_argument(ind++, sizeof(cl_mem), Shear_inds.get_buf_add());
//	Shear_kernels[0].set_argument(ind++, sizeof(int), (void *)&vls.nBnodes);
//	Shear_kernels[0].set_argument(ind++, sizeof(double), (void *)&vlb.Fval);
//	Shear_kernels[0].set_argument(ind++, sizeof(cl_mem), vlb.Alpha_array.get_buf_add());
//
//	ind = 0;
//	Shear_kernels[1].set_argument(ind++, sizeof(cl_mem), vlb.Ro_array.get_buf_add());
//	Shear_kernels[1].set_argument(ind++, sizeof(cl_mem), vlb.J_array.get_buf_add());
//	Shear_kernels[1].set_argument(ind++, sizeof(cl_mem), Shear_coeffs.get_buf_add());
//	Shear_kernels[1].set_argument(ind++, sizeof(cl_mem), vlb.FB.get_buf_add());
//	Shear_kernels[1].set_argument(ind++, sizeof(cl_mem), Tau.get_buf_add());
//	Shear_kernels[1].set_argument(ind++, sizeof(cl_mem), Shear_inds.get_buf_add());
//	Shear_kernels[1].set_argument(ind++, sizeof(int), (void *)&vls.nBnodes);
//	Shear_kernels[1].set_argument(ind++, sizeof(double), (void *)&vlb.Fval);
//	Shear_kernels[1].set_argument(ind++, sizeof(cl_mem), vlb.Alpha_array.get_buf_add());
//
//
//	ind = 0;
//	Shear_kernels[2].set_argument(ind++, sizeof(cl_mem), BLindicies.get_buf_add());
//	Shear_kernels[2].set_argument(ind++, sizeof(cl_mem), Weights.get_buf_add());
//	Shear_kernels[2].set_argument(ind++, sizeof(cl_mem), Tau.get_buf_add());
//	Shear_kernels[2].set_argument(ind++, sizeof(cl_mem), Sind.get_buf_add());
//	Shear_kernels[2].set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
//	Shear_kernels[2].set_argument(ind++, sizeof(int), (void *)&numbl_bounds);
//
//
//	Get_Umax_kernel.set_argument(0, sizeof(cl_mem), vlb.J_array.get_buf_add());
//	Get_Umax_kernel.set_argument(1, sizeof(cl_mem), Umax_val.get_buf_add());
//	Get_Umax_kernel.set_argument(2, sizeof(cl_uint), (void*)&trP.Uvals_start);
//
//
//	Clump_Particle_kernels[0].set_argument(0, sizeof(cl_mem), P.get_buf_add());
//	Clump_Particle_kernels[0].set_argument(1, sizeof(cl_mem), vls.C0.get_buf_add());
//	Clump_Particle_kernels[0].set_argument(2, sizeof(Trparam), (void*)&trP);
//	Clump_Particle_kernels[0].set_argument(3, sizeof(cl_mem), RandList.get_buf_add());
//
//	Clump_Particle_kernels[1].set_argument(0, sizeof(cl_mem), P.get_buf_add());
//	Clump_Particle_kernels[1].set_argument(1, sizeof(cl_mem), BL.get_buf_add());
//	Clump_Particle_kernels[1].set_argument(2, sizeof(Trparam), (void*)&trP);
//	Clump_Particle_kernels[1].set_argument(3, sizeof(cl_uint), (void*)&vls.nBL);
//
//
//#ifdef WALL_SHEAR_STRESS
//	ind = 0;
//	Save_shear_kernel.set_argument(ind++, sizeof(cl_mem), BLindicies.get_buf_add());
//	Save_shear_kernel.set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
//	Save_shear_kernel.set_argument(ind++, sizeof(cl_mem), SS_output.get_buf_add());
//	Save_shear_kernel.set_argument(ind++, sizeof(int), (void *)&numbl_bounds);
//#endif
//
//
//#ifdef SAVE_IO_DISTS
//	ind = 0;
//	Save_Dists_kernel.set_argument(ind++, sizeof(cl_mem), Ploc.get_buf_add());
//	Save_Dists_kernel.set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
//	Save_Dists_kernel.set_argument(ind++, sizeof(cl_mem), IO_inds.get_buf_add());
//	Save_Dists_kernel.set_argument(ind++, sizeof(cl_mem), Node_neigh.get_buf_add());
//	Save_Dists_kernel.set_argument(ind++, sizeof(cl_mem), IO_dists_save.get_buf_add());
//	Save_Dists_kernel.set_argument(ind++, sizeof(cl_uint2), (void*)&IO_inds_info);
//	Save_Dists_kernel.set_argument(ind++, sizeof(cl_int), (void*)&Save_IO_Loc);
//	Save_Dists_kernel.set_argument(ind++, sizeof(cl_int), (void*)&Nd);
//#endif
}

void clVariablesTR::Check_Neigh_Node(int i, int j, int *cur_ind, int nod_loc)
{
	if(Active_flag(i, j)  == 1)
	{
		Node_neigh(nod_loc, *cur_ind) = i*(p.nY-1)+j;
		(*cur_ind) += 1;
	}
}

void clVariablesTR::Find_node_neighs()
{
	int xstart = (int)floor(X_MIN_VAL);
	int xstop = (int)ceil(X_STOP_POS);
	Node_neigh.allocate(nActiveNodes, 9);
	Node_neigh.fill(-1);
	Winds.zeros(nActiveNodes);
	int inletx = (int)ceil(X_release);
	int outletx = (int)floor(X_STOP_POS - 1);
	Array1Di Outlet_inds;
	Outlet_inds.zeros(p.nY - 1);
	Array1Di Inlet_inds;
	Inlet_inds.zeros(p.nY - 1);
	int In_ind = 0;
	int Out_ind = 0;
	int count = 0;


	for (int ii = 0; ii < nActiveNodes; ii++)
	{
		int i = Active_indicies(ii).x;
		int j = Active_indicies(ii).y;

		if (i == inletx)
		{
			Inlet_inds(In_ind++) = ii;
		}
		if (i == outletx)
		{
			Outlet_inds(Out_ind++) = ii;
		}


		if (NodI(i, j).Wall_Flag != 0)
		{
			Winds(count) = ii;
			count++;
		}
		int cur_ind = 1;

		Node_neigh(ii,0) = i*(p.nY-1)+j;
		Check_Neigh_Node(i + 1, j, &cur_ind, ii);
		Check_Neigh_Node(i - 1, j, &cur_ind, ii);
		Check_Neigh_Node(i, j + 1, &cur_ind, ii);
		Check_Neigh_Node(i, j - 1, &cur_ind, ii);
		Check_Neigh_Node(i + 1, j + 1, &cur_ind, ii);
		Check_Neigh_Node(i - 1, j - 1, &cur_ind, ii);
		Check_Neigh_Node(i + 1, j - 1, &cur_ind, ii);
		Check_Neigh_Node(i - 1, j + 1, &cur_ind, ii);
	}

	IO_inds.zeros(In_ind + Out_ind);
	IO_inds_info.x = In_ind;
	IO_inds_info.y = In_ind + Out_ind;
	
	for (int i = 0; i < In_ind; i++)
	{
		IO_inds(i) = Inlet_inds(i);
	}
	for (int i = 0; i < Out_ind; i++)
	{
		IO_inds(i+IO_inds_info.x) = Outlet_inds(i);
	}
	
	Num_wall_nodes = count;
	Num_W_nodes.allocate(1);
	Num_W_nodes(0) = count;
	Num_W_nodes.allocate_buffer_w_copy();

	double gsize = (double)Num_wall_nodes / WORKGROUPSIZE_TR_WALL;
	int wallglobal = (int)(ceil(gsize)*WORKGROUPSIZE_TR_WALL);
	
	TR_Wall_Par_kernel[0].reset_global_size(Num_wall_nodes);
	TR_Wall_Node_kernel[0].reset_global_size(Num_wall_nodes);
	TR_Wall_Node_kernel[1].reset_global_size(Num_wall_nodes);

	

	Winds.reallocate_host_only(wallglobal);
	Winds.reset_sizes(Num_wall_nodes, wallglobal);
	Num_wall_nodes_max = wallglobal;

	//int gsize_update = (int)ceil((double)IO_inds_info.y / (double)WORKGROUPSIZE_TR_DISTS) * WORKGROUPSIZE_TR_DISTS;

	Save_Dists_kernel.set_size(IO_inds_info.y, WORKGROUPSIZE_TR_DISTS);
	IO_dists_save.zeros(OUTPUT_MAX_LINES_IO * 2 * Nd);
	IO_dists_save.allocate_buffer_w_copy();
	IO_inds.allocate_buffer_w_copy();
	IO_inds.FreeHost();
	Save_IO_Loc = 0;
}

void clVariablesTR::update_Shear_arrays()
{
	//if (vls.fullsize_Bnodes_old != vls.fullsize_Bnodes)
	//{
	//	vls.fullsize_Bnodes = vls.lengthBnodes * 2;
	//	
	//	Shear_array_len = vls.fullsize_Bnodes;
	//	Shear_inds.reallocate_device_only(Shear_array_len );
	//	Shear_coeffs.reallocate_device_only(2 * vls.fullsize_Bnodes);
	//	Tau.reallocate_device_only(Shear_array_len);
	//	vls.Bnodes.reallocate(vls.fullsize_Bnodes);

	//	clFinish(IOQUEUE);

	//	Shear_kernels[1].set_argument(2, sizeof(cl_mem), Shear_coeffs.get_buf_add());
	//	Shear_kernels[1].set_argument(4, sizeof(cl_mem), Tau.get_buf_add());
	//	Shear_kernels[1].set_argument(5, sizeof(cl_mem), Shear_inds.get_buf_add());

	//	Shear_kernels[0].set_argument(2, sizeof(cl_mem), Shear_coeffs.get_buf_add());
	//	Shear_kernels[0].set_argument(4, sizeof(cl_mem), Tau.get_buf_add());
	//	Shear_kernels[0].set_argument(5, sizeof(cl_mem), Shear_inds.get_buf_add());

	//	Shear_kernels[2].set_argument(2, sizeof(cl_mem), Tau.get_buf_add());
	//	
	//	Update_SS_kernel[0].set_argument(0, sizeof(cl_mem), vls.Bnodes.get_buf_add());
	//	Update_SS_kernel[0].set_argument(4, sizeof(cl_mem), Shear_coeffs.get_buf_add());
	//	Update_SS_kernel[0].set_argument(5, sizeof(cl_mem), Shear_inds.get_buf_add());

	//	Update_SS_kernel[1].set_argument(4, sizeof(cl_mem), vls.Bnodes.get_buf_add());
	//}

	//vls.Bnodes.copy_to_buffer(CL_FALSE, vls.nBnodes);
	//
	///*int new_size = (int)ceil((double)vls.nBnodes / (double)WORKGROUPSIZE_TR_SHEAR)*WORKGROUPSIZE_TR_SHEAR;*/
	//
	//Bnode_top_start = vls.bot_ind_end + 1;
	//Shear_kernels[0].reset_global_size(vls.nBnodes);
	//Shear_kernels[1].reset_global_size(vls.nBnodes);
	//Update_SS_kernel[0].reset_global_size(vls.nBnodes);

	//Update_SS_kernel[0].set_argument(7, sizeof(cl_int), (void*)&vls.num_el);
	//Update_SS_kernel[0].set_argument(9, sizeof(int), (void *)&vls.nBnodes);
	//Update_SS_kernel[0].set_argument(10, sizeof(int), (void*)&Bnode_top_start);

	//cl_int2 bindicies_end = { { Bnode_top_start, vls.nBnodes } };
	//Update_SS_kernel[1].set_argument(10, sizeof(cl_int2), (void*)&bindicies_end);

	//cl_event UpdateSSevt;

	//Update_SS_kernel[0].call_kernel();
	//Update_SS_kernel[1].call_kernel(&UpdateSSevt);
	//
	//Bindicies_loc.FillBuffer(TRQUEUE, { { -1, -1 } }, 1, &UpdateSSevt);

	//clReleaseEvent(UpdateSSevt);

	//Shear_kernels[0].set_argument(6, sizeof(int), (void*)&vls.nBnodes);
	//Shear_kernels[1].set_argument(6, sizeof(int), (void*)&vls.nBnodes);
}

double clVariablesTR::weight_kernel_func(double Sval)
{
	return (1. / sqrt(2.*p.Pi) * exp(-0.5*Sval*Sval));
	//if (Sval <= 1.)
	//{
	//	double tempval = 1. - Sval*Sval;
	//	return 35. / 32. * tempval * tempval * tempval;
	//}
	//return 0.;
}

par clVariablesTR::AddPars(par p1, par p2)
{
	p1.Num_rep += p2.Num_rep;
	p1.pos = Add2(p1.pos, p2.pos);
	int psum = p1.timer + p2.timer;
	p1.timer = MIN(psum, 20000);
	return p1;
}

bool myfunction(cl_int2 i, cl_int2 j) { return (i.x < j.x); }


void clVariablesTR::Clump_Particles()
{
//	sort_particles_for_clumping();
//	/////////////Fill vector with values corresponing to particle locations/////////
//
//	double size_multiplier = 15.;
//	int tr_size_y = (p.nY - 1) * (int)size_multiplier;
//	P.read_from_buffer();
//
//	std::vector<cl_int2> ParLocVec(nN, { { 0, 0 } });
//	int j;
//#pragma omp parallel for 
//	for (j = 0; j < nN; j++)
//	{
//		ParLocVec[j].y = j;
//		cl_double2 postemp = P(j).pos;
//		postemp = { { postemp.x * size_multiplier, postemp.y * size_multiplier } };
//		cl_int2 Posi = { { (int)floor(postemp.x), (int)floor(postemp.y) } };
//
//		int pcur_loc = Posi.x*tr_size_y + Posi.y;
//
//		if (P(j).loc == -2)
//		{
//			ParLocVec[j].x = -2;
//		}
//		else if (P(j).loc == -1)
//		{
//			ParLocVec[j].x = -1;
//		}
//		else
//		{
//			ParLocVec[j].x = pcur_loc;
//		}
//	}
//
//
//	//Sort vector
//	std::sort(ParLocVec.begin(), ParLocVec.end(), myfunction);
//
//
//	///Fill beginning of temporary particle array with particles about to be re-released
//	///and particles currently deposited on wall
//	int cur_ind = 0;
//
//	while (TRUE)
//	{
//		SortTmp(cur_ind) = P(ParLocVec[cur_ind].y);
//		cur_ind++;
//		if (ParLocVec[cur_ind].x > -1)
//			break;
//	}
//
//	////Combine same type particles located in each sub-box
//	int cur_red_ind = cur_ind;
//	int num_dep_particles = cur_ind;
//	while (TRUE)
//	{
//		int start_ind = cur_ind;
//		int cur_loc = ParLocVec[cur_ind].x;
//		int stop_ind = cur_ind;
//
//		//find index of last particle located in same box.
//		while (TRUE)
//		{
//			if (stop_ind == nN - 1)
//			{
//				break;
//			}
//			if (ParLocVec[stop_ind + 1].x != cur_loc)
//			{
//				break;
//			}
//			stop_ind++;
//		}
//
//		for (int i = 0; i < Nd; i++)
//		{
//			int count_pars = 0;
//			par Pcur;
//			Pcur.Dep_Flag = -1;
//			Pcur.Dep_timer = TIME_BEFORE_STICK;
//			Pcur.loc = P(ParLocVec[start_ind].y).loc;
//			Pcur.Num_rep = 0;
//			Pcur.pos = { { 0., 0. } };
//			Pcur.timer = 0;
//			Pcur.type = i;
//			for (int j = start_ind; j <= stop_ind; j++)
//			{
//				par Pcurtemp = P(ParLocVec[j].y);
//				if (Pcurtemp.type == i)
//				{
//					count_pars++;
//					Pcur = AddPars(Pcur, Pcurtemp);
//				}
//			}
//			if (count_pars)
//			{
//				//Pcur.Num_rep /= (double)count_pars;
//				Pcur.pos = Divide2(Pcur.pos, (double)count_pars);
//				Pcur.timer /= (double)count_pars;
//				cl_int2 Posi = { { (int)floor(Pcur.pos.x), (int)floor(Pcur.pos.y) } };
//				Pcur.loc = Posi.x*(p.nY-1) + Posi.y;
//				SortTmp(cur_red_ind++) = Pcur;
//			}
//		}
//
//		cur_ind = stop_ind + 1;
//		if (cur_ind == nN)
//			break;
//	}
//
//
//
//	P.write_array_to_buffer(CL_TRUE, SortTmp.get_array(), cur_red_ind-1);
//
//	Clump_Particle_kernels[1].set_argument(2, sizeof(Trparam), (void*)&trP);
//	Clump_Particle_kernels[1].set_argument(4, sizeof(int), (void*)&cur_red_ind);
//	Clump_Particle_kernels[1].set_argument(5, sizeof(int), (void*)&num_dep_particles);
//	
//	Clump_Particle_kernels[0].set_argument(2, sizeof(Trparam), (void*)&trP);
//	Clump_Particle_kernels[0].set_argument(4, sizeof(int), (void*)&cur_red_ind);
//
//
//
//	Clump_Particle_kernels[0].set_global_call_kernel(nN - cur_red_ind);
//
//	Clump_Particle_kernels[1].set_global_call_kernel(cur_red_ind - num_dep_particles);
//
//	clFinish(TRQUEUE);
//
//	initial_sort_particles();
//
}


void clVariablesTR::allocate_buffers()
{
	P.allocate_buffer_w_copy();
	Active_flag.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	NodI.allocate_buffer_w_copy();
	NodV.allocate_buffer_w_copy();
	NodC.allocate_buffer_w_copy();
	parP.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	RandList.allocate_buffer_w_copy();
	TR_indicies.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	Active_indicies.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	Tau.allocate_buffer_w_copy();
	Sind.allocate_buffer_w_copy();
	Weights.allocate_buffer_w_copy();
	Shear_coeffs.allocate_buffer_w_copy();
	Shear_inds.allocate_buffer_w_copy();
	BLindicies.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	BL.allocate_buffer_w_copy();
	Winds.allocate_buffer_w_copy();
	BL_dep.allocate_buffer_w_copy();
	Ploc.allocate_buffer_w_copy();
	Node_neigh.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	SS_output.zeros(OUTPUT_MAX_LINES_SS, numbl_bounds);
	SS_output.allocate_buffer_w_copy();
	BLind_ind.allocate_buffer_w_copy();
#ifdef USE_OPENGL
	ini_particle_colors();
	Par_Color_List.allocate_buffer_w_copy();
	Par_Color_List.FreeHost();
#endif


////////Allocation of Buffers only used for updates on Device////////////
	Bindicies_loc.zeros(p.nX);
	Bindicies_loc.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
	Bindicies_loc.FillBuffer({ { -1, -1 } });
	Bindicies_loc.FreeHost();

	Update_flag.zeros(nN);
	Update_flag.allocate_buffer_w_copy();
	Update_flag.FreeHost();

	
	SortTmp.allocate(nN);
	SortTmp.allocate_buffer();
	//SortTmp.FreeHost();

	PV.allocate(nN*2);
	PV.allocate_buffer();
	cl_double2 PVt = { { 0., 0. } };
	PV.FillBuffer(PVt);

	Reflect_info.zeros(nN * 2);
	Reflect_info.allocate_buffer_w_copy();
	Reflect_info.FreeHost();

	Reflect_inds.zeros(2);
	Reflect_inds.allocate_buffer_w_copy();
///////////////////////////////////////////////////////////////////////
}

void clVariablesTR::ini_shear_coeffs()
{
	//clFinish(IOQUEUE);

	//int ind = 0;
	//Bnode_top_start = vls.bot_ind_end + 1;
	//int Blinks_top_ind_t = BL.getSizeX() / 2;
	//double cutrad = CUTOFF_RADIUS;
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.Bnodes.get_buf_add());
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Stor.get_buf_add());
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.dir_array.get_buf_add());
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.ii0_array.get_buf_add());
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_mem), Shear_coeffs.get_buf_add());
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_mem), Shear_inds.get_buf_add());
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_mem), Bindicies_loc.get_buf_add());
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_int), (void*)&vls.num_el);
	//Update_SS_kernel[0].set_argument(ind++, sizeof(cl_double), (void*)&cutrad);
	//Update_SS_kernel[0].set_argument(ind++, sizeof(int), (void *)&vls.nBnodes);
	//Update_SS_kernel[0].set_argument(ind++, sizeof(int), (void*)&Bnode_top_start);

	//cutrad = 2.5;
	//int bind_rad = INDEX_RADIUS_SEARCH;
	//cl_int2 bindicies_end = { { Bnode_top_start, vls.nBnodes } };
	//ind = 0;
	//Update_SS_kernel[1].set_argument(ind++, sizeof(cl_mem), BLindicies.get_buf_add());
	//Update_SS_kernel[1].set_argument(ind++, sizeof(cl_mem), Weights.get_buf_add());
	//Update_SS_kernel[1].set_argument(ind++, sizeof(cl_mem), Sind.get_buf_add());
	//Update_SS_kernel[1].set_argument(ind++, sizeof(cl_mem), BL.get_buf_add());
	//Update_SS_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.Bnodes.get_buf_add());
	//Update_SS_kernel[1].set_argument(ind++, sizeof(cl_mem), Bindicies_loc.get_buf_add());
	//Update_SS_kernel[1].set_argument(ind++, sizeof(int), (void *)&numbl_bounds);
	//Update_SS_kernel[1].set_argument(ind++, sizeof(cl_double), (void*)&cutrad);
	//Update_SS_kernel[1].set_argument(ind++, sizeof(int), (void *)&Blinks_top_ind_t);
	//Update_SS_kernel[1].set_argument(ind++, sizeof(int), (void*)&bind_rad);
	//Update_SS_kernel[1].set_argument(ind++, sizeof(cl_int2), (void*)&bindicies_end);


	//Update_SS_kernel[0].call_kernel();

	//Update_SS_kernel[1].call_kernel();

	//Bindicies_loc.FillBuffer(TRQUEUE, { { -1, -1 } });

}

void clVariablesTR::freeHostMem()
{
	BLind_ind.FreeHost();
	//NodI.FreeHost();
	//NodC.FreeHost();
	//NodV.FreeHost();
	//RandList.FreeHost();
	//BLind_ind.FreeHost();
	//P.FreeHost();
	//BL_dep.FreeHost();
	//PV.FreeHost();
	//Tau.FreeHost();
	//SortTmp.FreeHost();
	//Umax_val.FreeHost();
}

void clVariablesTR::Re_Release_Par()
{
	Get_Umax_kernel.call_kernel();
	//Umax_val.read_from_buffer(IOQUEUE, CL_TRUE);
	//printf("Uval = %g\n", Umax_val(0));
	cl_uint maxval = Ploc(0).y;
	cl_uint par_in = MAX(Par_conc_inlet*Umean_Current / maxval,1);
	TR_ReRelease_kernel.set_argument(4, &maxval);
	TR_ReRelease_kernel.set_argument(5, &par_in);
	TR_ReRelease_kernel.set_global_call_kernel(maxval);
}

void clVariablesTR::Update_Par_Rem_Args()
{
	cl_uint offset = Ploc(1).x;
	cl_uint total_removal = Ploc(1).y - offset;


	TR_Shear_Removal_kernel.set_argument(4, &offset);
	TR_Shear_Removal_kernel.set_argument(5, &total_removal);

	TR_Shear_Removal_kernel.reset_global_size(total_removal);
}

void clVariablesTR::Save_SS()
{
	//Tau.save_from_device("Tau");
	//vls.dir_array.savetxt("dir_array");
	//vls.ii0_array.savetxt("ii0_array");
	//vls.Bnodes.savetxt("Bnodes");
	//BL.savetxt("BLtemp");
	//Sind.save_from_device("Sind");
	//Weights.save_from_device("Weights");
	//BLindicies.save_from_device("Blindicies");
	//Shear_coeffs.save_from_device("SScoeff");
	//Shear_inds.save_txt_from_device("Sinds");
	//vlb.J_array.save_txt_from_device("J_array");
	//vlb.Ro_array.save_txt_from_device("Ro_array");
	//vlb.FA.save_txt_from_device("FA");
	//vlb.FB.save_txt_from_device("FB");

}

void clVariablesTR::test_node(int i, int j, cl_double2 vL0, cl_double2 vL1, int blind, int bot_flag)
{
	cl_double2 vLd = Subtract2(vL1, vL0);
	int nind = i*(p.nY-1)+j;
	if (BLind_ind[nind] == MAX_BL_PER_NODE)
	{
		printf("Too many BL's in node %d\n", i);
		exit(0);

	}
	cl_double2 CC00 = { { (double)i - 0.05, (double)j - 0.05 } };
	cl_double2 CC10 = { { (double)i + 1.05, (double)j - 0.05 } };
	cl_double2 CC01 = { { (double)i - 0.05, (double)j + 1.05 } };
	cl_double2 CC11 = { { (double)i + 1.05, (double)j + 1.05 } };
	if (test_cross(vL0, vLd, CC00, CC10) || test_cross(vL0, vLd, CC00, CC01) || 
		test_cross(vL0, vLd, CC10, CC11) || test_cross(vL0, vLd, CC01, CC11))
	{
		Nod[nind].Wall_Flag = bot_flag;
		Nod[nind].BLind[BLind_ind[nind]] = blind;
		BLind_ind[nind]++;
		return;
	}
	
}

BOOL clVariablesTR::test_cross(cl_double2 vL0, cl_double2 vLd, cl_double2 vP0, cl_double2 vP1)
{
	cl_double2 vP10 = Subtract2(vP1,vP0);
	cl_double2 vNm = { { -1.*vP10.y, vP10.x } };
	double den = DOT_PROD(vNm, vLd);
	if (den == 0.)
		return FALSE;
	cl_double2 vPL0 = Subtract2(vP0, vL0);
	double dist = DOT_PROD(vNm, vPL0) / den;
	cl_double2 vC = { { vL0.x + vLd.x * dist, vL0.y + vLd.y * dist } };
	cl_double2 vd0 = Subtract2(vC, vP0), vd1 = Subtract2(vC,vP1);
	if (dist >= 0. && dist <= 1.)
	if (fabs(GETLEN(vP10) - GETLEN(vd0) - GETLEN(vd1)) < p.eps)
			return TRUE;
	return FALSE;
}

void clVariablesTR::ini_rand()
{
	for (int i = 0; i < nN; i++)
	{
		RandList[i].x = (cl_uint)rand();
		RandList[i].y = (cl_uint)rand();
	}
}

void clVariablesTR::ini_trp()
{
//	Umax_val.zeros(1);
//	Umax_val.allocate_buffer_w_copy();
//
//	for (int i = 0; i < vls.nBL; i++)
//	{
//		if (vls.C[vls.BL(i, 0)].x <= X_release && vls.C[vls.BL(i, 1)].x > X_release)
//		{
//			trP.BL_rel_bot = i;
//			break;
//		}
//	}
//
//	for (int i = (int)(vls.nBL / 2) + 1; i < vls.nBL; i++)
//	{
//		if (vls.C[vls.BL(i, 1)].x <= X_release && vls.C[vls.BL(i, 0)].x > X_release)
//		{
//			trP.BL_rel_top = i;
//			break;
//		}
//	}
//
//	trP.X_release = X_release;
//
//	int y0 = 0;
//	int x0 = (int)trP.X_release - 1;
//	while (vls.M(x0, y0) == LB_SOLID)
//		y0++;
//
//	int yred = vlb.Stor(x0, y0);
//	trP.Uvals_start = (vlb.Channel_Height + 1)*x0 + yred;
//
//	Num_nodes_at_release = 0;
//
//	while (vls.M(x0, y0) != LB_SOLID)
//	{
//		y0++;
//		Num_nodes_at_release++;
//	}
//#ifdef OPENCL_VERSION_1_2
//	Get_Umax_kernel.set_size_1D(1, 1);
//#else
//	Get_Umax_kernel.set_size_1D(Num_nodes_at_release, Num_nodes_at_release);
//#endif
//
//
//	trP.umax_val = 3. / 2.*RE_NUMBER*MU_NUMBER / p.Pipe_radius;
//	double Xtop_val = vls.C[vls.BL(trP.BL_rel_top, 1)].y + (X_release - vls.C[vls.BL(trP.BL_rel_top, 1)].x) *
//		(vls.C[vls.BL(trP.BL_rel_top, 0)].y - vls.C[vls.BL(trP.BL_rel_top, 1)].y) /
//		(vls.C[vls.BL(trP.BL_rel_top, 0)].x - vls.C[vls.BL(trP.BL_rel_top, 1)].x);
//	trP.offset_y = vls.C[vls.BL(trP.BL_rel_bot, 1)].y + (X_release - vls.C[vls.BL(trP.BL_rel_bot, 1)].x) *
//		(vls.C[vls.BL(trP.BL_rel_bot, 0)].y - vls.C[vls.BL(trP.BL_rel_bot, 1)].y) /
//		(vls.C[vls.BL(trP.BL_rel_bot, 0)].x - vls.C[vls.BL(trP.BL_rel_bot, 1)].x);
//
//	trP.bval = Xtop_val - trP.offset_y;
//
//	trP.Bottom_location = ceil(trP.offset_y) - trP.offset_y;
//	trP.Top_location = floor(Xtop_val) - trP.offset_y;
//
//

}

void clVariablesTR::update_trp()
{
//	double Xtop_val = vls.C[vls.BL(trP.BL_rel_top, 1)].y + (X_release - vls.C[vls.BL(trP.BL_rel_top, 1)].x) *
//		(vls.C[vls.BL(trP.BL_rel_top, 0)].y - vls.C[vls.BL(trP.BL_rel_top, 1)].y) /
//		(vls.C[vls.BL(trP.BL_rel_top, 0)].x - vls.C[vls.BL(trP.BL_rel_top, 1)].x);
//	trP.offset_y = vls.C[vls.BL(trP.BL_rel_bot, 1)].y + (X_release - vls.C[vls.BL(trP.BL_rel_bot, 1)].x) *
//		(vls.C[vls.BL(trP.BL_rel_bot, 0)].y - vls.C[vls.BL(trP.BL_rel_bot, 1)].y) /
//		(vls.C[vls.BL(trP.BL_rel_bot, 0)].x - vls.C[vls.BL(trP.BL_rel_bot, 1)].x);
//	trP.bval = Xtop_val - trP.offset_y;
//	trP.Bottom_location = ceil(trP.offset_y) - trP.offset_y;
//	trP.Top_location = floor(Xtop_val) - trP.offset_y;
//	
//	int y0 = 0;
//	int x0 = (int)trP.X_release - 1;
//	while (vls.M(x0, y0) == LB_SOLID)
//		y0++;
//
//	int yred = vlb.Stor(x0, y0);
//	trP.Uvals_start = (vlb.Channel_Height + 1)*x0 + yred;
//
//	Num_nodes_at_release = 0;
//
//	while (vls.M(x0, y0) != LB_SOLID)
//	{
//		y0++;
//		Num_nodes_at_release++;
//	}
//
//#ifdef OPENCL_VERSION_1_2
//	Get_Umax_kernel.set_size_1D(1, 1);
//#else
//	Get_Umax_kernel.set_size_1D(Num_nodes_at_release, Num_nodes_at_release);
//#endif
//	Get_Umax_kernel.set_argument(2, sizeof(cl_uint), (void*)&trP.Uvals_start);
//		
//	TR_ReRelease_kernel.set_argument(2, sizeof(Trparam), (void *)&trP);
}

void clVariablesTR::save_box(double x1, double y1, double dx, double dy)
{
	double xlow = x1;
	double xhigh = x1 + dx;
	if (dx < 0.)
	{
		xlow = xhigh;
		xhigh = x1;
	}

	double ylow = y1;
	double yhigh = y1 + dy;
	if (dy < 0.)
	{
		ylow = yhigh;
		yhigh = y1;
	}

	cl_double2 Xbound = { { xlow, xhigh } };
	cl_double2 Ybound = { { ylow, yhigh } };

	int	psize = 100;
	Array1DP Pout;
	Array1Dv2d PVout;
	Pout.allocate(psize);
	PVout.allocate(psize);
	P.read_from_buffer();
	PV.read_from_buffer();
	int ind = 0;

	for (int i = 0; i < vtr.nN; i++)
	{
		cl_double2 pos = P(i).pos;
		if (pos.x >= Xbound.x && pos.x <= Xbound.y && pos.y >= Ybound.x && pos.y <= Ybound.y)
		{
			PVout(ind) = PV(i);
			Pout(ind++) = P(i);

			if (ind == psize)
			{
				psize *= 2;
				Pout.reallocate_host_only(psize);
				PV.reallocate_host_only(psize);
			}
		}
	}
	psize = ind;
	if (psize > 0)
	{
		Pout.reallocate_host_only(psize);
		PVout.reallocate_host_only(psize);
		Pout.savetxt_full("trp_out");
		PVout.savetxt("trv_out");
	}
	
	int blsize = (int)ceil(dx) * 1.5;

	P.FreeHost();
	PV.FreeHost();



	Array1DBL BLout;
	BLout.allocate(blsize);

	BL.read_from_buffer();
	ind = 0;
	for (int i = 0; i < vls.nBL; i++)
	{
		cl_double2 pos1 = BL(i).vP0, pos2 = BL(i).vP1;
		if (pos2.y >= Ybound.x && pos1.y <= Ybound.y && pos2.x >= Xbound.x && pos1.x <= Xbound.y)
		{
			BLout(ind++) = BL(i);
			if (ind == blsize)
			{
				blsize *= 1.5;
				BLout.reallocate_host_only(blsize);
			}
		}
	}

	blsize = ind;
	if (blsize)
	{
		BLout.reallocate_host_only(blsize);
		BLout.savetxt("trbl_out");
	}
	
	BL.FreeHost();
	Xbound.x = MAX(Xbound.x, 0.);
	Xbound.y = MIN(Xbound.y, (double)p.nX - 2.);
	Ybound.x = MAX(Ybound.x, 0.);
	Ybound.y = MIN(Ybound.y, (double)p.nY - 3.);
	int sizex = (int)ceil(Xbound.y) - (int)floor(Xbound.x) + 1;
	int sizey = (int)ceil(Ybound.y) - (int)floor(Ybound.x) + 1;

	Array2DNC NodCout;
	Array2DNI NodIout;
	Array2Di BLind_out;
	Array2DNV NodVout;

	NodCout.allocate(sizex, sizey);
	NodIout.allocate(sizex, sizey);
	NodVout.allocate(sizex, sizey);
	BLind_out.allocate(sizex, sizey);

	NodC.read_from_buffer();
	NodI.read_from_buffer();
	NodV.read_from_buffer();
	BLind_ind.read_from_buffer();

	int i_ind = 0;
	int j_ind = 0;
	for (int i = (int)floor(Xbound.x); i <= (int)ceil(Xbound.y); i++)
	{
		j_ind = 0;
		for (int j = (int)floor(Ybound.x); j <= (int)ceil(Ybound.y); j++)
		{
			NodCout(i_ind, j_ind) = NodC(i, j);
			NodIout(i_ind, j_ind) = NodI(i, j);
			NodVout(i_ind, j_ind) = NodV(i, j);
			BLind_out(i_ind, j_ind) = BLind_ind(i, j);
			j_ind++;
		}
		i_ind++;
	}

	NodIout.savetxt("trni_out");
	NodCout.savetxt("trnc_out");
	NodVout.savetxt("trnv_out");
	BLind_out.savetxt("trblind_out");

	//vfd.Temp.save_txt_from_device("fdt_out");

}

void clVariablesTR::ini_particle_colors()
{
	Par_Color_List.zeros(Nd * 3);

	if (Nd == 1)
	{
		Par_Color_List(0) = 1.f;
		return;
	}
	if (Nd == 2)
	{
		Par_Color_List(0) = 1.f;
		Par_Color_List(5) = 1.f;
		return;
	}
	if (Nd == 3)
	{
		Par_Color_List(0) = 1.f;
		Par_Color_List(5) = 1.f;
		Par_Color_List(7) = 1.f;
		return;
	}

	int P1 = Nd / 4;
	int P2 = Nd / 2;
	int P3 = Nd / 4 * 3;
	int P4 = Nd;

	float stepval = 1.f / ((float)P1);
	for (int i = 0; i < P1; i++)
	{
		Par_Color_List(i * 3 + 1) = stepval*(float)i;
		Par_Color_List(i * 3 + 2) = 1.f;
	}

	stepval = 1.f / (float)(P2 - P1);
	for (int i = P1; i < P2; i++)
	{
		Par_Color_List(i * 3 + 2) = 1.f - stepval*(float)(i - P1);
		Par_Color_List(i * 3 + 1) = 1.f;
	}

	stepval = 1.f / (float)(P3 - P2);
	for (int i = P2; i < P3; i++)
	{
		Par_Color_List(i * 3) = stepval*(float)(i - P2);
		Par_Color_List(i * 3 + 1) = 1.f;
	}

	stepval = 1.f / (float)(P4 - P3);
	for (int i = P3; i < P4; i++)
	{
		Par_Color_List(i * 3 + 1) = 1.f - stepval*(float)(i - P3);
		Par_Color_List(i * 3) = 1.f;
	}
}

// TODO: Figure out what DELTA_P represents and fix if incorrect

void clVariablesTR::ini_particle_properties()
{
	nN = TRC_NUM_TRACERS;
	Array1Dd Par_dist;

	if (Par_dist.loadtxtnew("load" SLASH "pdist") == FALSE)
	{
		printf("error loading particle distrbution\n");
		exit(0);
	}

	
	Nd = Par_dist.getSizeX() / 2;

	D_p_real.zeros(Nd);
	D_dists.zeros(Nd);

	for (int i = 0; i < Nd; i++)
	{
		D_p_real(i) = Par_dist(i * 2);
		D_dists(i) = Par_dist(i * 2 + 1);
	}



	double Num_per_m2 = SOOT_NUMBER_CONCENTRATION / p.DELTA_L / p.DELTA_L;
	double Num_per_m = Num_per_m2 * 2. * p.Pipe_radius;
	double Umean = vlb.Re * (vlb.tau - 0.5) / 3. / p.Pipe_radius;
	double Num_per_s = Num_per_m;

	Par_conc_inlet = (cl_uint)Num_per_s * (double)NUM_STEPS_BTW_SORT; //number of particles released after each sort step 
																		//must be multiplied by Umean to get true value
	
	WofA = 2.*sqrt(SURF_ENERGY_SURF * SURF_ENERGY_SOOT);
	Estar = 1. / ((1. - POISSON_SURF*POISSON_SURF) / Y_MOD_SURF + (1. - POISSON_SOOT*POISSON_SOOT) / Y_MOD_SOOT);
	Estar_s = Y_MOD_SOOT;
	WofA_s = SURF_ENERGY_SOOT;


	double MFP = MEAN_FREE_PATH_AIR;

	double ConstA = 1.2, ConstB = 0.4, ConstC = 1.1;
	double ConstCs = 1.14, ConstCm = 1.17, ConstCt = 2.18;

	double Kparticle = THERMAL_CONDUCTIVITY_PARTICLE;


	Kth_pars.zeros(Nd);
	D_p.zeros(Nd);
	Q_A_prime.zeros(Nd, 2);
	Q_A.zeros(Nd, 2);
	M_p.zeros(Nd);
	F_po.zeros(Nd);
	R_d.zeros(Nd, 2);
	Tau_crit.zeros(Nd, 2);
	V_par.zeros(Nd);
	Tau_crit_max = 0.;

	for (int i = 0; i < Nd; i++)
	{
		double Knud = MFP / D_p_real(i);

		double CCF = 1. + Knud * (ConstA + ConstB * exp(-ConstC / Knud));
		double k_ratio = THERMAL_CONDUCTIVITY_AIR / THERMAL_CONDUCTIVITY_PARTICLE;
		double Kth1 = 2. * ConstCs * CCF / (1. + 3.*ConstCm * Knud);
		double Kth2 = (k_ratio + ConstCt * Knud) / (1. + 2. * k_ratio + 2.*ConstCt*Knud);

		Kth_pars(i) = Kth1*Kth2;

		double a3 = 9.*PI_NUMBER*WofA*D_p_real(i)*D_p_real(i) / Estar / 4.;
		double Rdi = pow(a3, (1. / 3.));
		Rdi *= p.DELTA_L;
		R_d(i, 0) = Rdi;

		double a3_s = 9.*PI_NUMBER*WofA_s*D_p_real(i)*D_p_real(i) / Estar_s / 4.;
		double Rdi_s = pow(a3_s, (1. / 3.));
		Rdi_s *= p.DELTA_L;
		R_d(i, 1) = Rdi_s;

		double QAi = 2.*WofA*PI_NUMBER*R_d(i, 0)*R_d(i, 0);
		QAi *= (p.DELTA_F / p.DELTA_L);
		Q_A(i, 0) = QAi;

		double QAi_s = 2.*WofA_s*PI_NUMBER*R_d(i, 1)*R_d(i, 1);
		QAi *= (p.DELTA_F / p.DELTA_L);
		Q_A(i, 1) = QAi_s;

		double Q_A_pt = pow(WofA, 5.)*pow(D_p_real(i), 4.) / Estar / Estar / 16;
		double QA_pi = 7.09 * pow(Q_A_pt, (1. / 3.));
		QA_pi *= p.DELTA_P;
		Q_A_prime(i, 0) = QA_pi;

		double Q_A_pt_s = pow(WofA_s, 5.)*pow(D_p_real(i), 4.) / Estar_s / Estar_s / 16;
		double QA_pi_s = 7.09 * pow(Q_A_pt_s, (1. / 3.));
		QA_pi_s *= p.DELTA_P;
		Q_A_prime(i, 1) = QA_pi_s;

		double Mpi = PI_NUMBER*D_p_real(i)*D_p_real(i)*DENSITY_SOOT / 4.;
		Mpi *= p.DELTA_M;
		M_p(i) = Mpi;

		double F_poi = 625.*HAMAKER_CONST / 3. / D_p_real(i);
		F_poi *= p.DELTA_F;
		F_po(i) = F_poi;

		double Dpi = D_p_real(i)*p.DELTA_L;
		D_p(i) = Dpi;

		double Vpi = D_p(i) * D_p(i) * PI_NUMBER / 4. / (1. - DEP_POROSITY);
		V_par(i) = Vpi;

		double Tci = 4.*F_po(i)*R_d(i, 0) / (3. * PI_NUMBER * pow(D_p(i), 3.) * WALL_CORRECTION);
		
		Tau_crit(i, 0) = Tci;
		if (Tau_crit(i, 0) > Tau_crit_max)
			Tau_crit_max = Tau_crit(i, 0);

		double Tci_s = 4.*F_po(i)*R_d(i, 1) / (3. * PI_NUMBER * pow(D_p(i), 3.) * WALL_CORRECTION);
	

		Tau_crit(i, 1) = Tci_s;
		if (Tau_crit(i, 1) > Tau_crit_max)
			Tau_crit_max = Tau_crit(i, 1);
	}

	double Distsum = 0.;
	Mean_index = 0;
	parP.allocate(Nd);
	double dist_max = 0;
	int max_dist_num = 0;
	for (int i = 0; i < vtr.Nd; i++)
	{
		parP(i).Dp = D_p(i);
		parP(i).Q_A.x = Q_A(i, 0);
		parP(i).Q_A.y = Q_A(i, 1);
		parP(i).Q_A_prime.x = Q_A_prime(i, 0);
		parP(i).Q_A_prime.y = Q_A_prime(i,1);
		parP(i).tau_crit.x = Tau_crit(i, 0);
		parP(i).tau_crit.y = Tau_crit(i, 1);
		parP(i).Mp = M_p(i);
		parP(i).Kth = Kth_pars(i);
		parP(i).D_dist = Distsum;
		if (dist_max < D_dists(i))
			max_dist_num = i;
		Distsum += D_dists(i);
		parP(i).L_coeff = LIFT_COEFFICIENT * parP[i].Dp * parP[i].Dp * parP[i].Dp / MU_NUMBER / 8.;
		parP(i).D_coeff = 6 * PI_NUMBER * parP[i].Dp * parP[i].Dp * WALL_CORRECTION / 4.;
	}

	Tcrit_color.zeros(2);
	Tcrit_color(0) = parP(max_dist_num).tau_crit.x;
	Tcrit_color(1) = parP(max_dist_num).tau_crit.y;
	
	save_particle_arrays();

}

void clVariablesTR::save_particle_arrays()
{
	D_p.savetxt("D_p");
	D_p_real.savetxt("D_p_real");
	Q_A_prime.savetxt("Q_A_prime");
	Q_A.savetxt("Q_A");
	M_p.savetxt("M_p");
	F_po.savetxt("F_po");
	Tau_crit.savetxt("Tau_crit");
	V_par.savetxt("Vol_p");
	Kth_pars.savetxt("Kth");
	D_dists.savetxt("D_dist");
};

int clVariablesTR::get_type()
{
	double randval = rand1();
	for (int j = 0; j < vtr.Nd-1; j++)
	{
		if (randval >= parP[j].D_dist && randval < parP[j + 1].D_dist)
			return j;
	}
	return vtr.Nd - 1; 
}

double clVariablesTR::get_offset(double xval)
{
	int i = 0;
	while (TRUE)
	{
		if (vls.C0(i).x < xval && vls.C0(i + 1).x >= xval)
			break;
		i++; 
	}

	return vls.C0(i).y + (xval - vls.C0(i).x) * (vls.C0(i + 1).y - vls.C0(i).y) / (vls.C0(i + 1).x - vls.C0(i).x);

}

void clVariablesTR::ini_particles()
{
	for (int i = 0; i < nN; i++)
	{
		par Ptemp;

		Ptemp.pos.x = X_release + (rand1())*(X_STOP_POS - X_release);
		//Ptemp.pos.x = X_release + (0.25 + 0.75*rand1())*(X_STOP_POS - X_release);





		Ptemp.pos.y = get_offset(Ptemp.pos.x) + 7.*trP.bval / 16. + rand1()*trP.bval / 8.;
			
		Ptemp.type = get_type();
		Ptemp.timer = 1;
		Ptemp.Dep_Flag = -1;
		Ptemp.Num_rep = 0;
		Ptemp.Dep_timer = 1;
		cl_int2 Posi = { { (int)floor(Ptemp.pos.x), (int)floor(Ptemp.pos.y) } };

		Ptemp.loc = Posi.x*vtr.TrDomainSize.y + Posi.y;
		P[i] = Ptemp;
	}

}

int clVariablesTR::get_ind(cl_int2 ij)
{
	//int i = ij.x;
	//int j = ij.y;
	//int yred = vlb.Stor(i, j);
	//if (yred == -1)
	//	return i*(vlb.Channel_Height + 1) + vlb.Channel_Height;
	//else
	//	return i*(vlb.Channel_Height + 1) + vlb.Stor(i, j);
	return 0;
}

void clVariablesTR::ini_node()
{
	TR_indicies.allocate((p.nX - 1)*p.nY - 1);
	Active_indicies.allocate((p.nX - 1)*(p.nY - 1));
	Nod.allocate(p.nX - 1, p.nY - 1); //Temporary array storing lattice site info used to create NodI and NodV
	int ind = 0;
	int trind = 0;
	int MinXval = floor(X_RELEASE_POS - 2);
	int MaxXval = ceil(X_STOP_POS);

	//for (int i = 0; i < p.nX-1; i++)
	//{
	//	for (int j = 0; j < p.nY-1; j++)
	//	{
	//		cl_int2 ii00 = { { i, j } };
	//		cl_int2 ii10 = { { i + 1, j } };
	//		cl_int2 ii01 = { { i, j + 1 } };
	//		cl_int2 ii11 = { { i + 1, j + 1 } };

	//		Nact Neightemp;
	//		Neightemp.ii00 = ii00;
	//		Neightemp.ii10 = ii10;
	//		Neightemp.ii01 = ii01;
	//		Neightemp.ii11 = ii11;

	//		int jj00 = get_ind(ii00);
	//		int jj10 = get_ind(ii10);
	//		int jj01 = get_ind(ii01);
	//		int jj11 = get_ind(ii11);

	//		Nod(i, j).neigh.x = jj00;
	//		Nod(i, j).neigh.y = jj10;
	//		Nod(i, j).neigh.z = jj01;
	//		Nod(i, j).neigh.w = jj11;

	//		int tnum = 0;

	//		if (vls.M_o(ii00) != LB_SOLID)
	//			tnum += 1;
	//		if (vls.M_o(ii10) != LB_SOLID)
	//			tnum += 2;
	//		if (vls.M_o(ii11) != LB_SOLID)
	//			tnum += 4;
	//		if (vls.M_o(ii01) != LB_SOLID)
	//			tnum += 8;

	//		if (tnum > 0 && i > MinXval && i <= MaxXval)
	//		{
	//			Active_flag(i,j) = 1;
	//			cl_int2 trind_temp = { { i, j } };
	//			Active_indicies(trind) = trind_temp;
	//			TR_indicies(trind++) = i*(p.nY - 1) + j;
	//		}

	//		if (tnum == 0)
	//			tnum = -1;
	//		Nod(i, j).type = tnum;

	//		if (tnum == -1)
	//		{
	//			Nod(i, j).dX_cur.x = 1.;
	//			Nod(i, j).dX_cur.y = 1.;
	//		}
	//		else if (tnum == 0)
	//		{
	//			Nod(i, j).dX_cur.x = 1.;
	//			Nod(i, j).dX_cur.y = 1.;
	//		}
	//		else if (tnum == 1)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii00.x, ii00.y, 0);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii00.x, ii00.y, 2);
	//		}
	//		else if (tnum == 2)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii10.x, ii10.y, 1);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii10.x, ii10.y, 2);
	//		}
	//		else if (tnum == 3)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii00.x, ii00.y, 2);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii10.x, ii10.y, 2);
	//		}
	//		else if (tnum == 4)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii11.x, ii11.y, 1);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii11.x, ii11.y, 3);
	//		}
	//		else if (tnum == 5)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii00.x, ii00.y, 0);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii00.x, ii00.y, 2);
	//		}
	//		else if (tnum == 6)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii10.x, ii10.y, 1);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii11.x, ii11.y, 1);
	//		}
	//		else if (tnum == 7)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii00.x, ii00.y, 2);
	//			Nod(i, j).dX_cur.y = 1.;
	//		}
	//		else if (tnum == 8)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii01.x, ii01.y, 0);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii01.x, ii01.y, 3);
	//		}
	//		else if (tnum == 9)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii00.x, ii00.y, 0);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii01.x, ii01.y, 0);
	//		}
	//		else if (tnum == 10)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii10.x, ii10.y, 1);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii10.x, ii10.y, 2);
	//		}
	//		else if (tnum == 11)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii10.x, ii10.y, 2);
	//			Nod(i, j).dX_cur.y = 1.;
	//		}
	//		else if (tnum == 12)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii01.x, ii01.y, 3);
	//			Nod(i, j).dX_cur.y = vfd.dX_full(ii11.x, ii11.y, 3);
	//		}
	//		else if (tnum == 13)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii11.x, ii11.y, 3);
	//			Nod(i, j).dX_cur.y = 1.;
	//		}
	//		else if (tnum == 14)
	//		{
	//			Nod(i, j).dX_cur.x = vfd.dX_full(ii01.x, ii01.y, 3);
	//			Nod(i, j).dX_cur.y = 1.;
	//		}
	//		else if (tnum == 15)
	//		{
	//			Nod(i, j).dX_cur.x = 1.;
	//			Nod(i, j).dX_cur.y = 1.;
	//		}
	//		Nod(i, j).dX = Nod(i, j).dX_cur;
	//		ind++;
	//	}
	//}

	for (int i = 0; i < p.nX-1; i++)
	{
		for (int j = 0; j < p.nY-1; j++)
		{
			for (int k = 0; k < MAX_BL_PER_NODE; k++)
			{
				Nod(i, j).BLind[k] = -1;
			}
			Nod(i, j).Wall_Flag = 0;
		}
	}
	TR_indicies.reallocate_host_only(trind);
	Active_indicies.reallocate_host_only(trind);
	nActiveNodes = trind;
}

void clVariablesTR::saveParams()
{
	p.setParameter("Use Par Solver", trSolverFlag);

	if (p.Time > 0)
		p.setParameter("Restart Run", false);
	else
		p.setParameter("Restart Run", true);

	p.setParameter("Num Particles", nN);
	p.setParameter("Num Par Sizes", Nd);

	p.setParameter("Num Steps Btw Sort", numStepsBtwSort);
	p.setParameter("Soot Num Concentration", sootNumConc);

	p.setParameter("Par Themal Conductivity", parThermalCond);

	p.setParameter("Surface Energy Wall", surfEnergySurf);
	p.setParameter("Poisson Number Wall", poissonSurf);
	p.setParameter("Youngs Modulus Wall", yModSurf);

	p.setParameter("Surface Energy Soot", surfEnergySoot);
	p.setParameter("Poisson Number Soot", poissonSoot);
	p.setParameter("Youngs Modulus Soot", yModSoot);

	p.setParameter("MFP Air", mfpAir);

	p.setParameter("Hamaker Constant", hamakerConst);
	p.setParameter("Deposit Porosity", depPorosity);
	p.setParameter("Wall Correction", wallCorrection);
	p.setParameter("Lift Coefficient", liftCoeff);

	p.setParameter("Index Radius Search", indRadiusSearch);
	p.setParameter("Mass Flux Inlet", massFluxInlet);
	p.setParameter("Max BL Per Node", maxBLPerNode);
	p.setParameter("Particle Release Time", parReleaseTime);
	p.setParameter("X Stop Distance", stopDistX);
	p.setParameter("Time Before Re-release", timeBeforeReRelease*p.trSteps_wall);
	p.setParameter("Time Before Stick", timeBeforeStick*p.trSteps_wall);
	p.setParameter("Max Num BL Roll", maxNumBLRoll);
	p.setParameter("X Release Pos", xReleasePos);
	p.setParameter("X Stop Pos", xStopPos);
	p.setParameter("Reduce Dep Stop1", reduceDepStop1);
	p.setParameter("Reduce Dep Stop2", reduceDepStop2);
	p.setParameter("Start Thermo Pos", startThermoVel);
	p.setParameter("Amount Reduce Dep", amtReduceDep);
	p.setParameter("X Min Val", xMinVal);
	p.setParameter("Cutoff Radius", cutoffRadius);
	p.setParameter("Foul Size Switch Soot", foulSizeSwitchSoot);
	p.setParameter("Num Each Par", numEachPar);

	if (trSolverFlag)
	{
		*p.yamlOut << YAML::BeginMap;
		*p.yamlOut << YAML::Key << "Dp Dists";
		*p.yamlOut << YAML::Value << D_dists;
		*p.yamlOut << YAML::EndMap;

		*p.yamlOut << YAML::BeginMap;
		*p.yamlOut << YAML::Key << "Par Diameters";
		*p.yamlOut << YAML::Value << D_p_real;
		*p.yamlOut << YAML::EndMap;
	}
}
void clVariablesTR::Split_Node_arrays()
{
	double Ka = (THERMAL_CONDUCTIVITY_AIR * p.DELTA_F / p.DELTA_T);
	double Ks = (THERMAL_CONDUCTIVITY_FOUL * p.DELTA_F / p.DELTA_T);

	for (int i = 0; i < (p.nX-1)*(p.nY-1); i++)
	{
		NodC[i].neigh.x = Nod[i].neigh.x;
		NodC[i].neigh.y = Nod[i].neigh.y;
		NodC[i].neigh.z = Nod[i].neigh.z;
		NodC[i].neigh.w = Nod[i].neigh.w;
		double dXcx = Nod[i].dX_cur.x, dXcy = Nod[i].dX_cur.y;
		double dXx = Nod[i].dX.x, dXy = Nod[i].dX.y;
		int typet = Nod[i].type;
		switch (typet)
		{
		case -1:
		{
				   break;
		}
		case 1:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffT10 = { { (C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0. } };
				  NodC[i].CoeffT01 = { { (C1y - 1.) / dXcy + 1., 0., C2y / dXcy, 0. } };
				  NodC[i].CoeffT11 = { { (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., C2x / dXcx, C2y / dXcy, 0. } };

				  NodC[i].CoeffU00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffU10 = { { 1. - 1. / dXcx, 0., 0., 0. } };
				  NodC[i].CoeffU01 = { { 1. - 1. / dXcy, 0., 0., 0. } };
				  NodC[i].CoeffU11 = { { 1. - 1. / dXcy - 1. / dXcx, 0, 0, 0 } };
				  break;
		}
		case 2:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { C2x / dXcx, (C1x - 1.) / dXcx + 1., 0, 0 } };
				  NodC[i].CoeffT10 = { { 0., 1., 0., 0. } };
				  NodC[i].CoeffT01 = { { C2x / dXcx, (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., 0., C2y / dXcy } };
				  NodC[i].CoeffT11 = { { 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy } };

				  NodC[i].CoeffU00 = { { 0., 1. - 1. / dXcx, 0., 0. } };
				  NodC[i].CoeffU10 = { { 0., 1., 0., 0. } };
				  NodC[i].CoeffU01 = { { 0., 1. - 1. / dXcx - 1. / dXcy, 0., 0. } };
				  NodC[i].CoeffU11 = { { 0., 1. - 1. / dXcy, 0., 0. } };
				  break;
		}
		case 3:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffT10 = { { 0., 1., 0., 0. } };
				  NodC[i].CoeffT01 = { { (C1x - 1.) / dXcx + 1., 0, C2x / dXcx, 0. } };
				  NodC[i].CoeffT11 = { { 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy } };

				  NodC[i].CoeffU00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffU10 = { { 0., 1., 0., 0. } };
				  NodC[i].CoeffU01 = { { 1. - 1. / dXcx, 0., 0., 0. } };
				  NodC[i].CoeffU11 = { { 0., 1. - 1. / dXcy, 0., 0. } };
				  break;
		}
		case 4:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { 0., C2y / dXcy, C2x / dXcx, (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1. } };
				  NodC[i].CoeffT10 = { { 0., C2y / dXcy, 0., (C1y - 1.) / dXcy + 1. } };
				  NodC[i].CoeffT01 = { { 0., 0., C2x / dXcx, (C1x - 1.) / dXcx + 1. } };
				  NodC[i].CoeffT11 = { { 0., 0., 0., 1. } };

				  NodC[i].CoeffU00 = { { 0., 0., 0., 1. - 1. / dXcy - 1. / dXcx } };
				  NodC[i].CoeffU10 = { { 0., 0., 0., 1. - 1. / dXcy } };
				  NodC[i].CoeffU01 = { { 0., 0., 0., 1. - 1. / dXcx } };
				  NodC[i].CoeffU11 = { { 0., 0., 0., 1. } };
				  break;
		}
		case 5:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffT10 = { { (C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0. } };
				  NodC[i].CoeffT01 = { { (C1y - 1.) / dXcy + 1., 0., C2y / dXcy, 0. } };
				  NodC[i].CoeffT11 = { { 0., 0., 0., 1. } };

				  NodC[i].CoeffU00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffU10 = { { 1. - 1. / dXcx, 0., 0., 0. } };
				  NodC[i].CoeffU01 = { { 1. - 1. / dXcy, 0., 0., 0. } };
				  NodC[i].CoeffU11 = { { 0., 0., 0., 1. } };
				  break;
		}
		case 6:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { C2x / dXcx, (C1x - 1.) / dXcx + 1., 0., 0. } };
				  NodC[i].CoeffT10 = { { 0., 1., 0., 0. } };
				  NodC[i].CoeffT01 = { { 0., 0., C2y / dXcy, (C1y - 1.) / dXcy + 1. } };
				  NodC[i].CoeffT11 = { { 0., 0., 0., 1. } };

				  NodC[i].CoeffU00 = { { 0., 1. - 1. / dXcx, 0., 0. } };
				  NodC[i].CoeffU10 = { { 0., 1., 0., 0. } };
				  NodC[i].CoeffU01 = { { 0., 0., 0., 1. - 1. / dXcy } };
				  NodC[i].CoeffU11 = { { 0., 0., 0., 1. } };
				  break;
		}
		case 7:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffT10 = { { 0., 1., 0., 0. } };
				  NodC[i].CoeffT01 = { { 1. - (C1x - 1.) / dXcx, 0., -C2x / dXcx, 0. } };
				  NodC[i].CoeffT11 = { { 0., 0., 0., 1. } };

				  NodC[i].CoeffU00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffU10 = { { 0., 1., 0., 0. } };
				  NodC[i].CoeffU01 = { { 1. - 1. / dXcx, 0., 0., 0. } };
				  NodC[i].CoeffU11 = { { 0., 0., 0., 1. } };
				  break;
		}
		case 8:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { C2y / dXcy, 0., (C1y - 1.) / dXcy + 1., 0. } };
				  NodC[i].CoeffT10 = { { C2y / dXcy, 0., (C1x - 1.) / dXcx + (C1y - 1.) / dXcy + 1., C2x / dXcx } };
				  NodC[i].CoeffT01 = { { 0., 0., 1., 0. } };
				  NodC[i].CoeffT11 = { { 0., 0., (C1x - 1.) / dXcx + 1., C2x / dXcx } };

				  NodC[i].CoeffU00 = { { 0, 0, 1 - 1 / dXcy, 0 } };
				  NodC[i].CoeffU10 = { { 0, 0, 1. - 1. / dXcy - 1. / dXcx, 0 } };
				  NodC[i].CoeffU01 = { { 0., 0., 1., 0. } };
				  NodC[i].CoeffU11 = { { 0, 0, 1. - 1. / dXcx, 0 } };
				  break;
		}
		case 9:
		{
				  double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				  double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				  double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				  NodC[i].CoeffT00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffT10 = { { (C1x - 1.) / dXcx + 1., C2x / dXcx, 0., 0. } };
				  NodC[i].CoeffT01 = { { 0., 0., 1., 0. } };
				  NodC[i].CoeffT11 = { { 0., 0., (C1y - 1.) / dXcy + 1., C2y / dXcy } };

				  NodC[i].CoeffU00 = { { 1., 0., 0., 0. } };
				  NodC[i].CoeffU10 = { { 1. - 1. / dXcx, 0., 0., 0. } };
				  NodC[i].CoeffU01 = { { 0., 0., 1., 0. } };
				  NodC[i].CoeffU11 = { { 0., 0., 1. - 1. / dXcy, 0. } };
				  break;
		}
		case 10:
		{
				   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				   double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				   double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				   NodC[i].CoeffT00 = { { C2x / dXcx, (C1x - 1.) / dXcx + 1., 0, 0 } };
				   NodC[i].CoeffT10 = { { 0., 1., 0., 0. } };
				   NodC[i].CoeffT01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffT11 = { { 0., (C1y - 1.) / dXcy + 1., 0, C2y / dXcy } };

				   NodC[i].CoeffU00 = { { 0., 1. - 1. / dXcx, 0., 0. } };
				   NodC[i].CoeffU10 = { { 0., 1., 0., 0. } };
				   NodC[i].CoeffU01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffU11 = { { 0., 1. - 1. / dXcy, 0., 0. } };
				   break;
		}
		case 11:
		{
				   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));


				   NodC[i].CoeffT00 = { { 1., 0., 0., 0. } };
				   NodC[i].CoeffT10 = { { 0., 1., 0., 0. } };
				   NodC[i].CoeffT01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffT11 = { { 0., (C1x - 1.) / dXcx + 1., 0, C2x / dXcx } };

				   NodC[i].CoeffU00 = { { 1., 0., 0., 0. } };
				   NodC[i].CoeffU10 = { { 0., 1., 0., 0. } };
				   NodC[i].CoeffU01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffU11 = { { 0., 1. - 1. / dXcx, 0., 0. } };
				   break;
		}
		case 12:
		{
				   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				   double C1y = -(Ka*(dXcy - dXy)) / (Ks*dXcy - Ka*(dXcy - dXy));
				   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));
				   double C2y = (Ks*dXcy) / (Ks*dXcy - Ka*(dXcy - dXy));

				   NodC[i].CoeffT00 = { { C2x / dXcx, 0., (C1x - 1.) / dXcx + 1., 0 } };
				   NodC[i].CoeffT10 = { { 0., C2y / dXcy, 0., (C1y - 1.) / dXcy + 1. } };
				   NodC[i].CoeffT01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffT11 = { { 0., 0., 0., 1. } };

				   NodC[i].CoeffU00 = { { 0., 0., 1. - 1. / dXcx, 0. } };
				   NodC[i].CoeffU10 = { { 0., 0., 0., 1. - 1. / dXcy } };
				   NodC[i].CoeffU01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffU11 = { { 0., 0., 0., 1. } };
				   break;
		}
		case 13:
		{
				   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));

				   NodC[i].CoeffT00 = { { 1., 0., 0., 0 } };
				   NodC[i].CoeffT10 = { { 0., C2x / dXcx, 0., (C1x - 1.) / dXcy + 1. } };
				   NodC[i].CoeffT01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffT11 = { { 0., 0., 0., 1. } };

				   NodC[i].CoeffU00 = { { 1., 0., 0., 0. } };
				   NodC[i].CoeffU10 = { { 0., 0., 0., 1. - 1. / dXcx } };
				   NodC[i].CoeffU01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffU11 = { { 0., 0., 0., 1. } };
				   break;
		}
		case 14:
		{
				   double C1x = -(Ka*(dXcx - dXx)) / (Ks*dXcx - Ka*(dXcx - dXx));
				   double C2x = (Ks*dXcx) / (Ks*dXcx - Ka*(dXcx - dXx));

				   NodC[i].CoeffT00 = { { C2x / dXcx, 0., (C1x - 1.) / dXcx + 1., 0 } };
				   NodC[i].CoeffT10 = { { 0., 1., 0., 0. } };
				   NodC[i].CoeffT01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffT11 = { { 0., 0., 0., 1. } };

				   NodC[i].CoeffU00 = { { 0., 0., 1. - 1. / dXcx, 0. } };
				   NodC[i].CoeffU10 = { { 0., 1., 0., 0. } };
				   NodC[i].CoeffU01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffU11 = { { 0., 0., 0., 1. } };
				   break;
		}
		case 15:
		{
				   NodC[i].CoeffT00 = { { 1., 0., 0., 0. } };
				   NodC[i].CoeffT10 = { { 0., 1., 0., 0. } };
				   NodC[i].CoeffT01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffT11 = { { 0., 0., 0., 1. } };

				   NodC[i].CoeffU00 = { { 1., 0., 0., 0. } };
				   NodC[i].CoeffU10 = { { 0., 1., 0., 0. } };
				   NodC[i].CoeffU01 = { { 0., 0., 1., 0. } };
				   NodC[i].CoeffU11 = { { 0., 0., 0., 1. } };
				   break;
		}
		}

		NodV[i].Temps = { { 0., 0., 0., 0. } };
		NodV[i].U00 = { { 0., 0. } };
		NodV[i].U10 = { { 0., 0. } };
		NodV[i].U01 = { { 0., 0. } };
		NodV[i].U11 = { { 0., 0. } };

		for (int j = 0; j < MAX_BL_PER_NODE; j++)
		{
			NodI[i].BLind[j] = Nod[i].BLind[j];
		}
		NodI[i].Wall_Flag = Nod[i].Wall_Flag;

	}
	Nod.FreeHost();
}


double clVariablesTR::rand1()
{
	double randval = (double)rand();
	return randval / rmax;
}

void clVariablesTR::Release_Objects()
{// THis should be handled by destructor
//	TR_ReRelease_kernel.free_memory();
//	TR_Shear_Removal_kernel.free_memory();
//	TR_Update_Par_kernel.free_memory();
//	TR_Wall_Par_kernel[0].free_memory();
//	TR_Wall_Par_kernel[1].free_memory();
//	TR_Wall_Node_kernel[0].free_memory();
//	TR_Wall_Node_kernel[1].free_memory();
//	TR_Node_kernel[0].free_memory();
//	TR_Node_kernel[1].free_memory();
//	
//	Shear_kernels[0].free_memory();
//	Shear_kernels[1].free_memory();
//	Shear_kernels[2].free_memory();
//	Sort_kernel[0].free_memory();
//	Sort_kernel[1].free_memory();
//	Sort_kernel[2].free_memory();
//	Get_Umax_kernel.free_memory();
//
//	P.delete_array();
//	Nod.delete_array();
//	NodI.delete_array();
//	NodC.delete_array();
//	NodV.delete_array();
//	parP.delete_array();
//	
//	BL.delete_array();
//	RandList.delete_array();
//	BLind_ind.delete_array();
//	BL_dep.delete_array();
//	Sind.delete_array();
//	BLindicies.delete_array();
//	Weights.delete_array();
//	Shear_inds.delete_array();
//	Shear_coeffs.delete_array();
//	Ploc.delete_array();
//	PV.delete_array();
//	Tau.delete_array();
//	SortTmp.delete_array();
//	Node_neigh.delete_array();
//	Winds.delete_array();
//	Num_W_nodes.delete_array();
//	P_vbo.delete_array(IOQUEUE);
//	Umax_val.delete_array(IOQUEUE);
//
//#ifdef USE_OPENGL
//	TR_GL_kernel.free_memory();
//	Par_Color_List.delete_array();
//#endif
}



void clVariablesTR::ini_sort()
{
	Sort_kernel[0].set_argument(0, P.get_buf_add());
	Sort_kernel[0].set_argument(1, &nN);
	Sort_kernel[0].set_local_memory(2, WORKGROUPSIZE_SORT*sizeof(par));
	Sort_kernel[0].set_local_memory(3, WORKGROUPSIZE_SORT*sizeof(par));

	numMerges = 0;
	localRange = WORKGROUPSIZE_SORT;
	size_t log2BlockSize = nN / WORKGROUPSIZE_SORT;
	for (; log2BlockSize > 1; log2BlockSize >>= 1)
	{
		++numMerges;
	}

	if (numMerges & 1)
	{
		Sort_kernel[2].set_argument(0, SortTmp.get_buf_add());
	}
	else
	{
		Sort_kernel[2].set_argument(0, P.get_buf_add());
	}
	Sort_kernel[2].set_argument(1, Ploc.get_buf_add());


	size_t vecPow2 = (nN & (nN - 1));
	numMerges += vecPow2 ? 1 : 0;

	Sort_kernel[1].set_argument(2, &nN);
}

#ifdef USE_OPENGL

void clVariablesTR::sort_particles(cl_event *TR_prev_evt)
{
	cl_event kernelEvent;
	Sort_kernel[0].call_kernel();
	int n1 = -1;
	Ploc.FillBuffer(IOQUEUE, (void*)&Ploc_inds, sizeof(cl_int4), 1, 1, TR_prev_evt);
	Ploc.FillBuffer(IOQUEUE, (void*)&n1, sizeof(int), 2*vtr.TrDomainSize.x*vtr.TrDomainSize.y, 4, 0, NULL, &kernelEvent);
	for (size_t pass = 1; pass <= numMerges; ++pass)
	{
		unsigned srcLogicalBlockSize = static_cast< unsigned >(localRange << (pass - 1));
		Sort_kernel[1].set_argument(3, sizeof(int), (void*)&srcLogicalBlockSize);
		if (pass & 0x1)
		{
			Sort_kernel[1].set_argument(0, sizeof(cl_mem), P.get_buf_add());
			Sort_kernel[1].set_argument(1, sizeof(cl_mem), SortTmp.get_buf_add());
		}
		else
		{
			Sort_kernel[1].set_argument(1, sizeof(cl_mem), P.get_buf_add());
			Sort_kernel[1].set_argument(0, sizeof(cl_mem), SortTmp.get_buf_add());
		}
		Sort_kernel[1].call_kernel();
	}

	Sort_kernel[2].call_kernel(1, &kernelEvent);
	clReleaseEvent(kernelEvent);
	if (numMerges & 1)
	{
		P.enqueue_copy_to_buffer(TRQUEUE, SortTmp.get_buffer());
	}
	
	Ploc.read_from_buffer_size(TRQUEUE, CL_FALSE, 2);

}
#else
void clVariablesTR::sort_particles()
{
	//cl_event kernelEvent;
	//Sort_kernel[0].call_kernel();
	//int n1 = -1;
	//Ploc.FillBuffer(IOQUEUE, (void*)&Ploc_inds, sizeof(cl_int4), 1, 1);
	//Ploc.FillBuffer(IOQUEUE, (void*)&n1, sizeof(int), 2*vtr.TrDomainSize.x*vtr.TrDomainSize.y, 4, 0, NULL, &kernelEvent);
	//for (size_t pass = 1; pass <= numMerges; ++pass)
	//{
	//	unsigned srcLogicalBlockSize = static_cast< unsigned >(localRange << (pass - 1));
	//	Sort_kernel[1].set_argument(3, sizeof(int), (void*)&srcLogicalBlockSize);
	//	if (pass & 0x1)
	//	{
	//		Sort_kernel[1].set_argument(0, sizeof(cl_mem), P.get_buf_add());
	//		Sort_kernel[1].set_argument(1, sizeof(cl_mem), SortTmp.get_buf_add());
	//	}
	//	else
	//	{
	//		Sort_kernel[1].set_argument(1, sizeof(cl_mem), P.get_buf_add());
	//		Sort_kernel[1].set_argument(0, sizeof(cl_mem), SortTmp.get_buf_add());
	//	}
	//	Sort_kernel[1].call_kernel();
	//}

	//Sort_kernel[2].call_kernel(1, &kernelEvent);
	//clReleaseEvent(kernelEvent);
	//if (numMerges & 1)
	//{
	//	P.enqueue_copy_to_buffer(TRQUEUE, SortTmp.get_buffer());
	//}

	//Ploc.read_from_buffer_size(TRQUEUE, CL_FALSE, 2);
}


#endif

void clVariablesTR::initial_sort_particles()
{
	//cl_event kernelEvent;
	//Sort_kernel[0].call_kernel();
	//int n1 = -1;
	//Ploc.FillBuffer(IOQUEUE, (void*)&Ploc_inds, sizeof(cl_int4), 1);
	//Ploc.FillBuffer(IOQUEUE, (void*)&n1, sizeof(int), 2*vtr.TrDomainSize.x*vtr.TrDomainSize.y, 4, 0, NULL, &kernelEvent);
	//for (size_t pass = 1; pass <= numMerges; ++pass)
	//{
	//	unsigned srcLogicalBlockSize = static_cast< unsigned >(localRange << (pass - 1));
	//	Sort_kernel[1].set_argument(3, sizeof(int), (void*)&srcLogicalBlockSize);
	//	if (pass & 0x1)
	//	{
	//		Sort_kernel[1].set_argument(0, sizeof(cl_mem), P.get_buf_add());
	//		Sort_kernel[1].set_argument(1, sizeof(cl_mem), SortTmp.get_buf_add());
	//	}
	//	else
	//	{
	//		Sort_kernel[1].set_argument(1, sizeof(cl_mem), P.get_buf_add());
	//		Sort_kernel[1].set_argument(0, sizeof(cl_mem), SortTmp.get_buf_add());
	//	}
	//	Sort_kernel[1].call_kernel();
	//}

	//Sort_kernel[2].call_kernel(1, &kernelEvent);
	//clReleaseEvent(kernelEvent);
	//if (numMerges & 1)
	//{
	//	P.enqueue_copy_to_buffer(TRQUEUE, SortTmp.get_buffer());
	//}



	//Ploc.read_from_buffer_size(TRQUEUE, CL_TRUE, 2);
	//
	//if (vtr.Ploc(1).y > 0)
	//	vtr.Update_Par_Rem_Args();

	//if (vtr.Ploc(0).y > 0)
	//	vtr.Re_Release_Par();

	//vtr.Sort_timer = NUM_STEPS_BTW_SORT;
}


void clVariablesTR::sort_particles_for_clumping()
{
	//cl_event kernelEvent;
	//Sort_kernel[0].call_kernel();
	//int n1 = -1;
	//Ploc.FillBuffer(IOQUEUE, (void*)&Ploc_inds, sizeof(cl_int4), 1);
	//Ploc.FillBuffer(IOQUEUE, (void*)&n1, sizeof(int), 2*vtr.TrDomainSize.x*vtr.TrDomainSize.y, 4, 0, NULL, &kernelEvent);
	//for (size_t pass = 1; pass <= numMerges; ++pass)
	//{
	//	unsigned srcLogicalBlockSize = static_cast< unsigned >(localRange << (pass - 1));
	//	Sort_kernel[1].set_argument(3, sizeof(int), (void*)&srcLogicalBlockSize);
	//	if (pass & 0x1)
	//	{
	//		Sort_kernel[1].set_argument(0, sizeof(cl_mem), P.get_buf_add());
	//		Sort_kernel[1].set_argument(1, sizeof(cl_mem), SortTmp.get_buf_add());
	//	}
	//	else
	//	{
	//		Sort_kernel[1].set_argument(1, sizeof(cl_mem), P.get_buf_add());
	//		Sort_kernel[1].set_argument(0, sizeof(cl_mem), SortTmp.get_buf_add());
	//	}
	//	Sort_kernel[1].call_kernel();
	//}

	//Sort_kernel[2].call_kernel(1, &kernelEvent);
	//clReleaseEvent(kernelEvent);
	//if (numMerges & 1)
	//{
	//	P.enqueue_copy_to_buffer(TRQUEUE, SortTmp.get_buffer());
	//}

	//Ploc.read_from_buffer_size(TRQUEUE, CL_TRUE, 2);

	//vtr.Sort_timer = NUM_STEPS_BTW_SORT;
}

void clVariablesTR::iniGL()
{
	//numpargl = (int)floor((double)nN / NUM_PAR_GL_DIV);
	//
	//P_vbo.zeros(numpargl * 2);
	//P_vbo.AllocateColor(numpargl * 3);
	//P_vbo.set_color_array({ { 1.f, 0.f, 0.f } });
	//for (int i = 0; i < numpargl; i++)
	//{
	//	P_vbo(i * 2 + 0) = P(i).pos.x;
	//	P_vbo(i * 2 + 1) = P(i).pos.y;
	//}
	//
	//P_vbo.CreateVBO_position(POINT_VBO, , CL_MEM_WRITE_ONLY);
	//P_vbo.copy_to_buffer(IOQUEUE);
	//
	//P_vbo.CreateVBO_color(, CL_MEM_WRITE_ONLY);

	//TR_GL_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_opengl_par");
	//TR_GL_kernel.set_size_1D(TRC_NUM_TRACERS, WORKGROUPSIZE_TR_GL);
	//Tcrit_color.allocate_buffer_w_copy();
	//Shear_kernels[2].set_argument(6, sizeof(cl_mem), vls.LSb_vbo.get_col_buf_add());
	//Shear_kernels[2].set_argument(7, sizeof(cl_mem), vls.LSt_vbo.get_col_buf_add());
	//Shear_kernels[2].set_argument(8, sizeof(cl_mem), Tcrit_color.get_buf_add());

	//int ind = 0;
	//TR_GL_kernel.set_argument(ind++, sizeof(cl_mem), P.get_buf_add());
	//TR_GL_kernel.set_argument(ind++, sizeof(cl_mem), P_vbo.get_buf_add());
	//TR_GL_kernel.set_argument(ind++, sizeof(cl_mem), P_vbo.get_col_buf_add());
	//TR_GL_kernel.set_argument(ind++, sizeof(cl_mem), Par_Color_List.get_buf_add());

	//P_vbo.FreeHostColor();
	//P_vbo.FreeHost();
}
