// clVariables.cpp: implementation of the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////


// TODO: 
//	For all changes to classes make sure arrays spanning domain are allocated with padding in x direction
//		and make sure necessary changes are made to kernels, to account for these changes
//	1) Update vlb -------------------------------------------------------- DONE
//	2) Update vls -------------------------------------------------------- DONE
//	3) Test LB ----------------------------------------------------------- DONE
//	4) Update vfd to reflect changes in previous methods and to
//		use viennacl (just doing simple laminar HT)				---------- DONE
//	5) Test LB/FD
//	6) Check and make necessary changes to Reduction kernels and
//		output methods (esp. save states) for LB/FD
//	7) Implement turbulence model, be sure to include necessary 
//		save states for these
//	8) Test LB/Turbulence model
//	9) Add necessary terms to TFD for turbulence
//	10) Test turbulent HT
//	11) Update vtr to reflect changes in previous methods. Make 
//		sure to account for changes in viscosity and shape of arrays
//	12) Update vfl to reflect changes in previous methods. Make sure
//		to account for removed methods, changes in array shapes and 
//		add new methods necessary for turb model
//	13) Make sure reduction kernels and output methods for vtr and 
//		vfl account for changes, 
//	14) Test full implementation.

#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"
#include "clProblem.h"


void clVariablesLB::allocateArrays()
{
	FA.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 9, 9);
	FB.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 9, 9);
	Ux_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	Uy_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	Ro_array.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
	NodeType.zeros(p.nX, p.XsizeFull, p.nY, p.nY);
}

void clVariablesLB::allocateBuffers()
{
	FA.allocate_buffer_w_copy();
	FB.allocate_buffer_w_copy();

	double zer = 0.;
	Ro_array.allocate_buffer_size(p.FullSize);
	Ro_array.FillBuffer(0., p.FullSize);
	Ro_array.copy_to_buffer();

	Ux_array.allocate_buffer_size(p.FullSize);
	Ux_array.FillBuffer(0., p.FullSize);
	Ux_array.copy_to_buffer();

	Uy_array.allocate_buffer_size(p.FullSize);
	Uy_array.FillBuffer(0., p.FullSize);
	Uy_array.copy_to_buffer();
}



void clVariablesLB::calcFeqVals(double *Feq_vals, double rho_val, cl_double2 u_val)
{
	Feq_vals[0] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y -
		2. * u_val.x*u_val.x / 3 - 2 * u_val.y*u_val.y / 3. + 4. / 9.);
	
	Feq_vals[1] = rho_val*(-u_val.x*u_val.x * u_val.y*u_val.y +
		2. * u_val.x*u_val.x / 3 - u_val.x*u_val.y*u_val.y +
		2. * u_val.x / 3. - u_val.y*u_val.y / 3. + 2. / 9.) / 2.;
	
	Feq_vals[2] = rho_val*(-u_val.x*u_val.x * u_val.y*u_val.y +
		2. * u_val.x*u_val.x / 3 + u_val.x*u_val.y*u_val.y - 
		2. * u_val.x / 3. - u_val.y*u_val.y / 3. + 2. / 9.) / 2.;
	
	Feq_vals[3] = rho_val*(-u_val.x*u_val.x * u_val.y*u_val.y -
		u_val.x*u_val.x * u_val.y - u_val.x*u_val.x / 3. +
		2. * u_val.y*u_val.y / 3. + 2. * u_val.y / 3. + 2. / 9.) / 2.;
	
	Feq_vals[4] = rho_val*(-u_val.x*u_val.x * u_val.y*u_val.y +
		u_val.x*u_val.x * u_val.y - u_val.x*u_val.x / 3. +
		2. * u_val.y*u_val.y / 3. - 2. * u_val.y / 3. + 2. / 9.) / 2.;
	
	Feq_vals[5] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y +
		u_val.x*u_val.x * u_val.y + u_val.x*u_val.x / 3. +
		u_val.x*u_val.y*u_val.y + u_val.x*u_val.y + u_val.x / 3. +
		u_val.y*u_val.y / 3. + u_val.y / 3. + 1. / 9.) / 4.;
	
	Feq_vals[6] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y -
		u_val.x*u_val.x * u_val.y + u_val.x*u_val.x / 3. -
		u_val.x*u_val.y*u_val.y + u_val.x*u_val.y - u_val.x / 3. + 
		u_val.y*u_val.y / 3. - u_val.y / 3. + 1. / 9.) / 4.;
	
	Feq_vals[7] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y -
		u_val.x*u_val.x * u_val.y + u_val.x*u_val.x / 3. +
		u_val.x*u_val.y*u_val.y - u_val.x*u_val.y + u_val.x / 3. +
		u_val.y*u_val.y / 3. - u_val.y / 3. + 1. / 9.) / 4.;
	
	Feq_vals[8] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y +
		u_val.x*u_val.x * u_val.y + u_val.x*u_val.x / 3. - 
		u_val.x*u_val.y*u_val.y - u_val.x*u_val.y - u_val.x / 3. +
		u_val.y*u_val.y / 3. + u_val.y / 3. + 1. / 9.) / 4.;
}




double clVariablesLB::calcUmean()
{
	return sumUx.reduceSingle() / sumRo.reduceSingle();
}

void clVariablesLB::createKernels()
{
	if (kOmegaClass.kOmegaSolverFlag)
	{
		Collision_kernel.create_kernel(GetSourceProgram, LBQUEUE_REF,
			"LB_collision_SRT_Fluid_w_kOmega", "LB_collision_SRT_Fluid_w_kOmega");
		kOmegaClass.createKernels();
	}
	else
	{
		Collision_kernel.create_kernel(GetSourceProgram, LBQUEUE_REF,
			"LB_collision_SRT_Fluid", "LB_collision_SRT_Fluid");
	}

	int GlobalX = p.XsizeFull;
	int GlobalY = p.nY;

	Collision_kernel.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);
	
	calcRoJKernel.create_kernel(GetSourceProgram, LBQUEUE_REF, "Calc_Ro_J");
	calcRoJKernel.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);


	




	//IBB_kernel_Fluid[0].create_kernel(p.program, &p.LBqueue, "LB_IBB");
	//IBB_kernel_Fluid[1].create_kernel(p.program, &p.LBqueue, "LB_IBB");

	//IBB_kernel_Fluid[0].set_size_1D(vls.length_ibb, WORKGROUPSIZE_IBB);
	//IBB_kernel_Fluid[1].set_size_1D(vls.length_ibb, WORKGROUPSIZE_IBB);

	// Update_output_kernel[0].create_kernel(p.program, &p.LBqueue, "LB_reduce_Ro_J");
	// Update_output_kernel[0].set_size(p.FullSize, WORKGROUPSIZE_RED);

	// if(vfd.thermalSolverFlag)
	// {
	// 	Update_output_kernel[1].create_kernel(p.program, &p.LBqueue, "FD_reduce_Ro_J_T");
	// 	Update_output_kernel[1].set_size(5, 1);
	// }
	// else
	// {
	// 	Update_output_kernel[1].create_kernel(p.program, &p.LBqueue, "LB_reduce_Ro_J_2");
	// 	Update_output_kernel[1].set_size_1D(2, 1);
	// }

}

void clVariablesLB::freeHostArrays()
{
}


void clVariablesLB::getInitialFlow()
{
	int Time = 0;
	int option = (OPTION_SAVE_MACRO_FIELDS);
	Collision_kernel.setOption(&option);
	int NextDumpSteptime = tlbmIniDumpStep - 1;
	int dump_step_interval = NextDumpSteptime + 1;
	int stopTime = tlbmNumIniSteps_LB;

	clock_t start_time = clock();
	double sizedom = (double)p.nX * (double)p.Channel_Height * (double)CLOCKS_PER_SEC / 1.0e6;

	while (Time < stopTime && !fLoadedFlag)
	{
		Solve();
		clFinish(LBQUEUE);

		if (kOmegaClass.kOmegaSolverFlag)
			kOmegaClass.Solve();

		Time++;
		if (Time > NextDumpSteptime)
		{
			save2file();
			NextDumpSteptime += dump_step_interval;
		}

		if (MOD(Time, 1000) == 0)
		{
			clock_t elap = clock() - start_time;
			double mlup = sizedom*(double)Time / (double)elap;
			double redRho = sumRo.reduceSingle();
			double redUx = sumUx.reduceSingle() / redRho;

			double redKap = kOmegaClass.sumKappa.reduceSingle() / redRho;
			double redOme = kOmegaClass.sumOmega.reduceSingle() / redRho;
			double redKapMin = kOmegaClass.minKappa.reduceSingle();
			double redOmeMin = kOmegaClass.minOmega.reduceSingle();
			printf("Time = %d, Umean = %g, K = %g, Omega = %g, Re = %g MLUPS = %g\n",
				Time, redUx, redKap, redOme, redUx * p.Pipe_radius / MuVal, mlup);

			if (redOmeMin < 0.)
				printf("warning: omega contains negative values of %g\n", redOmeMin);
			if (redKapMin < 0.)
				printf("warning: omega contains negative values of %g\n", redKapMin);
		}
	}
	save2file();

#ifdef CREATE_BIN_FILES
	std::string NewDir = "load" SLASH "LBFD";
	p.MakeDir(NewDir);
	saveRestartFiles();
	p.RenameFile("load\\lbf.bin", "load\\LBFD\\lbf.bin");
	if (kOmegaClass.kOmegaSolverFlag)
	{
		p.RenameFile("load\\lbkappa.bin", "load\\LBFD\\lbkappa.bin");
		p.RenameFile("load\\lbomega.bin", "load\\LBFD\\lbomega.bin");
	}
#endif
}

void clVariablesLB::getIntialFlowAndTemp()
{
	vfd.solveSSThermal();
	
	int option = (OPTION_SAVE_MACRO_FIELDS);
	Collision_kernel.setOption(&option);
	int NextDumpSteptime = tlbmIniDumpStep - 1;
	int dump_step_interval = NextDumpSteptime + 1;
	
	int Time = 0;
	int stopTime = tlbmNumIniSteps_TLBM;

	clock_t start_time = clock();
	double sizedom = (double)p.nX * (double)p.Channel_Height * (double)CLOCKS_PER_SEC / 1.0e6;

	while (Time < stopTime)
	{
		Solve();
		clFinish(LBQUEUE);
		if (kOmegaClass.kOmegaSolverFlag)
			kOmegaClass.Solve();
		vfd.solveTemp();

		Time++;
		if (Time > NextDumpSteptime)
		{
			save2file();
			NextDumpSteptime += dump_step_interval;
		}

		if (MOD(Time, 1000) == 0)
		{
			clock_t elap = clock() - start_time;
			double mlup = sizedom*(double)Time / (double)elap;
			double redRho = sumRo.reduceSingle();
			double redUx = sumUx.reduceSingle() / redRho;
			double redTemp = vfd.sumTemp.reduceSingle() / redRho;
			double redKap = kOmegaClass.sumKappa.reduceSingle() / redRho;
			double redOme = kOmegaClass.sumOmega.reduceSingle() / redRho;
			double redKapMin = kOmegaClass.minKappa.reduceSingle();
			double redOmeMin = kOmegaClass.minOmega.reduceSingle();
			printf("Time = %d, Umean = %g, K = %g, Omega = %g, Re = %g, "\
				"Temp = %g, MLUPS = %g\n", Time, redUx, redKap, redOme, 
				redUx * p.Pipe_radius / MuVal, redTemp, mlup);

			if (redOmeMin < 0.)
				printf("warning: omega contains negative values of %g\n", redOmeMin);
			if (redKapMin < 0.)
				printf("warning: omega contains negative values of %g\n", redKapMin);
		}
	}

	save2file();
	vfd.save2file();

#ifdef CREATE_BIN_FILES
	saveRestartFiles();
	vfd.saveRestartFiles();
	p.RenameFile("load\\temp.bin", "load\\LBFD\\temp.bin");
	p.RenameFile("load\\lbf.bin", "load\\LBFD\\lbf.bin");
	if (kOmegaClass.kOmegaSolverFlag)
	{
		p.RenameFile("load\\lbkappa.bin", "load\\LBFD\\lbkappa.bin");
		p.RenameFile("load\\lbomega.bin", "load\\LBFD\\lbomega.bin");
	}
#endif
}


void clVariablesLB::ini()
{
	// Keep IBB arrays for now
	//IBB_array_len = vls.fullsize_ibb_arrays;
	//IBB_loc.allocate(IBB_array_len);
	//IBB_loc.fill({ { 0, 0 } });
	//IBB_coeff.allocate(IBB_array_len);
	//IBB_coeff.fill({ { 0., 0. } });

	alter = 0;
	Save_loc = 0;

	sourceGenerator::SourceInstance()->addFile2Kernel("lbKernels.cl");

	iniNodeType();

	if (!fLoadedFlag)
	{
		if (iniPoiseuille)
		{
			double Umaxxx = Re * MuVal / p.Pipe_radius * 3. / 2.;
			iniDists(Umaxxx);
		}
		else
		{
			iniDists();
		}
	}

	allocateBuffers();

	sumUx.ini(Ux_array, "redUx");
	sumUy.ini(Uy_array, "redUy");
	sumRo.ini(Ro_array, "redRo");

	setSourceDefines();
	
	if (kOmegaClass.kOmegaSolverFlag)
		kOmegaClass.ini();

	std::function<void(void)> createKerPtr = std::bind(&clVariablesLB::createKernels, this);
	std::function<void(void)> setArgsPtr = std::bind(&clVariablesLB::setKernelArgs, this);
	sourceGenerator::SourceInstance()->addIniFunction(createKerPtr, setArgsPtr);

	// For iterative flowrate method, doesnt hurt to just keep it here
	Next_Update_Time = p.Time + TIME_BETWEEN_INTERVALS;

	if(saveMacroStart)
		save2file();
}



void clVariablesLB::iniDists()
{
	double Rho_inlet = RhoVal;
	double Feq_temp[9];
	calcFeqVals(Feq_temp, Rho_inlet, { { 0, 0 } });
	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (vls.M(i, j) == LB_SOLID)
				continue;
			Ro_array(i, j) = RhoVal;
			for (int k = 0; k < 9; k++)
			{
				FA(i, j, k) = Feq_temp[k];
			}
		}
	}
}

void clVariablesLB::iniDists(double Umaxval)
{
	double kappa = 0.41;
	double Cp = 5.5;

	double Relam = Re;
	double Returb = sqrt(Re * 3);
	double Upara_max = Relam * MuVal / p.Pipe_radius * (double)1.5;

	double Utau_m = (double)Returb * (double)1.25 * (double)MuVal / (double)p.Pipe_radius;
	double Utau = (double)Returb * (double)MuVal / (double)p.Pipe_radius;

	double yplus_max = p.Pipe_radius * Utau / MuVal;
	double xplus_max = (double)p.nX * Utau / MuVal;
	double Ubar = Upara_max *2. / 3.;
	double duPlus = Ubar * kOmegaClass.perturbDUPlus / Utau;
	double betaPlus = (double)(2.*PI_NUMBER / (kOmegaClass.perturbAlpha));
	double alphaPlus = (double)(2.*PI_NUMBER / (kOmegaClass.perturbBeta));
	double epsilonPlus = Ubar / kOmegaClass.perturbEpsilon;


	for (int i = 0; i < p.nX; i++)
	{
		int j = 0;
		while (vls.M(i, j) == LB_SOLID)
			j++;

		double Yval = vls.dXArrCur(i, j, 3);
		double Feq_temp[9];

		while (vls.M(i, j) != LB_SOLID)
		{
			double deviation = boxMuller(1., 0.2);

			double ypara = (Yval < p.Pipe_radius) ? p.Pipe_radius - Yval : Yval - p.Pipe_radius;
			double yturb = (Yval < p.Pipe_radius) ? Yval : 2.*p.Pipe_radius - Yval;
			double yplus = yturb * Utau / MuVal;
			double xplus = (double)i * Utau / MuVal;
			double Uplus = (yplus < 11.44532166) ? yplus : log(yplus) / 0.41 + 5.5;

			double Utemp_loc;
			if (!kOmegaClass.iniTurbVelocity)
				Utemp_loc = Upara_max*(1. - ypara*ypara / p.Pipe_radius / p.Pipe_radius);
			else
				Utemp_loc = Uplus*Utau;

			double Vtemp_loc = 0.;

			if (kOmegaClass.perturbVelField)
			{
				Utemp_loc += (Utau * duPlus / 2.) * (yplus / 40.) * 
					exp(-kOmegaClass.perturbSigma*yplus*yplus + 0.5) *
					cos(betaPlus*xplus) * deviation;
				Vtemp_loc = yplus * epsilonPlus * exp(-kOmegaClass.perturbSigma*
					yplus*yplus) * sin(alphaPlus*xplus) * deviation;
			}

			calcFeqVals(Feq_temp, RhoVal, { { Utemp_loc, Vtemp_loc } });
			Ro_array(i, j) = RhoVal;
			Ux_array(i, j) = Utemp_loc;
			Uy_array(i, j) = Vtemp_loc;
			for (int k = 0; k < 9; k++)
			{
				FA(i, j, k) = Feq_temp[k];
			}
			Yval += 1.;
			j++;
		}
	}
}



void clVariablesLB::iniIBB()
{
	vfl.update_LB_kernel[1].reset_global_size(vls.length_ibb);
	vfl.update_LB_kernel[1].set_argument(7, &vls.num_el);
	vfl.update_LB_kernel[1].call_kernel();
}

void clVariablesLB::iniInletVels()
{
	Inlet_Vel.zeros(p.Channel_Height);
	double Mu = (1. / 3.)*(vlb.tau - 0.5);
	double Umaxval = (vlb.Re * Mu / p.Pipe_radius)*3. / 2.;
	double Y_b = vls.C[0].y;
	double Y_center = Y_b + p.Pipe_radius;
	double shift = ceil(Y_b);
	for (int j = 0; j < p.Channel_Height; j++)
	{
		int loc = j;
		double yloc = Y_center - ((double)j + shift);
		double Uvel = Umaxval * (1. - yloc * yloc / p.Pipe_radius / p.Pipe_radius);
		Inlet_Vel(j) = Uvel*6.;
	}
}

// TODO: look into seeing if this can replace M either in its entirety, 
// or atleast in certain parts
void clVariablesLB::iniNodeType()
{
	NodeType.fill((int)SOLID_NODE);

	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (vls.M(i, j) == LB_SOLID)
				continue;

			int norient = 0;
			for (int m = 1; m < 9; m++)
			{
				if (j == 0 || j == p.nY - 1)
				{
					norient |= boundsArr[m];
					norient |= boundsArrT1[m];
					continue;
				}
				int iis = MOD(i + vlb.Cx[m], p.nX);
				int jjs = j + Cy[m];
				if (vls.M(iis, jjs) == LB_SOLID)
				{
					norient |= boundsArr[m];
					if (vls.dXArrCur(i, j, m - 1) <= 0.5)
						norient |= boundsArrT1[m];
					else if (vls.dXArrCur(i, j, m - 1) > 0.5)
						norient |= boundsArrT2[m];
					else
						printf("dir at (%d,%d) points to solid node, but dx = 1\n", i, j);
				}
			}
			norient |= FLUID_NODE;
			NodeType(i, j) = norient;
		}
	}
	NodeType.allocate_buffer_w_copy();
	NodeType.savetxt("NodeType");
}


void clVariablesLB::loadParams()
{
	kOmegaClass.kOmegaSolverFlag = p.getParameter(
		"Use Turb Model", USE_TURBULENCE_MODEL);

	testRestartRun();
	
	int MuTauFl = p.getParameter("Mu", "Tau", MuVal, tau);
	if(MuTauFl == 1)
	{
		tauDefined = false;
		tau = 1. / (MuVal * 3. + 0.5);
	}
	else if(MuTauFl == 2)
	{
		tauDefined = true;
		MuVal = ((1./tau - 0.5)/3.);
	}
	else
	{
		tauDefined = true;
		tau = LB_TAU_NUMBER;
		MuVal = MuVal;
	}

	int ReFvalFl = p.getParameter("Re", "Fval", Re, Fval);
	if(ReFvalFl == 1)
	{
		reDefined = true;
		Fval = (Re * 3. * MuVal * MuVal / p.Pipe_radius / p.Pipe_radius / p.Pipe_radius);
	}
	else if(ReFvalFl == 2)
	{
		reDefined = false;
		Re = Fval * p.Pipe_radius * p.Pipe_radius * p.Pipe_radius / 3. / MuVal / MuVal;
	}
	else
	{
		reDefined = true;
		Re = RE_NUMBER;
		Fval = (Re * 3. * MuVal * MuVal / p.Pipe_radius / p.Pipe_radius / p.Pipe_radius);
	}

	RhoVal = p.getParameter("Rho", LB_RHO_NUMBER);
	Ux_inlet = p.getParameter("Ux inlet", UX_INLET);
	nuAir = p.getParameter("Nu air", VISC_AIR);
	rhoAir = p.getParameter("Rho air", DENSITY_AIR);

	p.DELTA_T = (nuAir * p.DELTA_L * p.DELTA_L / MuVal);
	p.DELTA_M = (RhoVal*p.DELTA_L*p.DELTA_L*p.DELTA_L / rhoAir);
	p.DELTA_F = (p.DELTA_M*p.DELTA_L / p.DELTA_T / p.DELTA_T);
	p.DELTA_P = (p.DELTA_F*p.DELTA_L);
	
	saveMacroStart = p.getParameter("Save Macros On Start", FLUID_SAVE_MACROS_ON_START);
	iniPoiseuille = p.getParameter("Ini Poiseuille", INI_POISEUILLE);
	runLBFDFirst = p.getParameter("Run LBFD First", RUN_LB_FD_FIRST);

	tlbmIniDumpStep = p.getParameter("Ini Dump Step", TLBM_INI_DUMP_STEP);
	tlbmNumIniSteps_LB = p.getParameter("Ini LB Steps", NUM_LB_INI_STEPS);
	tlbmNumIniSteps_TLBM = p.getParameter("Ini TLBM Steps", NUM_LBFD_INI_STEPS);

	numIntervalsPerAvg = p.getParameter("Num Intervals Per Avg", NUM_INTERVALS_PER_AVG);
	timeBetweenIntervals = p.getParameter("Time Between Intervals", TIME_BETWEEN_INTERVALS);
	pauseBtwAdj = p.getParameter("Pause Between Adjust", PAUSE_BETWEEN_ADJ);
	flowrateMaxPercentDiff = p.getParameter("Flowrate Max Percent Diff", MAX_FLOWRATE_PERCENT_DIFF);
	
	kOmegaClass.loadParams();

	// Used for interatively finding dP necessary for a set flowrate
	// The actual method has not been implemented yet
	Fval_prev = 0.;
	Um_keep = 0.;
	Um_prev = 0.;
	Fval_orig = 0.;
	Fval_Diff = 0.;
}

void clVariablesLB::save2file()
{
	Ux_array.save_txt_from_device("lbux");
	Uy_array.save_txt_from_device("lbuy");
	Ro_array.save_txt_from_device("lbro");
	if (kOmegaClass.kOmegaSolverFlag)
	{
		kOmegaClass.save2file();
	}
}



void clVariablesLB::saveDebug(int saveFl)
{
	kOmegaClass.saveDebug(saveFl);
}

void clVariablesLB::saveDistributions()
{
	if (alter)
	{
		FB.save_txt_from_device_as_multi2D("lbf");
	}
	else
	{
		FA.save_txt_from_device_as_multi2D("lbf");
	}
}

void clVariablesLB::saveParams()
{
	p.setParameter("ParamsName", fluidParamNum);

	if (p.Time > 0)
		p.setParameter("Restart Run", true);
	else
		p.setParameter("Restart Run", false);

	if (tauDefined)
		p.setParameter("Tau", tau);
	else
		p.setParameter("Mu", MuVal);

	if (reDefined)
		p.setParameter("Re", Re);
	else
		p.setParameter("Fval", Fval);

	p.setParameter("Rho", RhoVal);
	p.setParameter("Ux inlet", Ux_inlet);
	p.setParameter("Nu air", nuAir);
	p.setParameter("Rho air", rhoAir);

	p.setParameter("Num Intervals Per Avg", numIntervalsPerAvg);
	p.setParameter("Time Between Intervals", timeBetweenIntervals);
	p.setParameter("Pause Between Adjust", pauseBtwAdj);
	p.setParameter("Flowrate Max Percent Diff", flowrateMaxPercentDiff);


	// The values of these shouldnt matter during a restart
	//p.setParameter("Ini Turb Velocity", false);
	//p.setParameter("Save Macros On Start", false);
	//p.setParameter("Perturb Vel Field", false);
	//p.setParameter("Ini Poiseuille", false);
	//p.setParameter("Run LBFD First", false);
	//p.setParameter("Ini Dump Step", tlbmIniDumpStep);
	//p.setParameter("Ini LB Steps", tlbmNumIniSteps_LB);
	//p.setParameter("Ini TLBM Steps", tlbmNumIniSteps_TLBM);

	kOmegaClass.saveParams();
}



void clVariablesLB::saveRestartFiles()
{
	if (alter)
	{
		FB.save_bin_from_device("lbf");
	}
	else
	{
		FA.save_bin_from_device("lbf");
	}
	kOmegaClass.saveRestartFiles();
}


void clVariablesLB::saveTimeData()
{
	// TODO: move functions for saving time data for fluid
	// into this function from clProblem.

}


void clVariablesLB::setKernelArgs()
{
#ifdef DEBUG_TURBARR
	koDbgArr1.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 5, 5);
	koDbgArr1.allocate_buffer_w_copy();
	
	koDbgArr2.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 22, 22);
	koDbgArr2.allocate_buffer_w_copy();
#endif
	
	cl_int ind = 0;
	int zer = 0; 

	Collision_kernel.set_argument(ind++, Ro_array.get_buf_add());
	Collision_kernel.set_argument(ind++, Ux_array.get_buf_add());
	Collision_kernel.set_argument(ind++, Uy_array.get_buf_add());
	Collision_kernel.set_separate_arguments(ind++, FA.get_buf_add(), FB.get_buf_add());
	Collision_kernel.set_separate_arguments(ind++, FB.get_buf_add(), FA.get_buf_add());
	Collision_kernel.set_argument(ind++, NodeType.get_buf_add());
	Collision_kernel.set_argument(ind++, vls.dXArrCur.get_buf_add());
	Collision_kernel.set_argument(ind++, kOmegaClass.WallD.get_buf_add());
	Collision_kernel.set_argument(ind++, kOmegaClass.Nut_array.get_buf_add());
	Collision_kernel.set_argument(ind++, kOmegaClass.Kappa.get_add_Macro());
	Collision_kernel.set_argument(ind++, kOmegaClass.Omega.get_add_Macro());
	Collision_kernel.set_argument(ind++, kOmegaClass.Sxy_array.get_buf_add());
	Collision_kernel.setOptionInd(ind);
	Collision_kernel.set_argument(ind++, &zer);

	
	ind = 0;
	calcRoJKernel.set_argument(ind++, Ro_array.get_buf_add());
	calcRoJKernel.set_argument(ind++, Ux_array.get_buf_add());
	calcRoJKernel.set_argument(ind++, Uy_array.get_buf_add());
	calcRoJKernel.set_argument(ind++, FA.get_buf_add());
	calcRoJKernel.set_argument(ind++, NodeType.get_buf_add());

	if (fLoadedFlag)
	{
		calcRoJKernel.call_kernel();
		clFinish(LBQUEUE);
	}


	kOmegaClass.setKernelArgs();


	//ind = 0;
	//IBB_kernel_Fluid.set_argument(ind++, IBB_loc.get_buf_add());
	//IBB_kernel_Fluid.set_argument(ind++, IBB_coeff.get_buf_add());
	//IBB_kernel_Fluid.set_argument(ind++, FB.get_buf_add());
	//IBB_kernel_Fluid.set_argument(ind++, &vls.num_el);

//	ind = 0;
//	Update_output_kernel[0].set_argument(ind++, Ro_array.get_buf_add());
//	Update_output_kernel[0].set_argument(ind++, J_array.get_buf_add());
//	Update_output_kernel[0].set_argument(ind++, reduceBoth.get_buf_add());
//	Update_output_kernel[0].set_argument(ind++, 2 * WORKGROUPSIZE_RED * NULL);
}

#define SETSOURCEDEFINE SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr()

void clVariablesLB::setSourceDefines()
{
#ifdef INLET_OUTLET_BC
	SETSOURCEDEFINE, "INLET_OUTLET_BC");
#endif

	SETSOURCEDEFINE, "FTERM_VAL", Fval);
	SETSOURCEDEFINE, "INI_MASS", (double)(p.nX*p.Channel_Height));
	SETSOURCEDEFINE, "tau0", tau);
	SETSOURCEDEFINE, "UX_INLET_VAL", Ux_inlet);

}

#undef SETSOURCEDEFINE

void clVariablesLB::Solve()
{
	//cl_event col_evt;
	Collision_kernel.call_kernel(NULL, 0, NULL, NULL);
	clFinish(LBQUEUE);

//	clReleaseEvent(col_evt);
	alter ^= 1;
}


bool clVariablesLB::testRestartRun()
{
	allocateArrays();
	restartRunFlag = p.getParameter("Restart Run", false);

	bool fBool = FA.load("load" SLASH "lbf") && FB.load("load" SLASH "lbf");
	if (fBool)
	{
		fLoadedFlag = true;
	}
	else
	{
		fLoadedFlag = false;
		restartRunFlag = false;
	}

	restartRunFlag &= kOmegaClass.testRestartRun();

	return restartRunFlag;
}


// Not implemented
void clVariablesLB::updateIBBArrays()
{
	//int fullsize_ibb_old = vls.fullsize_ibb_arrays;
	//vls.fullsize_Bnodes_old = vls.fullsize_Bnodes;
	//vls.update_IBB_arrays();

	//if (fullsize_ibb_old <  vls.fullsize_ibb_arrays)
	//{
	//	int newsize = vls.length_ibb * 2;
	//	vls.fullsize_ibb_arrays = newsize;

	//	vls.ii0_array.reallocate(newsize, p.context, p.IOqueue);
	//	vls.iis1_array.reallocate(newsize, p.context, p.IOqueue);
	//	vls.dir_array.reallocate(newsize, p.context, p.IOqueue);
	//	vls.D_array.reallocate(newsize, p.context, p.IOqueue);

	//	IBB_loc.reallocate_device_only(newsize, p.context, p.IOqueue);
	//	IBB_coeff.reallocate_device_only(newsize, p.context, p.IOqueue);

	//	clFinish(p.IOqueue);

	//	vfl.update_LB_kernel[1].set_argument(0, IBB_loc.get_buf_add());
	//	vfl.update_LB_kernel[1].set_argument(1, IBB_coeff.get_buf_add());
	//	vfl.update_LB_kernel[1].set_argument(2, vls.ii0_array.get_buf_add());
	//	vfl.update_LB_kernel[1].set_argument(3, vls.iis1_array.get_buf_add());
	//	vfl.update_LB_kernel[1].set_argument(4, vls.dir_array.get_buf_add());
	//	vfl.update_LB_kernel[1].set_argument(6, vls.D_array.get_buf_add());

	//	vtr.Update_SS_kernel[0].set_argument(2, vls.dir_array.get_buf_add());
	//	vtr.Update_SS_kernel[0].set_argument(3, vls.ii0_array.get_buf_add());

	//	IBB_kernel_Fluid[0].set_argument(0, IBB_loc.get_buf_add());
	//	IBB_kernel_Fluid[1].set_argument(0, IBB_loc.get_buf_add());

	//	IBB_kernel_Fluid[0].set_argument(1, IBB_coeff.get_buf_add());
	//	IBB_kernel_Fluid[1].set_argument(1, IBB_coeff.get_buf_add());
	//}

	//vls.ii0_array.copy_to_buffer(p.IOqueue, CL_FALSE, vls.num_el);
	//vls.iis1_array.copy_to_buffer(p.IOqueue, CL_FALSE, vls.num_el);
	//vls.dir_array.copy_to_buffer(p.IOqueue, CL_FALSE, vls.num_el);
	//vls.D_array.copy_to_buffer(p.IOqueue, CL_FALSE, vls.num_el);

	//IBB_kernel_Fluid[0].reset_global_size(vls.length_ibb);
	//IBB_kernel_Fluid[1].reset_global_size(vls.length_ibb);
	//vfl.update_LB_kernel[1].reset_global_size(vls.length_ibb);

	//IBB_kernel_Fluid[0].set_argument(3, &vls.num_el);
	//IBB_kernel_Fluid[1].set_argument(3, &vls.num_el);

	//vfl.update_LB_kernel[1].set_argument(7, &vls.num_el);
	//clFinish(p.IOqueue);
	//vfl.update_LB_kernel[1].call_kernel();
}



void clVariablesLB::updateTimeData()
{
	//Update_output_kernel[1].set_argument(4, (void*)&Save_loc);
	//Update_output_kernel[0].call_kernel();
	//Update_output_kernel[1].call_kernel();
	//vls.Masses.read_from_buffer(p.LBqueue);
	//Save_loc++;
	//clFlush(p.LBqueue);
}



///////////////////////////////////////////////////
////      Static Variable Initialization       ////
///////////////////////////////////////////////////
const int clVariablesLB::rev[9] = { 0, 2, 1, 4, 3, 6, 5, 8, 7 };

const int clVariablesLB::Cx[9] = { 0, 1, -1, 0, 0, 1, -1, 1, -1 };
const int clVariablesLB::Cy[9] = { 0, 0, 0, 1, -1, 1, -1, -1, 1 };

const cl_int2 clVariablesLB::Cxy[9] = { { { 0, 0 } }, { { 1, 0 } }, { { -1, 0 } }, { { 0, 1 } }, { { 0, -1 } },
{ { 1, 1 } }, { { -1, -1 } }, { { 1, -1 } }, { { -1, 1 } } };

const cl_double2 clVariablesLB::Cxy_double[9] = { { { 0., 0. } }, { { 1., 0. } }, { { -1., 0. } }, { { 0., 1. } }, { { 0., -1. } }, { { 1., 1. } }, { { -1., -1. } }, { { 1., -1. } }, { { -1., 1. } } };

const double clVariablesLB::Weigh[9] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };
const double clVariablesLB::Cs2 = 1. / 3.;


const int clVariablesLB::boundsArr[LB_NUMBER_CONNECTIONS] = { C_BOUND, E_BOUND, W_BOUND, N_BOUND, S_BOUND, NE_BOUND, SW_BOUND, SE_BOUND, NW_BOUND };
const int clVariablesLB::boundsArrT1[LB_NUMBER_CONNECTIONS] = { C_BOUND, E_BOUND_T1, W_BOUND_T1, N_BOUND_T1, S_BOUND_T1, NE_BOUND_T1, SW_BOUND_T1, SE_BOUND_T1, NW_BOUND_T1 };
const int clVariablesLB::boundsArrT2[LB_NUMBER_CONNECTIONS] = { C_BOUND, E_BOUND_T2, W_BOUND_T2, N_BOUND_T2, S_BOUND_T2, NE_BOUND_T2, SW_BOUND_T2, SE_BOUND_T2, NW_BOUND_T2 };



