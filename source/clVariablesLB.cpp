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



// TODO: Check for bin files in load\\LBFD and make sure not to 
//		run getInitialFlow and/or getInitialFlowAndTemp


#include "HelperFuncs.h"
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
		collisionKernel.create_kernel(GetSourceProgram, LBQUEUE_REF,
			"LB_collision_SRT_Fluid_w_kOmega", "LB_collision_SRT_Fluid_w_kOmega");
		kOmegaClass.createKernels();
	}
	else
	{
		collisionKernel.create_kernel(GetSourceProgram, LBQUEUE_REF,
			"LB_collision_SRT_Fluid", "LB_collision_SRT_Fluid");
	}

	int GlobalX = p.XsizeFull;
	int GlobalY = p.nY;

	collisionKernel.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);
	
	calcRoJKernel.create_kernel(GetSourceProgram, LBQUEUE_REF, "Calc_Ro_J");
	calcRoJKernel.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);


	if (kOmegaClass.kOmegaSolverFlag)
		kOmegaClass.createKernels();

#ifndef IN_KERNEL_IBB
	ibbKernel.create_kernel(GetSourceProgram, LBQUEUE_REF,
		"lbIBB", "lbIBB");
	ibbKernel.set_size(vls.ibbArr.curSize(), WORKGROUPSIZE_IBB);
#endif
}

void clVariablesLB::freeHostArrays()
{
	FA.FreeHost();
	FB.FreeHost();

	if (kOmegaClass.kOmegaSolverFlag)
		kOmegaClass.freeHostArrays();
}


void clVariablesLB::getInitialFlow()
{
	int Time = 0;
	int option = (OPTION_SAVE_MACRO_FIELDS);
	collisionKernel.setOption(&option);
	int NextDumpSteptime = tlbmIniDumpStep - 1;
	int dump_step_interval = NextDumpSteptime + 1;
	int stopTime = tlbmNumIniSteps_LB;

	clock_t start_time = clock();
	double sizedom = (double)p.nX * (double)p.Channel_Height * (double)CLOCKS_PER_SEC / 1.0e6;

#ifdef _DEBUG
	//vls.saveDebug();
	//vlb.saveDebug();
#endif


	while (Time < stopTime && !fLoadedFlag)
	{
		Solve();

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
	MakeDir(NewDir);
	saveRestartFiles();
	RenameFile("load\\lbf.bin", "load\\LBFD\\lbf.bin");
	if (kOmegaClass.kOmegaSolverFlag)
	{
		RenameFile("load\\lbkappa.bin", "load\\LBFD\\lbkappa.bin");
		RenameFile("load\\lbomega.bin", "load\\LBFD\\lbomega.bin");
	}
#endif
}

void clVariablesLB::getIntialFlowAndTemp()
{
	// Not currently working
	//vfd.solveSSThermal();
	
	int option = (OPTION_SAVE_MACRO_FIELDS);
	collisionKernel.setOption(&option);
	int NextDumpSteptime = tlbmIniDumpStep - 1;
	int dump_step_interval = NextDumpSteptime + 1;
	
	int Time = 0;
	int stopTime = tlbmNumIniSteps_TLBM;

	clock_t start_time = clock();
	double sizedom = (double)p.nX * (double)p.Channel_Height * (double)CLOCKS_PER_SEC / 1.0e6;


	while (Time < stopTime)
	{
		Solve();
		
		vfd.Solve();

		Time++;
		if (Time > NextDumpSteptime)
		{
			save2file();
			vfd.save2file();
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
	std::string dirName = "load" SLASH "LBFD";
	MakeDir(dirName);

	saveRestartFiles();
	vfd.saveRestartFiles();
	RenameFile("load\\temp.bin","load\\LBFD\\temp.bin");
	RenameFile("load\\lbf.bin", "load\\LBFD\\lbf.bin");
	if (kOmegaClass.kOmegaSolverFlag)
	{
		RenameFile("load\\lbkappa.bin", "load\\LBFD\\lbkappa.bin");
		RenameFile("load\\lbomega.bin", "load\\LBFD\\lbomega.bin");
	}
#endif
}


void clVariablesLB::ini()
{
	alter = 0;

	sourceGenerator::SourceInstance()->addFile2Kernel("lbKernels.cl");

	if (!fLoadedFlag)
	{
		if (loadVelTxtFlag)
		{
			vLoadedFlag =	Ux_array.loadtxt("restart_files\\lbux") && 
							Uy_array.loadtxt("restart_files\\lbuy") && 
							Ro_array.loadtxt(("restart_files\\lbro"));
		}

		if (vLoadedFlag)
		{
			iniDistsFromVel();
		}
		else if (iniPoiseuille)
		{
			double Umaxxx = Re * MuVal / p.Pipe_radius * 3. / 2. * geomPenaltyFactor;
			iniDists(Umaxxx);
		}
		else
		{
			iniDists();
		}
	}
	else
	{
		iniRhoUFromDists();
	}

	allocateBuffers();

	sumUx.ini(Ux_array, restartRunFlag, "redUx");
	sumUy.ini(Uy_array, restartRunFlag, "redUy");
	sumRo.ini(Ro_array, restartRunFlag, "redRo");


	kOmegaClass.ini();


	setSourceDefines();
	

	std::function<void(void)> createKerPtr = std::bind(&clVariablesLB::createKernels, this);
	std::function<void(void)> setArgsPtr = std::bind(&clVariablesLB::setKernelArgs, this);
	sourceGenerator::SourceInstance()->addIniFunction(createKerPtr, setArgsPtr);

	// For iterative flowrate method, doesnt hurt to just keep it here
	Next_Update_Time = p.Time + timeBetweenIntervals;

	if(saveMacroStart && !restartRunFlag)
		save2file();

	LOGMESSAGE("vlb initialized");
}

void clVariablesLB::iniDistsFromVel()
{
	double Rho_inlet = RhoVal;
	double Feq_temp[9];

	for (int i = 0; i < p.XsizeFull; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (i >= p.nX)
			{
				for (int k = 0; k < 9; k++)
				{
					FA(i, j, k) = Weigh[k];
				}
				continue;
			}

			
			// just in case these arent 0.
			if (vls.nType(i, j) & M_SOLID_NODE)
			{
				Ro_array(i,j) = 0.;
				Ux_array(i, j) = 0.;
				Uy_array(i, j) = 0.;
			}

			cl_double2 uvals = { { Ux_array(i,j), Uy_array(i,j) } };
			double rhoval = Ro_array(i, j);
			
			calcFeqVals(Feq_temp, rhoval, uvals);

			for (int k = 0; k < 9; k++)
			{
				FA(i, j, k) = Feq_temp[k];
			}
		}
	}
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
			if (vls.nType(i, j) & M_SOLID_NODE)
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
		while (vls.nType(i, j) & M_SOLID_NODE)
			j++;

		double Yval = vls.dXArr(i, j, 3);
		double Feq_temp[9];

		while (vls.nType(i, j) & M_FLUID_NODE)
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

void clVariablesLB::iniRhoUFromDists()
{
	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (vls.nType(i, j) & M_SOLID_NODE)
			{
				Ux_array(i, j) = 0.;
				Uy_array(i, j) = 0.;
				Ro_array(i, j) = 0.;
				continue;
			}

			Ux_array(i, j) = FA(i, j, 1) - FA(i, j, 2) + FA(i, j, 5) -
				FA(i, j, 6) + FA(i, j, 7) - FA(i, j, 8);
			Uy_array(i, j) = FA(i, j, 3) - FA(i, j, 4) + FA(i, j, 5) -
				FA(i, j, 6) - FA(i, j, 7) + FA(i, j, 8);
			Ro_array(i, j) = FA(i,j,0) + FA(i, j, 1) + FA(i, j, 2) + FA(i, j, 3) +
				FA(i, j, 4) + FA(i, j, 5) + FA(i, j, 6) + FA(i, j, 7) + FA(i, j, 8);

			ERROR_CHECKING((Ro_array(i, j) < 0.9), "Rho value calculated from loaded "\
				"distibutions < 0.9, which\nmeans there either is error in saved distribution\n"\
				"or the file was read in incorrectly (possibly due to error\nin set domain size.", \
				ERROR_INITIALIZING_VLB);
		}

	}
}

void clVariablesLB::iniTimeData()
{
	// no need to initialize anything here, since
	// reduce kernels are used elsewhere and therefore
	// will be initialized in the ini function
}

void clVariablesLB::loadParams()
{
	kOmegaClass.kOmegaSolverFlag = p.getParameter(
		"Use Turb Model", USE_TURBULENCE_MODEL);

	saveAvgInfo = p.getParameter("Save Avg Info", LB_SAVE_AVG_INFO);

	lbOutSize = p.getParameter("lbOut Max Lines", OUTPUT_MAX_LINES);

	loadVelTxtFlag = p.getParameter("load Vel Text", LOAD_VEL_TEXT_NO_BIN);

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

	turbVelBCFlag = p.getParameter("Use Turb Vel BC", USE_TURB_VEL_BC);

	geomPenaltyFactor = p.getParameter("Geometric Penalty Factor", GEOMETRIC_PENALTY_FACTOR);

	kOmegaClass.loadParams();

	// Used for interatively finding dP necessary for a set flowrate
	// The actual method has not been implemented yet
	Fval_prev = 0.;
	Um_keep = 0.;
	Um_prev = 0.;
	Fval_orig = 0.;
	Fval_Diff = 0.;
}


void clVariablesLB::renameSaveFiles()
{
	return;
	Ux_array.RenameTxtFile();
	Uy_array.RenameTxtFile();
	Ro_array.RenameTxtFile();
	if (kOmegaClass.kOmegaSolverFlag)
	{
		kOmegaClass.renameSaveFiles();
	}
}

void clVariablesLB::save2file()
{
	Ux_array.save_txt_from_device();
	Uy_array.save_txt_from_device();
	Ro_array.save_txt_from_device();
	if (kOmegaClass.kOmegaSolverFlag)
	{
		kOmegaClass.save2file();
	}
}



void clVariablesLB::saveDebug(int saveFl)
{
	//saveDistributions();
	//IBB_coeff.save_txt_from_device();
	//IBB_loc.save_txt_from_device();
	save2file();	
	vls.nType.save_txt_from_device_short_to_int();
	vls.dXArr.save_txt_from_device_as_multi2D();
	if (kOmegaClass.kOmegaSolverFlag)
	{
		kOmegaClass.saveDebug(saveFl);
	}
}

void clVariablesLB::saveDistributions(bool saveOpposite)
{
	if (saveOpposite)
		alter ^= 1;
	if (alter)
	{
		FB.save_txt_from_device_as_multi2D("lbf");
	}
	else
	{
		FA.save_txt_from_device_as_multi2D("lbf");
	}
	if (saveOpposite)
		alter ^= 1;
}


void clVariablesLB::saveParams()
{
	p.setParameter("ParamsName", fluidParamNum);

	p.setParameter("Save Avg Info", saveAvgInfo);

	p.setParameter("lbOut Max Lines", lbOutSize);

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
	p.setParameter("Use Turb Vel BC", turbVelBCFlag);
	p.setParameter("load Vel Text", loadVelTxtFlag);

	p.setParameter("Geometric Penalty Factor", geomPenaltyFactor);

	if (kOmegaClass.kOmegaSolverFlag)
		kOmegaClass.saveParams();
	else
		p.setParameter("Use Turb Model", false);
	// The values of these shouldnt matter during a restart
	//p.setParameter("Ini Turb Velocity", false);
	//p.setParameter("Save Macros On Start", false);
	//p.setParameter("Perturb Vel Field", false);
	//p.setParameter("Ini Poiseuille", false);
	//p.setParameter("Run LBFD First", false);
	//p.setParameter("Ini Dump Step", tlbmIniDumpStep);
	//p.setParameter("Ini LB Steps", tlbmNumIniSteps_LB);
	//p.setParameter("Ini TLBM Steps", tlbmNumIniSteps_TLBM);


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

	if (kOmegaClass.kOmegaSolverFlag)
		kOmegaClass.saveRestartFiles();
}


void clVariablesLB::saveTimeData()
{
	// TODO: move functions for saving time data for fluid
	// into this function from clProblem.
	if (!saveAvgInfo)
		return;
	sumUx.appendToFileFromDevice();
	sumUy.appendToFileFromDevice();
	sumRo.appendToFileFromDevice();
	if (vfd.thermalSolverFlag)
	{
		vfd.sumTemp.appendToFileFromDevice();
	}
	if (kOmegaClass.kOmegaSolverFlag)
	{
		kOmegaClass.sumKappa.appendToFileFromDevice();
		kOmegaClass.sumOmega.appendToFileFromDevice();
	}
	// kOmegaClass has empty saveTimeData, so no need to call
}


void clVariablesLB::setKernelArgs()
{
#ifdef DEBUG_TURBARR
	kOmegaClass.koDbgArr1.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 5, 5);
	kOmegaClass.koDbgArr1.allocate_buffer_w_copy();
	
	kOmegaClass.koDbgArr2.zeros(p.nX, p.XsizeFull, p.nY, p.nY, 22, 22);
	kOmegaClass.koDbgArr2.allocate_buffer_w_copy();
#endif
	
	cl_int ind = 0;
	int zer = 0; 

	collisionKernel.set_argument(ind++, Ro_array.get_buf_add());
	collisionKernel.set_argument(ind++, Ux_array.get_buf_add());
	collisionKernel.set_argument(ind++, Uy_array.get_buf_add());
	collisionKernel.set_separate_arguments(ind++, FA.get_buf_add(), FB.get_buf_add());
	collisionKernel.set_separate_arguments(ind++, FB.get_buf_add(), FA.get_buf_add());
	collisionKernel.set_argument(ind++, vls.nType.get_buf_add());
#ifdef IN_KERNEL_IBB
	collisionKernel.set_argument(ind++, vls.dXArr.get_buf_add());
#endif
	if (kOmegaClass.kOmegaSolverFlag)
	{
		collisionKernel.set_argument(ind++, kOmegaClass.WallD.get_buf_add());
		collisionKernel.set_argument(ind++, kOmegaClass.Nut_array.get_buf_add());
		collisionKernel.set_argument(ind++, kOmegaClass.Kappa.get_add_Macro());
		collisionKernel.set_argument(ind++, kOmegaClass.Omega.get_add_Macro());
		collisionKernel.set_argument(ind++, kOmegaClass.Sxy_array.get_buf_add());
		if (turbVelBCFlag)
		{
			collisionKernel.set_argument(ind++, vls.ssArrIndMap.get_buf_add());
			collisionKernel.set_argument(ind++, vtr.wallShear.Tau.get_buf_add());
		}
	}
	collisionKernel.setOptionInd(ind);
	collisionKernel.set_argument(ind++, &zer);

	ind = 0;
	calcRoJKernel.set_argument(ind++, Ro_array.get_buf_add());
	calcRoJKernel.set_argument(ind++, Ux_array.get_buf_add());
	calcRoJKernel.set_argument(ind++, Uy_array.get_buf_add());
	calcRoJKernel.set_argument(ind++, FA.get_buf_add());
	calcRoJKernel.set_argument(ind++, vls.nType.get_buf_add());

	if (fLoadedFlag)
	{
		calcRoJKernel.call_kernel();
		clFinish(LBQUEUE);
	}

#ifndef IN_KERNEL_IBB
	ind = 0;
	ibbKernel.set_argument(0, vls.ibbArr.get_buf_add());
	ibbKernel.set_argument(1, vls.ibbDistArr.get_buf_add());
	ibbKernel.set_separate_arguments(2, FB.get_buf_add(),
		FA.get_buf_add());

	int ibbnumel_ = vls.ibbArrCurIndex(0);
	ibbKernel.set_argument(3, &ibbnumel_);
	ibbKernel.setOptionInd(3);
	ibbKernel.reset_global_size(ibbnumel_);
#endif

	if (kOmegaClass.kOmegaSolverFlag)
		kOmegaClass.setKernelArgs();
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
	if (turbVelBCFlag)
	{
		SETSOURCEDEFINE, "USE_TURB_VEL_BC");
	}
	// Must call kOmegaClass because source wont compile without
	// those defines set
	kOmegaClass.setSourceDefines();
}

#undef SETSOURCEDEFINE



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


void clVariablesLB::update()
{
	// vls updates IBB arrays, so only need to update wall Distances for kOmega

	if (kOmegaClass.kOmegaSolverFlag)
		kOmegaClass.updateWallDKernel.call_kernel(LBQUEUE_REF);


}

//void clVariablesLB::updateIBBArrays(bool reSizeFlag)
//{
//#ifndef IN_KERNEL_IBB
//	if (reSizeFlag)
//	{
//		ibbKernel.set_argument(0, vls.ibbArr.get_buf_add());
//	}
//	ibbKernel.setOption(vls.ibbArr.curSizeAdd());
//#endif
//}



void clVariablesLB::updateTimeData()
{
	if (!saveAvgInfo)
		return;
	sumUx.reduce(p.Time);
	sumUy.reduce(p.Time);
	sumRo.reduce(p.Time);
	if (vfd.thermalSolverFlag)
	{
		vfd.sumTemp.reduce(p.Time);
	}
	if (kOmegaClass.kOmegaSolverFlag)
	{
		kOmegaClass.sumKappa.reduce(p.Time);
		kOmegaClass.sumOmega.reduce(p.Time);
	}

	// Function empty in kOmegaClass so no need to call it
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
#ifdef IN_KERNEL_IBB
const int clVariablesLB::boundsArrT1[LB_NUMBER_CONNECTIONS] = { C_BOUND, E_BOUND_T1, W_BOUND_T1, N_BOUND_T1, S_BOUND_T1, NE_BOUND_T1, SW_BOUND_T1, SE_BOUND_T1, NW_BOUND_T1 };
const int clVariablesLB::boundsArrT2[LB_NUMBER_CONNECTIONS] = { C_BOUND, E_BOUND_T2, W_BOUND_T2, N_BOUND_T2, S_BOUND_T2, NE_BOUND_T2, SW_BOUND_T2, SE_BOUND_T2, NW_BOUND_T2 };
#endif

