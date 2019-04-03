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
	if (kOmegaSolverFlag)
	{
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

}

void clVariablesLB::calcDxTerms()
{
	/// Doing it this way will allow for more complicated corrections to be made
	/// without worrying about effects on speed.

	dXCoeffs.zeros(p.XsizeFull, p.nY, 12);

	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (vls.M(i, j) == LB_SOLID)
				continue;


			double dx_e = vls.dXArrCur(i, j, 0), dx_w = vls.dXArrCur(i, j, 1), dx = dx_e + dx_w;
			double dy_n = vls.dXArrCur(i, j, 2), dy_s = vls.dXArrCur(i, j, 3), dy = dy_n + dy_s;

			double Xe_coeff = dx_w / (dx_e*dx), Xw_coeff = -dx_e / (dx_w*dx), Xc_coeff = (dx_e - dx_w) / (dx_e*dx_w);
			double Yn_coeff = dy_s / (dy_n*dy), Ys_coeff = -dy_n / (dy_s*dy), Yc_coeff = (dy_n - dy_s) / (dy_n*dy_s);

			double Xe2_coeff = 2. / (dx_e*dx), Xw2_coeff = 2. / (dx_w*dx), Xc2_coeff = -2. / (dx_e*dx_w);
			double Yn2_coeff = 2. / (dy_n*dy), Ys2_coeff = 2. / (dy_s*dy), Yc2_coeff = -2. / (dy_n*dy_s);


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

	dXCoeffs.allocate_buffer_w_copy();
}

void clVariablesLB::calcFeqVals(double *Feq_vals, double rho_val, cl_double2 u_val)
{
	Feq_vals[0] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y - 2. * u_val.x*u_val.x / 3 - 2 * u_val.y*u_val.y / 3. + 4. / 9.);
	Feq_vals[1] = rho_val*(-u_val.x*u_val.x * u_val.y*u_val.y + 2. * u_val.x*u_val.x / 3 - u_val.x*u_val.y*u_val.y + 2. * u_val.x / 3. - u_val.y*u_val.y / 3. + 2. / 9.) / 2.;
	Feq_vals[2] = rho_val*(-u_val.x*u_val.x * u_val.y*u_val.y + 2. * u_val.x*u_val.x / 3 + u_val.x*u_val.y*u_val.y - 2. * u_val.x / 3. - u_val.y*u_val.y / 3. + 2. / 9.) / 2.;
	Feq_vals[3] = rho_val*(-u_val.x*u_val.x * u_val.y*u_val.y - u_val.x*u_val.x * u_val.y - u_val.x*u_val.x / 3. + 2. * u_val.y*u_val.y / 3. + 2. * u_val.y / 3. + 2. / 9.) / 2.;
	Feq_vals[4] = rho_val*(-u_val.x*u_val.x * u_val.y*u_val.y + u_val.x*u_val.x * u_val.y - u_val.x*u_val.x / 3. + 2. * u_val.y*u_val.y / 3. - 2. * u_val.y / 3. + 2. / 9.) / 2.;
	Feq_vals[5] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y + u_val.x*u_val.x * u_val.y + u_val.x*u_val.x / 3. + u_val.x*u_val.y*u_val.y + u_val.x*u_val.y + u_val.x / 3. + u_val.y*u_val.y / 3. + u_val.y / 3. + 1. / 9.) / 4.;
	Feq_vals[6] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y - u_val.x*u_val.x * u_val.y + u_val.x*u_val.x / 3. - u_val.x*u_val.y*u_val.y + u_val.x*u_val.y - u_val.x / 3. + u_val.y*u_val.y / 3. - u_val.y / 3. + 1. / 9.) / 4.;
	Feq_vals[7] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y - u_val.x*u_val.x * u_val.y + u_val.x*u_val.x / 3. + u_val.x*u_val.y*u_val.y - u_val.x*u_val.y + u_val.x / 3. + u_val.y*u_val.y / 3. - u_val.y / 3. + 1. / 9.) / 4.;
	Feq_vals[8] = rho_val*(u_val.x*u_val.x * u_val.y*u_val.y + u_val.x*u_val.x * u_val.y + u_val.x*u_val.x / 3. - u_val.x*u_val.y*u_val.y - u_val.x*u_val.y - u_val.x / 3. + u_val.y*u_val.y / 3. + u_val.y / 3. + 1. / 9.) / 4.;
}


// HERE
void clVariablesLB::calcNutArray()
{
	for (int i = 0; i < p.nX; i++)
	{
		for (int j = 0; j < p.nY; j++)
		{
			if (vls.M(i, j) != LB_SOLID)
				Nut_array(i, j) = Kappa(i, j) / Omega(i, j);
		}
	}
}

void clVariablesLB::createKernels()
{
	if (kOmegaSolverFlag)
	{
		Collision_kernel.create_kernel(GetSourceProgram, LBQUEUE_REF,
			"LB_collision_SRT_Fluid_w_kOmega", "LB_collision_SRT_Fluid_w_kOmega");
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

	if (vfd.thermalSolverFlag)
	{
		kOmegaUpdateDiffCoeffs.create_kernel(GetSourceProgram, LBQUEUE_REF,
			"Update_kOmega_Diffusivities_Alphat");
	}
	else
	{
		kOmegaUpdateDiffCoeffs.create_kernel(GetSourceProgram, LBQUEUE_REF,
			"Update_kOmega_Diffusivities");
	}

	kOmegaUpdateDiffCoeffs.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);
	
	kOmegaUpdateCoeffs.create_kernel(GetSourceProgram, LBQUEUE_REF, 
		"Update_kOmega_Coeffs_Implicit");
	kOmegaUpdateCoeffs.set_size(GlobalX, WORKGROUPSIZEX_LB, GlobalY, WORKGROUPSIZEY_LB);






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

void clVariablesLB::freeHostMem()
{
	// No need to worry about memory usage
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
		SolveLB();
		clFinish(LBQUEUE);

		if (kOmegaSolverFlag)
			SolveKOmega();

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

			double redKap = sumKappa.reduceSingle() / redRho;
			double redOme = sumOmega.reduceSingle() / redRho;
			double redKapMin = minKappa.reduceSingle();
			double redOmeMin = minOmega.reduceSingle();
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
	if (kOmegaSolverFlag)
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
		SolveLB();
		clFinish(LBQUEUE);
		if (kOmegaSolverFlag)
			SolveKOmega();
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
			double redKap = sumKappa.reduceSingle() / redRho;
			double redOme = sumOmega.reduceSingle() / redRho;
			double redKapMin = minKappa.reduceSingle();
			double redOmeMin = minOmega.reduceSingle();
			printf("Time = %d, Umean = %g, K = %g, Omega = %g, Re = %g, Temp = %g, MLUPS = %g\n",
				Time, redUx, redKap, redOme, redUx * p.Pipe_radius / MuVal, redTemp, mlup);

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
	if (kOmegaSolverFlag)
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

	allocateArrays();
	iniLBM();
	
	if (kOmegaSolverFlag)
		iniKOmega();

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
	double duPlus = Ubar * perturbDUPlus / Utau;
	double betaPlus = (double)(2.*PI_NUMBER / (perturbAlpha));
	double alphaPlus = (double)(2.*PI_NUMBER / (perturbBeta));
	double epsilonPlus = Ubar / perturbEpsilon;


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
			if (!iniTurbVelocity)
				Utemp_loc = Upara_max*(1. - ypara*ypara / p.Pipe_radius / p.Pipe_radius);
			else
				Utemp_loc = Uplus*Utau;

			double Vtemp_loc = 0.;

			if (perturbVelField)
			{
				Utemp_loc += (Utau * duPlus / 2.) * (yplus / 40.) * exp(-perturbSigma*
					yplus*yplus + 0.5) * cos(betaPlus*xplus) * deviation;
				Vtemp_loc = yplus * epsilonPlus * exp(-perturbSigma*
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

// TODO: implement wall function values for non-parallel channels
//			this will most likely require a calculation of wall function values
//			which can be done in the kernel
void clVariablesLB::iniKOmega()
{
	sourceGenerator::SourceInstance()->addFile2Kernel("kOmegaKernels.cl");
	/// Calc initial Re_dh
	
	// calculation of yplus_lam
	double ypl = 11.;
	for (int i = 0; i < 11; i++)
	{
		ypl = log(max(ROUGHNESS_FACTOR*ypl, 1)) / 0.41;
	}
	// specific to flow in straight channel
	// based on openfoam's implementation
	// see http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2016/FangqingLiu/openfoamFinal.pdf
	// and openfoam documentation
	
	ReTurbVal = sqrt(Re*3.);
	yPlusWallVal = ReTurbVal / p.Pipe_radius / 2.;
	UtauVal = ReTurbVal*MuVal / p.Pipe_radius;
	if (yPlusWallVal < ypl)
	{
		double Cf = (1. / ((yPlusWallVal + 11.)*(yPlusWallVal + 11.)) +
			2.*yPlusWallVal / 11. / 11. / 11. - 1. / 11. / 11.);
		kappaWallVal = UtauVal*UtauVal*(2400. / 1.9 / 1.9 * Cf);
	}
	else
	{
		// Nothing implemented for yPlus < yPlus_laminar
		exit(ERROR_INITIALIZING_KOMEGA);
	}
	nutWallVal = MuVal * (0.41 * yPlusWallVal / log(roughnessFactor*yPlusWallVal) - 1.);
	double omega_vis = 6.*MuVal / (0.075*yPlusWallVal*yPlusWallVal);
	double omega_log = UtauVal / (0.3*yPlusWallVal*0.41);
	omegaWallVal = sqrt(omega_vis*omega_vis + omega_log*omega_log);
	omegaWallVal = MAX(omegaWallVal, 0.);
	double umean_ini = Ux_array.reduce() / Ro_array.reduce();
	if (umean_ini == 0.)
	{
		umean_ini = Re*MuVal / p.Pipe_radius / 4.;
	}

	double Reini = p.Pipe_radius*Re / MuVal*4.;
	TurbIntensity = 0.16*pow(Reini,-1./8.);
	TurbLScale_inv = 1. / (0.036*p.Pipe_radius);
	
	double kIniVal = 1.5*(TurbIntensity*umean_ini)*(TurbIntensity*umean_ini);
	
	double omegaIniVal = sqrt(kIniVal) / KO_C_MU_0_25 * TurbLScale_inv;
	double turbViscIniVal = kIniVal / omegaIniVal;
	
	if (turbViscIniVal > 1.)
	{
		printf("Initial Turbulent viscosity too high\n");
		exit(ERROR_INITIALIZING_KOMEGA);
	}
	
	// Fill with initial values if not loaded from bin
	if (!kOmegaLoadedFlag)
	{
		Omega_array.fill(omegaIniVal);
		Kappa_array.fill(kIniVal);
	}

	KappaInds.ini(p.nX, p.XsizeFull, p.nY, p.nY, &vls.M);
	Kappa.CreateSolver(&Kappa_array, &KappaInds, LBQUEUE_REF, 
		kappaMaxIters, kappaMaxRelTol, kappaMaxAbsTol);
	
	OmegaInds.ini(p.nX, p.XsizeFull, p.nY, p.nY, &vls.M);
	Omega.CreateSolver(&Omega_array, &OmegaInds, LBQUEUE_REF,
		omegaMaxIters, omegaMaxRelTol, omegaMaxAbsTol);

	int i = 0;
	while (vls.M_o(0, i) == 0)
	{
		i++;
	}
	int wallindlow = i;
	while (vls.M_o(0, i) == 1)
	{
		i++;
	}
	
	std::vector<int> wallinds = { wallindlow, i-1 };
	Kappa.setInitialValueRows(kappaWallVal, wallinds); // set to inlet val
	Omega.setInitialValueRows(omegaWallVal, wallinds); // set to inlet val

	sumOmega.ini(*Omega.getMacroArray(), "redOmega");
	sumKappa.ini(*Kappa.getMacroArray(), "redKappa");
	minOmega.ini(*Omega.getMacroArray(), "minOmega");
	minKappa.ini(*Kappa.getMacroArray(), "minKappa");

	iniWallD();
	calcDxTerms();

	/////////////////////////  Parameters for kOmegaModel/////////////////////////
	/// Note, these do not have the necessary padding for reductions

	calcNutArray();
	Diff_Omega.fill(MuVal);
	Diff_K.fill(MuVal);

	Sxy_array.allocate_buffer_w_copy();
	Nut_array.allocate_buffer_w_copy();
	Diff_Omega.allocate_buffer_w_copy();
	Diff_K.allocate_buffer_w_copy();
	Fval_array.allocate_buffer_w_copy();
	dKdO_array.allocate_buffer_w_copy();

	setSourceDefinesKOmega();
}

void clVariablesLB::iniLBM()
{
	sourceGenerator::SourceInstance()->addFile2Kernel("lbKernels.cl");

	iniNodeType();

	if (!fLoadedFlag)
	{
		if (iniPoiseuille)
		{
			double Umaxxx = Re*MuVal / p.Pipe_radius * 3. / 2.;
			iniDists(Umaxxx);
		}
		else
		{
			iniDists();
		}
	}

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

	sumUx.ini(Ux_array, "redUx");
	sumUy.ini(Uy_array, "redUy");
	sumRo.ini(Ro_array, "redRo");

	setSourceDefinesLBM();
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


void clVariablesLB::iniWallD()
{
	double ystart = 0.;
	int j = 0;
	while (vls.M_o(0, j) == 0)
	{
		j++;
	}

	ystart = vls.dXArr(0, j, 3);
	int botind = j;
	double ycenter = (double)j + p.Pipe_radius - ystart;
	WallD.zeros(p.nX, p.XsizeFull, p.nY, p.nY);

	for (int i = 0; i < p.nX; i++)
	{
		for (j = 0; j < p.nY; j++)
		{
			if (vls.M_o(0, j) == 0)
				continue;

			double ytemp = (double)(j - botind) + ystart;

			if (ytemp > p.Pipe_radius)
			{
				ytemp = 2.*p.Pipe_radius - ytemp;
			}

			WallD(i, j) = ytemp;
		}
	}


	WallD.allocate_buffer_w_copy();
#ifdef DEBUG_TURBARR
	WallD.savetxt("WallD");
#endif
}


void clVariablesLB::loadParams()
{
	kOmegaSolverFlag = p.getParameter<bool>("Use Turb Model", USE_TURBULENCE_MODEL);

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

	int ReFvalFl = p.getParameter<double>("Re", "Fval", Re, Fval);
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

	RhoVal = p.getParameter<double>("Rho", LB_RHO_NUMBER);
	Ux_inlet = p.getParameter<double>("Ux inlet", UX_INLET);
	nuAir = p.getParameter<double>("Nu air", VISC_AIR);
	rhoAir = p.getParameter<double>("Rho air", DENSITY_AIR);

	p.DELTA_T = (nuAir * p.DELTA_L * p.DELTA_L / MuVal);
	p.DELTA_M = (RhoVal*p.DELTA_L*p.DELTA_L*p.DELTA_L / rhoAir);
	p.DELTA_F = (p.DELTA_M*p.DELTA_L / p.DELTA_T / p.DELTA_T);
	p.DELTA_P = (p.DELTA_F*p.DELTA_L);
	
	iniTurbVelocity = p.getParameter<bool>("Ini Turb Velocity", INI_TURB_PROFILE);
	saveMacroStart = p.getParameter<bool>("Save Macros On Start", FLUID_SAVE_MACROS_ON_START);
	perturbVelField = p.getParameter<bool>("Perturb Vel Field", PERTURB_VEL_FIELD);
	iniPoiseuille = p.getParameter<bool>("Ini Poiseuille", INI_POISEUILLE);
	runLBFDFirst = p.getParameter<bool>("Run LBFD First", RUN_LB_FD_FIRST);

	tlbmIniDumpStep = p.getParameter<unsigned int>("Ini Dump Step", TLBM_INI_DUMP_STEP);
	tlbmNumIniSteps_LB = p.getParameter<unsigned int>("Ini LB Steps", NUM_LB_INI_STEPS);
	tlbmNumIniSteps_TLBM = p.getParameter<unsigned int>("Ini TLBM Steps", NUM_LBFD_INI_STEPS);

	perturbDUPlus = p.getParameter<double>("DUplus Perturb Multiplier", DUPLUS_MULT);
	perturbEpsilon = p.getParameter<double>("Epsilon Perturb Multiplier", EPSILON_MULT);
	perturbBeta = p.getParameter<double>("Beta Perturb Multiplier", BETA_MULT);
	perturbAlpha = p.getParameter<double>("Alpha Perturb Multiplier", ALPHA_MULT);
	perturbSigma = p.getParameter<double>("Sigma Perturb Multiplier", SIGMA_MULT);
	
	roughnessFactor = p.getParameter<double>("Wall Roughness Factor", ROUGHNESS_FACTOR);
	
	numIntervalsPerAvg = p.getParameter<unsigned int>("Num Intervals Per Avg", NUM_INTERVALS_PER_AVG);
	timeBetweenIntervals = p.getParameter<unsigned int>("Time Between Intervals", TIME_BETWEEN_INTERVALS);
	pauseBtwAdj = p.getParameter<unsigned int>("Pause Between Adjust", PAUSE_BETWEEN_ADJ);
	flowrateMaxPercentDiff = p.getParameter<double>("Flowrate Max Percent Diff", MAX_FLOWRATE_PERCENT_DIFF);
	
	kappaMaxRelTol = p.getParameter<double>("Kappa Max Rel Tol", KOMEGA_MAX_REL_TOL);
	kappaMaxAbsTol = p.getParameter<double>("Kappa Max Abs Tol", KOMEGA_MAX_ABS_TOL);
	kappaMaxIters =  p.getParameter<int>("Kappa Max Iterations", KOMEGA_MAX_ITERS);

	omegaMaxRelTol = p.getParameter<double>("Omega Max Rel Tol", kappaMaxRelTol);
	omegaMaxAbsTol = p.getParameter<double>("Omega Max Abs Tol", kappaMaxAbsTol);
	omegaMaxIters = p.getParameter<int>("Omega Max Iterations", kappaMaxIters);


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
	if (kOmegaSolverFlag)
	{
		Kappa.savetxt_from_device();
		Omega.savetxt_from_device();
		//Nut_array.save_txt_from_device("lbnut");
		//Sxy_array.save_txt_from_device("lbtau");
	}
}



void clVariablesLB::saveDebug(int saveFl)
{
#ifdef DEBUG_TURBARR
	Array2Dd outtemp;
	outtemp.zeros(p.nX, p.nY);

	Ux_array.save_txt_from_device("lbux");
	Uy_array.save_txt_from_device("lbuy");
	Ro_array.save_txt_from_device("lbro");


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

	//Kappa.saveAxbCSR_from_device();
	//Omega.saveAxbCSR_from_device();

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
	p.setParameter("Use Turb Model", kOmegaSolverFlag);
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
	p.setParameter("Kappa Max Rel Tol", kappaMaxRelTol);
	p.setParameter("Kappa Max Abs Tol", kappaMaxAbsTol);
	p.setParameter("Kappa Max Iterations", kappaMaxIters);
	p.setParameter("Omega Max Rel Tol", omegaMaxRelTol);
	p.setParameter("Omega Max Abs Tol", omegaMaxAbsTol);
	p.setParameter("Omega Max Iterations", omegaMaxIters);



	// The values of these shouldnt matter, but
	// setting them to false just in case
	p.setParameter("Ini Turb Velocity", false);
	p.setParameter("Save Macros On Start", false);
	p.setParameter("Perturb Vel Field", false);
	p.setParameter("Ini Poiseuille", false);
	p.setParameter("Run LBFD First", false);
	p.setParameter("Ini Dump Step", tlbmIniDumpStep);
	p.setParameter("Ini LB Steps", tlbmNumIniSteps_LB);
	p.setParameter("Ini TLBM Steps", tlbmNumIniSteps_TLBM);
	p.setParameter("DUplus Perturb Multiplier", perturbDUPlus);
	p.setParameter("Epsilon Perturb Multiplier", perturbEpsilon);
	p.setParameter("Beta Perturb Multiplier", perturbBeta);
	p.setParameter("Alpha Perturb Multiplier", perturbAlpha);
	p.setParameter("Sigma Perturb Multiplier", perturbSigma);
	p.setParameter("Wall Roughness Factor", roughnessFactor);
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
	Kappa.saveCheckPoint();
	Omega.saveCheckPoint();
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
	Collision_kernel.set_argument(ind++, WallD.get_buf_add());
	Collision_kernel.set_argument(ind++, Nut_array.get_buf_add());
	Collision_kernel.set_argument(ind++, Kappa.get_add_Macro());
	Collision_kernel.set_argument(ind++, Omega.get_add_Macro());
	Collision_kernel.set_argument(ind++, Sxy_array.get_buf_add());
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



	ind = 0;
	kOmegaUpdateDiffCoeffs.set_argument(ind++, NodeType.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Kappa.get_add_Macro());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Omega.get_add_Macro());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Nut_array.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, WallD.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Diff_K.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Diff_Omega.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Fval_array.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, dKdO_array.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, Ro_array.get_buf_add());
	if (vfd.thermalSolverFlag)
		kOmegaUpdateDiffCoeffs.set_argument(ind++, vfd.Alphat.get_buf_add());
	kOmegaUpdateDiffCoeffs.set_argument(ind++, dXCoeffs.get_buf_add());
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
	kOmegaUpdateCoeffs.set_argument(ind++, Ro_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, Ux_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, Uy_array.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, Kappa.get_add_A());
	kOmegaUpdateCoeffs.set_argument(ind++, Omega.get_add_A());
	kOmegaUpdateCoeffs.set_argument(ind++, Kappa.get_add_b());
	kOmegaUpdateCoeffs.set_argument(ind++, Omega.get_add_b());
	kOmegaUpdateCoeffs.set_argument(ind++, dXCoeffs.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, WallD.get_buf_add());
	kOmegaUpdateCoeffs.set_argument(ind++, &p.dTlb);
	kOmegaUpdateCoeffs.set_argument(ind++, &TurbIntensity);
	kOmegaUpdateCoeffs.set_argument(ind++, &TurbLScale_inv);
#ifdef DEBUG_TURBARR
	kOmegaUpdateCoeffs.set_argument(ind++, koDbgArr2.get_buf_add());
#endif


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

void clVariablesLB::setSourceDefinesLBM()
{
#ifdef INLET_OUTLET_BC
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "INLET_OUTLET_BC");
#endif

	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "FTERM_VAL", Fval);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "INI_MASS", (double)(p.nX*p.Channel_Height));
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "tau0", tau);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "UX_INLET_VAL", Ux_inlet);
}

void clVariablesLB::setSourceDefinesKOmega()
{
#ifdef DEBUG_TURBARR
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "DEBUG_ARRAYS");
#endif
#ifdef DEBUG_TURBARR
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "DISABLE_TURBULENT_VISC");
#endif

	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "RE_TURBULENT", ReTurbVal);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "UTAU_VAL", UtauVal);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "YPLUS_WALL_NODE", yPlusWallVal);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "K_WALL_VALUE", kappaWallVal);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "NUT_WALL_VALUE", nutWallVal);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "OMEGA_WALL_VALUE", omegaWallVal);
}

void clVariablesLB::SolveLB()
{
	//cl_event col_evt;
	Collision_kernel.call_kernel(NULL, 0, NULL, NULL);
	clFinish(LBQUEUE);


//	clReleaseEvent(col_evt);
	alter ^= 1;
}


void clVariablesLB::SolveKOmega()
{
	kOmegaUpdateDiffCoeffs.call_kernel();
	kOmegaUpdateCoeffs.call_kernel();
	clFinish(LBQUEUE);
	Kappa.solve();
	Omega.solve();
}


void clVariablesLB::testRestartRun()
{
	allocateArrays();
	restartRunFlag = p.getParameter<bool>("Restart Run", false);

	BOOL fBool = FA.load("load" SLASH "lbf") && FB.load("load" SLASH "lbf");
	if (fBool)
	{
		fLoadedFlag = true;
	}
	else
	{
		fLoadedFlag = false;
		restartRunFlag = false;
	}

	BOOL koBool = Kappa_array.load("load" SLASH "lbkappa") && Omega_array.load("load" SLASH "lbomega");
	if (koBool)
	{
		kOmegaLoadedFlag = true;
	}
	else
	{
		kOmegaLoadedFlag = false;
		restartRunFlag = false;
	}
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



void clVariablesLB::updateOutputs()
{
	//Update_output_kernel[1].set_argument(4, (void*)&Save_loc);
	//Update_output_kernel[0].call_kernel();
	//Update_output_kernel[1].call_kernel();
	//vls.Masses.read_from_buffer(p.LBqueue);
	//Save_loc++;
	//clFlush(p.LBqueue);
}

//void functionPointerWrapper(void* pt2Object, FuncPtrType fptype_)
//{
//	clVariablesLB* mySelf = (clVariablesLB*)pt2Object;
//	if (fptype_ == ptr2CreateKernels)
//		mySelf->createKernels();
//	else if (fptype_ == ptr2SetKernelArgs)
//		mySelf->setKernelArgs();
//	else
//		mySelf->loadParams();
//}

///////////////////////////////////////////////////
////      Static Variable Initialization       ////
///////////////////////////////////////////////////
const int clVariablesLB::rev[9] = { 0, 2, 1, 4, 3, 6, 5, 8, 7 };
const int clVariablesLB::per[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

const int clVariablesLB::Cx[9] = { 0, 1, -1, 0, 0, 1, -1, 1, -1 };
const int clVariablesLB::Cy[9] = { 0, 0, 0, 1, -1, 1, -1, -1, 1 };

const cl_int2 clVariablesLB::Cxy[9] = { { { 0, 0 } }, { { 1, 0 } }, { { -1, 0 } }, { { 0, 1 } }, { { 0, -1 } },
{ { 1, 1 } }, { { -1, -1 } }, { { 1, -1 } }, { { -1, 1 } } };

const cl_double2 clVariablesLB::Cxy_double[9] = { { { 0., 0. } }, { { 1., 0. } }, { { -1., 0. } }, { { 0., 1. } }, { { 0., -1. } }, { { 1., 1. } }, { { -1., -1. } }, { { 1., -1. } }, { { -1., 1. } } };

const double clVariablesLB::Weigh[9] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };
const double clVariablesLB::Cs2 = 1. / 3.;

const int clVariablesLB::LF[LB_NUMBER_CONNECTIONS] = {	LB_NO_BOUNDARY_LINK, 
														LB_BOUNDARY_LINK_1, 
														LB_BOUNDARY_LINK_2,
														LB_BOUNDARY_LINK_3,
														LB_BOUNDARY_LINK_4,
														LB_BOUNDARY_LINK_5,
														LB_BOUNDARY_LINK_6,
														LB_BOUNDARY_LINK_7,
														LB_BOUNDARY_LINK_8};

const int clVariablesLB::boundsArr[LB_NUMBER_CONNECTIONS] = { C_BOUND, E_BOUND, W_BOUND, N_BOUND, S_BOUND, NE_BOUND, SW_BOUND, SE_BOUND, NW_BOUND };
const int clVariablesLB::boundsArrT1[LB_NUMBER_CONNECTIONS] = { C_BOUND, E_BOUND_T1, W_BOUND_T1, N_BOUND_T1, S_BOUND_T1, NE_BOUND_T1, SW_BOUND_T1, SE_BOUND_T1, NW_BOUND_T1 };
const int clVariablesLB::boundsArrT2[LB_NUMBER_CONNECTIONS] = { C_BOUND, E_BOUND_T2, W_BOUND_T2, N_BOUND_T2, S_BOUND_T2, NE_BOUND_T2, SW_BOUND_T2, SE_BOUND_T2, NW_BOUND_T2 };



