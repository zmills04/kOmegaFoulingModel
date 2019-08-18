// particleProperties.cpp: Implementation of methods for 
// particleProperties class declared in particleProperties.h.
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////

#include "particleProperties.h"
#include "clProblem.h"
#include "clVariablesTR.h"

// TODO: make sure that this is working after change to 
//		defining arrays as constants in opencl source

//void particleProperties::allocateArrays()
//{
//	this->allocateArrays();
//}


//void particleProperties::allocateBuffers()
//{
//}

void particleProperties::createKernels()
{
}

void particleProperties::createPropertyArrays()
{
	double ConstA = 1.2, ConstB = 0.4, ConstC = 1.1;
	double ConstCs = 1.14, ConstCm = 1.17, ConstCt = 2.18;
	

	double Distsum = 0.;
	Kth.zeros(Nd);
	F_po.zeros(Nd);
	R_d.zeros(Nd, 2);
	V_par.zeros(Nd);
	Tau_crit_max = 0.;

	for (int i = 0; i < Nd; i++)
	{
		double Knud = mfpAir / D_p_real(i);

		double CCF = 1. + Knud * (ConstA + ConstB * exp(-ConstC / Knud));
		double k_ratio = vfd.kAir / parThermalCond;
		double Kth1 = 2. * ConstCs * CCF / (1. + 3.*ConstCm * Knud);
		double Kth2 = (k_ratio + ConstCt * Knud) / (1. + 2. * k_ratio + 2.*ConstCt*Knud);

		Kth(i) = Kth1*Kth2;

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
		Q_A(i).x = QAi;

		double QAi_s = 2.*WofA_s*PI_NUMBER*R_d(i, 1)*R_d(i, 1);
		QAi *= (p.DELTA_F / p.DELTA_L);
		Q_A(i).y = QAi_s;

		double Q_A_pt = pow(WofA, 5.)*pow(D_p_real(i), 4.) / Estar / Estar / 16;
		double QA_pi = 7.09 * pow(Q_A_pt, (1. / 3.));
		QA_pi *= p.DELTA_P;
		Q_A_prime(i).x = QA_pi;

		double Q_A_pt_s = pow(WofA_s, 5.)*pow(D_p_real(i), 4.) / Estar_s / Estar_s / 16;
		double QA_pi_s = 7.09 * pow(Q_A_pt_s, (1. / 3.));
		QA_pi_s *= p.DELTA_P;
		Q_A_prime(i).y = QA_pi_s;

		double Mpi = PI_NUMBER*D_p_real(i)*D_p_real(i)*vfd.rhoSoot / 4.;
		Mpi *= p.DELTA_M;
		Mp(i) = Mpi;

		double F_poi = 625.*hamakerConst / 3. / D_p_real(i);
		F_poi *= p.DELTA_F;
		F_po(i) = F_poi;

		double Dpi = D_p_real(i)*p.DELTA_L;
		Dp(i) = Dpi;

		double Vpi = Dp(i) * Dp(i) * PI_NUMBER / 4. / (1. - depPorosity);
		V_par(i) = Vpi;

		double Tci = 4.*F_po(i)*R_d(i, 0) / (3. * PI_NUMBER * pow(Dp(i), 3.) * wallCorrection);

		tau_crit(i).x = Tci;
		if (tau_crit(i).x > Tau_crit_max)
			Tau_crit_max = tau_crit(i).x;

		double Tci_s = 4.*F_po(i)*R_d(i, 1) / (3. * PI_NUMBER * pow(Dp(i), 3.) * wallCorrection);


		tau_crit(i).y = Tci_s;
		if (tau_crit(i).y > Tau_crit_max)
			Tau_crit_max = tau_crit(i).y;

		D_dist(i) = Distsum;
		Distsum += D_dist(i);
		L_coeff(i) = liftCoeff * Dp(i) * Dp(i) * Dp(i) / vlb.MuVal / 8.;
		D_coeff(i) = 6.0 * PI_NUMBER * Dp(i) * Dp(i) * wallCorrection / 4.;
	}
}


// TODO: Figure out if any of these can be
//		freed now that there is no duplicated
//		temporary arrays
void particleProperties::freeHostArrays()
{
	// D_p_real and D_dists currently kept,
	// because they are needed for each
	// saveParameters step.
	//Kth_pars.FreeHost();
	//D_p.FreeHost();
	//Q_A_prime.FreeHost();
	//Q_A.FreeHost();
	//M_p.FreeHost();
	//F_po.FreeHost();
	//R_d.FreeHost();
	//Tau_crit.FreeHost();
	//V_par.FreeHost();
}

void particleProperties::ini()
{
	// fills Pparam base class with data contained in
	// this classes arrays
	//fillPropStructArray();

	// Writes arrays to file in a folder called particles
	saveParticleArrays();
	
	// Copy arrays to device
	//copyToDevice();

#if _DEBUG
	LOGMESSAGE("Initialized Particle Properties Class");
#endif

}

void particleProperties::loadParams()
{
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
	numEachPar = p.getParameter("Num Each Par", NUM_EACH_PAR);
	parVolMultiplier = numEachPar / (1. - depPorosity);

	this->setSizes(Nd, Nd, 1);

	D_p_real.zeros(Nd);
	p.yamlIn["Dp Dists"] >> D_dist;
	p.yamlIn["Par Diameters"] >> D_p_real;

	double Num_per_m2 = sootNumConc / p.DELTA_L / p.DELTA_L;
	inletConcPerDx = Num_per_m2 * 2. * p.Pipe_radius;

	//number of particles released after each sort step 
	//must be multiplied by Umean to get true value

	WofA = 2.*sqrt(surfEnergySurf * surfEnergySoot);
	Estar = 1. / ((1. - poissonSurf*poissonSurf) / yModSurf + (1. - poissonSoot*poissonSoot) / yModSoot);
	Estar_s = yModSoot;
	WofA_s = surfEnergySoot;

	allocateArrays();
	createPropertyArrays();
}




void particleProperties::save2file()
{
}


void particleProperties::saveDebug()
{
}

void particleProperties::saveParams()
{
	p.setParameter("Num Par Sizes", Nd);
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
	p.setParameter("Num Each Par", numEachPar);

	if (vtr.trSolverFlag)
	{
		*p.yamlOut << YAML::BeginMap;
		*p.yamlOut << YAML::Key << "Dp Dists";
		*p.yamlOut << YAML::Value << D_dist;
		*p.yamlOut << YAML::EndMap;

		*p.yamlOut << YAML::BeginMap;
		*p.yamlOut << YAML::Key << "Par Diameters";
		*p.yamlOut << YAML::Value << D_p_real;
		*p.yamlOut << YAML::EndMap;
	}
}

void particleProperties::saveParticleArrays()
{
	std::string NewDir = "particles";
	MakeDir(NewDir);
	Dp.savetxt("particles" SLASH "D_p");
	D_p_real.savetxt("particles" SLASH "D_p_real");
	Q_A_prime.savetxt("particles" SLASH "Q_A_prime");
	Q_A.savetxt("particles" SLASH "Q_A");
	Mp.savetxt("particles" SLASH "M_p");
	F_po.savetxt("particles" SLASH "F_po");
	tau_crit.savetxt("particles" SLASH "Tau_crit");
	V_par.savetxt("particles" SLASH "Vol_p");
	Kth.savetxt("particles" SLASH "Kth");
	D_dist.savetxt("particles" SLASH "D_dist");
}

void particleProperties::saveRestartFiles()
{
}

void particleProperties::saveTimeData()
{
}

void particleProperties::setKernelArgs()
{
}

#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(),
void particleProperties::setSourceDefines()
{
	setSrcDefinePrefix "NUM_PAR_SIZES", Nd);

// The defines are added even when trSolver is not being used just in case
// a kernel references one of these defines (due to methods occasionally
// overlapping in functions). This variable is not included in files outside
// of TR_kernels, so it will be fine to not include it in definition file.
	if (vtr.trSolverFlag)
	{
		char parMultVal[80];
		sprintf(parMultVal, "%20.16g", V_par[0] / vls.lsSpacing);
		std::string parmult = "constant double Par_multiplier[" + std::to_string(Nd) + "] = { " + parMultVal;

		sprintf(parMultVal, "%20.16g", Dp[0]);
		std::string dpString = "constant double paramDp[" + std::to_string(Nd) + "] = { " + parMultVal;

		sprintf(parMultVal, "%20.16g", Mp[0]);
		std::string mpString = "constant double paramMp[" + std::to_string(Nd) + "] = { " + parMultVal;

		sprintf(parMultVal, "%20.16g", Kth[0]);
		std::string kthString = "constant double paramKth[" + std::to_string(Nd) + "] = { " + parMultVal;

		sprintf(parMultVal, "%20.16g", L_coeff[0]);
		std::string lcoeffString = "constant double paramLcoeff[" + std::to_string(Nd) + "] = { " + parMultVal;

		sprintf(parMultVal, "%20.16g", D_coeff[0]);
		std::string dcoeffString = "constant double paramDcoeff[" + std::to_string(Nd) + "] = { " + parMultVal;

		sprintf(parMultVal, "%20.16g", D_dist[0]);
		std::string ddistString = "constant double paramDdist[" + std::to_string(Nd) + "] = { " + parMultVal;

		sprintf(parMultVal, "%20.16g", Q_A_prime(0).x);
		std::string qaprimeString = "constant double paramQaPrime[" + std::to_string(Nd * 2) + "] = { " + parMultVal;
		sprintf(parMultVal, ", %20.16g", Q_A_prime(0).y);
		qaprimeString.append(parMultVal);

		sprintf(parMultVal, "%20.16g", Q_A(0).x);
		std::string qaString = "constant double paramQa[" + std::to_string(Nd * 2) + "] = { " + parMultVal;
		sprintf(parMultVal, ", %20.16g", Q_A(0).y);
		qaString.append(parMultVal);

		sprintf(parMultVal, "%20.16g", tau_crit(0).x);
		std::string taucritString = "constant double paramTauCrit[" + std::to_string(Nd * 2) + "] = { " + parMultVal;
		sprintf(parMultVal, ", %20.16g", tau_crit(0).y);
		taucritString.append(parMultVal);




		for (int i = 1; i < Nd; i++)
		{
			sprintf(parMultVal, ", %20.16g", V_par[i] / vls.lsSpacing);
			parmult.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", Dp[i]);
			dpString.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", Mp[i]);
			mpString.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", Kth[i]);
			kthString.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", L_coeff[i]);
			lcoeffString.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", D_coeff[i]);
			dcoeffString.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", D_dist[i]);
			ddistString.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", Q_A(i).x);
			qaString.append(parMultVal);
			sprintf(parMultVal, ", %20.16g", Q_A(i).y);
			qaString.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", Q_A_prime(i).x);
			qaprimeString.append(parMultVal);
			sprintf(parMultVal, ", %20.16g", Q_A_prime(i).y);
			qaprimeString.append(parMultVal);

			sprintf(parMultVal, ", %20.16g", tau_crit(i).x);
			taucritString.append(parMultVal);
			sprintf(parMultVal, ", %20.16g", tau_crit(i).y);
			taucritString.append(parMultVal);
		}

		parmult.append("};");
		dpString.append("};");
		mpString.append("};");
		kthString.append("};");
		lcoeffString.append("};");
		dcoeffString.append("};");
		ddistString.append("};");
		qaString.append("};");
		qaprimeString.append("};");
		taucritString.append("};");

		SOURCEINSTANCE->addString2Defines(parmult);
		SOURCEINSTANCE->addString2Defines(dpString);
		SOURCEINSTANCE->addString2Defines(mpString);
		SOURCEINSTANCE->addString2Defines(kthString);
		SOURCEINSTANCE->addString2Defines(lcoeffString);
		SOURCEINSTANCE->addString2Defines(dcoeffString);
		SOURCEINSTANCE->addString2Defines(ddistString);
		SOURCEINSTANCE->addString2Defines(qaString);
		SOURCEINSTANCE->addString2Defines(qaprimeString);
		SOURCEINSTANCE->addString2Defines(taucritString);
	}
}
#undef setSrcDefinePrefix

bool particleProperties::testRestartRun()
{
	return true;
}

void particleProperties::updateTimeData()
{
}

