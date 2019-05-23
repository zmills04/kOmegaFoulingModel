#include "particleProperties.h"
#include "clProblem.h"
#include "clVariablesTR.h"



void particleProperties::createPropertyArrays()
{
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
}


void particleProperties::fillPropStructArray()
{
	double Distsum = 0.;
	this->setSizes(Nd, Nd, 1);
	for (int i = 0; i < Nd; i++)
	{
		legacyPParam struct_;
		struct_.Dp = D_p(i);
		struct_.Q_A.x = Q_A(i, 0);
		struct_.Q_A.y = Q_A(i, 1);
		struct_.Q_A_prime.x = Q_A_prime(i, 0);
		struct_.Q_A_prime.y = Q_A_prime(i, 1);
		struct_.tau_crit.x = Tau_crit(i, 0);
		struct_.tau_crit.y = Tau_crit(i, 1);
		struct_.Mp = M_p(i);
		struct_.Kth = Kth_pars(i);
		struct_.D_dist = Distsum;
		Distsum += D_dists(i);
		struct_.L_coeff = liftCoeff * struct_.Dp * struct_.Dp * struct_.Dp / vlb.MuVal / 8.;
		struct_.D_coeff = 6 * PI_NUMBER * struct_.Dp * struct_.Dp * wallCorrection / 4.;
		this->setStruct(struct_, i);
	}
}

void particleProperties::ini()
{
	fillPropStructArray();
	saveParticleArrays();
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

	D_p_real.zeros(Nd);
	D_dists.zeros(Nd);
	p.yamlIn["Dp Dists"] >> D_dists;
	p.yamlIn["Par Diameters"] >> D_p_real;

	double Num_per_m2 = sootNumConc / p.DELTA_L / p.DELTA_L;
	inletConcPerDx = Num_per_m2 * 2. * p.Pipe_radius;


	//number of particles released after each sort step 
	//must be multiplied by Umean to get true value

	WofA = 2.*sqrt(surfEnergySurf * surfEnergySoot);
	Estar = 1. / ((1. - poissonSurf*poissonSurf) / yModSurf + (1. - poissonSoot*poissonSoot) / yModSoot);
	Estar_s = yModSoot;
	WofA_s = surfEnergySoot;

	createPropertyArrays();
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
		*p.yamlOut << YAML::Value << D_dists;
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
	p.MakeDir(NewDir);
	D_p.savetxt("particles" SLASH "D_p");
	D_p_real.savetxt("particles" SLASH "D_p_real");
	Q_A_prime.savetxt("particles" SLASH "Q_A_prime");
	Q_A.savetxt("particles" SLASH "Q_A");
	M_p.savetxt("particles" SLASH "M_p");
	F_po.savetxt("particles" SLASH "F_po");
	Tau_crit.savetxt("particles" SLASH "Tau_crit");
	V_par.savetxt("particles" SLASH "Vol_p");
	Kth_pars.savetxt("particles" SLASH "Kth");
	D_dists.savetxt("particles" SLASH "D_dist");
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
		std::string parmult = "const double Par_multiplier[" + std::to_string(Nd) + "] = { " + parMultVal;

		for (int i = 1; i < Nd; i++)
		{
			sprintf(parMultVal, ", %20.16g, ", V_par[i] / vls.lsSpacing);
			parmult.append(parMultVal);
		}

		parmult.append("};");
		SOURCEINSTANCE->addString2Defines(parmult);
	}
}
#undef setSrcDefinePrefix

bool particleProperties::testRestartRun()
{
	return true;
}


