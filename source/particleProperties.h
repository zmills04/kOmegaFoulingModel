// particleProperties.h: Class storing particle property information
// and function used to initialize properties, write to file, etc.
// Contains both wall and particle material properties
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PARTICLEPROPERTIES_H__INCLUDED_)
#define AFX_PARTICLEPROPERTIES_H__INCLUDED_

#include "StdAfx.h"
#include "particleStructs.h"

#pragma once

// Will be able to treat class as though its an array of pParams
// (specifically, will be able to easily access elements of array and
// pass buffer to kernels)
class particleProperties : public PParam
{
public:
	particleProperties() : PParam("parP"),
		D_p("Dp"), D_p_real("Dp_real"), M_p("Mp"), 
		V_par("Vp"), D_dists("parDists"), F_po("Fpo"),
		Kth_pars("KthVal"), Q_A_prime("Qa_prime"),
		Q_A("Qa"), Tau_crit("tauCrit"), R_d("Rd")
	{}

	~particleProperties()
	{}

	int Nd;
	double WofA, Estar, WofA_s, Estar_s;
	double Tau_crit_max;
	double sootNumConc, parThermalCond, mfpAir;
	double surfEnergySurf, poissonSurf, yModSurf;
	double surfEnergySoot, poissonSoot, yModSoot;
	double hamakerConst, depPorosity, wallCorrection;
	double liftCoeff, parVolMultiplier, numEachPar;
	double inletConcPerDx; // concentration at inlet in #particles/dx
	// Arrays containing particle properties (D_p is particle
	// diameters in LB units, D_p_real is diameters in meters,
	// all others are in LB units)

	Array1Dd D_p, D_p_real, M_p, V_par, D_dists, F_po, Kth_pars;
	Array2Dd Q_A_prime, Q_A, Tau_crit, R_d;
	

	////////////////////////////
	//    Basic Functions     //
	////////////////////////////
	//void allocateArrays();
	//void createKernels();
	void ini();
	void loadParams();
	//void save2file();
	//void saveDebug();
	void saveParams();
	//void saveRestartFiles();
	//void setKernelArgs();
	void setSourceDefines();
	bool testRestartRun();

	////////////////////////////
	//  Class Specific Funcs  //
	////////////////////////////
	
	// creates Property arrays
	void createPropertyArrays();

	// fills Array1DPP array
	void fillPropStructArray();

	// Saves particle property infomation
	// for post-processing
	// writes files to particles directory
	void saveParticleArrays();

};


















#endif 