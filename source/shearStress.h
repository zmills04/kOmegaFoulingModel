// particleProperties.h: Class storing particle property information
// and function used to initialize properties, write to file, etc.
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SHEARSTRESS_H__INCLUDED_)
#define AFX_SHEARSTRESS_H__INCLUDED_

#include "StdAfx.h"


#pragma once

class shearStress
{
public:
	
	shearStress() : Sind("ssSind"), BLindicies("ssBLindicies"),
		Weights("ssWeights"), Shear_inds("ssShear_inds"),
		Shear_coeffs("ssShear_coeffs"), Bindicies_loc("ssBindicies_loc"),
		Tau("ssTau"), SS_output("SS_output")
	{}

	~shearStress()
	{}

	// Shear calculation kernels
	// [0] and [1] calculate at boundary
	// nodes ([0] using vlb.fB and [1] using vlb.fA)
	// [2] interpolates to wall boundary locations
	Kernel shearKernels[3];

	// Shear removal of deposited particles
	Kernel trShearRemovalKernel;
	
	// Updates Shear Stress coefficients (called in vfl)
	// two step process
	Kernel updateSSKernel[2];

	// Prepares shear values to be written to file
	Kernel saveShearKernel;


	int bNodeTopStart;
	int Save_loc_SS;
	int Shear_array_len;
	bool calcSSFlag;
	int shearSize;
	int numbl_bounds;
	int maxOutLinesSS;

	//Index of first node used in averaging to get shear at wall
	Array1Dv4i Sind;

	//Index of BL array each index in SS arrays corresponds to	
	Array1Di BLindicies;

	// Weights for averaging shear when interpolating from SS at nodes to SS at wall
	Array1Dv4d Weights;

	// Used to update SS coefficients
	Array1Dv2i Bindicies_loc;

	//locations of nodes where shear is calculated
	Array1Dv2i Shear_inds;

	//coefficients used in calc. of shear at nodes
	Array1Dv8d Shear_coeffs;

	//Shear stresses calculated at boundary nodes
	Array1Dd Tau;

	//Array to store Shear Stresses for output
	Array2Dd SS_output;


	////////////////////////////
	//    Basic Functions     //
	////////////////////////////
	void allocateArrays();
	void allocateBuffers();
	void createKernels();
	void freeHostArrays();
	void ini();
	void loadParams();
	void save2file();
	void saveDebug();
	void saveParams();
	void saveRestartFiles();
	void setKernelArgs();
	void setSourceDefines();
	bool testRestartRun();
	


	////////////////////////////
	//  Class Specific Funcs  //
	////////////////////////////

	void iniShearCoeffs();
	void updateShearArrays();

	// Functions for writing Shear Stresses to file
	void updateSS();
	void resetSSOut();
	void saveSS();
	void updateParRemArgs();
};



#endif