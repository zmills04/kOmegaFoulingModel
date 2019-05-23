// clVariables.cpp: implementation of the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////


//TODO:	save necessary variables, finish organizing functions and variables,
//		setup all load from bin files to allow a restart
//TODO: Add ability to use timer for delaying release of particles
//			when Dep_Flag == -2. This will allow for stagger release
//			of particles at beginning and after clumping particles
//			without having to track empty particles throughout domain
//			(which is currently the method being used)
//TODO: Add functions to trStructs for releasing host memory

#include "StdAfx.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"
#include "clProblem.h"
#include <random>




/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////////////////////                                /////////////////////////
////////////////////     Model Intitialization      /////////////////////////
////////////////////                                /////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


void clVariablesTR::allocateArrays()
{
	P.allocate(nN);
	NodV.allocate(p.nX - 1, p.nY - 1);//Stores Velocities and Temps used for calculating particle vels
	NodI.allocate(p.nX - 1, p.nY - 1);//Contains boundary links near lattice site as well as flag for wall 
	NodC.allocate(p.nX - 1, p.nY - 1);//Coefficients used to calculate NodV each time step
	BL.allocate(vls.nBL);//Boundary link info
	Active_flag.zeros(p.nX - 1, p.nY - 1);
	BL_dep.zeros(vls.nBL, parP.Nd);
	RandList.allocate(nN);//rng array
	BLind_ind.zeros(p.nX - 1, p.nY - 1);
	TR_indicies.allocate((p.nX - 1)*p.nY - 1);
	Active_indicies.allocate((p.nX - 1)*(p.nY - 1));
	Nod = new Array2DNT(p.nX - 1, p.nY - 1, "trNod"); //Temporary array storing lattice site info used to create NodI and NodV



	wallShear.allocateArrays();
	parSort.allocateArrays();
	glParticles.allocateArrays();
}

void clVariablesTR::ini()
{
	if (!trSolverFlag)
	{ // makes sure that there are no undefined variables when not using particle solver
		setSourceDefines();
		return;
	}
	rmax = (double)RAND_MAX + 1.;
	
	iniRand();

	if (restartRunFlag == false)
	{
		iniParticles();
		//call necessary restart functions
		parSort.sortTimer = 0;
	}
	
	
	iniNode();
	iniBLinks();
	splitNodeArrays();
		
	//parP class is initialized 
	wallShear.ini();
	parSort.ini();

	findNodeNeighs();

	


	setSourceDefines();
	std::function<void(void)> createKerPtr = std::bind(&clVariablesTR::createKernels, this);
	std::function<void(void)> setArgsPtr = std::bind(&clVariablesTR::setKernelArgs, this);

	sourceGenerator::SourceInstance()->addIniFunction(createKerPtr, setArgsPtr);


	if (saveMacroStart)
		save2file();

#ifdef SAVE_BIN_IN_LOAD
	saveRestartFilesIni();
	saveRestartFiles();
#endif
	
	Winds.save_txt_from_device("Winds");
	Node_neigh.save_txt_from_device("Nneigh");


	if (p.useOpenGL)
		vls.iniGL();

}



void clVariablesTR::loadParams()
{
	trSolverFlag = p.getParameter("Use Par Solver",
		USE_PARTICLE_SOLVER);

	if (!trSolverFlag && p.useOpenGL)
		p.useOpenGL = false;

	TrDomainSize.x = p.nX - 1;
	TrDomainSize.y = p.nY - 1;
	
	saveMacroStart = p.getParameter("Save Macros On Start", TR_SAVE_MACROS_ON_START);
	nN = p.getParameter("Num Particles", TRC_NUM_TRACERS);
	parP.Nd = p.getParameter<int>("Num Par Sizes");
	
	testRestartRun();
	
	indRadiusSearch = p.getParameter("Index Radius Search", INDEX_RADIUS_SEARCH);
	massFluxInlet = p.getParameter("Mass Flux Inlet", MASS_FLUX_INLET);
	maxBLPerNode = p.getParameter("Max BL Per Node", MAX_BL_PER_NODE);
	parReleaseTime = p.getParameter("Particle Release Time", PARTICLE_RELEASE_TIME);
	stopDistX = p.getParameter("X Stop Distance", STOP_DIST_X);
	timeBeforeReRelease = p.getParameter("Time Before Re-release", TIME_BEFORE_RERELEASE) / p.trSteps_wall;
	timeBeforeStick = p.getParameter("Time Before Stick", TIME_BEFORE_STICK) / p.trSteps_wall;
	maxNumBLRoll = p.getParameter("Max Num BL Roll", MAX_NUM_BL_ROLL);
	X_release = p.getParameter("X Release Pos", X_RELEASE_POS);
	xReleasePos = (int)X_release;
	xStopPos = p.getParameter("X Stop Pos", X_STOP_POS);
	reduceDepStop1 = p.getParameter("Reduce Dep Stop1", REDUCE_DEP_STOP1);
	reduceDepStop2 = p.getParameter("Reduce Dep Stop2", REDUCE_DEP_STOP2);
	startThermoVel = p.getParameter("Start Thermo Pos", START_THERMO_VEL);
	amtReduceDep = p.getParameter("Amount Reduce Dep", AMT_REDUCE_DEP);
	xMinVal = p.getParameter("X Min Val", X_MIN_VAL);
	cutoffRadius = p.getParameter("Cutoff Radius", CUTOFF_RADIUS);
	foulSizeSwitchSoot = p.getParameter("Foul Size Switch Soot", FOUL_SIZE_SWITCH_SOOT2);

	if (trSolverFlag == false)
		return;

	parP.loadParams();
	wallShear.loadParams();
	parSort.loadParams();
}

void clVariablesTR::testRestartRun()
{
	allocateArrays();

	restartRunFlag = p.getParameter("Restart Run", false);
	if (!restartRunFlag)
		return;
	restartRunFlag &= p.getParameter("Num Active Nodes", nActiveNodes, 0);
	restartRunFlag &= p.getParameter("Min BL Bottom", Bounds.MIN_BL_BOT, 0);
	restartRunFlag &= p.getParameter("Max BL Bottom", Bounds.MAX_BL_BOT, 0);
	restartRunFlag &= p.getParameter("Min BL Top", Bounds.MIN_BL_TOP, 0);
	restartRunFlag &= p.getParameter("Max BL Top", Bounds.MAX_BL_TOP, 0);
	restartRunFlag &= BL_dep.load("load" SLASH "BLdep_temp");
	restartRunFlag &= P.load("load" SLASH "trp");
	restartRunFlag &= parSort.testRestartRun();
}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////////////////////                                /////////////////////////
////////////////////     Kernel Intitialization     /////////////////////////
////////////////////                                /////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


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

	TR_Update_Par_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_no_wall");
	TR_Update_Par_kernel.set_size(nActiveNodes, WORKGROUPSIZE_TR);

	TR_Wall_Par_kernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_along_wall");
	TR_Wall_Par_kernel[0].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);
	//Global size is set after num_wall_nodes is calculated

	TR_Wall_Par_kernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_reflect_particles");
	TR_Wall_Par_kernel[1].set_size(2 * WORKGROUPSIZE_TR_WALL_REFLECT, WORKGROUPSIZE_TR_WALL_REFLECT);

	TR_Wall_Par_kernel[2].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_contact_wall");
	TR_Wall_Par_kernel[2].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	TR_Wall_Par_kernel[3].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_contact_wall2");
	TR_Wall_Par_kernel[3].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	TR_Wall_Par_kernel[4].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_deposit_particles_on_wall");
	TR_Wall_Par_kernel[4].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	Save_Dists_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_save_dists");
}



#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(),
void clVariablesTR::setSourceDefines()
{
	setSrcDefinePrefix "START_THERMO_VEL", startThermoVel);
	setSrcDefinePrefix "MAX_BL_PER_NODE", maxBLPerNode);
	setSrcDefinePrefix "TRC_NUM_TRACERS", nN);
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
	setSrcDefinePrefix "ALPHA_FLUID", vfd.Alpha_fluid);
	setSrcDefinePrefix "ALPHA_FOUL", vfd.Alpha_foul);
	if (MAX_NUM_BL_ROLL != -1)
	{
		setSrcDefinePrefix "MAX_NUM_BL_ROLL", maxNumBLRoll);
	}

	parP.setSourceDefines();
	parSort.setSourceDefines();
	wallShear.setSourceDefines();
}
#undef setSrcDefinePrefix


void clVariablesTR::updateIODists()
{
	Save_Dists_kernel.set_argument(6, &Save_IO_Loc);
	Save_Dists_kernel.call_kernel();
	clFinish(TRQUEUE);
	Save_IO_Loc++;
}

void clVariablesTR::resetIODists()
{
	Save_IO_Loc = 0;
	IO_dists_save.FillBuffer(0);
}


void clVariablesTR::saveRestartFilesIni()
{
	Active_indicies.save_bin_from_device("Active_indicies");
	Active_flag.save_bin_from_device("Active_flag");
	TR_indicies.save_bin_from_device("TR_indicies");
}

// TODO: add necessary calls to subclasses and set those functions
// up in the subclasses
void clVariablesTR::saveRestartFiles()
{
	P.save_bin_from_device("trP");
	
}

// TODO: decide what to save and set it up
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
	Node_neigh.save_txt_from_device("Node_neigh");
	Winds.save_txt_from_device("Winds");
	TR_indicies.save_txt_from_device("TR_indicies");
	Active_indicies.save_txt_from_device("Active_indicies");
	Active_flag.save_txt_from_device("Active_flag");
}


// TODO: split this into multiple functions to increase readability
void clVariablesTR::iniBLinks()
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
	
			Nod->operator[](nind).Wall_Flag = 1;
			Nod->operator[](nind).BLind[BLind_ind[nind]] = i;
			BLind_ind[nind]++;
		}
		else
		{
			for (int i0 = Cmin.x; i0 < Cmax.x; i0++)
			{
				for (int j0 = Cmin.y; j0 < Cmax.y; j0++)
				{
					testNode(i0, j0, c0, c1, i, 1);
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
			Nod->operator[](nind).Wall_Flag = 2;
			Nod->operator[](nind).BLind[BLind_ind[nind]] = i;
			BLind_ind[nind]++;
		}
		else
		{
			for (int i0 = Cmin.x; i0 < Cmax.x; i0++)
			{
				for (int j0 = Cmin.y; j0 < Cmax.y; j0++)
				{
					testNode(i0, j0, c0, c1, i, 2);
				}
			}
		}
	}
}

void clVariablesTR::callUpdateWallParticles(cl_event *fill_evt)
{
	int i = 0;
	vtr.TR_Wall_Par_kernel[0].set_argument(13, &i);
	vtr.TR_Wall_Par_kernel[0].call_kernel(nullptr, 1, fill_evt);
	for (i = 1; i < 9; i++)
	{
		vtr.TR_Wall_Par_kernel[0].set_argument(13, &i);
		vtr.TR_Wall_Par_kernel[0].call_kernel();
	}

	Reflect_inds.read_from_buffer(TRQUEUE_REF, CL_TRUE);

	if (Reflect_inds(0) > 0)
	{
		TR_Wall_Par_kernel[1].set_argument(4, &Reflect_inds(0));
		TR_Wall_Par_kernel[1].set_global_call_kernel(Reflect_inds(0));
		//TODO: Figure out what this is doing and make necessary changes
		Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 0);
	}

	if (Reflect_inds(1) == 0)
		return;

	while (TRUE)
	{
		TR_Wall_Par_kernel[2].set_argument(8, &Reflect_inds(1));
		Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 1);
		TR_Wall_Par_kernel[2].set_global_call_kernel(Reflect_inds(1));
		Reflect_inds.read_from_buffer(TRQUEUE_REF, CL_TRUE);
		if (Reflect_inds(1) == 0)
			break;

		TR_Wall_Par_kernel[3].set_argument(6, &Reflect_inds(1));
		Reflect_inds.FillBuffer(TRQUEUE, 0, 1, 1);
		TR_Wall_Par_kernel[3].set_global_call_kernel(Reflect_inds(1));;
		Reflect_inds.read_from_buffer(TRQUEUE_REF, CL_TRUE);
		if (Reflect_inds(1) == 0)
			break;
	}

	if (Reflect_inds(0) > 0)
	{
		vtr.TR_Wall_Par_kernel[4].set_argument(3, &Reflect_inds(0));
		vtr.TR_Wall_Par_kernel[4].set_global_call_kernel(Reflect_inds(0));

	}
	Reflect_inds.FillBuffer(0, TRQUEUE_REF);
}

void clVariablesTR::setKernelArgs()
{
	cl_int ind = 0;
	TR_Node_kernel[0].set_argument(ind++, vfd.Temp.get_add_Macro());
	TR_Node_kernel[0].set_argument(ind++, NodV.get_buf_add());
	TR_Node_kernel[0].set_argument(ind++, NodI.get_buf_add());
	TR_Node_kernel[0].set_argument(ind++, NodC.get_buf_add());
	TR_Node_kernel[0].set_argument(ind++, TR_indicies.get_buf_add());
	TR_Node_kernel[0].set_argument(ind++, &nActiveNodes);

	ind = 0;
	TR_Wall_Node_kernel[1].set_argument(ind++, vlb.Ux_array.get_buf_add());
	TR_Wall_Node_kernel[1].set_argument(ind++, vlb.Uy_array.get_buf_add());
	TR_Node_kernel[1].set_argument(ind++, NodV.get_buf_add());
	TR_Node_kernel[1].set_argument(ind++, NodI.get_buf_add());
	TR_Node_kernel[1].set_argument(ind++, NodC.get_buf_add());
	TR_Node_kernel[1].set_argument(ind++, TR_indicies.get_buf_add());
	TR_Node_kernel[1].set_argument(ind++, &nActiveNodes);

	ind = 0;
	TR_Wall_Node_kernel[0].set_argument(ind++, vfd.Temp.get_add_Macro());
	TR_Wall_Node_kernel[0].set_argument(ind++, NodC.get_buf_add());
	TR_Wall_Node_kernel[0].set_argument(ind++, NodV.get_buf_add());
	TR_Wall_Node_kernel[0].set_argument(ind++, TR_indicies.get_buf_add());
	TR_Wall_Node_kernel[0].set_argument(ind++, Winds.get_buf_add());
	TR_Wall_Node_kernel[0].set_argument(ind++, &Num_wall_nodes);

	ind = 0;
	TR_Wall_Node_kernel[1].set_argument(ind++, vlb.Ux_array.get_buf_add());
	TR_Wall_Node_kernel[1].set_argument(ind++, vlb.Uy_array.get_buf_add());
	TR_Wall_Node_kernel[1].set_argument(ind++, NodC.get_buf_add());
	TR_Wall_Node_kernel[1].set_argument(ind++, NodV.get_buf_add());
	TR_Wall_Node_kernel[1].set_argument(ind++, TR_indicies.get_buf_add());
	TR_Wall_Node_kernel[1].set_argument(ind++, Winds.get_buf_add());
	TR_Wall_Node_kernel[1].set_argument(ind++, &Num_wall_nodes);

	ind = 0; //Update Particles
	TR_Update_Par_kernel.set_argument(ind++, NodV.get_buf_add());
	TR_Update_Par_kernel.set_argument(ind++, parSort.Ploc.get_buf_add());
	TR_Update_Par_kernel.set_argument(ind++, P.get_buf_add());
	TR_Update_Par_kernel.set_argument(ind++, parP.get_buf_add());
	TR_Update_Par_kernel.set_argument(ind++, Node_neigh.get_buf_add());
	TR_Update_Par_kernel.set_argument(ind++, Update_flag.get_buf_add());
	TR_Update_Par_kernel.set_argument(ind++, NodI.get_buf_add());
	TR_Update_Par_kernel.set_argument(ind++, &nActiveNodes);


	ind = 0; //Update Particles
	TR_Wall_Par_kernel[0].set_argument(ind++, NodV.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, parSort.Ploc.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, P.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, PV.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, parP.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, Node_neigh.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, Update_flag.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, Winds.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, &Num_wall_nodes);
	TR_Wall_Par_kernel[0].set_argument(ind++, NodI.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, BL.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, Reflect_inds.get_buf_add());
	TR_Wall_Par_kernel[0].set_argument(ind++, Reflect_info.get_buf_add());
	//index 13 for neigh_ind


	ind = 0; //update Walls
	TR_Wall_Par_kernel[1].set_argument(ind++, BL.get_buf_add());
	TR_Wall_Par_kernel[1].set_argument(ind++, P.get_buf_add());
	TR_Wall_Par_kernel[1].set_argument(ind++, PV.get_buf_add());
	TR_Wall_Par_kernel[1].set_argument(ind++, Reflect_info.get_buf_add());
	//index 4 for number of pars

	ind = 0; //update Walls
	TR_Wall_Par_kernel[2].set_argument(ind++, BL.get_buf_add());
	TR_Wall_Par_kernel[2].set_argument(ind++, P.get_buf_add());
	TR_Wall_Par_kernel[2].set_argument(ind++, PV.get_buf_add());
	TR_Wall_Par_kernel[2].set_argument(ind++, parP.get_buf_add());
	TR_Wall_Par_kernel[2].set_argument(ind++, RandList.get_buf_add());
	TR_Wall_Par_kernel[2].set_argument(ind++, Reflect_info.get_buf_add());
	TR_Wall_Par_kernel[2].set_argument(ind++, Reflect_inds.get_buf_add());
	TR_Wall_Par_kernel[2].set_argument(ind++, &nN);
	//index 8 for num particles

	ind = 0; //update Walls
	TR_Wall_Par_kernel[3].set_argument(ind++, BL.get_buf_add());
	TR_Wall_Par_kernel[3].set_argument(ind++, P.get_buf_add());
	TR_Wall_Par_kernel[3].set_argument(ind++, PV.get_buf_add());
	TR_Wall_Par_kernel[3].set_argument(ind++, Reflect_info.get_buf_add());
	TR_Wall_Par_kernel[3].set_argument(ind++, Reflect_inds.get_buf_add());
	TR_Wall_Par_kernel[3].set_argument(ind++, &nN);
	//index 6 for num particles

	ind = 0; //update Walls
	TR_Wall_Par_kernel[4].set_argument(ind++, BL.get_buf_add());
	TR_Wall_Par_kernel[4].set_argument(ind++, P.get_buf_add());
	TR_Wall_Par_kernel[4].set_argument(ind++, Reflect_info.get_buf_add());
	//index 3 for num particles

	if (calcIOFlag)
	{
		ind = 0;
		Save_Dists_kernel.set_argument(ind++, parSort.Ploc.get_buf_add());
		Save_Dists_kernel.set_argument(ind++, P.get_buf_add());
		Save_Dists_kernel.set_argument(ind++, IO_inds.get_buf_add());
		Save_Dists_kernel.set_argument(ind++, Node_neigh.get_buf_add());
		Save_Dists_kernel.set_argument(ind++, IO_dists_save.get_buf_add());
		Save_Dists_kernel.set_argument(ind++, &IO_inds_info);
		Save_Dists_kernel.set_argument(ind++, &Save_IO_Loc);
		Save_Dists_kernel.set_argument(ind++, &parP.Nd);
	}
}

void clVariablesTR::checkNeighNode(int i, int j, int *cur_ind, int nod_loc)
{
	if(Active_flag(i, j)  == 1)
	{
		Node_neigh(nod_loc, *cur_ind) = i*(p.nY-1)+j;
		(*cur_ind) += 1;
	}
}

void clVariablesTR::findNodeNeighs()
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
		checkNeighNode(i + 1, j, &cur_ind, ii);
		checkNeighNode(i - 1, j, &cur_ind, ii);
		checkNeighNode(i, j + 1, &cur_ind, ii);
		checkNeighNode(i, j - 1, &cur_ind, ii);
		checkNeighNode(i + 1, j + 1, &cur_ind, ii);
		checkNeighNode(i - 1, j - 1, &cur_ind, ii);
		checkNeighNode(i + 1, j - 1, &cur_ind, ii);
		checkNeighNode(i - 1, j + 1, &cur_ind, ii);
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
	IO_dists_save.zeros(OUTPUT_MAX_LINES_IO * 2 * parP.Nd);
	IO_dists_save.allocate_buffer_w_copy();
	IO_inds.allocate_buffer_w_copy();
	IO_inds.FreeHost();
	Save_IO_Loc = 0;
}



double clVariablesTR::weightKernelFunc(double Sval)
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


void clVariablesTR::allocateBuffers()
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
	BL.allocate_buffer_w_copy();
	Winds.allocate_buffer_w_copy();
	BL_dep.allocate_buffer_w_copy();
	Node_neigh.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	BLind_ind.allocate_buffer_w_copy();

////////Allocation of Buffers only used for updates on Device////////////
	Update_flag.zeros(nN);
	Update_flag.allocate_buffer_w_copy();
	Update_flag.FreeHost();

	
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





void clVariablesTR::testNode(int i, int j, cl_double2 vL0, cl_double2 vL1, int blind, int bot_flag)
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
	if (testCross(vL0, vLd, CC00, CC10) || testCross(vL0, vLd, CC00, CC01) || 
		testCross(vL0, vLd, CC10, CC11) || testCross(vL0, vLd, CC01, CC11))
	{
		Nod->operator[](nind).Wall_Flag = bot_flag;
		Nod->operator[](nind).BLind[BLind_ind[nind]] = blind;
		BLind_ind[nind]++;
		return;
	}
	
}

bool clVariablesTR::testCross(cl_double2 vL0, cl_double2 vLd, cl_double2 vP0, cl_double2 vP1)
{
	cl_double2 vP10 = Subtract2(vP1,vP0);
	cl_double2 vNm = { { -1.*vP10.y, vP10.x } };
	double den = DOT_PROD(vNm, vLd);
	if (den == 0.)
		return false;
	cl_double2 vPL0 = Subtract2(vP0, vL0);
	double dist = DOT_PROD(vNm, vPL0) / den;
	cl_double2 vC = { { vL0.x + vLd.x * dist, vL0.y + vLd.y * dist } };
	cl_double2 vd0 = Subtract2(vC, vP0), vd1 = Subtract2(vC,vP1);
	if (dist >= 0. && dist <= 1.)
	if (fabs(GETLEN(vP10) - GETLEN(vd0) - GETLEN(vd1)) < p.eps)
			return true;
	return false;
}

void clVariablesTR::iniRand()
{
	if (RandList.load("load" SLASH "RandList") == true)
		return;
	for (int i = 0; i < nN; i++)
	{
		RandList[i].x = (cl_uint)rand();
		RandList[i].y = (cl_uint)rand();
	}
}




void clVariablesTR::saveBox(double x1, double y1, double dx, double dy)
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


int clVariablesTR::getType()
{
	double randval = rand1();
	for (int j = 0; j < parP.Nd-1; j++)
	{
		if (randval >= parP[j].D_dist && randval < parP[j + 1].D_dist)
			return j;
	}
	return parP.Nd - 1; 
}

double clVariablesTR::getOffset(double xval)
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

void clVariablesTR::iniParticles()
{
	if (restartRunFlag)
		return;
	for (int i = 0; i < nN; i++)
	{
		par Ptemp;
		// Particles will be released a portion at a time using reReleasePar kernel 
		Ptemp.pos.x = 0.;
		Ptemp.pos.y = 0.;
		
			
		Ptemp.type = getType();
		Ptemp.timer = 1;
		Ptemp.Dep_Flag = -2;
		Ptemp.Num_rep = 0;
		Ptemp.Dep_timer = 1;
		cl_int2 Posi = { { (int)floor(Ptemp.pos.x), (int)floor(Ptemp.pos.y) } };

		Ptemp.loc = Posi.x*vtr.TrDomainSize.y + Posi.y;
		P[i] = Ptemp;
	}

}


// TODO: figure out changes necessary for 
// the old reduced domain implementation and
// and switching to i+j*nX ordering
int clVariablesTR::getInd(cl_int2 ij)
{
	return ij.x + ij.y*p.nX;
}

void clVariablesTR::iniNode()
{
	int ind = 0;
	int trind = 0;
	int MinXval = floor(X_release - 2.);
	int MaxXval = ceil(xStopPos);

	for (int i = 0; i < p.nX-1; i++)
	{
		for (int j = 0; j < p.nY-1; j++)
		{
			cl_int2 ii00 = { { i, j } };
			cl_int2 ii10 = { { i + 1, j } };
			cl_int2 ii01 = { { i, j + 1 } };
			cl_int2 ii11 = { { i + 1, j + 1 } };

			Nact Neightemp;
			Neightemp.ii00 = ii00;
			Neightemp.ii10 = ii10;
			Neightemp.ii01 = ii01;
			Neightemp.ii11 = ii11;

			int jj00 = getInd(ii00);
			int jj10 = getInd(ii10);
			int jj01 = getInd(ii01);
			int jj11 = getInd(ii11);

			Nod->operator()(i, j).neigh.x = jj00;
			Nod->operator()(i, j).neigh.y = jj10;
			Nod->operator()(i, j).neigh.z = jj01;
			Nod->operator()(i, j).neigh.w = jj11;

			int tnum = 0;

			if (vls.M_o(ii00) != LB_SOLID)
				tnum += 1;
			if (vls.M_o(ii10) != LB_SOLID)
				tnum += 2;
			if (vls.M_o(ii11) != LB_SOLID)
				tnum += 4;
			if (vls.M_o(ii01) != LB_SOLID)
				tnum += 8;

			if (tnum > 0 && i > MinXval && i <= MaxXval)
			{
				Active_flag(i,j) = 1;
				cl_int2 trind_temp = { { i, j } };
				Active_indicies(trind) = trind_temp;
				TR_indicies(trind++) = i*(p.nY - 1) + j;
			}

			if (tnum == 0)
				tnum = -1;
			Nod->operator()(i, j).type = tnum;

			if (tnum == -1)
			{
				Nod->operator()(i, j).dX_cur.x = 1.;
				Nod->operator()(i, j).dX_cur.y = 1.;
			}
			else if (tnum == 0)
			{
				Nod->operator()(i, j).dX_cur.x = 1.;
				Nod->operator()(i, j).dX_cur.y = 1.;
			}
			else if (tnum == 1)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii00.x, ii00.y, 0);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii00.x, ii00.y, 2);
			}
			else if (tnum == 2)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii10.x, ii10.y, 1);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii10.x, ii10.y, 2);
			}
			else if (tnum == 3)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii00.x, ii00.y, 2);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii10.x, ii10.y, 2);
			}
			else if (tnum == 4)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii11.x, ii11.y, 1);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii11.x, ii11.y, 3);
			}
			else if (tnum == 5)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii00.x, ii00.y, 0);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii00.x, ii00.y, 2);
			}
			else if (tnum == 6)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii10.x, ii10.y, 1);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii11.x, ii11.y, 1);
			}
			else if (tnum == 7)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii00.x, ii00.y, 2);
				Nod->operator()(i, j).dX_cur.y = 1.;
			}
			else if (tnum == 8)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii01.x, ii01.y, 0);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii01.x, ii01.y, 3);
			}
			else if (tnum == 9)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii00.x, ii00.y, 0);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii01.x, ii01.y, 0);
			}
			else if (tnum == 10)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii10.x, ii10.y, 1);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii10.x, ii10.y, 2);
			}
			else if (tnum == 11)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii10.x, ii10.y, 2);
				Nod->operator()(i, j).dX_cur.y = 1.;
			}
			else if (tnum == 12)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii01.x, ii01.y, 3);
				Nod->operator()(i, j).dX_cur.y = vls.dXArrCur(ii11.x, ii11.y, 3);
			}
			else if (tnum == 13)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii11.x, ii11.y, 3);
				Nod->operator()(i, j).dX_cur.y = 1.;
			}
			else if (tnum == 14)
			{
				Nod->operator()(i, j).dX_cur.x = vls.dXArrCur(ii01.x, ii01.y, 3);
				Nod->operator()(i, j).dX_cur.y = 1.;
			}
			else if (tnum == 15)
			{
				Nod->operator()(i, j).dX_cur.x = 1.;
				Nod->operator()(i, j).dX_cur.y = 1.;
			}
			Nod->operator()(i, j).dX = Nod->operator()(i, j).dX_cur;
			ind++;
		}
	}

	for (int i = 0; i < p.nX-1; i++)
	{
		for (int j = 0; j < p.nY-1; j++)
		{
			for (int k = 0; k < MAX_BL_PER_NODE; k++)
			{
				Nod->operator()(i, j).BLind[k] = -1;
			}
			Nod->operator()(i, j).Wall_Flag = 0;
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
	p.setParameter("Index Radius Search", indRadiusSearch);
	p.setParameter("Mass Flux Inlet", massFluxInlet);
	p.setParameter("Max BL Per Node", maxBLPerNode);
	p.setParameter("Particle Release Time", parReleaseTime);
	p.setParameter("X Stop Distance", stopDistX);
	p.setParameter("Time Before Re-release", timeBeforeReRelease*p.trSteps_wall);
	p.setParameter("Time Before Stick", timeBeforeStick*p.trSteps_wall);
	p.setParameter("Max Num BL Roll", maxNumBLRoll);
	p.setParameter("X Release Pos", X_release);
	p.setParameter("X Stop Pos", xStopPos);
	p.setParameter("Reduce Dep Stop1", reduceDepStop1);
	p.setParameter("Reduce Dep Stop2", reduceDepStop2);
	p.setParameter("Start Thermo Pos", startThermoVel);
	p.setParameter("Amount Reduce Dep", amtReduceDep);
	p.setParameter("X Min Val", xMinVal);
	p.setParameter("Cutoff Radius", cutoffRadius);
	p.setParameter("Foul Size Switch Soot", foulSizeSwitchSoot);

	parSort.saveParams();
	wallShear.saveParams();
	parP.saveParams();
	glParticles.saveParams();
}

void clVariablesTR::splitNodeArrays()
{
	double Ka = (vfd.kAir * p.DELTA_F / p.DELTA_T);
	double Ks = (vfd.kSoot * p.DELTA_F / p.DELTA_T);

	for (int i = 0; i < (p.nX-1)*(p.nY-1); i++)
	{
		NodC[i].neigh.x = Nod->operator[](i).neigh.x;
		NodC[i].neigh.y = Nod->operator[](i).neigh.y;
		NodC[i].neigh.z = Nod->operator[](i).neigh.z;
		NodC[i].neigh.w = Nod->operator[](i).neigh.w;
		double dXcx = Nod->operator[](i).dX_cur.x, dXcy = Nod->operator[](i).dX_cur.y;
		double dXx = Nod->operator[](i).dX.x, dXy = Nod->operator[](i).dX.y;
		int typet = Nod->operator[](i).type;
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
			NodI[i].BLind[j] = Nod->operator[](i).BLind[j];
		}
		NodI[i].Wall_Flag = Nod->operator[](i).Wall_Flag;

	}
	delete Nod;
}


double clVariablesTR::rand1()
{
	double randval = (double)rand();
	return randval / rmax;
}

void clVariablesTR::releaseObjects()
{// All Kernel and Array objects are handled by their destructor
	delete[] Nod; // this should be deleted after initialization, so this is just in case
}
