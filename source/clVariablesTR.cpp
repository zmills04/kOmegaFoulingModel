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
//TODO: Need to handle conversion of vNvec to vTvec differently for top and bottom
//		BLinks: on bottom vTvec = {{-vNvec.y, vNvec.x}}, on top vTvec = {{vNvec.y, -vNvec.x}}
//		i.e. vTvec must be multiplied by -1 for top BLinks
//TODO: Implement storing distribution info of particles leaving outlet
//		Particles leaving outlet will be assigned depFlag = -3, so that
//		they can be counted before being rereleased. Released particles will
//		also be stored in this array. Each row will store the number of particles
//		entering and leaving domain for each diameter for the time period defined by the
//		save step. THere will not be a kernel to calculate anything, but the update function
//		will increment the row index.


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


// TODO: reduce TrDomainSize.x to the size of xMinVal to xStopVal
//		this will significantly reduce memory usage

void clVariablesTR::loadParams()
{
	trSolverFlag = p.getParameter("Use Par Solver",
		USE_PARTICLE_SOLVER);

	if (!trSolverFlag && p.useOpenGL)
		p.useOpenGL = false;

	saveMacroStart = p.getParameter("Save Macros On Start", TR_SAVE_MACROS_ON_START);
	nN = p.getParameter("Num Particles", TRC_NUM_TRACERS);
	parP.Nd = p.getParameter<int>("Num Par Sizes");
	calcIOFlag = p.getParameter("Calc IO Dists", SAVE_IO_DISTS);
	maxOutputLinesIO = p.getParameter("Max Out Lines IO", OUTPUT_MAX_LINES_SS);

	indRadiusSearch = p.getParameter("Index Radius Search", INDEX_RADIUS_SEARCH);
	massFluxInlet = p.getParameter("Mass Flux Inlet", MASS_FLUX_INLET);
	maxBLPerNode = p.getParameter("Max BL Per Node", MAX_BL_PER_NODE);
	parReleaseTime = p.getParameter("Particle Release Time", PARTICLE_RELEASE_TIME);
	
	timeBeforeReRelease = p.getParameter("Time Before Re-release", TIME_BEFORE_RERELEASE) / p.trSteps_wall;
	timeBeforeStick = p.getParameter("Time Before Stick", TIME_BEFORE_STICK) / p.trSteps_wall;
	maxNumBLRoll = p.getParameter("Max Num BL Roll", MAX_NUM_BL_ROLL);

	lsSpacing = p.getParameter("LS Spacing", LS_SPACING);

	// When initializing various x location variables, code checks to
	// make sure positions make sense and throws error if not.
	stopDistX = p.getParameter("X Stop Distance", STOP_DIST_X);
	xReleaseVal = p.getParameter("X Release Pos", X_RELEASE_POS);
	xReleasePos = (int)floor(xReleaseVal);
	xMinVal = p.getParameter("X Min Val", xReleaseVal - 1.25); 

	double xRelDiff = xReleaseVal - xMinVal;
	ERROR_CHECKING(((xReleaseVal - xMinVal) <= 0.), "xMinVal >= xReleaseVal", ERROR_INITIALIZING_VTR);
	ERROR_CHECKING(((xReleaseVal - xMinVal) > 10.), "xMinVal - xReleaseVal > 10, which is an "
		"unnecessarily large difference", ERROR_INITIALIZING_VTR);

	xStopVal = p.getParameter("X Stop Pos", X_STOP_POS);
	xStopPos = (int)ceil(xStopVal);
	ERROR_CHECKING(((xStopVal - xReleaseVal) < 100.), "xStopVal - reduceDepStop1 < 100",
		ERROR_INITIALIZING_VTR);

	reduceDepStop1 = p.getParameter("Reduce Dep Stop1", REDUCE_DEP_STOP1);
	ERROR_CHECKING(((reduceDepStop1 - xReleaseVal) < 0.), "xReleaseVal > reduceDepStop1",
		ERROR_INITIALIZING_VTR);

	reduceDepStop2 = p.getParameter("Reduce Dep Stop2", REDUCE_DEP_STOP2);
	ERROR_CHECKING(((reduceDepStop2 - reduceDepStop1) < 0.), "reduceDepStop1 > reduceDepStop2",
		ERROR_INITIALIZING_VTR);

	startThermoVel = p.getParameter("Start Thermo Pos", START_THERMO_VEL);
	ERROR_CHECKING(((startThermoVel - xReleaseVal) < 0.), "xReleaseVal > startThermoVel",
		ERROR_INITIALIZING_VTR);

	amtReduceDep = p.getParameter("Amount Reduce Dep", AMT_REDUCE_DEP);
	ERROR_CHECKING((amtReduceDep <= 0.) || (amtReduceDep > 1.), "amtReduceDep must be (0,1]",
		ERROR_INITIALIZING_VTR);

	cutoffRadius = p.getParameter("Cutoff Radius", CUTOFF_RADIUS);

	blSearchRad = p.getParameter("BL_SEARCH_RAD", BL_SEARCH_RAD);

	foulSizeSwitchSoot = p.getParameter("Foul Size Switch Soot", FOUL_SIZE_SWITCH_SOOT2);

	trDomainXBounds = { {(int)floor(xMinVal), (int)floor(xStopVal)} };

	trDomainSize.x = (trDomainXBounds.y - trDomainXBounds.x);
	trDomainSize.y = p.nY - 1;
	
	testRestartRun();
}

// TODO: make sure to check that all files needed for a restart are loaded here
bool clVariablesTR::testRestartRun()
{
	allocateArrays();

	restartRunFlag = p.getParameter("Restart Run", false);
	
	if (!restartRunFlag)
		return;
	
	// Check that parameters necessary for restart are found in the parameter file
	restartRunFlag &= p.getParameter("Num Active Nodes", nActiveNodes, 0);
	restartRunFlag &= p.getParameter("Min BL Bottom", Bounds.MIN_BL_BOT, 0);
	restartRunFlag &= p.getParameter("Max BL Bottom", Bounds.MAX_BL_BOT, 0);
	restartRunFlag &= p.getParameter("Min BL Top", Bounds.MIN_BL_TOP, 0);
	restartRunFlag &= p.getParameter("Max BL Top", Bounds.MAX_BL_TOP, 0);
	
	// Check that bin/text files associated with arrays and particleStructs which 
	// are necessary for the restart are found and loaded
	restartRunFlag &= blDep.load("load" SLASH "BLdep_temp");
	restartRunFlag &= P.load();

	// Check that subclasses return true from their testRestartRun functions
	restartRunFlag &= parSort.testRestartRun();
	restartRunFlag &= wallShear.testRestartRun();
	restartRunFlag &= parP.testRestartRun();
	restartRunFlag &= glParticles.testRestartRun();

	return restartRunFlag;
}



void clVariablesTR::saveParams()
{
	p.setParameter("Use Par Solver", trSolverFlag);

	if (p.Time > 0)
		p.setParameter("Restart Run", false);
	else
		p.setParameter("Restart Run", true);

	p.setParameter("Num Particles", nN);
	p.setParameter("Max Out Lines IO", maxOutputLinesIO);
	p.setParameter("Calc IO Dists", calcIOFlag);
	p.setParameter("Index Radius Search", indRadiusSearch);
	p.setParameter("Mass Flux Inlet", massFluxInlet);
	p.setParameter("Max BL Per Node", maxBLPerNode);
	p.setParameter("Particle Release Time", parReleaseTime);
	p.setParameter("Time Before Re-release", timeBeforeReRelease * p.trSteps_wall);
	p.setParameter("Time Before Stick", timeBeforeStick * p.trSteps_wall);
	p.setParameter("Max Num BL Roll", maxNumBLRoll);
	p.setParameter("X Release Pos", xReleaseVal);
	p.setParameter("X Stop Pos", xStopVal);
	p.setParameter("X Stop Distance", stopDistX);
	p.setParameter("Reduce Dep Stop1", reduceDepStop1);
	p.setParameter("Reduce Dep Stop2", reduceDepStop2);
	p.setParameter("Start Thermo Pos", startThermoVel);
	p.setParameter("Amount Reduce Dep", amtReduceDep);
	p.setParameter("X Min Val", xMinVal);
	p.setParameter("Cutoff Radius", cutoffRadius);
	p.setParameter("Foul Size Switch Soot", foulSizeSwitchSoot);
	p.setParameter("LS Spacing", lsSpacing);
	p.setParameter("BL_SEARCH_RAD", blSearchRad);

	parSort.saveParams();
	wallShear.saveParams();
	parP.saveParams();
	glParticles.saveParams();
}




void clVariablesTR::allocateArrays()
{
	trDomainFullSizeX = getGlobalSizeMacro(trDomainSize.x, WORKGROUPSIZE_TR);

	// Allocate particleStructs
	P.setSizes(nN, nN, 1);
	P.allocateArrays();

	BL.setSizes(vls.nBL, vls.nBL, 1);
	BL.allocateArrays();

	NodC.setSizes(trDomainSize.x, trDomainFullSizeX, trDomainSize.y);
	NodC.allocateArrays();

	NodI.setSizes(trDomainSize.x, trDomainFullSizeX, trDomainSize.y);
	NodI.allocateArrays();
	
	NodV.setSizes(trDomainSize.x, trDomainFullSizeX, trDomainSize.y);
	NodV.allocateArrays();

	// Allocate Arrays 
	activeInds.allocateDynamic(trDomainSize.x * trDomainSize.y);
	
	wallInds.allocateDynamic(trDomainSize.x * trDomainSize.y / 10);
	
	/*activeFlag.zeros(trDomainSize.x, trDomainFullSizeX, trDomainSize.y, trDomainSize.y);*/
	
	blDep.zeros(vls.nBL, parP.Nd);
	
	PV.allocate(nN);
	PV.fill({ {0., 0.} });

	RandList.allocate(nN);
	
	reflectInfo.zeros(nN);

	reflectInds.zeros(3);
	
	updateFlag.zeros(nN);

	if (calcIOFlag)
	{
		ioDistsSave.zeros(parP.Nd*2, maxOutputLinesIO);
		ioDistsSave.createTimeArray();
	}


	// Not sure if these will be used, so keeping commented for now
	// blIndInd;
	// nodeNeigh;
	// numWallNodes;
	// trIndicies;

	parP.allocateArrays();
	wallShear.allocateArrays();
	parSort.allocateArrays();
	glParticles.allocateArrays();
}

void clVariablesTR::iniBLinksBottom()
{
	for (int i = 0; i < vls.nBL / 2; i++)
	{
		legacyBLinks bltemp;

		bltemp.P01ind = { {static_cast<cl_ushort>(vls.BL(i, 0)),
						  static_cast<cl_ushort>(vls.BL(i, 1))} };

		cl_double2 vP0, vP1;

		vP0 = { {vls.C[vls.BL(i, 0)].x, vls.C[vls.BL(i, 0)].y} };
		vP1 = { {vls.C[vls.BL(i, 1)].x, vls.C[vls.BL(i, 1)].y} };

		cl_double2 tanvec = Subtract2(vP1, vP0);
		cl_double2 normvec = { { -tanvec.y, tanvec.x } };

		bltemp.blLen = GETLEN(tanvec);
		bltemp.vNvec = Divide2(normvec, bltemp.blLen);
		tanvec = Divide2(tanvec, bltemp.blLen);
		bltemp.Tau = 0.;

		cl_double2 vP0not = { { vls.C0[vls.BL(i, 0)].x, vls.C0[vls.BL(i, 0)].y } };
		cl_double2 vP1not = { { vls.C0[vls.BL(i, 1)].x, vls.C0[vls.BL(i, 1)].y } };
		cl_double2 tanvecnot = Subtract2(vP1not, vP0not);
		tanvecnot = Divide2(tanvecnot, 2.);
		cl_double2 centernot = Subtract2(vP1not, tanvecnot);

		cl_double2 centpos = { { vP0.x + tanvec.x * bltemp.blLen / 2.,
			vP0.y + tanvec.y * bltemp.blLen / 2. } };

		cl_double2 centerdist = Subtract2(centpos, centernot);

		if (GETLEN(centerdist) <= FOUL_SIZE_SWITCH_SOOT2)
		{
			bltemp.int_type = 0;
		}
		else
		{
			bltemp.int_type = 1;
		}
		int xind = (int)floor(centpos.x);
		int yind = (int)floor(centpos.y);


		// Since this will most often be used to refer to indicies in particle structs,
		// Node_loc will refer to linear index of those arrays rather than the full
		// domain 
		if (xind >= trDomainXBounds.x && xind < trDomainXBounds.y)
			bltemp.Node_loc = (xind - trDomainXBounds.x) + (trDomainFullSizeX) * yind;
		else
			bltemp.Node_loc = -1;

		BL.setStruct(bltemp, i);
	}
}

void clVariablesTR::iniBLinksTop()
{
	for (int ii = vls.nBL - 1; ii >= vls.nBL / 2; ii--)
	{
		int i = (vls.nBL - 1 - ii) + vls.nBL / 2;

		legacyBLinks bltemp;

		bltemp.P01ind = { {static_cast<cl_ushort>(vls.BL(i, 1)),
						  static_cast<cl_ushort>(vls.BL(i, 0))} };

		cl_double2 vP0, vP1;

		vP0 = { {vls.C[vls.BL(i, 1)].x, vls.C[vls.BL(i, 1)].y} };
		vP1 = { {vls.C[vls.BL(i, 0)].x, vls.C[vls.BL(i, 0)].y} };

		cl_double2 tanvec = Subtract2(vP1, vP0);
		cl_double2 normvec = { { tanvec.y, -tanvec.x } };

		bltemp.blLen = GETLEN(tanvec);
		bltemp.vNvec = Divide2(normvec, bltemp.blLen);
		tanvec = Divide2(tanvec, bltemp.blLen);
		bltemp.Tau = 0.;

		cl_double2 centpos = { { vP0.x + tanvec.x * bltemp.blLen / 2.,
			vP0.y + tanvec.y * bltemp.blLen / 2. } };

		cl_double2 vP0not = { { vls.C0[vls.BL(i, 1)].x, vls.C0[vls.BL(i, 1)].y } };
		cl_double2 vP1not = { { vls.C0[vls.BL(i, 0)].x, vls.C0[vls.BL(i, 0)].y } };
		cl_double2 tanvecnot = Subtract2(vP1not, vP0not);
		tanvecnot = Divide2(tanvecnot, 2.);
		cl_double2 centernot = Subtract2(vP1not, tanvecnot);

		cl_double2 centerdist = Subtract2(centpos, centernot);

		if (GETLEN(centerdist) <= FOUL_SIZE_SWITCH_SOOT2)
		{
			bltemp.int_type = 0;
		}
		else
		{
			bltemp.int_type = 1;
		}

		int xind = (int)floor(centpos.x);
		int yind = (int)floor(centpos.y);

		// Since this will most often be used to refer to indicies in particle structs,
		// Node_loc will refer to linear index of those arrays rather than the full
		// domain 
		if (xind >= trDomainXBounds.x && xind < trDomainXBounds.y)
			bltemp.Node_loc = (xind - trDomainXBounds.x) + (trDomainFullSizeX) * yind;
		else
			bltemp.Node_loc = -1;

		BL.setStruct(bltemp, i);
	}
}


void clVariablesTR::iniNodI()
{
	// Fill out wallFlag info, which contains info about 
	// node type, neighbors, and if BL is located within node
	NodI.wallFlag.fill(WF_EMPTY);

	for (int i = 0; i < trDomainSize.x; i++)
	{
		int ii = i + trDomainXBounds.x;
		for (int j = 0; i < trDomainSize.y; i++)
		{
			// test all 4 LB nodes to see if they are
			// solid or fluid and set appropriate
			// flags
			if (vls.M(ii, j) & M_SOLID_NODE)
				NodI.wallFlag(i, j) |= WF_00_SOLID;
			if (vls.M(ii + 1, j) & M_SOLID_NODE)
				NodI.wallFlag(i, j) |= WF_10_SOLID;
			if (vls.M(ii, j + 1) & M_SOLID_NODE)
				NodI.wallFlag(i, j) |= WF_01_SOLID;
			if (vls.M(ii + 1, j + 1) & M_SOLID_NODE)
				NodI.wallFlag(i, j) |= WF_11_SOLID;

			// if all are solid, set node as solid
			if (NodI.wallFlag(i, j) == WF_TEST_ALL_SOLID)
			{
				NodI.wallFlag(i, j) |= WF_SOLID;
			}

			// if any are not solid, set as fluid node
			if (NodI.wallFlag(i, j) ^= WF_TEST_ALL_SOLID)
			{
				NodI.wallFlag(i, j) |= WF_FLUID;
			}
		}
	}
}


void clVariablesTR::iniBLandNodI()
{
	// BLinks contains all elements of vls.BL, so all nBL elements will
	// be initialized, since it doesnt really matter
	iniBLinksBottom();
	iniBLinksTop();

	

	// the remaining particleStructs have sizes = trDomainSize, so
	// need to account for this when initializing other structs.
	int xstart = trDomainXBounds.x;
	int xstop = trDomainXBounds.y;



	// going through BL list to find the first BLs
	// located at the tracer release location.
	// since the top and bottom BLs share identical starting
	// x locations, the top index is just the bottom start
	// index + nBL/2
	for (int i = 0; i < vls.nBL / 2; i++)
	{
		cl_double c0t = vls.C[vls.BL(i, 0)].x;
		cl_double c1t = vls.C[vls.BL(i, 1)].x;
		cl_double pMiddle = static_cast<int>((c0t + c1t)/2.);
			   
		if (pMiddle >= xstart)
		{
			Bounds.MIN_BL_BOT = i;
			Bounds.MIN_BL_TOP = vls.nBL / 2 + i;
			break;
		}
	}

	// getting end indicies for BL array. As with start, they
	// are separated by nBL/2 on top and bottom.
	for (int i = Bounds.MIN_BL_BOT + 1; i < vls.nBL / 2; i++)
	{
		cl_double c0t = vls.C[vls.BL(i, 0)].x;
		cl_double c1t = vls.C[vls.BL(i, 1)].x;
		cl_double pMiddle = static_cast<int>((c0t + c1t) / 2.);

		// breaks once the middle of a BL passes beyond xstop
		// leaving i = MAX_BL_BOT, which allows indexing to 
		// remaining similar to zero based array indexing
		// i.e. will use i < MAX_BL_BOT in for loops, 
		// length = MAX_BL_BOT-MIN_BL_BOT, etc
		if (pMiddle > xstop)
		{
			Bounds.MAX_BL_BOT = i;
			Bounds.MAX_BL_TOP = vls.nBL / 2 + i;
			break;
		}
		
	}

	Array2Dd dist2Center(trDomainSize.x, trDomainSize.y, "Dist2Center");
	dist2Center.fill(100.);
	NodI.BLind.fill(-1);

	//Fill NodI array with BL info from bottom wall
	for (int i = Bounds.MIN_BL_BOT; i < Bounds.MIN_BL_BOT; i++)
	{
		fillNodeI(i, dist2Center);
		fillNodeI(i + vls.nBL / 2, dist2Center);
	}
}

void clVariablesTR::fillNodeI(int i, Array2Dd &dist2Center)
{
	int n0 = static_cast<int>(BL.P01ind(i).x);
	int n1 = static_cast<int>(BL.P01ind(i).y);

	// get the corners defining a region surrounding the BLink
	// with a small amount of padding on either side (reason of 0.05*vTvecTemp)
	cl_double2 c0t = vls.C[n0], c1t = vls.C[n1];

	cl_double2 vTvecTemp = Subtract2(c1t, c0t);
	vTvecTemp = VNORM(vTvecTemp);

	cl_double2 c0 = { { c0t.x - vTvecTemp.x * 0.05, c0t.y - vTvecTemp.y * 0.05 } };
	cl_double2 c1 = { { c1t.x + vTvecTemp.x * 0.05, c1t.y + vTvecTemp.y * 0.05 } };

	cl_int2 Cmin = vls.min2(c0, c1);
	cl_int2 Cmax = vls.max2(c0, c1);

	cl_int2 trCmin = { { Cmin.x - trDomainXBounds.x, Cmin.y} };
	cl_int2 trCmax = { { Cmax.x - trDomainXBounds.x, Cmax.y} };

	// Make sure to avoid accesing OOB indicies
	if (trCmin.x < 0)
	{
		trCmin.x = 0;
		Cmin.x = trDomainXBounds.x;
	}
	if (trCmax.x >= trDomainXBounds.y)
	{
		trCmin.x = trDomainXBounds.y - 1;
		Cmin.x = trDomainXBounds.y - 1 + trDomainXBounds.x;
	}

	if (Cmin.y < 0)
	{
		Cmin.y = 0;
		trCmin.y = 0;
	}

	if (Cmax.y >= trDomainSize.y)
	{
		Cmax.y = trDomainSize.y - 1;
		trCmax.y = Cmax.y;
	}

	for (int i0 = trCmin.x; i0 < trCmax.x; i0++)
	{
		for (int j0 = trCmin.y; j0 < trCmax.y; j0++)
		{
			testNode(i0, j0, c0t, c1t, i, dist2Center);
		}
	}
}

void clVariablesTR::testNode(int i0, int j0, cl_double2 c0t, cl_double2 c1t, 
	int blind, Array2Dd& dist2Center)
{
	int nind = i0 + trDomainFullSizeX * j0;
	cl_double2 nCenter = { {(double)(i0 +
			trDomainXBounds.x) + 0.5, (double)(j0)+0.5} };
	cl_double2 pCenter = Add2(c0t, c1t);
	pCenter = Divide2(pCenter, 2.);
	cl_double2 distVec = Subtract2(nCenter, pCenter);
	double distTemp = GETLEN(distVec);



	if (distTemp < dist2Center(i0, j0))
	{
		NodI.BLind(i0, j0) = blind;
		NodI.wallFlag(i0, j0) |= (blind > vls.nBL/2) ? WF_TOP_WALL : WF_BOT_WALL;
	}
	
}

#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(),
void clVariablesTR::setSourceDefines()
{
	setSrcDefinePrefix "BL_SEARCH_RAD", blSearchRad);
	setSrcDefinePrefix "TR_X_IND_START", trDomainXBounds.x);
	setSrcDefinePrefix "START_THERMO_VEL", startThermoVel);
	setSrcDefinePrefix "MAX_BL_PER_NODE", maxBLPerNode);
	setSrcDefinePrefix "TRC_NUM_TRACERS", nN);
	setSrcDefinePrefix "FULLSIZEX_TR", trDomainSize.x);
	setSrcDefinePrefix "FULLSIZEY_TR", trDomainSize.y);
	setSrcDefinePrefix "FULLSIZEX_TR_PADDED", trDomainFullSizeX);
	setSrcDefinePrefix "MU_NUMBER", vlb.MuVal);
	setSrcDefinePrefix "DTTR", p.dTtr);
	setSrcDefinePrefix "DTTR_WALL", p.dTtr_wall);
	setSrcDefinePrefix "FOUL_SIZE_SWITCH_SOOT2", foulSizeSwitchSoot);
	setSrcDefinePrefix "KAIR", vfd.kAir* p.DELTA_F / p.DELTA_T);
	setSrcDefinePrefix "KSOOT", vfd.kSoot* p.DELTA_F / p.DELTA_T);
	setSrcDefinePrefix "ALPHA_DEN_AIR", vfd.Alpha_den_air);
	setSrcDefinePrefix "ALPHA_DEN_SOOT", vfd.Alpha_den_soot);
	setSrcDefinePrefix "ALPHA_TEST_AIR", ((vfd.kAir* p.DELTA_F / p.DELTA_T) / vfd.Alpha_den_air));
	setSrcDefinePrefix "ALPHA_TEST_SOOT", ((vfd.kSoot* p.DELTA_F / p.DELTA_T) / vfd.Alpha_den_soot));
	setSrcDefinePrefix "X_START_VAL_INT", xReleasePos);
	setSrcDefinePrefix "X_MAX_VAL_INT", xStopPos);
	setSrcDefinePrefix "X_START_VAL", (double)xReleasePos);
	setSrcDefinePrefix "X_MIN_VAL", (double)xStopPos - 0.5);
	setSrcDefinePrefix "X_MAX_VAL", xStopPos);
	setSrcDefinePrefix "Y_MIN_VAL", 0.5);
	setSrcDefinePrefix "Y_MAX_VAL", p.nY - 0.5);
	setSrcDefinePrefix "XVALS_HEIGHT", p.Channel_Height + 2);
	setSrcDefinePrefix "ALPHA_FLUID", vfd.Alpha_fluid);
	setSrcDefinePrefix "ALPHA_FOUL", vfd.Alpha_foul);

	int wallTimer = p.trSteps_wall;
	int trTimer = p.trSteps;
	int divTrSteps = trTimer / wallTimer;
	int reReleaseTimeTemp = MIN(timeBeforeReRelease, 65000); // dont want to go above max ushort value
	reReleaseTimeTemp /= wallTimer;
	int stickTimeTemp = MIN(timeBeforeStick, 65000); // dont want to go above max ushort value
	stickTimeTemp /= wallTimer;

	setSrcDefinePrefix "TIMER_DECREMENT", divTrSteps);
	// TIMER_DECREMENT for wall will be 1
	setSrcDefinePrefix "PAR_TIMER_START", reReleaseTimeTemp);
	setSrcDefinePrefix "DEP_TIMER_START", stickTimeTemp);

	setSrcDefinePrefix "MIN_BL_BOT", Bounds.MIN_BL_BOT);
	setSrcDefinePrefix "MIN_BL_TOP", Bounds.MIN_BL_TOP);
	setSrcDefinePrefix "MAX_BL_BOT", Bounds.MAX_BL_BOT);
	setSrcDefinePrefix "MAX_BL_TOP", Bounds.MAX_BL_TOP);

	setSrcDefinePrefix "REFLECT_INFO_OFFSET", nN/2);





	if (MAX_NUM_BL_ROLL != -1)
	{
		setSrcDefinePrefix "MAX_NUM_BL_ROLL", maxNumBLRoll);
	}

	parP.setSourceDefines();
	parSort.setSourceDefines();
	wallShear.setSourceDefines();
	glParticles.setSourceDefines();
}
#undef setSrcDefinePrefix

void clVariablesTR::createKernels()
{
	//trNodeKernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_nodes_temp");
	//trNodeKernel[0].set_size(trDomainSize.x, WORKGROUPSIZE_TR, trDomainSize.y, 1);

	//trNodeKernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_nodes_vel");
	//trNodeKernel[1].set_size(trDomainSize.x, WORKGROUPSIZE_TR, trDomainSize.y, 1);


	// The actual global size of any kernels below that are set with dummy 2*LocalSize
	// will be set before launching, (calling set_size w/ dummy arg to ensure dim is set)
	//trWallNodeKernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_nodes_along_wall_temp");
	//trWallNodeKernel[0].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	//trWallNodeKernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_nodes_along_wall_vel");
	//trWallNodeKernel[1].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	trUpdateParKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_no_wall");
	trUpdateParKernel.set_size(nActiveNodes, WORKGROUPSIZE_TR);

	trWallParKernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_along_wall");
	trWallParKernel[0].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);
	//Global size is set after num_wall_nodes is calculated

	trWallParKernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_reflect_particles");
	trWallParKernel[1].set_size(2 * WORKGROUPSIZE_TR_WALL_REFLECT, WORKGROUPSIZE_TR_WALL_REFLECT);

	trWallParKernel[2].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_contact_wall");
	trWallParKernel[2].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	trWallParKernel[3].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_update_par_contact_wall2");
	trWallParKernel[3].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

	trWallParKernel[4].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_deposit_particles_on_wall");
	trWallParKernel[4].set_size(2 * WORKGROUPSIZE_TR_WALL, WORKGROUPSIZE_TR_WALL);

//	saveDistsKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_save_dists");
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
	iniBLandNodI();
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

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////////////////////                                /////////////////////////
////////////////////     Kernel Intitialization     /////////////////////////
////////////////////                                /////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////






// Dont think this is necessary as all of these can be reconstructed
void clVariablesTR::saveRestartFilesIni()
{
	//Active_indicies.save_bin_from_device("Active_indicies");
	//Active_flag.save_bin_from_device("Active_flag");
	//TR_indicies.save_bin_from_device("TR_indicies");
}

// TODO: add necessary calls to subclasses and set those functions
// up in the subclasses
void clVariablesTR::saveRestartFiles()
{
	P.saveFromDevice(true, P.saveBinFl);

	parSort.saveRestartFiles();
	wallShear.saveRestartFiles();
	glParticles.saveRestartFiles();
	parP.saveRestartFiles();
}

// TODO: decide what to save and set it up
void clVariablesTR::save2file()
{

}

void clVariablesTR::saveDebug()
{
	//NodI.save_txt_from_device("NodI");
	//NodC.save_txt_from_device("NodC");
	//NodV.save_txt_from_device("NodV");
	//BL.save_txt_from_device("BL");
	//BLind_ind.save_txt_from_device("BLind_ind");
	//Node_neigh.save_txt_from_device("Node_neigh");
	//Winds.save_txt_from_device("Winds");
	//TR_indicies.save_txt_from_device("TR_indicies");
	//Active_indicies.save_txt_from_device("Active_indicies");
	//Active_flag.save_txt_from_device("Active_flag");
}


// TODO: add check to make sure that total particles reflecting reflectInds[0]
// is not more than 1/8 of total particles
void clVariablesTR::updateWallParticles()
{
	// Update locations of tracers near walls
	trWallParKernel[0].call_kernel(TRQUEUE_REF);

	// read counts of number of particles crossing a wall
	reflectInds.read_from_buffer(TRQUEUE_REF, CL_TRUE);

	// if particles cross wall behind release location, 
	// this function will reflect without allowing for deposition
	if (reflectInds(0) > 0)
	{
		// set number of particles to reflect argument and global size
		// before calling kernel
		trWallParKernel[1].setOptionGlobalCallKernel(reflectInds(0));
		
		// Fill index 0 with 0
		reflectInds.FillBuffer(0, sizeof(int), 1, 0, TRQUEUE_REF);
	}

	// if no particles cross wall, return
	if (reflectInds(1) == 0)
		return;

	// continue until no tracers have crossed a wall
	while (TRUE)
	{
		// set number of particles to reflect argument and global size
		// before calling kernel
		trWallParKernel[2].setOptionGlobalCallKernel(reflectInds(1));

		// Read reflectInds to see if any more particles need reflecting,
		// if not exit loop
		reflectInds.read_from_buffer(TRQUEUE_REF, CL_TRUE);
		if (reflectInds(0) == 0)
			break;

		// fill reflectInds(1) with 0, since this will be storing
		// the count of particles reflecting from wall in 
		// trWallParKernel[3]
		reflectInds.FillBuffer(0, sizeof(int), 1, 1, TRQUEUE_REF);

		// set number of particles to reflect argument and global size
		// before calling kernel
		trWallParKernel[3].setOptionGlobalCallKernel(reflectInds(0));

		// check to see if there are still particles to reflect
		reflectInds.read_from_buffer(TRQUEUE_REF, CL_TRUE);
		if (reflectInds(1) == 0)
			break;

		// Fill index 0 with 0, since it will be used to count
		// number of particles reflecting off of walls
		reflectInds.FillBuffer(0, sizeof(int), 1, 0, TRQUEUE_REF);
	}

	// If there are particles that deposited, call trWallParKernel[4]
	if (reflectInds(2) > 0)
	{
		// set number of particles to reflect argument and global size
		// before calling kernel
		trWallParKernel[2].setOptionGlobalCallKernel(reflectInds(2));
	}
	// Go ahead and fill reflectInds with 0s for next time.
	reflectInds.FillBuffer(0, TRQUEUE_REF);
}



void clVariablesTR::setKernelArgs()
{
	cl_int ind = 0;

	NodC.setTempBuffers(trWallParKernel[0], ind);
	NodC.setVelBuffers(trWallParKernel[0], ind);
	
	Par::arrName arrList_[] = { Par::locArr, Par::posArr, 
		Par::typeArr, Par::depFlagArr, Par::timerArr };
	P.setBuffers(trWallParKernel[0], ind, arrList_, 5);

	trWallParKernel[0].set_argument(ind++, BL.P01ind.get_buf_add());
	trWallParKernel[0].set_argument(ind++, BL.vNvec.get_buf_add());
	trWallParKernel[0].set_argument(ind++, BL.blLen.get_buf_add());

	trWallParKernel[0].set_argument(ind++, NodI.wallFlag.get_buf_add());
	trWallParKernel[0].set_argument(ind++, NodI.BLind.get_buf_add());

	trWallParKernel[0].set_argument(ind++, wallInds.get_buf_add());
	trWallParKernel[0].set_argument(ind++, vlb.Ux_array.get_buf_add());
	trWallParKernel[0].set_argument(ind++, vlb.Uy_array.get_buf_add());
	trWallParKernel[0].set_argument(ind++, vfd.Temp.get_add_Macro());
	trWallParKernel[0].set_argument(ind++, vls.C.get_buf_add());
	trWallParKernel[0].set_argument(ind++, parSort.Ploc.get_buf_add());
	trWallParKernel[0].set_argument(ind++, updateFlag.get_buf_add());
	trWallParKernel[0].set_argument(ind++, reflectInfo.get_buf_add());
	trWallParKernel[0].set_argument(ind++, PV.get_buf_add());
	trWallParKernel[0].set_argument(ind++, reflectInds.get_buf_add());
	trWallParKernel[0].setOptionInd(ind);
	trWallParKernel[0].setOption(wallInds.curSizeAdd());
	

	ind = 0;
	trWallParKernel[1].set_argument(ind++, P.loc.get_buf_add());
	trWallParKernel[1].set_argument(ind++, P.pos.get_buf_add());
	trWallParKernel[1].set_argument(ind++, P.Dep_Flag.get_buf_add());
	trWallParKernel[1].set_argument(ind++, BL.vNvec.get_buf_add());
	trWallParKernel[1].set_argument(ind++, PV.get_buf_add());
	trWallParKernel[1].set_argument(ind++, reflectInfo.get_buf_add());
	trWallParKernel[1].setOptionInd(ind);

	// starting location for indicies stored in reflectInfo and PV arrays
	// and corresponding to reflectInds(1)
	int offsetVal = 

	ind = 0;
	trWallParKernel[2].set_argument(ind++, P.pos.get_buf_add());
	trWallParKernel[2].set_argument(ind++, P.type.get_buf_add());
	trWallParKernel[2].set_argument(ind++, P.Dep_Flag.get_buf_add());
	trWallParKernel[2].set_argument(ind++, BL.vNvec.get_buf_add());
	trWallParKernel[2].set_argument(ind++, BL.Tau.get_buf_add());
	trWallParKernel[2].set_argument(ind++, BL.int_type.get_buf_add());
	trWallParKernel[2].set_argument(ind++, PV.get_buf_add());
	trWallParKernel[2].set_argument(ind++, RandList.get_buf_add());
	trWallParKernel[2].set_argument(ind++, reflectInfo.get_buf_add());
	trWallParKernel[2].set_argument(ind++, reflectInds.get_buf_add());
	trWallParKernel[2].setOptionInd(ind);

	ind = 0;
	trWallParKernel[3].set_argument(ind++, P.loc.get_buf_add());
	trWallParKernel[3].set_argument(ind++, P.pos.get_buf_add());
	trWallParKernel[3].set_argument(ind++, P.Dep_Flag.get_buf_add());
	trWallParKernel[3].set_argument(ind++, BL.vNvec.get_buf_add());
	trWallParKernel[3].set_argument(ind++, BL.P01ind.get_buf_add());
	trWallParKernel[3].set_argument(ind++, BL.blLen.get_buf_add());
	trWallParKernel[3].set_argument(ind++, vls.C.get_buf_add());
	trWallParKernel[3].set_argument(ind++, PV.get_buf_add());
	trWallParKernel[3].set_argument(ind++, reflectInfo.get_buf_add());
	trWallParKernel[3].set_argument(ind++, reflectInds.get_buf_add());
	trWallParKernel[3].setOptionInd(ind);

	ind = 0;
	trWallParKernel[4].set_argument(ind++, P.loc.get_buf_add());
	trWallParKernel[4].set_argument(ind++, P.pos.get_buf_add());
	trWallParKernel[4].set_argument(ind++, P.Dep_Flag.get_buf_add());
	trWallParKernel[4].set_argument(ind++, P.Dep_timer.get_buf_add());
	trWallParKernel[4].set_argument(ind++, BL.P01ind.get_buf_add());
	trWallParKernel[4].set_argument(ind++, BL.Tau.get_buf_add());
	trWallParKernel[4].set_argument(ind++, vls.C.get_buf_add());
	trWallParKernel[4].set_argument(ind++, reflectInfo.get_buf_add());
	trWallParKernel[4].setOptionInd(ind);

	ind = 0;
	P.setBuffers(trUpdateParKernel, ind, arrList_, 5);
	trUpdateParKernel.set_argument(ind++, NodI.wallFlag.get_buf_add());
	trUpdateParKernel.set_argument(ind++, vlb.Ux_array.get_buf_add());
	trUpdateParKernel.set_argument(ind++, vlb.Uy_array.get_buf_add());
	trUpdateParKernel.set_argument(ind++, vfd.Temp.get_add_Macro());
	trUpdateParKernel.set_argument(ind++, vls.C.get_buf_add());
	trUpdateParKernel.set_argument(ind++, parSort.Ploc.get_buf_add());
	trUpdateParKernel.set_argument(ind++, updateFlag.get_buf_add());
	trUpdateParKernel.set_argument(ind++, activeInds.get_buf_add());
	trUpdateParKernel.setOptionInd(ind);
	trUpdateParKernel.setOption(activeInds.curSizeAdd());
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
double clVariablesTR::rand1()
{
	double randval = (double)rand();
	return randval / rmax;
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
		if (randval >= parP.D_dist[j] && randval < parP.D_dist[j + 1])
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
		legacyPar Ptemp;
		// Particles will be released a portion at a time using reReleasePar kernel 
		Ptemp.pos.x = 0.;
		Ptemp.pos.y = 0.;
		
			
		Ptemp.type = getType();
		Ptemp.timer = 1;
		Ptemp.Dep_Flag = -2;
		Ptemp.Num_rep = 0;
		Ptemp.Dep_timer = 1;
		cl_int2 Posi = { { (int)floor(Ptemp.pos.x), (int)floor(Ptemp.pos.y) } };

		Ptemp.loc = Posi.x + trDomainFullSizeX * Posi.y;
		P.setStruct(Ptemp, i);
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

	for (int i = trDomainXBounds.x; i < trDomainXBounds.y; i++)
	{
		for (int j = 0; j < p.nY-1; j++)
		{
			cl_int2 ii00 = { { i, j } };
			cl_int2 ii10 = { { i + 1, j } };
			cl_int2 ii01 = { { i, j + 1 } };
			cl_int2 ii11 = { { i + 1, j + 1 } };

			int tnum = 0;

			if (vls.M(ii00.x, ii00.y) & M_SOLID_NODE)
				tnum += 1;
			if (vls.M(ii10.x, ii10.y) & M_SOLID_NODE)
				tnum += 2;
			if (vls.M(ii11.x, ii11.y) & M_SOLID_NODE)
				tnum += 4;
			if (vls.M(ii01.x, ii01.y) & M_SOLID_NODE)
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


void clVariablesTR::iniNodC()
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


