// clVariablesFL.cpp: class for fouling layer
//
// (c) Zach Mills, 2015 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clVariablesFL.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"

void clVariablesFL::loadParams()
{
	flSolverFlag = p.getParameter("Use FL Solver", USE_FL_SOLVER);
	saveOnStartFlag = p.getParameter("FL Save On Start", FL_SAVE_ON_START);
	Num_IO_indicies = p.getParameter("Number IO Indicies", NUM_INLET_OUTLET_NODES);
	smoothingPct = p.getParameter("Smoothing Percent", PERCENT_USED_IN_SMOOTHING);
	neighsPerSideSmoothing = p.getParameter("Neighs Per Side Smoothing", NEIGHS_PER_SIDE_SMOOTHING);

	testRestartRun();
}

#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(),
void clVariablesFL::setSourceDefines()
{
	setSrcDefinePrefix "PERCENT_USED_IN_SMOOTHING", vfl.smoothingPct);
	setSrcDefinePrefix "PERCENT_NOT_USED_IN_SMOOTHING", (1. - vfl.smoothingPct));
	setSrcDefinePrefix "NEIGHS_PER_SIDE_SMOOTHING", vfl.neighsPerSideSmoothing);
}
#undef setSrcDefinePrefix

void clVariablesFL::saveParams()
{
	if (p.Time > 0)
		p.setParameter("Restart Run", false);
	else
		p.setParameter("Restart Run", true);

	p.setParameter("Use FL Solver", flSolverFlag);
	p.setParameter("Number Active Nodes", Num_active_nodes);
	p.getParameter("Steps Btw Smoothing", flTimePerSmooth);
	p.getParameter("FL Update Time", flTimePerUpdate);
	p.setParameter("Cur Smooth Timer", Smooth_timer);
	p.setParameter("Cur Fl Update Timer", FL_timer);
	p.setParameter("Neighs Per Side Smoothing", neighsPerSideSmoothing);
	*p.yamlOut << YAML::BeginMap;
	*p.yamlOut << YAML::Key << "IO End" << YAML::Value << IO_end << YAML::EndMap;
	p.setParameter("Number IO Indicies", Num_IO_indicies);
	p.setParameter("Smoothing Percent", smoothingPct);

}

void clVariablesFL::allocateArrays()
{
	FI.allocate(Num_active_nodes);
	BLdep_tot.zeros(Num_active_nodes, vtr.Nd);
	BLdep_tot_temp.zeros(Num_active_nodes, vtr.Nd);
}

void clVariablesFL::ini()
{
	setSourceDefines();
	if (!flSolverFlag)
		return;
	if (restartRunFlag == false)
	{
		ini_foulI();
	}

	ini_IO_vars();
	allocateBuffers();

	std::function<void(void)> createKerPtr = std::bind(&clVariablesFL::createKernels, this);
	std::function<void(void)> setArgsPtr = std::bind(&clVariablesFL::setKernelArgs, this);

	sourceGenerator::SourceInstance()->addIniFunction(createKerPtr, setArgsPtr);
	if (saveOnStartFlag)
		save2file();
}

void clVariablesFL::testRestartRun()
{	
	restartRunFlag = p.getParameter("Restart Run", false);
	if (restartRunFlag == false)
	{
		flTimePerUpdate = p.getParameter("FL Update Time", UPDATE_TIME);
		FL_timer = flTimePerUpdate;
		flTimePerSmooth = p.getParameter("Steps Btw Smoothing", NUM_UPDATES_BTW_SMOOTHING);
		Smooth_timer = flTimePerSmooth;

		if (flSolverFlag == false)
		{
			IO_end.x = 0;
			IO_end.y = 0;
			IO_end.z = 0;
			IO_end.w = 0;
			Num_active_nodes = 0;
			return;
		}
		int ind = 0;
		while (vls.C0[ind].x < vtr.xReleasePos)
		{
			ind++;
		}
		IO_end.x = ind;

		while (vls.C0[ind].x < vtr.xStopPos)
		{
			ind++;
		}

		IO_end.y = ind;

		ind = vls.nN - 1;

		while (vls.C0[ind].x < vtr.xReleasePos)
		{
			ind--;
		}
		IO_end.z = ind;

		while (vls.C0[ind].x < vtr.xStopPos)
		{
			ind--;
		}

		IO_end.w = ind;
		
		Num_active_nodes = (IO_end.y - IO_end.x) + (IO_end.z - IO_end.w) + 2;
		allocateArrays();
	}
	else
	{
		Num_active_nodes = p.getParameter<int>("Number Active Nodes");
		flTimePerSmooth = p.getParameter<int>("Steps Btw Smoothing");
		flTimePerUpdate = p.getParameter<int>("FL Update Time");
		Smooth_timer = p.getParameter<int>("Cur Smooth Timer");
		FL_timer = p.getParameter<int>("Cur Fl Update Timer");
		p.yamlIn["IO End"] >> IO_end;

		allocateArrays();

		ERROR_CHECKING(BLdep_tot.load("load" SLASH "BLdep_tot") == false,
			"BLdep_tot must be provided in load file to restart from checkpoint",
			ERROR_INITIALIZING_VFL);

		ERROR_CHECKING(FI.load("load" SLASH "fli") == false,
			"fli must be provided in load file to restart from checkpoint",
			ERROR_INITIALIZING_VFL);
	}
}

void clVariablesFL::createKernels()
{
	int gsize_FL0 = WORKGROUPSIZE_FL_SHIFT;
	while (gsize_FL0 < vls.nBL)
		gsize_FL0 += WORKGROUPSIZE_FL_SHIFT;

	int gsize_FD0 = WORKGROUPSIZE_UPDATEFD;
	while (gsize_FD0 < vls.nBL)
		gsize_FD0 += WORKGROUPSIZE_UPDATEFD;

	int gsize_LS0 = WORKGROUPSIZE_UPDATEM;
	while (gsize_LS0 < vls.nN - 1)
		gsize_LS0 += WORKGROUPSIZE_UPDATEM;
	
	int gsizex_LS1 = (int)ceil(double(p.nX) / WORKGROUPSIZEX_LB)*WORKGROUPSIZEX_LB;
	int gsizey_LS1 = (int)ceil(double(p.nY) / WORKGROUPSIZEY_LB)*WORKGROUPSIZEY_LB;

	int gsizex_LS4_LB = (int)ceil((double)p.nX / WORKGROUPSIZEX_LB) * WORKGROUPSIZEX_LB;
	int gsizey_LS4_LB = (int)ceil((double)p.Channel_Height / WORKGROUPSIZEY_LB) * WORKGROUPSIZEY_LB;

	int gsizex_FD1 = (int)ceil((double)p.nX / WORKGROUPSIZEX_FD) * WORKGROUPSIZEX_FD;
	int gsizey_FD1 = (int)ceil((double)p.Channel_Height / WORKGROUPSIZEY_FD) * WORKGROUPSIZEY_FD;

	

	int gsizex_TR0 = (int)ceil(double(vtr.nActiveNodes) / (double)WORKGROUPSIZE_TR) *
		WORKGROUPSIZE_TR;

	int gsize_TR1 = WORKGROUPSIZE_TR_WALL;
	while (gsize_TR1 < vls.nBL)
		gsize_TR1 += WORKGROUPSIZE_TR_WALL;

	int gsizex_TR3 = (int)ceil(double(vtr.nActiveNodes) / (double)WORKGROUPSIZE_RERELEASE) *
		WORKGROUPSIZE_RERELEASE;

	int gsizex_TR4 = (int)ceil(double(vtr.nActiveNodes) / (double)WORKGROUPSIZE_UPDATEWALL) *
		WORKGROUPSIZE_UPDATEWALL;

	update_FL_kernel[0].create_kernel(GetSourceProgram, LBQUEUE_REF, "Shift_walls");
	update_FL_kernel[0].set_size(gsize_FL0, WORKGROUPSIZE_FL_SHIFT);

	update_FL_kernel[1].create_kernel(GetSourceProgram, LBQUEUE_REF, "Ramp_ends");
	update_FL_kernel[1].set_size(4 * Num_IO_indicies, Num_IO_indicies);

	update_FL_kernel[2].create_kernel(GetSourceProgram, LBQUEUE_REF, "Smooth_walls1");
	update_FL_kernel[2].set_size(gsize_FL0, WORKGROUPSIZE_FL_SHIFT);

	update_FL_kernel[3].create_kernel(GetSourceProgram, LBQUEUE_REF, "Smooth_walls2");
	update_FL_kernel[3].set_size(gsize_FL0, WORKGROUPSIZE_FL_SHIFT);

	update_LS_kernel[0].create_kernel(GetSourceProgram, LBQUEUE_REF, "Update_M");
	update_LS_kernel[0].set_size(gsize_LS0, WORKGROUPSIZE_UPDATEM);

	update_LS_kernel[1].create_kernel(GetSourceProgram, LBQUEUE_REF, "Fill_Sf_kernel");
	update_LS_kernel[1].set_size(gsizex_LS1, WORKGROUPSIZEX_LB, gsizey_LS1, WORKGROUPSIZEY_LB);

	update_LS_kernel[2].create_kernel(GetSourceProgram, LBQUEUE_REF, "Sum_M_arrays1");
	update_LS_kernel[2].set_size(p.FullSize / 2, WORKGROUPSIZE_RED);

	update_LS_kernel[3].create_kernel(GetSourceProgram, LBQUEUE_REF, "Sum_M_arrays2");
	update_LS_kernel[3].set_size(2, 1);
	
	update_LS_kernel[4].create_kernel(GetSourceProgram, LBQUEUE_REF, "Reduce_M_array");
	update_LS_kernel[4].set_size(gsizex_LS4_LB, WORKGROUPSIZEX_LB, gsizey_LS4_LB, WORKGROUPSIZEY_LB);

	update_LB_kernel[0].create_kernel(GetSourceProgram, LBQUEUE_REF, "update_LB_pos");
	update_LB_kernel[0].set_size(gsizex_LS4_LB, WORKGROUPSIZEX_LB, gsizey_LS4_LB, WORKGROUPSIZEY_LB);

	update_LB_kernel[1].create_kernel(GetSourceProgram, LBQUEUE_REF, "LB_Update_IBB");
	update_LB_kernel[1].set_local_size(WORKGROUPSIZE_IBB);

	update_FD_kernel.create_kernel(GetSourceProgram, FDQUEUE_REF, "Update_FD_dir_arrays");
	update_FD_kernel.set_size(gsize_FD0, WORKGROUPSIZE_UPDATEFD);

	update_TR_kernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "Update_BD_Nodes");
	update_TR_kernel[0].set_size(gsizex_TR0, WORKGROUPSIZE_TR);
	
	update_TR_kernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "Update_BD_Wall_Nodes");
	update_TR_kernel[1].set_size(gsize_TR1, WORKGROUPSIZE_TR_WALL);
	
	update_TR_kernel[2].create_kernel(GetSourceProgram, TRQUEUE_REF, "Shift_particles");
	update_TR_kernel[2].set_size(TRC_NUM_TRACERS, WORKGROUPSIZE_SHIFT_PAR);

	update_TR_kernel[3].create_kernel(GetSourceProgram, TRQUEUE_REF, "Shift_deposited_particles");
	update_TR_kernel[3].set_size(gsizex_TR3, WORKGROUPSIZE_RERELEASE);
	
	update_TR_kernel[4].create_kernel(GetSourceProgram, TRQUEUE_REF, "find_wall_nodes");
	/*update_TR_kernel[4].set_size(gsizex_TR4, WORKGROUPSIZE_UPDATEWALL);*/
	update_TR_kernel[4].set_size(1,1);


	if (p.useOpenGL)
	{
		update_GL_kernel.create_kernel(GetSourceProgram, IOQUEUE_REF, "update_GL_wall");
		update_GL_kernel.set_size((double)vls.nN / 2., WORKGROUPSIZE_UPDATE_GL);
	}


}

double clVariablesFL::Gaussian_Kernel(double input)
{
	return (1. / sqrt(2.*p.Pi) * exp(-0.5*input*input));
}

void clVariablesFL::ini_foulI()
{
	int fiind = 0;
	int blind = 0;

	cl_double8 ioWeight;
	
	for (int kk = 0; kk < 2; kk++)
	{
		ioWeight.s[kk] = 0.;
	}
	for (int kk = 2; kk < 8; kk++)
	{
		ioWeight.s[kk] = 1. / 6. * AMT_REDUCE_DEP;
	}


	double cutoff_rad = 4.5*LS_SPACING;

	while (TRUE)
	{
		blind++;
		if (vls.BL(blind, 0) == IO_end.x)
			break;
	}

	blind--;
	cl_uint bl_inlet = blind;
	for (int ii = 0; ii < 4; ii++)
	{
		foulI FItemp;
		int i = IO_end.x + ii;
		FItemp.C_ind = i;
		FItemp.BL_ind = bl_inlet;
		cl_double2 LP = { { vls.C0[i - 1].x, vls.C0[i - 1].y } };
		cl_double2 CP = { { vls.C0[i].x, vls.C0[i].y } };
		cl_double2 RP = { { vls.C0[i + 1].x, vls.C0[i + 1].y } };
		cl_double2 tanvec0 = Subtract2(CP, LP);
		cl_double2 normvec0 = { { -tanvec0.y, tanvec0.x } };
		cl_double2 tanvec1 = Subtract2(RP, CP);
		cl_double2 normvec1 = { { -tanvec1.y, tanvec1.x } };
		cl_double len0 = GETLEN(tanvec0);
		cl_double len1 = GETLEN(tanvec1);
		normvec0 = Divide2(normvec0, len0);
		normvec1 = Divide2(normvec1, len1);

		FItemp.vN = { { (normvec0.x + normvec1.x) / 2., (normvec0.y + normvec1.y) / 2. } };
		FItemp.disp = 0.;
		FItemp.WeightsL = ioWeight.lo;
		FItemp.WeightsR = ioWeight.hi;
		blind++;
#pragma warning(suppress: 6386)
		FI[fiind++] = FItemp;

	}

	int fiind_inlet = fiind;

	for (cl_uint i = IO_end.x + 4; i <= IO_end.y - 5; i++)
	{
		foulI FItemp;
		cl_uint bl_start = blind - 3;
		FItemp.C_ind = i;
		FItemp.BL_ind = bl_start;
		cl_double2 LP = { { vls.C0[i - 1].x, vls.C0[i - 1].y } };
		cl_double2 CP = { { vls.C0[i].x, vls.C0[i].y } };
		cl_double2 RP = { { vls.C0[i + 1].x, vls.C0[i + 1].y } };
		cl_double2 tanvec0 = Subtract2(CP,LP);
		cl_double2 normvec0 = { { -tanvec0.y, tanvec0.x } };
		cl_double2 tanvec1 = Subtract2(RP, CP);
		cl_double2 normvec1 = { { -tanvec1.y, tanvec1.x } };
		cl_double len0 = GETLEN(tanvec0);
		cl_double len1 = GETLEN(tanvec1);
		normvec0 = Divide2(normvec0, len0);
		normvec1 = Divide2(normvec1, len1);

		FItemp.vN = { { (normvec0.x + normvec1.x) / 2., (normvec0.y + normvec1.y) / 2. } };
		FItemp.disp = 0.;

		cl_double8 Wtemp;
		double sumW = 0.;
		for (int kk = 0; kk < 8; kk++)
		{
			cl_double2 neigh_pos = get_center(vls.C0[vls.BL(bl_start + kk, 0)], vls.C0[vls.BL(bl_start+kk, 1)]);
			cl_double2 dVec = Subtract2(CP, neigh_pos);
			double dist = GETLEN(dVec);
			double cutoff_temp = cutoff_rad;
			//if (dist >= cutoff_rad)
			//	cutoff_temp *= 1.25;
			double Si = dist / cutoff_temp;
			double wtt = Gaussian_Kernel(Si);

			double red_amt;
			if(neigh_pos.x < REDUCE_DEP_STOP1)
				red_amt = AMT_REDUCE_DEP;
			else if (neigh_pos.x < REDUCE_DEP_STOP2)
			{
				red_amt = AMT_REDUCE_DEP + (neigh_pos.x - REDUCE_DEP_STOP1) * (1. - AMT_REDUCE_DEP) / (REDUCE_DEP_STOP2 - REDUCE_DEP_STOP1);
			}
			else
				red_amt = 1.;

			Wtemp.s[kk] = wtt * red_amt;
			sumW += wtt;
		}
		for (int kk = 0; kk < 8; kk++)
		{
			Wtemp.s[kk] /= sumW;
		}
		FItemp.WeightsL = Wtemp.lo;
		FItemp.WeightsR = Wtemp.hi;
		
		blind++;
		FI[fiind++] = FItemp;
		if (fiind == fiind_inlet)
		{
			for (int zz = 0; zz < fiind_inlet; zz++)
			{
				FI[zz].WeightsL = FI[fiind].WeightsL;
				FI[zz].WeightsR = FI[fiind].WeightsR;
			}
		}
	}
	
	cl_uint bl_outlet = blind-3;
	for (cl_uint i = IO_end.y - 4; i <= IO_end.y; i++)
	{
		foulI FItemp;
		FItemp.C_ind = i;
		FItemp.BL_ind = bl_outlet;
		cl_double2 LP = { { vls.C0[i - 1].x, vls.C0[i - 1].y } };
		cl_double2 CP = { { vls.C0[i].x, vls.C0[i].y } };
		cl_double2 RP = { { vls.C0[i + 1].x, vls.C0[i + 1].y } };
		cl_double2 tanvec0 = Subtract2(CP, LP);
		cl_double2 normvec0 = { { -tanvec0.y, tanvec0.x } };
		cl_double2 tanvec1 = Subtract2(RP, CP);
		cl_double2 normvec1 = { { -tanvec1.y, tanvec1.x } };
		cl_double len0 = GETLEN(tanvec0);
		cl_double len1 = GETLEN(tanvec1);
		normvec0 = Divide2(normvec0, len0);
		normvec1 = Divide2(normvec1, len1);

		FItemp.vN = { { (normvec0.x + normvec1.x) / 2., (normvec0.y + normvec1.y) / 2. } };
		FItemp.disp = 0.;
		FItemp.WeightsL = { { 1. / 7., 1. / 7., 1. / 7., 1. / 7. } };
		FItemp.WeightsR = { { 1. / 7., 1. / 7., 1. / 7., 0. } };
		blind++;
		FI[fiind++] = FItemp;
	}

	
	blind = vls.nBL / 2;
	while (TRUE)
	{
		blind++;
		if (vtr.BL[blind].P1ind == IO_end.z)
			break;
	}

	bl_inlet = blind;

	int fiind_inlet_start = fiind;
	for (int ii = 0; ii < 4; ii++)
	{
		foulI FItemp;
		int i = IO_end.z - ii;
		FItemp.C_ind = i;
		FItemp.BL_ind = bl_inlet;
		cl_double2 RP = { { vls.C0[i - 1].x, vls.C0[i - 1].y } };
		cl_double2 CP = { { vls.C0[i].x, vls.C0[i].y } };
		cl_double2 LP = { { vls.C0[i + 1].x, vls.C0[i + 1].y } };
		cl_double2 tanvec0 = Subtract2(LP, CP);
		cl_double2 normvec0 = { { -tanvec0.y, tanvec0.x } };
		cl_double2 tanvec1 = Subtract2(CP, RP);
		cl_double2 normvec1 = { { -tanvec1.y, tanvec1.x } };
		cl_double len0 = GETLEN(tanvec0);
		cl_double len1 = GETLEN(tanvec1);
		normvec0 = Divide2(normvec0, len0);
		normvec1 = Divide2(normvec1, len1);

		FItemp.vN = { { (normvec0.x + normvec1.x) / 2., (normvec0.y + normvec1.y) / 2. } };
		FItemp.disp = 0.;
		FItemp.WeightsL = ioWeight.lo;
		FItemp.WeightsR = ioWeight.hi;
		blind++;
		FI[fiind++] = FItemp;

	}

	fiind_inlet = fiind;
	
	for (cl_uint i = IO_end.z - 4; i >= IO_end.w + 5; i--)
	{
		cl_uint bl_start = blind - 3;
		foulI FItemp;
		FItemp.C_ind = i;
		FItemp.BL_ind = bl_start;
		cl_double2 RP = { { vls.C0[i - 1].x, vls.C0[i - 1].y } };
		cl_double2 CP = { { vls.C0[i].x, vls.C0[i].y } };
		cl_double2 LP = { { vls.C0[i + 1].x, vls.C0[i + 1].y } };
		cl_double2 tanvec0 = Subtract2(LP, CP);
		cl_double2 normvec0 = { { -tanvec0.y, tanvec0.x } };
		cl_double2 tanvec1 = Subtract2(CP, RP);
		cl_double2 normvec1 = { { -tanvec1.y, tanvec1.x } };
		cl_double len0 = GETLEN(tanvec0);
		cl_double len1 = GETLEN(tanvec1);
		normvec0 = Divide2(normvec0, len0);
		normvec1 = Divide2(normvec1, len1);

		FItemp.vN = { { (normvec0.x + normvec1.x) / 2., (normvec0.y + normvec1.y) / 2. } };
		FItemp.disp = 0.;

		

		cl_double8 Wtemp;
		double sumW = 0.;
		for (int kk = 0; kk < 8; kk++)
		{
			cl_double2 neigh_pos = get_center(vtr.BL[bl_start + kk].vP0, vtr.BL[bl_start + kk].vP1);
			cl_double2 dVec = Subtract2(CP, neigh_pos);
			double dist = GETLEN(dVec);
			double cutoff_temp = cutoff_rad;
			//if (dist >= cutoff_rad)
			//	cutoff_temp *= 1.25;
			double Si = dist / cutoff_temp;
			double wtt = Gaussian_Kernel(Si);

			double red_amt;
			if (neigh_pos.x < REDUCE_DEP_STOP1)
				red_amt = AMT_REDUCE_DEP;
			else if (neigh_pos.x < REDUCE_DEP_STOP2)
			{
				red_amt = AMT_REDUCE_DEP + (neigh_pos.x - REDUCE_DEP_STOP1) * (1. - AMT_REDUCE_DEP) / (REDUCE_DEP_STOP2 - REDUCE_DEP_STOP1);
			}
			else
				red_amt = 1.;

			Wtemp.s[kk] = wtt * red_amt;
			sumW += wtt;
		}
		for (int kk = 0; kk < 8; kk++)
		{
			Wtemp.s[kk] /= sumW;
		}
		FItemp.WeightsL = Wtemp.lo;
		FItemp.WeightsR = Wtemp.hi;

		FI[fiind++] = FItemp;
		if (fiind == fiind_inlet)
		{
			for (int zz = fiind_inlet_start; zz < fiind_inlet; zz++)
			{
				FI[zz].WeightsL = FI[fiind].WeightsL;
				FI[zz].WeightsR = FI[fiind].WeightsR;
			}
		}
		blind++;
	}

	bl_outlet = blind-3;
	for (cl_uint i = IO_end.w + 4; i >= IO_end.w; i--)
	{
		foulI FItemp;
		FItemp.C_ind = i;
		FItemp.BL_ind = bl_outlet;
		cl_double2 RP = { { vls.C0[i - 1].x, vls.C0[i - 1].y } };
		cl_double2 CP = { { vls.C0[i].x, vls.C0[i].y } };
		cl_double2 LP = { { vls.C0[i + 1].x, vls.C0[i + 1].y } };
		cl_double2 tanvec0 = Subtract2(LP, CP);
		cl_double2 normvec0 = { { -tanvec0.y, tanvec0.x } };
		cl_double2 tanvec1 = Subtract2(CP, RP);
		cl_double2 normvec1 = { { -tanvec1.y, tanvec1.x } };
		cl_double len0 = GETLEN(tanvec0);
		cl_double len1 = GETLEN(tanvec1);
		normvec0 = Divide2(normvec0, len0);
		normvec1 = Divide2(normvec1, len1);

		FItemp.vN = { { (normvec0.x + normvec1.x) / 2., (normvec0.y + normvec1.y) / 2. } };
		FItemp.disp = 0.;
		FItemp.WeightsL = { { 1. / 7., 1. / 7., 1. / 7., 1. / 7. } };
		FItemp.WeightsR = { { 1. / 7., 1. / 7., 1. / 7., 0. } };
		FI[fiind++] = FItemp;
	}
}

cl_double2 clVariablesFL::get_center(cl_double2 P0, cl_double2 P1)
{
	cl_double2 temp = Add2(P0, P1);
	return Divide2(temp,2.);
}

void clVariablesFL::ini_IO_vars()
{
	Num_IO_indicies = NUM_INLET_OUTLET_NODES;
	IO_ind_dist.allocate(Num_IO_indicies * 4);
	RI.allocate(Num_IO_indicies*4);
	cl_uint4 cur_ind = IO_end;

	cl_double4 IO_endx;

	IO_endx.x = vls.C0[cur_ind.x].x;
	IO_endx.y = vls.C0[cur_ind.y].x;
	IO_endx.z = vls.C0[cur_ind.z].x;
	IO_endx.w = vls.C0[cur_ind.w].x;

	cl_uint4 Ind_dir;

	Ind_dir.x = -1;
	Ind_dir.y = 1;
	Ind_dir.z = 1;
	Ind_dir.w = -1;

	for (cl_uint i = 1; i <= Num_IO_indicies; i++)
	{
		cur_ind.x--;
		cur_ind.y++;
		cur_ind.z++;
		cur_ind.w--;
		
		RI(i-1).Ybegin = vls.C0[IO_end.x].y;
		RI(Num_IO_indicies + i - 1).Ybegin = vls.C0[IO_end.y].y;
		RI(Num_IO_indicies*2 + i - 1).Ybegin = vls.C0[IO_end.z].y;
		RI(Num_IO_indicies*3 + i - 1).Ybegin = vls.C0[IO_end.w].y;

		RI(i - 1).IOind = IO_end.x;
		RI(Num_IO_indicies + i - 1).IOind = IO_end.y;
		RI(Num_IO_indicies * 2 + i - 1).IOind = IO_end.z;
		RI(Num_IO_indicies * 3 + i - 1).IOind = IO_end.w;

		RI(i - 1).Coeff = vls.C0[cur_ind.x].x - IO_endx.x;
		RI[Num_IO_indicies + i - 1].Coeff = IO_endx.y - vls.C0[cur_ind.y].x;
		RI[Num_IO_indicies * 2 + i - 1].Coeff = IO_endx.z - vls.C0[cur_ind.z].x;
		RI[Num_IO_indicies * 3 + i - 1].Coeff = vls.C0[cur_ind.w].x - IO_endx.w;

		RI(i - 1).Cind = cur_ind.x;
		RI[Num_IO_indicies + i - 1].Cind = cur_ind.y;
		RI[Num_IO_indicies * 2 + i - 1].Cind = cur_ind.z;
		RI[Num_IO_indicies * 3 + i - 1].Cind = cur_ind.w;
	}

	cur_ind.x--;
	cur_ind.y++;
	cur_ind.z++;
	cur_ind.w--;

	cl_double4 tot_dist;
	tot_dist.x = IO_endx.x - vls.C0[cur_ind.x].x;
	tot_dist.y = IO_endx.y - vls.C0[cur_ind.y].x;
	tot_dist.z = IO_endx.z - vls.C0[cur_ind.z].x;
	tot_dist.w = IO_endx.w - vls.C0[cur_ind.w].x;

	for (cl_uint i = 0; i < Num_IO_indicies; i++)
	{
		RI(i).Coeff /= -tot_dist.x;
		RI[Num_IO_indicies + i].Coeff /= tot_dist.y;
		RI[Num_IO_indicies * 2 + i].Coeff /= tot_dist.z;
		RI[Num_IO_indicies * 3 + i].Coeff /= -tot_dist.w;
	}
}

void clVariablesFL::update_LS(cl_event *evt)
{
	update_LS_kernel[0].call_kernel();
	update_LS_kernel[1].call_kernel();
	update_LS_kernel[2].call_kernel();
	update_LS_kernel[3].call_kernel();
	update_LS_kernel[4].call_kernel(NULL, 0, NULL, evt);
}

void clVariablesFL::update_FL()
{
	//cl_event Evt1;
	//update_FL_kernel[0].call_kernel(&Evt1);
	////Smooth_timer--;
	////if (Smooth_timer == 0)
	////{
	////	update_FL_kernel[2].call_kernel();
	////	update_FL_kernel[3].call_kernel();
	////}
	//update_FL_kernel[1].call_kernel();
	//cl_uint zer = 0;
	//vtr.BL_dep.FillBuffer(IOQUEUE, 0, 1, &Evt1);
	//
	////if (Smooth_timer == 0)
	////{
	////	BLdep_tot.enqueue_copy_to_buffer(IOQUEUE, BLdep_tot_temp.get_buffer());
	////	Smooth_timer = NUM_UPDATES_BTW_SMOOTHING;
	////}

	//clReleaseEvent(Evt1);
}

void clVariablesFL::update_LB()
{
	update_LB_kernel[0].call_kernel();
}

void clVariablesFL::update()
{//TODO: test if case is handled when dX becomes > 1.
//	cl_event Red_T_Evt;
//	//vfd.FD_reduce_T.call_kernels(&Red_T_Evt);
//	//vfd.Sum_Temp.read_from_buffer(IOQUEUE, 1, &Red_T_Evt);
//	clReleaseEvent(Red_T_Evt);
//
//
//	cl_event LS_Evt, Fill_Evt;
//	
//	//vfd.dX_cur.enqueue_copy_to_buffer(FDQUEUE, vfd.dX.get_buffer());
//	
//	update_FL();
//	update_LS(&LS_Evt);
//		
//	vls.C.read_from_buffer(IOQUEUE, 1, &LS_Evt);
//	vls.M.read_from_buffer(IOQUEUE);
//
//	update_LB();
//
//#ifdef USE_OPENGL
//	update_GL();
//#endif
//
//	update_FD(&LS_Evt);
//
//	clFlush(FDQUEUE);
//	
//	nodeI Ntemp;
//	for (int i = 0; i < MAX_BL_PER_NODE; i++)
//		Ntemp.BLind[i] = -1;
//	Ntemp.Wall_Flag = 0;
//	vtr.NodI.FillBuffer(IOQUEUE, Ntemp);
//	vtr.BLind_ind.FillBuffer(IOQUEUE, 0);
//
//	int Num_Wnodes_temp = vtr.Num_wall_nodes_max;
//	vtr.Num_W_nodes.FillBuffer(IOQUEUE, 0, &Fill_Evt);
//
//	clFlush(IOQUEUE);
//
//	update_TR(&Fill_Evt, Num_Wnodes_temp);
//
//	p.flushQueues();
//
//	clReleaseEvent(LS_Evt);
//	clReleaseEvent(Fill_Evt);
}


bool clVariablesFL::test_bounds()
{
	//bool test = false;
	//test = vlb.Act.Test_Bounds_Buffer(IOQUEUE, 0, p.nY); 
	//test = vlb.Stor.Test_Bounds_Buffer(IOQUEUE, -1, p.Channel_Height);
	////test = vlb.Prop_loc.Test_Bounds_Buffer(IOQUEUE, 0, p.nX*p.Channel_Height*8);
	//

	//test = FI.Test_Bounds_Buffer(IOQUEUE, 0, vls.nN-1, 0, vls.nBL-1);
	//if (test)
	//{
	//	printf("FI error\n");
	//	FI.save_txt_from_device("FI", IOQUEUE);
	//}

	//test = vtr.P.Test_Bounds_Buffer(IOQUEUE, -2, (p.nX - 1)*(p.nY - 1)-1);
	//if (test)
	//{
	//	printf("P error\n");
	//	vtr.P.save_txt_from_device_full("Par", IOQUEUE);
	//}

	//test = vtr.NodI.Test_Bounds_Buffer(IOQUEUE, -1, vls.nBL-1);
	//if (test)
	//{
	//	printf("NodI error\n");
	//	vtr.NodI.save_txt_from_device("NodI", IOQUEUE);
	//}
	//
	//test = vtr.NodC.Test_Bounds_Buffer(IOQUEUE, 0, (p.nX - 1)*(p.nY - 1)-1);
	//if (test)
	//{
	//	printf("NodC error\n");
	//	vtr.NodC.save_txt_from_device("NodC", IOQUEUE);
	//}

	//test = vtr.Shear_inds.Test_Bounds_Buffer(IOQUEUE, { { 0, 0 } }, { { p.nX*p.Channel_Height - 1, p.nX*p.Channel_Height - 1 } });
	//if (test)
	//{
	//	printf("Shear_inds error %d\n", vls.nBnodes);

	//	vtr.Shear_inds.save_txt_from_device("Shear_inds", IOQUEUE);
	//}

	//test = vtr.Node_neigh.Test_Bounds_Buffer(IOQUEUE, -1, (p.nX-1)*(p.nY-1)-1);
	//if (test)
	//{
	//	printf("Node_neigh error\n");
	//	vtr.Node_neigh.save_txt_from_device("Node_neigh", IOQUEUE);
	//}

	//test = vtr.Winds.Test_Bounds_Buffer(IOQUEUE, 0, (p.nX - 1)*(p.nY - 1)-1);
	//if (test)
	//{
	//	printf("Winds error\n");
	//	vtr.Winds.save_from_device("Winds", IOQUEUE);
	//}

//	return test;
	return 0;
}


void clVariablesFL::freeHostMem()
{

}
void clVariablesFL::save_variables()
{
	vls.save2file();
	vlb.save2file();
	vfd.save2file();
	vtr.save2file();
	save2file();

	vtr.saveDebug();
	vls.saveDebug();
	vlb.saveDebug();
	vfd.saveDebug();
	vtr.P.save_txt_from_device_full("trc");
}

void clVariablesFL::CallRename(char *file, const char *fol)
{
	char Buf[80];
	char Buf2[80];
	sprintf(Buf2, "%s.txt", file);
	sprintf(Buf, "%s" SLASH "%s.txt", fol, file);
	p.RenameFile(Buf2, Buf);
}

void clVariablesFL::RenameDebug_Files(int dirnumber)
{
	std::string NewDir = "debugFiles_" + std::to_string(dirnumber);

	p.MakeDir(NewDir);
	CallRename("lsc", NewDir.c_str());
	CallRename("lbm", NewDir.c_str());
	CallRename("lbms", NewDir.c_str());
	CallRename("lbm_red", NewDir.c_str());
	CallRename("lbmsf", NewDir.c_str());
	CallRename("ii0_array", NewDir.c_str());
	CallRename("iis1_array", NewDir.c_str());
	CallRename("dir_array", NewDir.c_str());

	CallRename("D_array", NewDir.c_str());
	CallRename("fdt", NewDir.c_str());
	CallRename("dXcur", NewDir.c_str());
	CallRename("dX", NewDir.c_str());
	CallRename("Alpha", NewDir.c_str());
	CallRename("Acoeff", NewDir.c_str());
	CallRename("Bcoeff", NewDir.c_str());
	CallRename("Ccoeff", NewDir.c_str());
	CallRename("BLdep_tot", NewDir.c_str());

	CallRename("Prop_loc", NewDir.c_str());
	CallRename("IBB_loc", NewDir.c_str());
	CallRename("IBB_coeff", NewDir.c_str());

	CallRename("NodC", NewDir.c_str());
	CallRename("NodI", NewDir.c_str());
	CallRename("NodV", NewDir.c_str());
	CallRename("BL", NewDir.c_str());
	CallRename("BLind_ind", NewDir.c_str());
	CallRename("BLindicies", NewDir.c_str());
	CallRename("Weights", NewDir.c_str());
	CallRename("Shear_inds", NewDir.c_str());
	CallRename("Shear_coeffs", NewDir.c_str());
	CallRename("Ploc", NewDir.c_str());
	CallRename("Node_neigh", NewDir.c_str());
	CallRename("Winds", NewDir.c_str());
	CallRename("Active_indicies", NewDir.c_str());
	CallRename("TR_indicies", NewDir.c_str());
	CallRename("Active_flag", NewDir.c_str());
	CallRename("trc", NewDir.c_str());

}
void clVariablesFL::update_TR(cl_event *wait_fill, int Num_Wnodes_temp)
{

	cl_uint offset = vtr.Ploc(1).x;
	cl_uint total_removal = vtr.Ploc(1).y - offset;
	update_TR_kernel[3].set_argument(2, &offset);
	update_TR_kernel[3].set_argument(3, &total_removal);


	update_TR_kernel[0].call_kernel(NULL, 1, wait_fill);  ///Update Nodes
	update_TR_kernel[1].call_kernel();  ///Update Wall Nodes
	update_TR_kernel[2].call_kernel();	///Update Particles
	update_TR_kernel[3].call_kernel();	///Update Wall Particles

	vtr.Num_W_nodes.read_from_buffer();
	clFlush(TRQUEUE);

	clFinish(FDQUEUE);

	// May or may not need to keep this (will require re-implementing method in LBsolver)
	//vlb.update_IBB_arrays();
	vtr.update_Shear_arrays();
	vtr.update_trp();

	clFinish(TRQUEUE);
	vtr.Num_wall_nodes = vtr.Num_W_nodes(0);
	//cl_event Fill_Evt;
	//vtr.Num_W_nodes.FillBuffer(IOQUEUE, 0, &Fill_Evt);

	if (Num_Wnodes_temp < vtr.Num_W_nodes(0))
	{
		double gsize = (double)vtr.Num_W_nodes(0) / WORKGROUPSIZE_TR_WALL;
		vtr.Num_wall_nodes_max = (int)(ceil(gsize)*WORKGROUPSIZE_TR_WALL);

		vtr.TR_Wall_Par_kernel[0].reset_global_size(vtr.Num_W_nodes(0));
		vtr.TR_Wall_Node_kernel[0].reset_global_size(vtr.Num_W_nodes(0));
		vtr.TR_Wall_Node_kernel[1].reset_global_size(vtr.Num_W_nodes(0));


		vtr.Winds.reallocate(vtr.Num_wall_nodes_max);
		clFinish(IOQUEUE);

		update_TR_kernel[4].set_argument(1, vtr.Winds.get_buf_add());


		vtr.TR_Wall_Par_kernel[0].set_argument(7, vtr.Winds.get_buf_add());
		vtr.TR_Wall_Node_kernel[0].set_argument(4, vtr.Winds.get_buf_add());
		vtr.TR_Wall_Node_kernel[1].set_argument(4, vtr.Winds.get_buf_add());
	}

	vtr.TR_Wall_Par_kernel[0].set_argument(8, &vtr.Num_wall_nodes);
	vtr.TR_Wall_Node_kernel[0].set_argument(5, &vtr.Num_wall_nodes);
	vtr.TR_Wall_Node_kernel[1].set_argument(5, &vtr.Num_wall_nodes);
	
	////update_TR_kernel[4].call_kernel(1, &Fill_Evt);
	update_TR_kernel[4].call_kernel();

	clFlush(TRQUEUE);
	//clReleaseEvent(Fill_Evt);
}



void clVariablesFL::update_GL()
{
	update_GL_kernel.call_kernel();
}

void clVariablesFL::update_FD(cl_event *wait)
{
	update_FD_kernel.call_kernel(NULL, 1, wait);

	int offset = (p.nX - 1)*p.Channel_Height* 5;
}

void clVariablesFL::setKernelArgs()
{

	//cl_int ind = 0;
	//update_FL_kernel[0].set_argument(ind++, sizeof(cl_mem), vtr.BL_dep.get_buf_add());
	//update_FL_kernel[0].set_argument(ind++, sizeof(cl_mem), BLdep_tot.get_buf_add());
	//update_FL_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.C.get_buf_add());
	//update_FL_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.C0.get_buf_add());
	//update_FL_kernel[0].set_argument(ind++, sizeof(cl_mem), FI.get_buf_add());
	//update_FL_kernel[0].set_argument(ind++, sizeof(cl_uint), (void *)&Num_active_nodes);

	//ind = 0;
	//update_FL_kernel[1].set_argument(ind++, sizeof(cl_mem), RI.get_buf_add());
	//update_FL_kernel[1].set_argument(ind++, sizeof(cl_mem), IO_ind_dist.get_buf_add());
	//update_FL_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.C.get_buf_add());

	//cl_uint2 SSbot = { { 0, Num_active_nodes / 2 } };
	//cl_uint2 SStop = { { Num_active_nodes / 2, Num_active_nodes } };
	//
	//ind = 0;
	//update_FL_kernel[2].set_argument(ind++, sizeof(cl_mem), BLdep_tot.get_buf_add());
	//update_FL_kernel[2].set_argument(ind++, sizeof(cl_mem), BLdep_tot_temp.get_buf_add());
	//update_FL_kernel[2].set_argument(ind++, sizeof(cl_mem), FI.get_buf_add());
	//update_FL_kernel[2].set_argument(ind++, sizeof(cl_uint), (void *)&Num_active_nodes);
	//update_FL_kernel[2].set_argument(ind++, sizeof(cl_int2), (void *)&SSbot);
	//update_FL_kernel[2].set_argument(ind++, sizeof(cl_int2), (void *)&SStop);

	//ind = 0;
	//update_FL_kernel[3].set_argument(ind++, sizeof(cl_mem), BLdep_tot_temp.get_buf_add());
	//update_FL_kernel[3].set_argument(ind++, sizeof(cl_mem), vls.C.get_buf_add());
	//update_FL_kernel[3].set_argument(ind++, sizeof(cl_mem), vls.C0.get_buf_add());
	//update_FL_kernel[3].set_argument(ind++, sizeof(cl_mem), FI.get_buf_add());
	//update_FL_kernel[3].set_argument(ind++, sizeof(cl_uint), (void *)&Num_active_nodes);


	//cl_int2 Dmax = { { p.nX - 1, p.nY - 1 } };
	//int nCnodes = vls.nN - 1;
	//int SkipNode = nCnodes / 2 - 1;
	//ind = 0;
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.M.get_buf_add());
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.M_o.get_buf_add());
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Stor.get_buf_add());
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.C.get_buf_add());
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.C0.get_buf_add());
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_int2), (void *)&Dmax);
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_int), (void *)&nCnodes);
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_int), (void *)&SkipNode);
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.J_array.get_buf_add());
	//update_LS_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Ro_array.get_buf_add());

	//ind = 0;
	//update_LS_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.M.get_buf_add());
	//update_LS_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.M_o.get_buf_add());
	//update_LS_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.s.get_buf_add());
	//update_LS_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.sf.get_buf_add());

	//ind = 0;
	//update_LS_kernel[2].set_argument(ind++, sizeof(cl_mem), vls.M.get_buf_add());
	//update_LS_kernel[2].set_argument(ind++, sizeof(cl_mem), vls.sf.get_buf_add());
	//update_LS_kernel[2].set_argument(ind++, sizeof(cl_mem), Sum_M_temp.get_buf_add());
	//update_LS_kernel[2].set_argument(ind++, 2 * WORKGROUPSIZE_RED * sizeof(double), NULL);

	//ind = 0;
	//update_LS_kernel[3].set_argument(ind++, sizeof(cl_mem), Sum_M_temp.get_buf_add());
	//update_LS_kernel[3].set_argument(ind++, sizeof(cl_mem), vls.Masses.get_buf_add());
	//update_LS_kernel[3].set_argument(ind++, sizeof(cl_int), (void *)&vls.numBlocksM);


	//ind = 0;
	//update_LS_kernel[4].set_argument(ind++, sizeof(cl_mem), vls.M.get_buf_add());
	//update_LS_kernel[4].set_argument(ind++, sizeof(cl_mem), vls.M_red.get_buf_add());
	//update_LS_kernel[4].set_argument(ind++, sizeof(cl_mem), vlb.Act.get_buf_add());

	//ind = 0;
	//update_LB_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.M.get_buf_add());
	//update_LB_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Prop_loc.get_buf_add());
	//update_LB_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Act.get_buf_add());
	//update_LB_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Stor.get_buf_add());

	//ind = 0;
	//update_LB_kernel[1].set_argument(ind++, sizeof(cl_mem), vlb.IBB_loc.get_buf_add());
	//update_LB_kernel[1].set_argument(ind++, sizeof(cl_mem), vlb.IBB_coeff.get_buf_add());
	//update_LB_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.ii0_array.get_buf_add());
	//update_LB_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.iis1_array.get_buf_add());
	//update_LB_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.dir_array.get_buf_add());
	//update_LB_kernel[1].set_argument(ind++, sizeof(cl_mem), vlb.Stor.get_buf_add());
	//update_LB_kernel[1].set_argument(ind++, sizeof(cl_mem), vls.D_array.get_buf_add());


	//ind = 0;
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vtr.BL.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vls.C.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vls.M.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vls.s.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vls.sf.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vls.M_o.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vfd.dX.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vfd.dX_cur.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vlb.Stor.get_buf_add());
	//update_FD_kernel.set_argument(ind++, sizeof(cl_int2), (void *)&Dmax);
	//update_FD_kernel.set_argument(ind++, sizeof(cl_int), (void *)&vls.nBL);
	//update_FD_kernel.set_argument(ind++, sizeof(cl_mem), vls.C0.get_buf_add());

	//ind = 0;
	//int Xnum = p.TrDomainSize.x, Ynum = p.TrDomainSize.y;

	//update_TR_kernel[0].set_argument(ind++, sizeof(cl_mem), vls.M.get_buf_add());
	//update_TR_kernel[0].set_argument(ind++, sizeof(cl_mem), vtr.NodC.get_buf_add());
	//update_TR_kernel[0].set_argument(ind++, sizeof(cl_mem), vfd.dX_cur.get_buf_add());
	//update_TR_kernel[0].set_argument(ind++, sizeof(cl_mem), vfd.dX.get_buf_add());
	//update_TR_kernel[0].set_argument(ind++, sizeof(cl_mem), vlb.Stor.get_buf_add());
	//update_TR_kernel[0].set_argument(ind++, sizeof(cl_mem), vtr.Active_indicies.get_buf_add());
	//update_TR_kernel[0].set_argument(ind++, sizeof(cl_uint), (void *)&vtr.nActiveNodes);
	//

	//int xstart = (int)floor(X_MIN_VAL);
	//int xstop = (int)ceil(X_STOP_POS) - 1;
	//ind = 0;
	//update_TR_kernel[1].set_argument(ind++, sizeof(cl_mem), vtr.BL.get_buf_add());
	//update_TR_kernel[1].set_argument(ind++, sizeof(cl_mem), vtr.NodI.get_buf_add());
	//update_TR_kernel[1].set_argument(ind++, sizeof(cl_mem), vtr.BLind_ind.get_buf_add());
	//update_TR_kernel[1].set_argument(ind++, sizeof(cl_mem), vtr.Num_W_nodes.get_buf_add());
	//update_TR_kernel[1].set_argument(ind++, sizeof(BLbound), (void *)&vtr.Bounds);
	//update_TR_kernel[1].set_argument(ind++, sizeof(cl_uint), (void *)&xstart);
	//update_TR_kernel[1].set_argument(ind++, sizeof(cl_uint), (void *)&xstop);
	//update_TR_kernel[1].set_argument(ind++, sizeof(cl_mem), vtr.Active_flag.get_buf_add());

	//double Ymin = vls.C0[0].y;
	//double Ymax = vls.C0[vls.nN - 1].y;
	//ind = 0;
	//update_TR_kernel[2].set_argument(ind++, sizeof(cl_mem), vtr.P.get_buf_add());
	//update_TR_kernel[2].set_argument(ind++, sizeof(cl_mem), vtr.NodI.get_buf_add());
	//update_TR_kernel[2].set_argument(ind++, sizeof(cl_mem), vtr.BL.get_buf_add());
	//update_TR_kernel[2].set_argument(ind++, sizeof(double), (void*)&vtr.trP.X_release);
	//update_TR_kernel[2].set_argument(ind++, sizeof(double), (void*)&Ymin);
	//update_TR_kernel[2].set_argument(ind++, sizeof(double), (void*)&Ymax);

	//ind = 0;
	//update_TR_kernel[3].set_argument(ind++, sizeof(cl_mem), vtr.P.get_buf_add());
	//update_TR_kernel[3].set_argument(ind++, sizeof(cl_mem), vtr.BL.get_buf_add());


	//ind = 0;
	//update_TR_kernel[4].set_argument(ind++, sizeof(cl_mem), vtr.NodI.get_buf_add());
	//update_TR_kernel[4].set_argument(ind++, sizeof(cl_mem), vtr.Winds.get_buf_add());
	//update_TR_kernel[4].set_argument(ind++, sizeof(cl_mem), vtr.Active_indicies.get_buf_add());
	//update_TR_kernel[4].set_argument(ind++, sizeof(cl_uint), (void *)&vtr.nActiveNodes);

	//if (p.useOpenGL)
	//{
	//	ind = 0;
	//	int num_nodes = vls.nN / 2;
	//	vfl.update_GL_kernel.set_argument(ind++, sizeof(cl_mem), vls.C.get_buf_add());
	//	vfl.update_GL_kernel.set_argument(ind++, sizeof(cl_mem), vls.LSb_vbo.get_buf_add());
	//	vfl.update_GL_kernel.set_argument(ind++, sizeof(cl_mem), vls.LSt_vbo.get_buf_add());
	//	vfl.update_GL_kernel.set_argument(ind++, sizeof(int), (void*)& num_nodes);
	//}
}

void clVariablesFL::Release_Objects()
{// Should be handled by destructors
	//update_FL_kernel[0].free_memory();
	//update_FL_kernel[1].free_memory();

	//update_LS_kernel[0].free_memory();
	//update_LS_kernel[1].free_memory();
	//update_LS_kernel[2].free_memory();
	//update_LS_kernel[3].free_memory();
	//update_LS_kernel[4].free_memory();
	//update_LB_kernel[0].free_memory();
	//update_LB_kernel[1].free_memory();

	//update_FD_kernel.free_memory();

	//update_TR_kernel[0].free_memory();
	//update_TR_kernel[1].free_memory();
	//update_TR_kernel[2].free_memory();
	//update_TR_kernel[3].free_memory();
	//update_TR_kernel[4].free_memory();


	//update_GL_kernel.free_memory();

	//FI.delete_array(IOQUEUE);
	//RI.delete_array(IOQUEUE);
	//BLdep_tot.delete_array(IOQUEUE);
	//BLdep_tot_temp.delete_array(IOQUEUE);
	//IO_ind_dist.delete_array(IOQUEUE);
	//Sum_M_temp.delete_array(IOQUEUE);

}

void clVariablesFL::allocateBuffers()
{
	FI.allocate_buffer_w_copy();
	RI.allocate_buffer_w_copy();
	IO_ind_dist.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	BLdep_tot.allocate_buffer_w_copy();
	BLdep_tot_temp.allocate_buffer_w_copy();

	//FIgure out whats going on here, and how to replace vls.numBlocksM
	//Sum_M_temp.zeros(2 * vls.numBlocksM);
	//Sum_M_temp.allocate_buffer_w_copy();
}

void clVariablesFL::save2file()
{
	BLdep_tot.save_txt_from_device("BLdep_tot");
}

void clVariablesFL::saveRestartFiles()
{
	FI.save_from_device("fli");
	BLdep_tot.save_bin_from_device("BLdep_tot");
}


void clVariablesFL::UpdateRestart()
{
//
//	cl_event LS_Evt, Fill_Evt;
//
//	//vfd.dX_cur.enqueue_copy_to_buffer(FDQUEUE, vfd.dX.get_buffer());
//
//	update_LS(&LS_Evt);
//
//	vls.C.read_from_buffer(IOQUEUE, 1, &LS_Evt);
//	vls.M.read_from_buffer(IOQUEUE);
//
//	update_LB();
//
//#ifdef USE_OPENGL
//	update_GL();
//#endif
//
//	update_FD(&LS_Evt);
//
//	clFlush(FDQUEUE);
//
//	nodeI Ntemp;
//	for (int i = 0; i < MAX_BL_PER_NODE; i++)
//		Ntemp.BLind[i] = -1;
//	Ntemp.Wall_Flag = 0;
//	vtr.NodI.FillBuffer(IOQUEUE, Ntemp);
//	vtr.BLind_ind.FillBuffer(IOQUEUE, 0);
//
//	int Num_Wnodes_temp = vtr.Num_wall_nodes_max;
//	vtr.Num_W_nodes.FillBuffer(IOQUEUE, 0, &Fill_Evt);
//
//	clFlush(IOQUEUE);
//
//	update_TR(&Fill_Evt, Num_Wnodes_temp);
//
//	p.flushQueues();
//
//	clReleaseEvent(LS_Evt);
//	clReleaseEvent(Fill_Evt);
//

}

