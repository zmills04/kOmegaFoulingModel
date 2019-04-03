// LBLS.cpp : Defines the entry point for the console application.
// (c) Alexander Alexeev, 2006 
//



// TODO: remove old files from source folder that are no longer needed.
//			specifically clOutput, FileReader/Writer and Structures
// TODO: consolidate error printing functions
// TODO: get rid of BOOL typedef and replace with regular bools

#include "StdAfx.h"
#include "oclEnvironment.h"
#include "Reducer.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"
#include "clProblem.h"
#include "clDisplay.h"


clDisplay d;
clProblem p;
clVariablesLS vls;
clVariablesLB vlb;
clVariablesFD vfd;
clVariablesTR vtr;
clVariablesFL vfl;

///////////////////////////////////////////////////////////////////////////////
///////////////// Make sure to test for output before increasing Time.

void ini();
void run();
void finish(void);
void Solve_LB_only();


int main(int argc, char* argv[])
{

	ini();
	run();
	finish();
	return 0;
}

void ini()
{
	omp_set_num_threads(OPENMP_NUMBER_THREADS);

	p.ini();
	vls.ini();
	//vls.ini_IBB_arrays();
	vlb.ini();
	vfd.ini();
	vtr.ini();
	vfl.ini();


	//vlb.ini_IBB();

	if(vfd.calcNuFlag)
		vfd.initialize_Nu_kernel();

#ifdef USE_OPENGL
	vls.iniGL();
	vtr.iniGL();
#endif

	sourceGenerator::SourceInstance()->buildSource();
	ReduceGenerator::ReduceInstance()->buildSource();
	
	// Meta data for all CSR matricies have been generated,
	// clSparse objects can be deleted
	BiCGStabGenerator::BiCGStabInstance()->deleteClSparseObjects(); 

	if (!vlb.restartRunFlag)
	{
		vlb.getInitialFlow();
	}

	if (vfd.thermalSolverFlag && !vfd.tempLoadedFlag)
	{
		vlb.getIntialFlowAndTemp();
	}


	if (!vtr.restartRunFlag)
		vtr.initial_sort_particles();

	clEnv::instance()->finishQueues();


	// Shouldnt be necessary since CPU mem >> GPU mem,
	// but may be useful when utilizing 2 GPUs on same computer
	vls.freeHostMem();
	vlb.freeHostMem();
	vfd.freeHostMem();
	vtr.freeHostMem();
	vfl.freeHostMem();

	vtr.TR_Node_kernel[0].call_kernel();
	vtr.TR_Wall_Node_kernel[0].call_kernel();


	if (vtr.restartRunFlag)
	{
		if (vtr.Ploc(1).y > 0)
			vtr.Update_Par_Rem_Args();
	}


	if ((vfl.restartRunFlag && vlb.restartRunFlag && vtr.restartRunFlag && vfd.restartRunFlag))
	{
		vfl.UpdateRestart();
	}

	clEnv::instance()->finishQueues();
	//vlb.Calculate_U_mean();
	d.ini();
	d.start();
}

void step()
{
	p.step();
	d.step();
#ifdef PROFILING
	vlb.count_profile++;
#endif
};

BOOL testFinish()
{
#ifdef PROFILING
	if (vlb.count_profile == 1001)
		return TRUE;
#endif

	if (p.Time >= p.StopTime)
		return TRUE;

	return FALSE;
};

void Solve_LB_only()
{
	cl_event Col_Fluid;
	vlb.Collision_kernel.call_kernel(LBQUEUE_REF, 0, NULL, &Col_Fluid);
	clFlush(LBQUEUE);
	
	//need to implement class for calculating mass losses to make necessary adjustments
	//vlb.Reduce_Ro.reduce();//LBQUEUE);

	vlb.alter ^= 1;
	clFlush(LBQUEUE);
	clFlush(FDQUEUE);
	clReleaseEvent(Col_Fluid);
}


void Solve_all()
{
//	cl_event TR_Node_Evt, Update_Flag_Evt, LB_Evt, TR_Shear_Evt;
//
//
//	vtr.Update_flag.FillBuffer(p.IOqueue,0, * &Update_Flag_Evt);
//	clFlush(p.IOqueue);
//
//	vtr.TR_Node_kernel[1].call_kernel();
//	vtr.TR_Wall_Node_kernel[1].call_kernel(&TR_Node_Evt);
//
//
//	vtr.TR_Update_Par_kernel.call_kernel(1, &Update_Flag_Evt);
//
//	vtr.call_update_wall_particles();
//
//#ifdef USE_OPENGL
//	cl_event TR_Par_Evt;
//	vtr.TR_GL_kernel.call_kernel(&TR_Par_Evt);
//#endif
//
//	vlb.Collision_kernel.call_kernel(1, &TR_Node_Evt, &LB_Evt);	//updates vel
//	
//	vtr.Shear_kernels[vlb.alter].call_kernel();
//
//	vtr.Shear_kernels[2].call_kernel(&TR_Shear_Evt);
//
//	//vlb.Reduce_Ro.call_kernels(LBQUEUE);
//
//	vlb.alter ^= 1;
//
//	if (vtr.Ploc(1).y > 0) //particle removal
//	{
//		vtr.TR_Shear_Removal_kernel.call_kernel(1, &TR_Shear_Evt);
//	}
//
//	if (vtr.Sort_timer <= 0)
//	{
//#ifdef USE_OPENGL
//		vtr.sort_particles(&TR_Par_Evt);
//#else
//		vtr.sort_particles();
//#endif
//
//	}
//
//	p.flushQueues();
//
//	vtr.TR_Node_kernel[0].call_kernel();
//	vtr.TR_Wall_Node_kernel[0].call_kernel();
//
//	if (vtr.Sort_timer <= 0)
//	{
//		clFinish(p.TRqueue);
//
//		if (vtr.Ploc(1).y > 0)
//			vtr.Update_Par_Rem_Args();
//
//		if (vtr.Ploc(0).y > 0)
//			vtr.Re_Release_Par();
//
//		vtr.Sort_timer = NUM_STEPS_BTW_SORT;
//		p.sort_called = 1;
//	}
//	else
//	{
//		p.sort_called = 0;
//		vtr.Sort_timer -= c.trSteps;
//	}
//
//	if (vfl.FL_timer <= 0)
//	{
//		p.update_walls = 1;
//	}
//	else
//	{
//		vfl.FL_timer -= c.trSteps;
//	}
//
//	if (p.update_walls * p.sort_called)
//	{
//		p.finishQueues();
//		vfl.update();
//		p.update_walls = 0;
//		vfl.FL_timer = UPDATE_TIME;
//
//	}
//
//	p.finishQueues();
//
//#ifdef USE_OPENGL
//	d.Render_Domain();
//	clReleaseEvent(TR_Par_Evt);
//#endif
//
//	clReleaseEvent(TR_Shear_Evt);
//	clReleaseEvent(LB_Evt);
//	clReleaseEvent(TR_Node_Evt);
//	clReleaseEvent(Update_Flag_Evt);
}

void Solve_LB_TRwalls()
{
//	cl_event TR_Node_Evt, Update_Flag_Evt, TR_Shear_Evt, LB_Evt;
//
//	vtr.Update_flag.FillBuffer(p.IOqueue, 0, &Update_Flag_Evt);
//	clFlush(p.IOqueue);
//
//	vtr.TR_Wall_Node_kernel[1].call_kernel(&TR_Node_Evt);
//
//	vtr.call_update_wall_particles(&Update_Flag_Evt);
//
//	vlb.Collision_kernel.call_kernel(1, &TR_Node_Evt, &LB_Evt);	//updates vel
//	
//	vtr.Shear_kernels[vlb.alter].call_kernel();
//
//	vtr.Shear_kernels[2].call_kernel(&TR_Shear_Evt);
//
//	//vlb.Reduce_Ro.call_kernels(LBQUEUE);
//	
//	vlb.alter ^= 1;
//
//	if (vtr.Ploc(1).y > 0) //particle removal
//	{
//		vtr.TR_Shear_Removal_kernel.call_kernel(1, &TR_Shear_Evt);
//	}
//
//	p.finishQueues();
//
//	clReleaseEvent(TR_Shear_Evt);
//	clReleaseEvent(TR_Node_Evt);
//	clReleaseEvent(Update_Flag_Evt);
//	clReleaseEvent(LB_Evt);
}
//
//
void run()
{
	//while (1)
	//{
	//	if (c.TimeN % c.trSteps_wall)
	//	{
	//		Solve_LB_only();
	//	}
	//	else
	//	{

	//		if (c.TimeN % c.trSteps == 0)
	//		{
	//			Solve_all();
	//		}
	//		else
	//			Solve_LB_TRwalls();
	//	}

	//	step();

	//	if (testFinish() == TRUE)
	//	{
	//		p.finishQueues();
	//		break;
	//	}
	//}

}

void finish(void)
{
	d.finish();
	//vls.Release_Objects();
	//vlb.releaseObjects();
	//vfd.Release_Objects();
	//vfl.Release_Objects();
	//vtr.Release_Objects();
	p.finish();
}




