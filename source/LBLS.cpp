// LBLS.cpp : Defines the entry point for the console application.
// (c) Alexander Alexeev, 2006 
//


// TODO: split the saving of parameters into static and dynamic variables,
//		since there is no reason to keep writing the save variables out as
//		it wastes cycles.
// TODO: Make sure that all variables are same between host and device.
//		i.e. make sure int is int, uint is uint, etc in kernels
// TODO: remove old files from source folder that are no longer needed.
//			specifically clOutput, FileReader/Writer and Structures
// TODO: consolidate error printing functions
// TODO: get rid of bool typedef and replace with regular bools
// TODO: Make sure all classes include setSourceDefines function,
//			and that they are called before generating source code
// TODO: Implement checks to ensure that places where short is used instead of 
//		int do not require values greater than 32768 for short and 65535 for ushort
// TODO: Make sure that all uses of nType in kernels use NTYPE_TYPE rather than int
//		or short. This is crucial because it will not throw any build errors,
//		but will screw everything up if NTYPE_TYPE isnt used.
// TODO: make sure that since Ro remains 1 throughout simulation, that this
//		wont create any issues by assuming that Ro = 0 in the fouling layer.
// TODO: set something up to update vls.lsMap every so often. Doesnt need to be
//		frequent since the nodes wont move much during the simulation (maybe every
//		clump step)
// TODO: Go ahead and implement own types for variables utilizing short and ushort
//		in order to reduce cost if domains require values larger than what they can represent.
//		specifically anything associated with BL and vls.C as this is where it most often was
//		used
// TODO: Make sure that dynamic array updates are coupled with resetting of kernel global sizes

#include "StdAfx.h"
#include "oclEnvironment.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"
#include "clProblem.h"
#include "clDisplay.h"
#include "Reducer.h"

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
	staticBaseVar::setArrayContext(CLCONTEXT_REF);

	p.ini();

	omp_set_num_threads(OPENMP_NUMBER_THREADS);
	LOGMESSAGE("Setting number of OpenMP threads to " + std::to_string(OPENMP_NUMBER_THREADS));

	vls.ini();
	vlb.ini();
	vfd.ini();
	vtr.ini();
	vfl.ini();

	if (vfd.calcNuFlag)
		vfd.iniNuCoeffs();

	sourceGenerator::SourceInstance()->buildSource();
	LOGMESSAGE("OpenCL source built");
	ReduceGenerator::ReduceInstance()->buildSource();
	LOGMESSAGE("Reducer source built");

	// Meta data for all CSR matricies have been generated,
	// clSparse objects can be deleted
	BiCGStabGenerator::BiCGStabInstance()->deleteClSparseObjects(); 

	// If f distributions are not loaded, and kappa and omega are not
	// loaded when using turbulence model, the initial flow must be
	// solved for.
	//if (!vlb.fLoadedFlag && (vlb.kOmegaClass.kOmegaSolverFlag &&
	//	!vlb.kOmegaClass.kOmegaLoadedFlag))
	//{
	//	LOGMESSAGE("Solving for initial flow field");
	//	vlb.getInitialFlow();
	//	LOGMESSAGE("Initial flow field obtained");
	//}

	// if using temperature solver, and no initial temperature is
	// loaded, must get initial temperature
	if (vfd.thermalSolverFlag && !vfd.tempLoadedFlag)
	{
		LOGMESSAGE("Solving for initial temperature distribution");
		vlb.getIntialFlowAndTemp();
		LOGMESSAGE("Solved for initial temperature distribution");
	}
	   
	if (!vtr.restartRunFlag)
		vtr.parSort.initialSort();

	clEnv::instance()->finishQueues();


	// Shouldnt be necessary since CPU mem >> GPU mem,
	// but may be useful when utilizing 2 GPUs on same computer
	vls.freeHostArrays();
	vlb.freeHostArrays();
	vfd.freeHostArrays();
	vtr.freeHostArrays();
	vfl.freeHostArrays();

	//vtr.TR_Node_kernel[0].call_kernel();
	//vtr.TR_Wall_Node_kernel[0].call_kernel();


	if (vtr.restartRunFlag)
	{
		//if (vtr.parSort.Ploc(1).y > 0)
		//	vtr.Update_Par_Rem_Args();
	}

	//// Not sure if need to do this after re-factoring of code.
	//if ((vfl.restartRunFlag && vlb.restartRunFlag && vtr.restartRunFlag && vfd.restartRunFlag))
	//{
	//	vfl.UpdateRestart();
	//}

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

bool testFinish()
{
#ifdef PROFILING
	if (vlb.count_profile == 1001)
		return true;
#endif

	if (p.Time >= p.StopTime)
		return true;

	return false;
};

void Solve_LB_only()
{
	//cl_event Col_Fluid;
	//vlb.Collision_kernel.call_kernel(LBQUEUE_REF, 0, NULL, &Col_Fluid);
	clFlush(LBQUEUE);
	
	//need to implement class for calculating mass losses to make necessary adjustments
	//vlb.Reduce_Ro.reduce();//LBQUEUE);

	vlb.alter ^= 1;
	clFlush(LBQUEUE);
	clFlush(FDQUEUE);
//	clReleaseEvent(Col_Fluid);
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

	//	if (testFinish() == true)
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




