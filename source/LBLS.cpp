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

// TODO: Make sure that when anything is updated on GPU based on vtr.BL, it will not be an
//		issue if vls.BL is used to initialize data (since top row indicies increase from left to right
//		for vtr.BL, and decreases for vls.BL 


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

#ifdef _DEBUG
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "GPU_DEBUG");
#endif

	vls.ini();
	vlb.ini();
	vfd.ini();
	vtr.ini();
	vfl.ini();

	if (vfd.calcNuFlag)
	{
		vfd.iniNuCoeffs();
	}

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
	if (!(vlb.fLoadedFlag || vlb.vLoadedFlag) || (vlb.kOmegaClass.kOmegaSolverFlag &&
		!vlb.kOmegaClass.kOmegaLoadedFlag))
	{
		LOGMESSAGE("Solving for initial flow field");
		vlb.getInitialFlow();
		LOGMESSAGE("Initial flow field obtained");
	}

	// if using temperature solver, and no initial temperature is
	// loaded (or temp is loaded from text file), must get initial temperature
	if (vfd.thermalSolverFlag && (!vfd.tempLoadedFlag))// || vlb.vLoadedFlag))
	{
		LOGMESSAGE("Solving for initial temperature distribution");
		vlb.getIntialFlowAndTemp();
		LOGMESSAGE("Solved for initial temperature distribution");
	}
	  
	vtr.callRand.call_kernel();
	clFinish(TRQUEUE);
	if (!vtr.restartRunFlag)
	{
		vtr.parSort.initialSort();
	}

	clEnv::instance()->finishQueues();

	p.freeHostArrays();

	if (vtr.restartRunFlag)
	{
		if (vtr.parSort.Ploc(1).y > 0)
			vtr.wallShear.updateParRemArgs();
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

	//vfl.testFlUpdate();


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

void Solve_no_TR()
{
	cl_event LB_Evt;
#ifndef IN_KERNEL_IBB
	vlb.collisionKernel.call_kernel(LBQUEUE_REF);
	vlb.ibbKernel.call_kernel(LBQUEUE_REF, 0, nullptr, &LB_Evt);
#else
	vlb.collisionKernel.call_kernel(LBQUEUE_REF, 0, nullptr, &LB_Evt);
#endif
	if (vlb.kOmegaClass.kOmegaSolverFlag)
		vlb.kOmegaClass.Solve();

	if (p.TimeN % p.tfdSteps == 0)
	{
		vfd.Solve(nullptr, 1, &LB_Evt);
	}

	clEnv::instance()->finishQueues();

	clReleaseEvent(LB_Evt);
}

void Solve_all()
{
	cl_event LB_Evt, TR_Par_Evt;

	vtr.updateFlag.FillBuffer(0, TRQUEUE_REF);
	clFlush(TRQUEUE);
	vtr.trUpdateParKernel.call_kernel();

	vtr.updateWallParticles();
	clFinish(TRQUEUE);

	if (p.useOpenGL)
		vtr.glParticles.TR_GL_kernel.call_kernel(nullptr, 0, nullptr, &TR_Par_Evt);

#ifndef IN_KERNEL_IBB
	vlb.collisionKernel.call_kernel(LBQUEUE_REF);
	vlb.ibbKernel.call_kernel(LBQUEUE_REF, 0, nullptr, &LB_Evt);
#else
	vlb.collisionKernel.call_kernel(LBQUEUE_REF, 0, nullptr, &LB_Evt);
#endif
	if (vlb.kOmegaClass.kOmegaSolverFlag)
		vlb.kOmegaClass.Solve();

	if (p.TimeN % p.tfdSteps == 0)
	{
		vfd.Solve(nullptr, 1, &LB_Evt);
	}

	vlb.alter ^= 1;

	vtr.wallShear.nodeShearKernel.call_kernel(
		static_cast<DualKernel::kernelID>(vlb.alter),
		TRQUEUE_REF, 1, &LB_Evt, nullptr);

	vtr.wallShear.wallShearKernel.call_kernel(TRQUEUE_REF);

	if (vtr.parSort.Ploc(1).y > 0) //particle removal
	{
		vtr.wallShear.trShearRemovalKernel.call_kernel(TRQUEUE_REF);
	}

	if (vtr.parSort.sortTimer <= 0)
	{
		
		p.sort_called = 1;
		if(p.useOpenGL)
			vtr.parSort.sortParticles(&TR_Par_Evt, 1);
		else
			vtr.parSort.sortParticles();
	}
	else
	{
		p.sort_called = 0;
		vtr.parSort.sortTimer -= p.trSteps;
	}

	clEnv::instance()->flushQueues();
		


	if (vfl.FL_timer <= 0)
	{
		if (p.sort_called == 0)
		{
			if (p.useOpenGL)
				vtr.parSort.sortParticles(&TR_Par_Evt, 1);
			else
				vtr.parSort.sortParticles();
		}
		clEnv::instance()->finishQueues();
		vfl.update();
		vfl.FL_timer = vfl.flTimePerUpdate;
	}
	else
	{
		clEnv::instance()->finishQueues();
		vfl.FL_timer -= p.trSteps;
	}


	if (p.useOpenGL)
	{
		FINISH_QUEUES;
		clEnv::instance()->renderDomain();
		clReleaseEvent(TR_Par_Evt);
	}

	clReleaseEvent(LB_Evt);
}

void Solve_LB_TRwalls()
{
	cl_event LB_Evt;

	vtr.updateFlag.FillBuffer(0, TRQUEUE_REF);
	clFlush(TRQUEUE);
	vtr.updateWallParticles();
	clFinish(TRQUEUE);

#ifndef IN_KERNEL_IBB
	vlb.collisionKernel.call_kernel(LBQUEUE_REF);
	vlb.ibbKernel.call_kernel(LBQUEUE_REF, 0, nullptr, &LB_Evt);
#else
	vlb.collisionKernel.call_kernel(LBQUEUE_REF, 0, nullptr, &LB_Evt);
#endif
	if (vlb.kOmegaClass.kOmegaSolverFlag)
		vlb.kOmegaClass.Solve();

	if (p.TimeN % p.tfdSteps == 0)
	{
		vfd.Solve(nullptr, 1, &LB_Evt);
	}
	
	vlb.alter ^= 1;

	vtr.wallShear.nodeShearKernel.call_kernel(
		static_cast<DualKernel::kernelID>(vlb.alter),
		TRQUEUE_REF, 1, &LB_Evt, nullptr);

	vtr.wallShear.wallShearKernel.call_kernel(TRQUEUE_REF);

	if (vtr.parSort.Ploc(1).y > 0) //particle removal
	{
		vtr.wallShear.trShearRemovalKernel.call_kernel(TRQUEUE_REF);
	}

	clEnv::instance()->finishQueues();

	clReleaseEvent(LB_Evt);
}


void run()
{
	if (p.useOpenGL)
	{
		vtr.glParticles.TR_GL_kernel.call_kernel();
		FINISH_QUEUES;
		clEnv::instance()->renderDomain();
	}

	int option = (OPTION_SAVE_MACRO_FIELDS);
	vlb.collisionKernel.setOption(&option);


	while (1)
	{
		// Every trStep timestep, solve all (temp only solved if 
		// p.TimeN % p.tfdSteps == 0
		if(p.TimeN % p.trSteps == 0)
		{
			Solve_all();
		}
		// If not a full trStep, but is a trSteps_wall,
		// solve lb, tracers near wall, and if necessary, temp
		else if (p.TimeN % p.trSteps_wall == 0)
		{
			Solve_LB_TRwalls();
		}
		// if not a tr or tr_wall step, solve lb and if necessary
		// temp
		else
		{
			Solve_no_TR();
		}

		step();

		if (testFinish() == true)
		{
			clEnv::instance()->finishQueues();
			break;
		}
	}
}

void finish(void)
{
	d.finish();
	p.finish();
}




