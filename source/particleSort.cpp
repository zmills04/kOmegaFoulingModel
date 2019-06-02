#include "particleSort.h"
#include "clVariablesTR.h"


void particleSort::allocateArrays()
{
	Ploc.allocate(vtr.TrDomainSize.x*vtr.TrDomainSize.y + 2);
	Ptemp.setSizes(vtr.nN, vtr.nN, 1);
}

void particleSort::allocateBuffers()
{
	Ploc.allocate_buffer_w_copy();
	int status;
	
	sortLocs1 = clCreateBuffer(CLCONTEXT, CL_MEM_READ_WRITE,
		sizeof(int) * vtr.nN, NULL, &status);
	ERROR_CHECKING_OCL(status, "Error creating sortLocs1 "\
		"in particlesSort class");
	sortLocs1 = clCreateBuffer(CLCONTEXT, CL_MEM_READ_WRITE,
		sizeof(int) * vtr.nN, NULL, &status);
	ERROR_CHECKING_OCL(status, "Error creating sortLocs2 "\
		"in particlesSort class");
}

bool myfunction(cl_int2 i, cl_int2 j) { return (i.x < j.x); }

void particleSort::clumpParticles()
{
	sortParticlesForClumping();
	/////////////Fill vector with values corresponing to particle locations/////////

	// domain split into cells with size 1/size_multiplier and
	// all particles within a given cell are lumped together.
	double size_multiplier = 15.;
	int tr_size_y = (p.nY - 1) * (int)size_multiplier;
	vtr.P.copyToHost();


	// stores [location in domain, location in Par array] for each particle
	std::vector<cl_int2> ParLocVec(vtr.nN, { { 0, 0 } });
	int j;
#pragma omp parallel for 
	for (j = 0; j < vtr.nN; j++)
	{
		ParLocVec[j].y = j;
		cl_double2 postemp = vtr.P(j).pos;
		postemp = { { postemp.x * size_multiplier, postemp.y * size_multiplier } };
		cl_int2 Posi = { { (int)floor(postemp.x), (int)floor(postemp.y) } };

		int pcur_loc = Posi.x*tr_size_y + Posi.y;

		if (vtr.P.Dep_Flag(j) == -2)
		{
			ParLocVec[j].x = -2;
		}
		else if (vtr.P.Dep_Flag(j) > -1)
		{
			ParLocVec[j].x = -1;
		}
		else
		{
			ParLocVec[j].x = pcur_loc;
		}
	}
	
	//Sort vector
	std::sort(ParLocVec.begin(), ParLocVec.end(), myfunction);


	///Fill beginning of temporary particle array with particles about to be re-released
	///and particles currently deposited on wall
	int cur_ind = 0;


	while (true)
	{
		Ptemp.setStruct(vtr.P.getStruct(ParLocVec[cur_ind].y), cur_ind);
		cur_ind++;
		if (ParLocVec[cur_ind].x > -1)
			break;
	}

	////Combine same type particles located in each sub-box
	int cur_red_ind = cur_ind;
	int num_dep_particles = cur_ind;
	while (true)
	{
		int start_ind = cur_ind;
		int cur_loc = ParLocVec[cur_ind].x;
		int stop_ind = cur_ind;

		//find index of last particle located in same box.
		while (true)
		{
			if (stop_ind == vtr.nN - 1)
			{
				break;
			}
			if (ParLocVec[stop_ind + 1].x != cur_loc)
			{
				break;
			}
			stop_ind++;
		}

		for (int i = 0; i < vtr.parP.Nd; i++)
		{
			int count_pars = 0;
			legacyPar Pcur;
			int intTimer = 0;

			Pcur.Dep_Flag = -1;
			Pcur.Dep_timer = TIME_BEFORE_STICK;
			Pcur.loc = vtr.P.loc(ParLocVec[start_ind].y);
			Pcur.Num_rep = 0;
			Pcur.pos = { { 0., 0. } };
			Pcur.timer = 0;
			Pcur.type = i;
			for (int j = start_ind; j <= stop_ind; j++)
			{
				legacyPar Pcurtemp = vtr.P.getStruct(ParLocVec[j].y);
				if (Pcurtemp.type == i)
				{
					count_pars++;
					Pcur.addPar(Pcurtemp);
					intTimer += Pcurtemp.timer;
				}
			}
			if (count_pars)
			{
				//Pcur.Num_rep /= (double)count_pars;
				Pcur.pos = Divide2(Pcur.pos, (double)count_pars);
				Pcur.timer = (cl_ushort)(intTimer/count_pars);
				cl_int2 Posi = { { (int)floor(Pcur.pos.x), (int)floor(Pcur.pos.y) } };
				Pcur.loc = Posi.x + vtr.TrDomainSize.x * Posi.y;
				Ptemp.setStruct(Pcur, cur_red_ind++);
			}
		}

		cur_ind = stop_ind + 1;
		if (cur_ind == vtr.nN)
			break;
	}



	// Adding removed particles as deposited, but with Release
	// timer (using timer array) set to release approx
	// avgParPerRelease particles per release step.
	int numReRelease = (vtr.nN - cur_red_ind);
	int numRelSteps = MAX(numReRelease / avgParPerRelease,1);
	
	for (int i = cur_red_ind; i < vtr.nN; i++)
	{
		legacyPar Pcur;
		Pcur.pos = { {0., 0.} };
		Pcur.type = 0;
		Pcur.Dep_Flag = -2;
		Pcur.Dep_timer = 0;
		Pcur.timer = (int)floor((double)(i - cur_red_ind) / (double)numRelSteps);
		Pcur.loc = -2;
		Ptemp.setStruct(Pcur, i);
	}

	// Writing Ptemp arrays on host to P device buffers
	vtr.P.writeParToBuffer(Ptemp, cur_red_ind - 1);
	
	// Update indexing arguments for clump kernel
	clumpParticleKernels.set_argument(5, &cur_red_ind);
	clumpParticleKernels.set_argument(6, &num_dep_particles);

	clumpParticleKernels.set_global_call_kernel(cur_red_ind - num_dep_particles);
	clFinish(TRQUEUE);

	initialSort();
}

void particleSort::createKernels()
{
	trReReleaseKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_release_par");
	trReReleaseKernel.set_size(2 * WORKGROUPSIZE_RERELEASE, WORKGROUPSIZE_RERELEASE);
	//Global size is set before enqueing kernel (calling set size w/ dummy arg to ensure dim is set)

	sortKernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "Sort_merge_local");
	sortKernel[0].set_size(vtr.nN, WORKGROUPSIZE_SORT);

	sortKernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "Sort_merge_global");
	sortKernel[1].set_size(vtr.nN, WORKGROUPSIZE_SORT);

	sortKernel[2].create_kernel(GetSourceProgram, TRQUEUE_REF, "Sort_merge_global");
	sortKernel[2].set_size(vtr.nN, WORKGROUPSIZE_SORT);

	sortKernel[3].create_kernel(GetSourceProgram, TRQUEUE_REF, "Sort_update_loc");
	sortKernel[3].set_size(vtr.nN, WORKGROUPSIZE_PLOC);

	getUmaxKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "find_umax");

	// Getting work size. For OpenCL version 1.2 the implementation
	// is basically the same as the final reduce step for generic reduce
	// kernels. For OpenCL 2.0, implementation uses built-in max reduce
	// function. Current implementation only allows for 1 work group,
	// so size of channel height must be <= max WG size for opencl 2.0 or
	// twice the max WG size for opencl 1.2 (throws error if not).
	// (currently, max WG size hardcoded as 256).

	int globalLocalSize = 2;
	while (globalLocalSize < p.Channel_Height)
		globalLocalSize *= 2;

#ifdef OPENCL_VERSION_1_2
	globalLocalSize /= 2;
#endif

	ERROR_CHECKING(globalLocalSize > 256, "Must update implementation "\
		"before channel height > 256 can be used", ERROR_INITIALIZING_VTR);
	getUmaxKernel.set_size(globalLocalSize, globalLocalSize);
	   	 
	clumpParticleKernels.create_kernel(GetSourceProgram, TRQUEUE_REF, "Test_Bounds_Clumped_Particles");
	clumpParticleKernels.set_size(2 * WORKGROUPSIZE_TR, WORKGROUPSIZE_TR);
}

void particleSort::setKernelArgs()
{
	int	ind = 0;
	Par::arrName arrList_[] = { Par::posArr, Par::depFlagArr, Par::locArr };
	vtr.P.setBuffers(clumpParticleKernels, ind, arrList_, 3);
	clumpParticleKernels.set_argument(ind++, trP.fVals.get_buf_add());
	clumpParticleKernels.set_argument(ind++, vls.C.get_buf_add());
	// Remaining arguments set before call

	ind = 0;
	getUmaxKernel.set_argument(ind++, vlb.Ux_array.get_buf_add());
	trP.setBuffers(getUmaxKernel, ind);
#ifdef OPENCL_VERSION_1_2
	getUmaxKernel.set_local_memory(ind, getUmaxKernel.getLocalSize() * sizeof(double));
#endif

	Par::arrName arrList2_[] = { Par::posArr, Par::numRepArr, Par::typeArr,
		Par::depFlagArr, Par::timerArr, Par::locArr };
	ind = 0;
	vtr.P.setBuffers(trReReleaseKernel, ind, arrList2_, 6);
	trReReleaseKernel.set_argument(ind++, vtr.RandList.get_buf_add());
	trP.setBuffers(trReReleaseKernel, ind);
	trReReleaseKernel.set_argument(ind++, vlb.Ux_array.get_buf_add());
	trReReleaseKernel.set_argument(ind++, vtr.parP.D_dist.get_buf_add());
	// These are set before the kernel is called
	//trReReleaseKernel.set_argument(10, maxel);
	//trReReleaseKernel.set_argument(11, Concentration Number);


	ind = 0;
	sortKernel[0].set_argument(ind++, vtr.P.loc.get_buf_add());
	sortKernel[0].set_argument(ind++, &sortLocs1);
	sortKernel[0].set_local_memory(ind++, WORKGROUPSIZE_SORT * sizeof(int));
	sortKernel[0].set_local_memory(ind++, WORKGROUPSIZE_SORT * sizeof(int));
	sortKernel[0].set_local_memory(ind++, WORKGROUPSIZE_SORT * sizeof(int));
	sortKernel[0].set_local_memory(ind++, WORKGROUPSIZE_SORT * sizeof(int));

	ind = 0;
	sortKernel[1].set_argument(ind++, vtr.P.loc.get_buf_add());
	sortKernel[1].set_argument(ind++, &sortLocs1);
	sortKernel[1].set_argument(ind++, Ptemp.loc.get_buf_add());
	sortKernel[1].set_argument(ind++, &sortLocs2);
	sortKernel[1].setOptionInd(ind); // option is srcLogicalBlockSize

	ind = 0;
	sortKernel[2].set_argument(ind++, Ptemp.loc.get_buf_add());
	sortKernel[2].set_argument(ind++, &sortLocs2);
	sortKernel[2].set_argument(ind++, vtr.P.loc.get_buf_add());
	sortKernel[2].set_argument(ind++, &sortLocs1);
	sortKernel[2].setOptionInd(ind); // option is srcLogicalBlockSize


	// All P arrays will be copied to Ptemp buffers while sort kernels
	// running (on a different queue), so that there is some overlap in
	// operation, then the elements of Ptemp are sorted into their correct
	// location of the P arrays.
	Par::arrName arrList3_[] = { Par::posArr, Par::numRepArr, Par::typeArr,
		Par::depFlagArr, Par::depTimerArr, Par::timerArr };
	ind = 0;
	Ptemp.setBuffers(sortKernel[3], ind, arrList3_, 6);
	vtr.P.setBuffers(sortKernel[3], ind, arrList3_, 6);
	if (oddMergesFlag)
	{ // if odd number, then Ptemp has correct loc array info, so need
	  // to copy loc array info over from Ptemp in kernel
		sortKernel[3].set_argument(ind++, Ptemp.loc.get_buf_add());
		sortKernel[3].set_argument(ind++, vtr.P.loc.get_buf_add());
		sortKernel[3].set_argument(ind++, &sortLocs2);
	}
	else
	{
		sortKernel[3].set_argument(ind++, vtr.P.loc.get_buf_add());
		sortKernel[3].set_argument(ind++, &sortLocs1);
	}
	sortKernel[3].set_argument(ind, Ploc.get_buf_add());
}


void particleSort::ini()
{
	iniTrp();

	// TODO: figure out what this is doing and make sure
	// it still works after re-organization
	Ploc_inds = { { -1, -1, -1, 0 } };
	
	if (vtr.restartRunFlag == false)
	{
		Ploc(0) = Ploc_inds.lo;
		Ploc(1) = Ploc_inds.hi;
	}

	numMerges = 0;
	localRange = WORKGROUPSIZE_SORT;
	size_t log2BlockSize = vtr.nN / WORKGROUPSIZE_SORT;
	for (; log2BlockSize > 1; log2BlockSize >>= 1)
	{
		++numMerges;
	}

	if (numMerges & 1)
		oddMergesFlag = true;
	else
		oddMergesFlag = false;
	

	size_t vecPow2 = (vtr.nN & (vtr.nN - 1));
	numMerges += vecPow2 ? 1 : 0;
}


// TODO: Check that the implementations using offset_y, top_location and bot_location
//			still work correctly
void particleSort::iniTrp()
{
	// find start inds
	for (int i = 0; i < vls.nBL; i++)
	{
		if ((vls.C[vls.BL(i, 0)].x <= vtr.X_release) &&
			(vls.C[vls.BL(i, 1)].x > vtr.X_release))
		{
			BL_rel_bot = i;
			LSC_rel_bot = vls.BL(i, 0);
			break;
		}
	}
	
	for (int i = (int)(vls.nBL / 2) + 1; i < vls.nBL; i++)
	{
		if ((vls.C[vls.BL(i, 1)].x <= vtr.X_release) &&
			(vls.C[vls.BL(i, 0)].x > vtr.X_release))
		{
			BL_rel_top = i;
			LSC_rel_top = vls.BL(i, 1);
			break;
		}
	}

	for (int i = 0; i < vls.nBL; i++)
	{
		if ((vls.C[vls.BL(i, 0)].x <= (double)vtr.xStopPos) &&
			(vls.C[vls.BL(i, 1)].x > (double)vtr.xStopPos))
		{
			BL_stop_bot = i;
			LSC_stop_bot = vls.BL(i, 1);
			break;
		}
	}

	for (int i = (int)(vls.nBL / 2) + 1; i < vls.nBL; i++)
	{
		if ((vls.C[vls.BL(i, 1)].x <= (double)vtr.xStopPos) && 
			(vls.C[vls.BL(i, 0)].x > (double)vtr.xStopPos))
		{
			BL_stop_top = i;
			LSC_stop_top = vls.BL(i, 0);
			break;
		}
	}

	
	
	int y0 = 0;
	int x0 = (int)vtr.X_release - 1;
	while (vls.M(x0, y0) == LB_SOLID)
		y0++;
	
	while (vls.M(x0, y0) != LB_SOLID)
		y0++;
	

	
	trP[TRP_UMAX_VAL_IND] = 3. / 2. * vlb.Re * vlb.MuVal / p.Pipe_radius;
	double Xtop_val = vls.C[vls.BL(BL_rel_top, 1)].y + (vtr.X_release - vls.C[vls.BL(BL_rel_top, 1)].x) *
		(vls.C[vls.BL(BL_rel_top, 0)].y - vls.C[vls.BL(BL_rel_top, 1)].y) /
		(vls.C[vls.BL(BL_rel_top, 0)].x - vls.C[vls.BL(BL_rel_top, 1)].x);
	trP[TRP_OFFSET_Y_IND] = vls.C[vls.BL(BL_rel_bot, 1)].y + (vtr.X_release - vls.C[vls.BL(BL_rel_bot, 1)].x) *
		(vls.C[vls.BL(BL_rel_bot, 0)].y - vls.C[vls.BL(BL_rel_bot, 1)].y) /
		(vls.C[vls.BL(BL_rel_bot, 0)].x - vls.C[vls.BL(BL_rel_bot, 1)].x);
	
	trP[TRP_BVAL_IND] = Xtop_val - trP[TRP_OFFSET_Y_IND];
	
	trP[TRP_BOT_LOC_IND] = ceil(trP[TRP_OFFSET_Y_IND]) - trP[TRP_OFFSET_Y_IND];
	trP[TRP_TOP_LOC_IND] = floor(Xtop_val) - trP[TRP_OFFSET_Y_IND];
}


void particleSort::initialSort()
{
	this->operator()();

	if (Ploc(1).y > 0)
		vtr.wallShear.updateParRemArgs();

	if (Ploc(0).y > 0)
		reReleasePar();

	sortTimer = numStepsBtwSort;
}

void particleSort::loadParams()
{
	avgParPerRelease = p.getParameter("Avg Par Per Release", AVG_PAR_PER_RELEASE);
	numStepsBtwSort = p.getParameter("Num Steps Btw Sort", NUM_STEPS_BTW_SORT);
	inletConcDtDivDx = vtr.parP.inletConcPerDx * (double)numStepsBtwSort;
	sortTimer = p.getParameter("Sort Timer", NUM_STEPS_BTW_SORT);
}


void particleSort::operator()(cl_event *TR_prev_evt, int numevt)
{
	cl_event ioEvt;
	int n1 = -1;
	Ploc.FillBuffer((void*)& Ploc_inds, sizeof(cl_int4), 1,
		0, IOQUEUE_REF, numevt, TR_prev_evt);
	Ploc.FillBuffer((void*)& n1, sizeof(int), 2 * vtr.TrDomainSize.x *
		vtr.TrDomainSize.y, 4, IOQUEUE_REF);

	Ptemp.copyToParOnDevice(vtr.P, false, vtr.nN, IOQUEUE_REF, 0,
		nullptr, &ioEvt);
	   	
	sortKernel[0].call_kernel();
	
	for (size_t pass = 1; pass <= numMerges; ++pass)
	{
		unsigned int srcLogicalBlockSize = static_cast< unsigned int>(localRange << (pass - 1));
		if (pass & 0x1)
		{
			sortKernel[1].setOptionCallKernel(&srcLogicalBlockSize);
		}
		else
		{
			sortKernel[2].setOptionCallKernel(&srcLogicalBlockSize);
		}
	}

	sortKernel[3].call_kernel(IOQUEUE_REF, 1, &ioEvt);
	clReleaseEvent(ioEvt);

	Ploc.read_from_buffer_size(2, TRQUEUE_REF, CL_FALSE);
	sortTimer = numStepsBtwSort;
}

void particleSort::reReleasePar()
{
	getUmaxKernel.call_kernel();
	cl_uint maxval = Ploc(0).y;
	double Umean = vlb.calcUmean();
	cl_uint par_in = MAX((cl_uint)(inletConcDtDivDx*Umean) / maxval, 1);
	trReReleaseKernel.set_argument(4, &maxval);
	trReReleaseKernel.set_argument(5, &par_in);
	trReReleaseKernel.set_global_call_kernel(maxval);
}

void particleSort::saveDebug()
{
	Ploc.save_txt_from_device();
}


void particleSort::saveParams()
{
	p.setParameter("Avg Par Per Rel", avgParPerRelease);
	p.setParameter("Num Steps Btw Sort", numStepsBtwSort);
	p.setParameter("Sort Timer", sortTimer);
}

void particleSort::saveRestartFiles()
{
	Ploc.save_bin_from_device("Ploc");
}


#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr()
void particleSort::setSourceDefines()
{
	int count = 1;
	int tempvar = 2;
	while (tempvar < WORKGROUPSIZE_SORT)
	{
		count++;
		tempvar *= 2;
	}

	setSrcDefinePrefix, "SORT_NUM_MERGES", count);
	
	// Trp indicies
	setSrcDefinePrefix, "TRP_X_RELEASE", vtr.X_release);
	setSrcDefinePrefix, "TRP_TOP_LOC_IND", TRP_TOP_LOC_IND);
	setSrcDefinePrefix, "TRP_BOT_LOC_IND", TRP_BOT_LOC_IND);
	setSrcDefinePrefix, "TRP_UMAX_VAL_IND", TRP_UMAX_VAL_IND);
	setSrcDefinePrefix, "TRP_BVAL_IND", TRP_BVAL_IND);
	setSrcDefinePrefix, "TRP_OFFSET_Y_IND", TRP_OFFSET_Y_IND);

	// BL and LSC info
	setSrcDefinePrefix, "BL_REL_BOT", BL_rel_bot);
	setSrcDefinePrefix, "BL_REL_TOP", BL_rel_top);
	setSrcDefinePrefix, "BL_STOP_BOT", BL_stop_bot);
	setSrcDefinePrefix, "BL_STOP_TOP", BL_stop_top);
	setSrcDefinePrefix, "LSC_REL_BOT", LSC_rel_bot);
	setSrcDefinePrefix, "LSC_REL_TOP", LSC_rel_top);
	setSrcDefinePrefix, "LSC_STOP_BOT", LSC_stop_bot);
	setSrcDefinePrefix, "LSC_STOP_TOP", LSC_stop_top);

	int y0 = 0;
	int x0 = (int)vtr.X_release - 1;
	while (vls.M(x0, y0) == LB_SOLID)
		y0++;

	int Uvals_start = x0 + y0 * p.XsizeFull;
	setSrcDefinePrefix, "UVALS_START_INDEX", Uvals_start);

	if (oddMergesFlag)
	{
		setSrcDefinePrefix, "ODD_NUM_MERGES");
	}

}
#undef setSrcDefinePrefix




void particleSort::sortParticlesForClumping()
{
	this->operator()();
	sortTimer = numStepsBtwSort;	
}



bool particleSort::testRestartRun()
{
	bool ret = Ploc.load("load" SLASH "Ploc");
	ret &= p.getParameter("Sort Timer", sortTimer, NUM_STEPS_BTW_SORT);
	return ret;
}




void particleSort::updateTrp()
{
	double Xtop_val = vls.C[vls.BL(BL_rel_top, 1)].y + (vtr.X_release - vls.C[vls.BL(BL_rel_top, 1)].x) *
		(vls.C[vls.BL(BL_rel_top, 0)].y - vls.C[vls.BL(BL_rel_top, 1)].y) /
		(vls.C[vls.BL(BL_rel_top, 0)].x - vls.C[vls.BL(BL_rel_top, 1)].x);
	trP[TRP_OFFSET_Y_IND] = vls.C[vls.BL(BL_rel_bot, 1)].y + (vtr.X_release - vls.C[vls.BL(BL_rel_bot, 1)].x) *
		(vls.C[vls.BL(BL_rel_bot, 0)].y - vls.C[vls.BL(BL_rel_bot, 1)].y) /
		(vls.C[vls.BL(BL_rel_bot, 0)].x - vls.C[vls.BL(BL_rel_bot, 1)].x);
	trP[TRP_BVAL_IND] = Xtop_val - trP[TRP_OFFSET_Y_IND];
	trP[TRP_BOT_LOC_IND] = ceil(trP[TRP_OFFSET_Y_IND]) - trP[TRP_OFFSET_Y_IND];
	trP[TRP_TOP_LOC_IND] = floor(Xtop_val) - trP[TRP_OFFSET_Y_IND];
		
	int y0 = 0;
	int x0 = (int)vtr.X_release - 1;
	while (vls.M(x0, y0) == LB_SOLID)
		y0++;
	
	while (vls.M(x0, y0) != LB_SOLID)
	{
		y0++;
	}
	trP.copyToDevice();
}

