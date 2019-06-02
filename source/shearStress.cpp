#include "shearStress.h"
#include "clVariablesTR.h"

void shearStress::allocateArrays()
{
	shearNodeSize = vls.ssArr.curSize();
	shearNodeDynSize = vls.ssArr.curFullSize();
	
	shearBLStartBot = vtr.Bounds.MIN_BL_BOT;
	shearBLStartTop = vtr.Bounds.MIN_BL_TOP;
	shearBLStopBot = vtr.Bounds.MAX_BL_BOT;
	shearBLStopTop = vtr.Bounds.MAX_BL_TOP;
	shearBLSizeTop = shearBLStartTop - shearBLStopTop + 1;
	shearBLSizeBot = shearBLStartTop - shearBLStopTop + 1;
	shearBLSizeTotal = shearBLSizeTop + shearBLSizeBot;

	sInds.allocate(shearBLSizeTotal);
	blInds.allocate(shearBLSizeTotal);
	ssWeights.allocate(shearBLSizeTotal);

	shearCoeffs.allocate(shearNodeSize*2, shearNodeDynSize * 2);
	Tau.allocate(shearNodeSize, shearNodeDynSize);

	// Not necessary for restart, so no need to place this in 
	if (Tau.load("load" SLASH "Tau") == false)
		Tau.fill(0);

	ssOutput.zeros(shearBLSizeBot, maxOutLinesSS);
	
	blIndsLoc.zeros(p.nX);
}

void shearStress::allocateBuffers()
{
	Tau.allocate_buffer_w_copy();
	
	sInds.allocate_buffer_w_copy();
	
	ssWeights.allocate_buffer_w_copy();
	
	shearCoeffs.allocate_buffer_w_copy();

	blInds.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	
	ssOutput.allocate_buffer_w_copy();
	
	blIndsLoc.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
	blIndsLoc.FillBuffer({ { -1, -1 } });
}

// TODO: set kernel sizes now that variables are created before initializing kernels
void shearStress::createKernels()
{
	// Worksize is just a dummy. Actual size will be set before each kernel call.
	trShearRemovalKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_shear_removal");
	trShearRemovalKernel.set_size(2 * WORKGROUPSIZE_RERELEASE, WORKGROUPSIZE_RERELEASE);
	
	nodeShearKernel.create_kernel(GetSourceProgram, LBQUEUE_REF, "TR_shear_1", "TR_shear_1");
	nodeShearKernel.set_size(shearNodeSize, WORKGROUPSIZE_TR_SHEAR);

	wallShearKernel.create_kernel(GetSourceProgram, LBQUEUE_REF, "TR_shear_2");
	wallShearKernel.set_size(shearBLSizeTotal, WORKGROUPSIZE_TR_SHEAR);

	saveShearKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_shear_save_out");
	saveShearKernel.set_size(shearBLSizeTotal, WORKGROUPSIZE_TR_SHEAR);

	updateSSKernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_Find_Shear_Coeffs1");
	updateSSKernel[0].set_size(shearNodeSize, WORKGROUPSIZE_TR_SHEAR);

	updateSSKernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_Find_Shear_Coeffs2");
	updateSSKernel[1].set_size(shearBLSizeTotal, WORKGROUPSIZE_TR_SHEAR);
}

void shearStress::freeHostArrays()
{
	blIndsLoc.FreeHost();
}

void shearStress::ini()
{
	Save_loc_SS = 0;

	// Sizes of arrays and kernels are function of number of BL's
	// and therefore constant (will not need reallocation at update step)
	numbl_bounds = (vtr.Bounds.MAX_BL_BOT - vtr.Bounds.MIN_BL_BOT) + (vtr.Bounds.MAX_BL_TOP - vtr.Bounds.MIN_BL_TOP) + 2;
	shearSize = (int)ceil((double)numbl_bounds / WORKGROUPSIZE_TR_SHEAR) * WORKGROUPSIZE_TR_SHEAR;
	
	////////////////////////////////////////////////////////////////////////////////////////

	// Sizes of arrays and kernels are not constant
	// and may need reallocation at update step
	Shear_array_len = vls.fullsize_Bnodes;
}


void shearStress::iniShearCoeffs()
{
	updateSSKernel[0].call_kernel();
	updateSSKernel[1].call_kernel();
	cl_int2 negOnes = {{ -1, -1 }};
	Bindicies_loc.FillBuffer(negOnes,TRQUEUE_REF);
}


void shearStress::loadParams()
{
	bool saveDefault = vtr.trSolverFlag & SAVE_WALL_SHEAR;
	calcSSFlag = p.getParameter("Save Shear Stress", saveDefault);
	maxOutLinesSS = p.getParameter("Max Out Lines SS", OUTPUT_MAX_LINES_SS);
	
	std::string saveloc_ = p.getParameter<std::string>("Save Shear Loc", SAVE_SHEAR_LOC);
	if (saveloc_.compare("bothWalls") == 0)
		ssLocSave = bothWalls;
	else if (saveloc_.compare("topWall") == 0)
		ssLocSave = topWall;
	else if (saveloc_.compare("botWall") == 0)
		ssLocSave = botWall;
	else
	{
		ERROR_CHECKING(true, "Save Shear Loc parameter must be either bothWalls, "\
			"topWall or botWall\n", ERROR_INITIALIZING_VTR)
	}
}

void shearStress::saveDebug()
{
	sInds.save_txt_from_device();
	blInds.save_txt_from_device();
	ssWeights.save_txt_from_device();
	shearInds.save_txt_from_device();
	shearCoeffs.save_txt_from_device();
}

void shearStress::saveParams()
{
	p.setParameter("Save Shear Stress", calcSSFlag);
	
	// Not reading in SS arrays, so will start as 0 always
	//p.setParameter("Save_loc_SS", Save_loc_SS);
	p.setParameter("Max Out Lines SS", maxOutLinesSS);
}


void shearStress::saveSS()
{ // For debugging, not used in production runs.
	//Tau.save_from_device("Tau");
	//vls.dir_array.savetxt("dir_array");
	//vls.ii0_array.savetxt("ii0_array");
	//vls.Bnodes.savetxt("Bnodes");
	//BL.savetxt("BLtemp");
	//Sind.save_from_device("Sind");
	//Weights.save_from_device("Weights");
	//BLindicies.save_from_device("Blindicies");
	//Shear_coeffs.save_from_device("SScoeff");
	//Shear_inds.save_txt_from_device("Sinds");
	//vlb.J_array.save_txt_from_device("J_array");
	//vlb.Ro_array.save_txt_from_device("Ro_array");
	//vlb.FA.save_txt_from_device("FA");
	//vlb.FB.save_txt_from_device("FB");
}


void shearStress::setKernelArgs()
{
	// Kernel to calcululate shear at boundary nodes
	int ind = 0;
	nodeShearKernel.set_argument(ind++, vlb.Ro_array.get_buf_add());
	nodeShearKernel.set_argument(ind++, vlb.Ux_array.get_buf_add());
	nodeShearKernel.set_argument(ind++, vlb.Uy_array.get_buf_add());
	nodeShearKernel.set_argument(ind++, shearCoeffs.get_buf_add());
	nodeShearKernel.set_separate_arguments(ind++, vlb.FB.get_buf_add(), vlb.FA.get_buf_add());
	nodeShearKernel.set_argument(ind++, Tau.get_buf_add());
	nodeShearKernel.set_argument(ind++, vls.ssArr.get_buf_add());
	nodeShearKernel.set_argument(ind++, vlb.kOmegaClass.Nut_array.get_buf_add());
	nodeShearKernel.setOptionInd(ind);
	nodeShearKernel.set_argument(ind++, vls.ssArr.curSizeAdd());
	nodeShearKernel.set_local_memory(ind++, 8 * WORKGROUPSIZE_TR_SHEAR * sizeof(double));

	// Kernel to calcululate shear at boundary links
	ind = 0;
	wallShearKernel.set_argument(ind++, ssWeights.get_buf_add());
	wallShearKernel.set_argument(ind++, Tau.get_buf_add());
	wallShearKernel.set_argument(ind++, sInds.get_buf_add());
	wallShearKernel.set_argument(ind++, vtr.BL.Tau.get_buf_add());



	// Kernel to determine if shear is sufficient to remove particles, and removes them
	// if it is.
	BLinks::arrName BLArrList[] = { BLinks::P01Arr, BLinks::vNArr,
		BLinks::lenArr, BLinks::tauArr, BLinks::typeArr };

	Par::arrName ParArrList[] = { Par::posArr, Par::numRepArr,
		Par::typeArr, Par::depFlagArr, Par::depTimerArr, Par::locArr };

	PParam::arrName PParamArrList[] = { PParam::tauCritArr, PParam::dCoeffArr,
		PParam::lCoeffArr, PParam::mpArr };

	ind = 0;
	cl_int2 zer2 = { {0,0} };
	trShearRemovalKernel.set_argument(ind++, vls.C.get_buf_add());
	vtr.BL.setBuffers(trShearRemovalKernel, ind, BLArrList, 5);
	vtr.P.setBuffers(trShearRemovalKernel, ind, ParArrList, 6);
	vtr.parP.setBuffers(trShearRemovalKernel, ind, PParamArrList, 4);
	trShearRemovalKernel.set_argument(ind++, vtr.BL_dep.get_buf_add());
	trShearRemovalKernel.setOptionInd(ind);
	trShearRemovalKernel.set_argument(ind++, &zer2);

	updateSSKernel[0].set_argument(ind++, vls.ssArr.get_buf_add());
	updateSSKernel[0].set_argument(ind++, vlb.NodeType.get_buf_add());
	updateSSKernel[0].set_argument(ind++, shearCoeffs.get_buf_add());
	updateSSKernel[0].setOptionInd(ind);
	updateSSKernel[0].set_argument(ind++, vls.ssArr.curSizeAdd());

	ind = 0;
	updateSSKernel[1].set_argument(ind++, vls.C.get_buf_add());
	vtr.BL.setBuffers(updateSSKernel[1], ind, BLArrList, 3);
	updateSSKernel[1].set_argument(ind++, vls.ssArr.get_buf_add());
	updateSSKernel[1].set_argument(ind++, vls.ssArrInds.get_buf_add());
	updateSSKernel[1].set_argument(ind++, ssWeights.get_buf_add());
	updateSSKernel[1].set_argument(ind++, sInds.get_buf_add());
	
	if (calcSSFlag)
	{
		ind = 0;
		saveShearKernel.set_argument(ind++, BLindicies.get_buf_add());
		saveShearKernel.set_argument(ind++, vtr.BL.get_buf_add());
		saveShearKernel.set_argument(ind++, SS_output.get_buf_add());
		saveShearKernel.set_argument(ind++, &vtr.numbl_bounds);
		saveShearKernel.setOptionInd(4);
	}

	iniShearCoeffs();
}

#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr()
void shearStress::setSourceDefines()
{
	setSrcDefinePrefix, "BL_BOT_START", shearBLStartBot);
	setSrcDefinePrefix, "BL_TOP_START", shearBLStartTop);
	setSrcDefinePrefix, "NUM_BL_BOT", shearBLSizeBot);
	setSrcDefinePrefix, "NUM_BL_TOP", shearBLSizeTop);
	setSrcDefinePrefix, "NUM_BL_TOTAL", shearBLSizeTotal);
	setSrcDefinePrefix, "SHEAR_CUTOFF_RADIUS", vtr.cutoffRadius);
	setSrcDefinePrefix, "INDEX_RADIUS_SEARCH", vtr.indRadiusSearch);
	
	switch (ssLocSave)
	{
	case bothWalls:
	{
		setSrcDefinePrefix, "SAVE_FULL_SHEAR");
		break;
	}
	case topWall:
	{
		setSrcDefinePrefix, "SAVE_SHEAR_TOP");
		break;
	}
	case botWall:
	{
		setSrcDefinePrefix, "SAVE_SHEAR_BOT");
		break;
	}
	}



}
#undef setSrcDefinePrefix
bool shearStress::testRestartRun()
{
	return true;
}

void shearStress::updateSS()
{
	saveShearKernel.setOptionCallKernel(&Save_loc_SS);
	clFlush(TRQUEUE);
	Save_loc_SS++;
}

void shearStress::resetSSOut()
{
	Save_loc_SS = 0;
}


void shearStress::updateParRemArgs()
{
	cl_uint offset = vtr.parSort.Ploc(1).x;
	cl_uint total_removal = vtr.parSort.Ploc(1).y - offset;


	trShearRemovalKernel.set_argument(4, &offset);
	trShearRemovalKernel.set_argument(5, &total_removal);
	trShearRemovalKernel.reset_global_size(total_removal);
}



void shearStress::updateShearArrays()
{
	if (vls.fullsize_Bnodes_old != vls.fullsize_Bnodes)
	{
		vls.fullsize_Bnodes = vls.lengthBnodes * 2;
		
		Shear_array_len = vls.fullsize_Bnodes;
		Shear_inds.reallocate_device_only(Shear_array_len );
		Shear_coeffs.reallocate_device_only(2 * vls.fullsize_Bnodes);
		Tau.reallocate_device_only(Shear_array_len);
		vls.Bnodes.reallocate(vls.fullsize_Bnodes);

		clFinish(IOQUEUE);

		shearKernels[1].set_argument(2, Shear_coeffs.get_buf_add());
		shearKernels[1].set_argument(4, Tau.get_buf_add());
		shearKernels[1].set_argument(5, Shear_inds.get_buf_add());

		shearKernels[0].set_argument(2, Shear_coeffs.get_buf_add());
		shearKernels[0].set_argument(4, Tau.get_buf_add());
		shearKernels[0].set_argument(5, Shear_inds.get_buf_add());

		shearKernels[2].set_argument(2, Tau.get_buf_add());
		
		updateSSKernel[0].set_argument(0, vls.Bnodes.get_buf_add());
		updateSSKernel[0].set_argument(4, Shear_coeffs.get_buf_add());
		updateSSKernel[0].set_argument(5, Shear_inds.get_buf_add());

		updateSSKernel[1].set_argument(4, vls.Bnodes.get_buf_add());
	}

	vls.Bnodes.copy_to_buffer(CL_FALSE, vls.nBnodes);
	
	/*int new_size = (int)ceil((double)vls.nBnodes / (double)WORKGROUPSIZE_TR_SHEAR)*WORKGROUPSIZE_TR_SHEAR;*/
	
	Bnode_top_start = vls.bot_ind_end + 1;
	shearKernels[0].reset_global_size(vls.nBnodes);
	shearKernels[1].reset_global_size(vls.nBnodes);
	updateSSKernel[0].reset_global_size(vls.nBnodes);

	updateSSKernel[0].set_argument(7, (void*)&vls.num_el);
	updateSSKernel[0].set_argument(9, (void *)&vls.nBnodes);
	updateSSKernel[0].set_argument(10, (void*)&Bnode_top_start);

	cl_int2 bindicies_end = { { Bnode_top_start, vls.nBnodes } };
	updateSSKernel[1].set_argument(10, (void*)&bindicies_end);

	cl_event UpdateSSevt;

	updateSSKernel[0].call_kernel();
	updateSSKernel[1].call_kernel(&UpdateSSevt);
	
	Bindicies_loc.FillBuffer(TRQUEUE, { { -1, -1 } }, 1, &UpdateSSevt);

	clReleaseEvent(UpdateSSevt);

	shearKernels[0].set_argument(6, (void*)&vls.nBnodes);
	shearKernels[1].set_argument(6, (void*)&vls.nBnodes);
}

















