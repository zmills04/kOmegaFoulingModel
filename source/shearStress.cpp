// shearStress.cpp: Implementation of methods defined in shearStress.h 
// which are used to calculate the shear stress in the channel, and
// remove deposited particles based on the magnitude of the local
// shear.
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////


#include "shearStress.h"
#include "clVariablesTR.h"
#include "logger.h"

void shearStress::allocateArrays()
{
	shearNodeSize = vls.ssArr.curSize();
	shearNodeDynSize = vls.ssArr.curFullSize();
	
	shearBLStartBot = vtr.Bounds.MIN_BL_BOT;
	shearBLStartTop = vtr.Bounds.MIN_BL_TOP;
	shearBLStopBot = vtr.Bounds.MAX_BL_BOT;
	shearBLStopTop = vtr.Bounds.MAX_BL_TOP;
	shearBLSizeTop = shearBLStopTop - shearBLStartTop + 1;
	shearBLSizeBot = shearBLStopBot - shearBLStartBot + 1;
	shearBLSizeTotal = shearBLSizeTop + shearBLSizeBot;

	sInds.allocate(shearBLSizeTotal);
	ssWeights.allocate(shearBLSizeTotal);

	shearCoeffs.allocate(shearNodeDynSize * 2);
	Tau.allocate(shearNodeDynSize);

	// Not necessary for restart, so no need to place this in
	// testRestartRun
	if (Tau.load("load" SLASH "Tau") == false)
		Tau.fill(0);

	if (calcSSFlag)
	{
		ssOutput.iniFile(vtr.restartRunFlag);

		if (ssLocSave == bothWalls) { ssOutputSize = shearBLSizeTotal; }
		else if (ssLocSave == topWall) { ssOutputSize = shearBLSizeTop; }
		else { ssOutputSize = shearBLSizeBot; }
	
		ssOutput.zeros(ssOutputSize, maxOutLinesSS);
		ssOutput.createTimeArray();
	}

	blIndsLoc.zeros(p.nX);
}

void shearStress::allocateBuffers()
{
	Tau.allocate_buffer_w_copy();
	
	sInds.allocate_buffer_w_copy();
	
	ssWeights.allocate_buffer_w_copy();
	
	shearCoeffs.allocate_buffer_w_copy();
	
	if(calcSSFlag)
		ssOutput.allocate_buffer_w_copy();
	
	blIndsLoc.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
	blIndsLoc.FillBuffer({ { -1, -1 } });
}

void shearStress::calculateShear(cl_command_queue* que_, cl_event* waitevt,
	int numwait, cl_event* evt_)
{
	// since these need to be called sequentually in the same queue, 
	// the wait event will only be passed to nodeShearKernel, and 
	// only wallShearKernel will be added to evt_.

	// This is not necessarily called everytime collisionKernel is called,
	// so we need to get the current value of alter in collisionKernel.
	nodeShearKernel.call_kernel(vlb.collisionKernel.getAlterID(), que_,
		numwait, waitevt);
	
	wallShearKernel.call_kernel(que_, 0, nullptr, evt_);
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

	if (calcSSFlag)
	{
		saveShearKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_shear_save_out");
		saveShearKernel.set_size(ssOutputSize, WORKGROUPSIZE_TR_SHEAR);
	}

	updateSSKernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_Find_Shear_Coeffs1");
	updateSSKernel[0].set_size(shearNodeSize, WORKGROUPSIZE_TR_SHEAR);

	updateSSKernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_Find_Shear_Coeffs2");
	updateSSKernel[1].set_size(shearBLSizeTotal, WORKGROUPSIZE_TR_SHEAR);
}

// TODO: make sure that all these can be freed
void shearStress::freeHostArrays()
{
	blIndsLoc.FreeHost();
	Tau.FreeHost();
	sInds.FreeHost();
	ssWeights.FreeHost();
	shearCoeffs.FreeHost();
}

void shearStress::ini()
{
	allocateArrays();

	allocateBuffers();

	// Add kernels to kernel source string for compilation
	sourceGenerator::SourceInstance()->addFile2Kernel("trKernelsShear.cl");


	LOGMESSAGE("wall shear class within vtr initialized");
}


void shearStress::iniShearCoeffs()
{
	updateSSKernel[0].call_kernel();
	updateSSKernel[1].call_kernel();
}


void shearStress::loadParams()
{
	calcSSFlag = p.getParameter("Save Shear Stress", SAVE_WALL_SHEAR) &
		vtr.trSolverFlag;

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

void shearStress::save2file()
{
	// nothing gets saved here since only output is 
	// ssOutput, which is a time based data set (saved in saveTimeData)
}

void shearStress::saveDebug()
{
	Tau.save_txt_from_device();
	sInds.save_txt_from_device();
	ssWeights.save_txt_from_device();
	shearCoeffs.save_txt_from_device();
	blIndsLoc.save_txt_from_device();
	//vls.ssArrIndMap.savetxt();
}

void shearStress::saveParams()
{
	p.setParameter("Save Shear Stress", calcSSFlag);
	std::string ssLocSaveString;
	if (ssLocSave == bothWalls) {
		ssLocSaveString = "bothWalls";
	}
	else if (ssLocSave == topWall) {
		ssLocSaveString = "topWall";
	}
	else {
		ssLocSaveString = "botWall";
		
	}
	p.setParameter<std::string>("Save Shear Loc", ssLocSaveString);
	p.setParameter("Max Out Lines SS", maxOutLinesSS);
}

void shearStress::saveRestartFiles()
{
	//nothing to needed for restart
}

void shearStress::saveTimeData()
{
	if (calcSSFlag)
	{
		ssOutput.appendData_from_device(TRQUEUE_REF);
	}
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

	// Add additional kernel arguments when using opengl (first 4 arguments
	// are the same regardless of whether or not opengl is used)
	if (p.useOpenGL)
	{
		wallShearKernel.set_argument(ind++, vtr.BL.int_type.get_buf_add());
		wallShearKernel.set_argument(ind++, vtr.BL.colorInds.get_buf_add());
		wallShearKernel.set_argument(ind++, vtr.glParticles.LSb_vbo.get_col_buf_add());
		wallShearKernel.set_argument(ind++, vtr.glParticles.LSt_vbo.get_col_buf_add());
		wallShearKernel.set_argument(ind++, vtr.glParticles.Tcrit_color.get_buf_add());
	}

	// Kernel to determine if shear is sufficient to remove particles, and removes them
	// if it is.
	BLinks::arrName BLArrList[] = { BLinks::P01Arr, BLinks::vNArr,
		BLinks::lenArr, BLinks::tauArr, BLinks::typeArr };

	Par::arrName ParArrList[] = { Par::posArr, Par::numRepArr,
		Par::typeArr, Par::depFlagArr, Par::depTimerArr, Par::locArr, Par::timerArr };

	ind = 0;


	cl_int2 zer2 = { {0,0} };
	trShearRemovalKernel.set_argument(ind++, vls.C.get_buf_add());
	vtr.BL.setBuffers(trShearRemovalKernel, ind, BLArrList, 5);
	vtr.P.setBuffers(trShearRemovalKernel, ind, ParArrList, 7);
	trShearRemovalKernel.set_argument(ind++, vtr.blDep.get_buf_add());
	trShearRemovalKernel.setOptionInd(ind);
	trShearRemovalKernel.set_argument(ind++, &zer2);


	ind = 0;
	updateSSKernel[0].set_argument(ind++, vls.ssArr.get_buf_add());
	updateSSKernel[0].set_argument(ind++, vls.nType.get_buf_add());
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
		saveShearKernel.set_argument(ind++, vtr.BL.Tau.get_buf_add());
		saveShearKernel.set_argument(ind++, ssOutput.get_buf_add());
		saveShearKernel.setOptionInd(ind);
		//save index is set before calling so no need to set here
		//saveShearKernel.set_argument(ind++, &Save_loc_SS); 
	}

	// shear coefficients are initialized in a kernel, so those
	// kernels are called now that arguments have been set
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
	
	if(ssLocSave == bothWalls) { setSrcDefinePrefix, "SAVE_FULL_SHEAR"); }
	else if (ssLocSave == topWall) { setSrcDefinePrefix, "SAVE_SHEAR_TOP"); }
	else { setSrcDefinePrefix, "SAVE_SHEAR_BOT"); }

}
#undef setSrcDefinePrefix


void shearStress::shearRemoval(cl_command_queue *que_,
	cl_event* waitevt, int numwait, cl_event *evt_)
{
	if ((vtr.parSort.Ploc(1).y - vtr.parSort.Ploc(1).x) > 0)
		trShearRemovalKernel.call_kernel(que_, numwait, waitevt, evt_);
}

bool shearStress::testRestartRun()
{
	// This will be called after vls.iniShearArrays is called.
	//allocateArrays();
	return true;
}

void shearStress::updateParRemArgs()
{
	// {offset, number of particles temporarily deposited}
	cl_int2 removalOptionVal = { { vtr.parSort.Ploc(1).x,
		vtr.parSort.Ploc(1).y - vtr.parSort.Ploc(1).x} };
	   
	trShearRemovalKernel.setOption(&removalOptionVal);
	trShearRemovalKernel.reset_global_size(removalOptionVal.y);
}

// Moved to clVariablesFL to overcome issues associated with 
// using std::thread to update arrays. 

//void shearStress::updateShearArrays()
//{
//	if (!vtr.wallShear.calcSSFlag)
//		return;
//	bool reSizeFlag = vls.updateShearArrays();
//	if (reSizeFlag)
//	{
//		shearNodeDynSize = vls.ssArr.curFullSize();
//		shearCoeffs.reallocate_device_only(2 * shearNodeDynSize);
//		Tau.reallocate_device_only(shearNodeDynSize);
//
//		clFinish(IOQUEUE);
//
//		updateSSKernel[0].set_argument(0, vls.ssArr.get_buf_add());
//
//		updateSSKernel[0].set_argument(2, shearCoeffs.get_buf_add());
//		updateSSKernel[1].set_argument(4, vls.ssArr.get_buf_add());
//		
//		nodeShearKernel.set_argument(3, shearCoeffs.get_buf_add());
//		nodeShearKernel.set_argument(5, Tau.get_buf_add());
//		nodeShearKernel.set_argument(6, vls.ssArr.get_buf_add());
//		
//		wallShearKernel.set_argument(1, Tau.get_buf_add());
//	}
//
//	shearNodeSize = vls.ssArr.curSize();
//
//	updateSSKernel[0].setOption(&shearNodeSize);
//	nodeShearKernel.setOption(&shearNodeSize);
//	
//	nodeShearKernel.reset_global_size(shearNodeSize);
//	updateSSKernel[0].reset_global_size(shearNodeSize);
//
//	updateSSKernel[0].call_kernel();
//	updateSSKernel[1].call_kernel();
//}




void shearStress::updateTimeData()
{
	if (calcSSFlag)
	{
		saveShearKernel.setOptionCallKernel(ssOutput.getCurIndAdd());
		ssOutput.setTimeAndIncrement(p.Time, TRQUEUE_REF);
	}
}












