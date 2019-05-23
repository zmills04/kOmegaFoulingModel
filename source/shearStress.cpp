#include "shearStress.h"
#include "clVariablesTR.h"

void shearStress::allocateArrays()
{
	Sind.allocate(numbl_bounds);
	BLindicies.allocate(numbl_bounds);
	Weights.allocate(numbl_bounds);

	Shear_inds.allocate(Shear_array_len);
	Shear_coeffs.allocate(Shear_array_len * 2);
	Tau.allocate(Shear_array_len);
	if (Tau.load("load" SLASH "Tau") == false)
		Tau.fill(0);

	SS_output.zeros(maxOutLinesSS, numbl_bounds);

	Bindicies_loc.zeros(p.nX);

}

void shearStress::allocateBuffers()
{
	Tau.allocate_buffer_w_copy();
	Sind.allocate_buffer_w_copy();
	Weights.allocate_buffer_w_copy();
	Shear_coeffs.allocate_buffer_w_copy();
	Shear_inds.allocate_buffer_w_copy();
	BLindicies.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	
	SS_output.allocate_buffer_w_copy();

	
	Bindicies_loc.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
	Bindicies_loc.FillBuffer({ { -1, -1 } });
}

// TODO: set kernel sizes now that variables are created before initializing kernels
void shearStress::createKernels()
{
	trShearRemovalKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_shear_removal");
	trShearRemovalKernel.set_size(2 * WORKGROUPSIZE_RERELEASE, WORKGROUPSIZE_RERELEASE);
	

	int shearnodesize = (int)ceil((double)vls.nBnodes / WORKGROUPSIZE_TR_SHEAR) * WORKGROUPSIZE_TR_SHEAR;
		
	shearKernels[0].create_kernel(GetSourceProgram, LBQUEUE_REF, "TR_shear_1");
	shearKernels[0].set_size(shearnodesize, WORKGROUPSIZE_TR_SHEAR);

	shearKernels[1].create_kernel(GetSourceProgram, LBQUEUE_REF, "TR_shear_1");
	shearKernels[1].set_size(shearnodesize, WORKGROUPSIZE_TR_SHEAR);

	shearKernels[2].create_kernel(GetSourceProgram, LBQUEUE_REF, "TR_shear_2");
	shearKernels[2].set_size(shearSize, WORKGROUPSIZE_TR_SHEAR);

	saveShearKernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_shear_save_out");
	saveShearKernel.set_size(shearSize, WORKGROUPSIZE_TR_SHEAR);

	updateSSKernel[0].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_Find_Shear_Coeffs1");
	updateSSKernel[0].set_size(shearnodesize, WORKGROUPSIZE_TR_SHEAR);

	updateSSKernel[1].create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_Find_Shear_Coeffs2");
	updateSSKernel[1].set_size(shearSize, WORKGROUPSIZE_TR_SHEAR);
}

void shearStress::freeHostArrays()
{
	Bindicies_loc.FreeHost();
}

void shearStress::ini()
{
	Save_loc_SS = 0;

	// Sizes of arrays and kernels are function of number of BL's
	// and therefore constant (will not need reallocation at update step)
	numbl_bounds = (vtr.Bounds.MAX_BL_BOT - vtr.Bounds.MIN_BL_BOT) + (vtr.Bounds.MAX_BL_TOP - vtr.Bounds.MIN_BL_TOP) + 2;
	shearSize = (int)ceil((double)numbl_bounds / WORKGROUPSIZE_TR_SHEAR) * WORKGROUPSIZE_TR_SHEAR;
	
	int curind = 0;

	//TODO: rather than using BLindicies, pass start and number of top
	//      and bot indicies to kernel. Then 
	//		ind = (i < num bot indicies) ? i+bot_start : i + top_start; (i = get_global_id(0))
	for (int i = vtr.Bounds.MIN_BL_BOT; i <= vtr.Bounds.MAX_BL_BOT; i++)
	{
		BLindicies[curind++] = i;
	}

	for (int i = vtr.Bounds.MIN_BL_TOP; i <= vtr.Bounds.MAX_BL_TOP; i++)
	{
		BLindicies[curind++] = i;
	}
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
}

void shearStress::saveDebug()
{
	Sind.save_txt_from_device();
	BLindicies.save_txt_from_device();
	Weights.save_txt_from_device();
	Shear_inds.save_txt_from_device();
	Shear_coeffs.save_txt_from_device();
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
	int ind = 0;
	shearKernels[0].set_argument(ind++, vlb.Ro_array.get_buf_add());
	shearKernels[0].set_argument(ind++, vlb.U_array.get_buf_add());
	shearKernels[0].set_argument(ind++, Shear_coeffs.get_buf_add());
	shearKernels[0].set_argument(ind++, vlb.FA.get_buf_add());
	shearKernels[0].set_argument(ind++, Tau.get_buf_add());
	shearKernels[0].set_argument(ind++, Shear_inds.get_buf_add());
	shearKernels[0].set_argument(ind++, &vls.nBnodes);
	shearKernels[0].set_argument(ind++, &vlb.Fval);
	shearKernels[0].set_argument(ind++, vlb.Alpha_array.get_buf_add());

	ind = 0;
	shearKernels[1].set_argument(ind++, vlb.Ro_array.get_buf_add());
	shearKernels[1].set_argument(ind++, vlb.J_array.get_buf_add());
	shearKernels[1].set_argument(ind++, Shear_coeffs.get_buf_add());
	shearKernels[1].set_argument(ind++, vlb.FB.get_buf_add());
	shearKernels[1].set_argument(ind++, Tau.get_buf_add());
	shearKernels[1].set_argument(ind++, Shear_inds.get_buf_add());
	shearKernels[1].set_argument(ind++, &vls.nBnodes);
	shearKernels[1].set_argument(ind++, &vlb.Fval);
	shearKernels[1].set_argument(ind++, vlb.Alpha_array.get_buf_add());


	ind = 0;
	shearKernels[2].set_argument(ind++, BLindicies.get_buf_add());
	shearKernels[2].set_argument(ind++, Weights.get_buf_add());
	shearKernels[2].set_argument(ind++, Tau.get_buf_add());
	shearKernels[2].set_argument(ind++, Sind.get_buf_add());
	shearKernels[2].set_argument(ind++, vtr.BL.get_buf_add());
	shearKernels[2].set_argument(ind++, &vtr.numbl_bounds);



	ind = 0;//removal
	trShearRemovalKernel.set_argument(ind++, vtr.BL.get_buf_add());
	trShearRemovalKernel.set_argument(ind++, vtr.P.get_buf_add());
	trShearRemovalKernel.set_argument(ind++, vtr.parP.get_buf_add());
	trShearRemovalKernel.set_argument(ind++, vtr.BL_dep.get_buf_add());


	clFinish(IOQUEUE);

	int ind = 0;
	bNodeTopStart = vls.bot_ind_end + 1;
	int Blinks_top_ind_t = vtr.BL.getSizeX() / 2;
	double cutrad = CUTOFF_RADIUS;
	updateSSKernel[0].set_argument(ind++, vls.Bnodes.get_buf_add());
	updateSSKernel[0].set_argument(ind++, vlb.Stor.get_buf_add());
	updateSSKernel[0].set_argument(ind++, vls.dir_array.get_buf_add());
	updateSSKernel[0].set_argument(ind++, vls.ii0_array.get_buf_add());
	updateSSKernel[0].set_argument(ind++, Shear_coeffs.get_buf_add());
	updateSSKernel[0].set_argument(ind++, Shear_inds.get_buf_add());
	updateSSKernel[0].set_argument(ind++, Bindicies_loc.get_buf_add());
	updateSSKernel[0].set_argument(ind++, &vls.num_el);
	updateSSKernel[0].set_argument(ind++, &cutrad);
	updateSSKernel[0].set_argument(ind++, &vls.nBnodes);
	updateSSKernel[0].set_argument(ind++, &bNodeTopStart);

	cutrad = 2.5;
	int bind_rad = INDEX_RADIUS_SEARCH;
	cl_int2 bindicies_end = { { bNodeTopStart, vls.nBnodes } };
	ind = 0;
	updateSSKernel[1].set_argument(ind++, BLindicies.get_buf_add());
	updateSSKernel[1].set_argument(ind++, Weights.get_buf_add());
	updateSSKernel[1].set_argument(ind++, Sind.get_buf_add());
	updateSSKernel[1].set_argument(ind++, vtr.BL.get_buf_add());
	updateSSKernel[1].set_argument(ind++, vls.Bnodes.get_buf_add());
	updateSSKernel[1].set_argument(ind++, Bindicies_loc.get_buf_add());
	updateSSKernel[1].set_argument(ind++, &vtr.numbl_bounds);
	updateSSKernel[1].set_argument(ind++, &cutrad);
	updateSSKernel[1].set_argument(ind++, &Blinks_top_ind_t);
	updateSSKernel[1].set_argument(ind++, &bind_rad);
	updateSSKernel[1].set_argument(ind++, &bindicies_end);

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

void shearStress::setSourceDefines()
{

}

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

















