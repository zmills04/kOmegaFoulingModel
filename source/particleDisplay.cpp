// particleDisplay.cpp: Implementation of particleDisplay class
// methods declared in particleDisplay.h
//
// (c) Zachary Mills, 2019 
//////////////////////////////////////////////////////////////////////


#include "particleDisplay.h"
#include "clVariablesTR.h"

void particleDisplay::allocateArrays()
{
	if (!p.useOpenGL)
		return;

	Tcrit_color.zeros(2);
	
	// Particle array
	P_vbo.zeros(numParGL * 2);
	P_vbo.AllocateColor(numParGL * 3);

	// Top wall array
	LSt_vbo.allocate(vls.nN);
	LSt_vbo.AllocateColor(vls.nN / 2 * 3);
	LSt_vbo.set_color_array({ { 0.f, 0.f, 0.f } });

	// Bottom wall array
	LSb_vbo.allocate(vls.nN);
	LSb_vbo.AllocateColor(vls.nN / 2 * 3);
	LSb_vbo.set_color_array({ { 0.f, 0.f, 0.f } });

	// Original wall position arrays
	LSt0_vbo.allocate(vls.nN);
	LSb0_vbo.allocate(vls.nN);


	// If gridline flag is set, create gridline arrays
	if (glGridFlag)
	{
		numLinesH = p.nX + 1;
		numLinesV = p.nY + 1;
		LinesV.allocate(numLinesV * 4);
		LinesH.allocate(numLinesH * 4);
	}

	// create Color inds in BL class
	vtr.BL.createColorInds();
}

void particleDisplay::allocateBuffers()
{
	if (!p.useOpenGL)
		return;

	iniParticleColors();
	parColorList.allocate_buffer_w_copy();


	P_vbo.CreateVBO_position(POINT_VBO, CL_MEM_WRITE_ONLY);
	P_vbo.copy_to_buffer();
	P_vbo.CreateVBO_color(CL_MEM_WRITE_ONLY);
	Tcrit_color.allocate_buffer_w_copy(CL_MEM_READ_ONLY);

	LSb_vbo.CreateVBO_position(LINE_VBO);
	LSb_vbo.CreateVBO_color(CL_MEM_READ_WRITE);

	LSt_vbo.CreateVBO_position(LINE_VBO, CL_MEM_READ_WRITE);
	LSt_vbo.CreateVBO_color(CL_MEM_READ_WRITE);

	LSb0_vbo.CreateVBO_position(LINE_VBO);
	LSb0_vbo.set_color_vector({ { 0.f, 0.f, 0.f } });

	LSt0_vbo.CreateVBO_position(LINE_VBO);
	LSt0_vbo.set_color_vector({ { 0.f, 0.f, 0.f } });

	if (glGridFlag)
	{
		LinesH.set_color_vector({ { 0.f, 0.f, 0.f } });
		LinesV.set_color_vector({ { 0.f, 0.f, 0.f } });
		LinesH.CreateVBO_position(LINE_VBO);
		LinesV.CreateVBO_position(LINE_VBO);
	}

	vtr.BL.colorInds.copy_to_buffer();

}

void particleDisplay::createKernels()
{
	if (!p.useOpenGL)
		return;
	TR_GL_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_opengl_par");
	TR_GL_kernel.set_size(TRC_NUM_TRACERS, WORKGROUPSIZE_TR_GL);
	
	int updateGLGlobalSize = getGlobalSizeMacro(vls.nN / 2, WORKGROUPSIZE_UPDATE_GL);
	updateGLKernel.create_kernel(GetSourceProgram, IOQUEUE_REF, "update_GL_wall");
	updateGLKernel.set_size(updateGLGlobalSize, WORKGROUPSIZE_UPDATE_GL);
}

void particleDisplay::freeHostArrays()
{
	if (!p.useOpenGL)
		return;
	parColorList.FreeHost();
	Tcrit_color.FreeHost();

	P_vbo.FreeHostColor();
	P_vbo.FreeHost();
	
	LSt_vbo.FreeHostColor();
	LSt_vbo.FreeHost();
	LSb_vbo.FreeHostColor();	
	LSb_vbo.FreeHost();
	LSt0_vbo.FreeHost();
	LSb0_vbo.FreeHost();

	if (glGridFlag)
	{
		LinesH.FreeHost();
		LinesV.FreeHost();
	}
}


void particleDisplay::ini()
{
	if (!p.useOpenGL)
		return;
	Tcrit_color(0) = vtr.parP.tau_crit(vtr.parP.Nd / 2).x;
	Tcrit_color(1) = vtr.parP.tau_crit(vtr.parP.Nd / 2).y;

	// Add kernels to kernel source string for compilation
	sourceGenerator::SourceInstance()->addFile2Kernel("trKernelsDisplay.cl");

	iniBLArrays();
	iniGridLines();
	iniParticleColors();

	allocateBuffers();
	freeHostArrays();
}

void particleDisplay::iniBLArrays()
{
	for (int i = 0; i < vls.nN / 2; i++)
	{
		LSb_vbo[i * 2] = vls.C[i].x;
		LSb_vbo[i * 2 + 1] = vls.C[i].y;
		LSb0_vbo[i * 2] = vls.C0[i].x;
		LSb0_vbo[i * 2 + 1] = vls.C0[i].y;

		LSt_vbo[i * 2] = vls.C[vls.nN - 1 - i].x;
		LSt_vbo[i * 2 + 1] = vls.C[vls.nN - 1 - i].y;
		LSt0_vbo[i * 2] = vls.C0[vls.nN - 1 - i].x;
		LSt0_vbo[i * 2 + 1] = vls.C0[vls.nN - 1 - i].y;
	}


	for (int i = 0; i < vls.nBL / 2; i++)
	{
		vtr.BL.colorInds(i) = i;
		vtr.BL.colorInds(vls.nBL / 2 + i) = -i;
	}
}


void particleDisplay::iniGridLines()
{
	if (!glGridFlag)
		return;

	for (int i = 0; i < numLinesV; i++)
	{
		LinesV(i * 4) = 0.f + (float)i;
		LinesV(i * 4 + 1) = 0.f;
		LinesV(i * 4 + 2) = 0.f + (float)i;
		LinesV(i * 4 + 3) = (float)p.nY;
	}
	for (int i = 0; i < numLinesH; i++)
	{
		LinesH(i * 4) = 0.;
		LinesH(i * 4 + 1) = (float)i;
		LinesH(i * 4 + 2) = (float)p.nX;
		LinesH(i * 4 + 3) = (float)i;
	}
}



void particleDisplay::iniParticleColors()
{
	parColorList.zeros(vtr.parP.Nd * 3);

	P_vbo.set_color_array({ { 1.f, 0.f, 0.f } });

	for (int i = 0; i < numParGL; i++)
	{
		P_vbo(i * 2 + 0) = vtr.P(i).pos.x;
		P_vbo(i * 2 + 1) = vtr.P(i).pos.y;
	}

	if (vtr.parP.Nd == 1)
	{
		parColorList(0) = 1.f;
		return;
	}
	if (vtr.parP.Nd == 2)
	{
		parColorList(0) = 1.f;
		parColorList(5) = 1.f;
		return;
	}
	if (vtr.parP.Nd == 3)
	{
		parColorList(0) = 1.f;
		parColorList(5) = 1.f;
		parColorList(7) = 1.f;
		return;
	}

	int P1 = vtr.parP.Nd / 4;
	int P2 = vtr.parP.Nd / 2;
	int P3 = vtr.parP.Nd / 4 * 3;
	int P4 = vtr.parP.Nd;

	float stepval = 1.f / ((float)P1);
	for (int i = 0; i < P1; i++)
	{
		parColorList(i * 3 + 1) = stepval*(float)i;
		parColorList(i * 3 + 2) = 1.f;
	}

	stepval = 1.f / (float)(P2 - P1);
	for (int i = P1; i < P2; i++)
	{
		parColorList(i * 3 + 2) = 1.f - stepval*(float)(i - P1);
		parColorList(i * 3 + 1) = 1.f;
	}

	stepval = 1.f / (float)(P3 - P2);
	for (int i = P2; i < P3; i++)
	{
		parColorList(i * 3) = stepval*(float)(i - P2);
		parColorList(i * 3 + 1) = 1.f;
	}

	stepval = 1.f / (float)(P4 - P3);
	for (int i = P3; i < P4; i++)
	{
		parColorList(i * 3 + 1) = 1.f - stepval*(float)(i - P3);
		parColorList(i * 3) = 1.f;
	}
}


void particleDisplay::loadParams()
{
	numParPerGLObj = p.getParameter("Paricles Per GL Obj", NUM_PAR_GL_DIV);
	numParGL = (int)floor((double)vtr.nN / numParPerGLObj);
	glGridFlag = p.getParameter<bool>("GL Use Gridlines", OPENGL_GRIDLINES);
	pointSizes = p.getParameter <double> ("GL Point Sizes", POINT_SIZES);
	lineSizes = p.getParameter <double> ("GL Line Sizes", LINE_SIZES);
	glScreenWidth = p.getParameter<int>("GL Screen Width", screenWidth);
	glScreenHeight = p.getParameter<int>("GL Screen Height", screenHeight);
}

void particleDisplay::saveParams()
{
	p.setParameter("Paricles Per GL Obj", numParPerGLObj);
	p.setParameter("GL Use Gridlines", glGridFlag);
	p.setParameter("GL Point Sizes", pointSizes);
	p.setParameter("GL Line Sizes", lineSizes);
	p.setParameter("GL Screen Width", glScreenWidth);
	p.setParameter("GL Screen Height", glScreenHeight);
}

void particleDisplay::setKernelArgs()
{
	if (!p.useOpenGL)
		return;
	
	int ind = 0;
	TR_GL_kernel.set_argument(ind++, vtr.P.pos.get_buf_add());
	TR_GL_kernel.set_argument(ind++, vtr.P.Dep_Flag.get_buf_add());
	TR_GL_kernel.set_argument(ind++, vtr.P.type.get_buf_add());
	TR_GL_kernel.set_argument(ind++, P_vbo.get_buf_add());
	TR_GL_kernel.set_argument(ind++, P_vbo.get_col_buf_add());
	TR_GL_kernel.set_argument(ind++, parColorList.get_buf_add());

	ind = 0;
	int num_nodes = vls.nN / 2;
	updateGLKernel.set_argument(ind++, vls.C.get_buf_add());
	updateGLKernel.set_argument(ind++, LSb_vbo.get_buf_add());
	updateGLKernel.set_argument(ind++, LSt_vbo.get_buf_add());
}

#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(),
void particleDisplay::setSourceDefines()
{
	setSrcDefinePrefix "NUM_PAR_GL_DIV", numParPerGLObj);
}
#undef setSrcDefinePrefix


bool particleDisplay::testRestartRun() 
{
	allocateArrays();
	return true; 
};