#include "particleDisplay.h"
#include "clVariablesTR.h"

void particleDisplay::allocateArrays()
{
	if (!p.useOpenGL)
		return;

	Tcrit_color.zeros(2);
	P_vbo.zeros(numParGL * 2);
	P_vbo.AllocateColor(numParGL * 3);
}

void particleDisplay::allocateBuffers()
{
	if (!p.useOpenGL)
		return;

	iniParticleColors();
	parColorList.allocate_buffer_w_copy();
	parColorList.FreeHost();

	P_vbo.CreateVBO_position(POINT_VBO, CL_MEM_WRITE_ONLY);
	P_vbo.copy_to_buffer();
	P_vbo.CreateVBO_color(CL_MEM_WRITE_ONLY);
	Tcrit_color.allocate_buffer_w_copy();

	// These arrays will never be read back and therefore can be freed
	P_vbo.FreeHostColor();
	P_vbo.FreeHost();
}

void particleDisplay::createKernels()
{
	if (!p.useOpenGL)
		return;
	TR_GL_kernel.create_kernel(GetSourceProgram, TRQUEUE_REF, "TR_opengl_par");
	TR_GL_kernel.set_size(TRC_NUM_TRACERS, WORKGROUPSIZE_TR_GL);
}


void particleDisplay::ini()
{
	if (!p.useOpenGL)
		return;
	Tcrit_color(0) = vtr.parP(vtr.parP.Nd / 2).tau_crit.x;
	Tcrit_color(1) = vtr.parP(vtr.parP.Nd / 2).tau_crit.y;
		
	P_vbo.set_color_array({ { 1.f, 0.f, 0.f } });
	
	for (int i = 0; i < numParGL; i++)
	{
		P_vbo(i * 2 + 0) = vtr.P(i).pos.x;
		P_vbo(i * 2 + 1) = vtr.P(i).pos.y;
	}
}
void particleDisplay::iniParticleColors()
{
	parColorList.zeros(vtr.parP.Nd * 3);

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
}

void particleDisplay::setKernelArgs()
{
	if (!p.useOpenGL)
		return;

	vtr.wallShear.shearKernels[2].set_argument(6, vls.LSb_vbo.get_col_buf_add());
	vtr.wallShear.shearKernels[2].set_argument(7, vls.LSt_vbo.get_col_buf_add());
	vtr.wallShear.shearKernels[2].set_argument(8, Tcrit_color.get_buf_add());

	int ind = 0;
	TR_GL_kernel.set_argument(ind++, vtr.P.get_buf_add());
	TR_GL_kernel.set_argument(ind++, P_vbo.get_buf_add());
	TR_GL_kernel.set_argument(ind++, P_vbo.get_col_buf_add());
	TR_GL_kernel.set_argument(ind++, parColorList.get_buf_add());
}

#define setSrcDefinePrefix		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(),
void particleDisplay::setSourceDefines()
{
	setSrcDefinePrefix "NUM_PAR_GL_DIV", numParPerGLObj);
}
#undef setSrcDefinePrefix