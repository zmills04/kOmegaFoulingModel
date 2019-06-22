
// Prepares CL-GL buffers for openGL calls to display data on screen
__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_TR_GL, 1, 1)))
void TR_opengl_par(__global double2 *__restrict__ pPos,
	__global int* __restrict__ pDepFlag,
	__global short* __restrict__ pType,
	__global float2* __restrict__ gl_par,
	__global float* __restrict__ gl_color,
	__global float* __restrict__ gl_color_array)
{
	int i = get_global_id(0);


	if (i >= TRC_NUM_TRACERS)
		return;

	gl_par[i] = convert_float2(pPos[i]);
	float Dep_mult = (pDepFlag[i] > -2) ? (1.f) : (0.f);
	i *= 3;
	int pTypeTemp = pType[i] * 3;
	gl_color[i] = gl_color_array[pTypeTemp] * Dep_mult;
	gl_color[i + 1] = gl_color_array[pTypeTemp + 1] * Dep_mult;
	gl_color[i + 2] = gl_color_array[pTypeTemp + 2] * Dep_mult;
}
