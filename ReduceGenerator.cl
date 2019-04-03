#define WG_SIZE		(256)
#define FTYPE		double

#define SIZE_TYPE	uint

__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernelBase_Sum_Size256_double ( 
__global const FTYPE *__restrict__ input,
__global FTYPE *__restrict__ output,
__local FTYPE *sdata)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE bid = get_group_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;

	sdata[tid] = input[stride] + input[stride+1];

	barrier(CLK_LOCAL_MEM_FENCE);
	for (SIZE_TYPE s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			//sdata[tid] += sdata[(tid + s)];
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0) { output[bid] = sdata[0]; }
}

#undef SIZE_TYPE
#undef WG_SIZE
#undef FTYPE


#define WG_SIZE		(1)
#define FTYPE		double

#define SIZE_TYPE	uint

__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernelFinal_Sum_Size1_double (
__global const FTYPE *__restrict__ input,
__global FTYPE *__restrict__ output,
__local FTYPE *sdata,
const int Offset
)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;

	sdata[tid] = input[stride] + input[stride+1];

	barrier(CLK_LOCAL_MEM_FENCE);
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		output[Offset] = sdata[0];
	}
}
#undef SIZE_TYPE
#undef WG_SIZE
#undef FTYPE


#define WG_SIZE		(256)
#define FTYPE		double

#define SIZE_TYPE	uint

__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernelBase__fmin_Size256_double ( 
__global const FTYPE *__restrict__ input,
__global FTYPE *__restrict__ output,
__local FTYPE *sdata)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE bid = get_group_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;

	sdata[tid] = fmin((input[stride]), (input[stride + 1]));

	barrier(CLK_LOCAL_MEM_FENCE);
	for (SIZE_TYPE s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			//sdata[tid] += sdata[(tid + s)];
			sdata[tid] = fmin(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0) { output[bid] = sdata[0]; }
}

#undef SIZE_TYPE
#undef WG_SIZE
#undef FTYPE


#define WG_SIZE		(1)
#define FTYPE		double

#define SIZE_TYPE	uint

__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernelFinal__fmin_Size1_double (
__global const FTYPE *__restrict__ input,
__global FTYPE *__restrict__ output,
__local FTYPE *sdata,
const int Offset
)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;

	sdata[tid] = fmin(input[stride], input[stride + 1]);

	barrier(CLK_LOCAL_MEM_FENCE);
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] = fmin(sdata[tid], sdata[(tid + s)]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		output[Offset] = sdata[0];
	}
}
#undef SIZE_TYPE
#undef WG_SIZE
#undef FTYPE


