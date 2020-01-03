#include "BiCGStabGenerator.h"
#include "BiCGStabSolver.h"
#include "GMRESGenerator.h"
#include "ReduceGenerator.h"
#include "GMRESSolver.h"


const cl_uint GMRESGenerator::reduce_group_size = WORKGROUPSIZE_RED;
GMRESGenerator* GMRESGenerator::b_instance = nullptr;


// Reduce kernel generators are more or less duplicated, but i didnt
// want to waste time trying to get them to work together.
void GMRESGenerator::ini(const int xsize_, const int xsizefull_, const int ysize_)
{
	// all domain sized arrays will be initialized with x size rounded to 
	// nearest lbm blocksize (usually 128). For reduce kernels, the fact that
	// the x dimension is padded is irrelevant, so we can just pass the padded
	// value to the ReduceGenerator initializer. For the BiCGStabSolver, we 
	// will be initializing variable arrays, so we need the actual x and padded
	// x to make sure we keep the padding set to zero. (important for reductions)


	// Just in case ReduceGenerator hasnt been initialized yet
	// we call ini. It will return immediatly if already done
	int wgsize_reduce_base = reduce_group_size;
	reduceNextToLastSize = reduce_group_size;

	ReduceGenerator::ReduceInstance()->ini(xsizefull_, ysize_);
	int wgsize_reduce_final = ReduceGenerator::ReduceInstance()->getFinalWGSize();
	if (wgsize_reduce_final < 128)
	{
		cl_uint last2Size = reduceNextToLastSize * wgsize_reduce_final;
		switch (wgsize_reduce_final)
		{
		case 2:
		{
			wgsize_reduce_final = 16;
			reduceNextToLastSize = 32;
			break;
		}
		case 4:
		{
			wgsize_reduce_final = 32;
			reduceNextToLastSize = 32;
			break;
		}
		case 8:
		{
			wgsize_reduce_final = 32;
			reduceNextToLastSize = 64;
			break;
		}
		case 16:
		{
			wgsize_reduce_final = 64;
			reduceNextToLastSize = 64;
			break;
		}
		case 32:
		{
			wgsize_reduce_final = 64;
			reduceNextToLastSize = 128;
			break;
		}
		case 64:
		{
			wgsize_reduce_final = 128;
			reduceNextToLastSize = 128;
			break;
		}
		} 
	}

	// initializes static variables used in BiCGStab Solver 
	GMRESSolver::ini(xsize_, xsizefull_, ysize_);


	std::ifstream in_kernels("source\\Kernels\\gmresKernels.cl");
	std::string result_kernels((std::istreambuf_iterator<char>(in_kernels)), std::istreambuf_iterator<char>());
	programSrc.append("#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n");
	programSrc.append(result_kernels);

	std::string finalNorm = genFinalNorm2Kernel(wgsize_reduce_final);
	programSrc.append(finalNorm);

	std::string intermedNormFull = genIntermediateNorm2Kernel(reduce_group_size);
	programSrc.append(intermedNormFull);

	
	std::string intermedNormN2L = genIntermediateNorm2KernelNext2Last(wgsize_reduce_final);
	programSrc.append(intermedNormN2L);

	programSrc.append(addVectors);
	programSrc.append(NormalizeResidualKernel);
	programSrc.append(Norm2BaseKernel);
	programSrc.append(NormalizeResidualKernelWithNorm2Base);
	programSrc.append(ReduceRVkKernel);

	std::string finalAddker = addIniGuess;
	findAndReplace(finalAddker, "(REPLACE_WORKGROUP_SIZE)", std::to_string(WORKGROUPSIZE_GMRESPROD));
	programSrc.append(finalAddker);

	// All kernels are generated above, so no need to have a separate call
	// to build the source.
	buildSource();

	LOGMESSAGE("GMRESGenerator singleton initialized");


}




std::string GMRESGenerator::genIntermediateNorm2Kernel(int wgsize)
{
	std::string redker = "";
	addDefine(redker, "WG_SIZE", wgsize);
	redker.append(Norm2IntermediateKernel);
	findAndReplace(redker, "(REPLACE_NAME_MODIFIER)", "intermediateNorm2Full");

	return redker;
}

std::string GMRESGenerator::genIntermediateNorm2KernelNext2Last(int wgsize)
{
	std::string redker = "";
	addDefine(redker, "WG_SIZE", wgsize);
	redker.append(Norm2IntermediateKernel);
	findAndReplace(redker, "(REPLACE_NAME_MODIFIER)", "intermediateNorm2");

	return redker;
}


std::string GMRESGenerator::genFinalNorm2Kernel(int wgsize)
{
	std::string redker = "";
	addDefine(redker, "WG_SIZE", wgsize);
	redker.append(Norm2FinalKernel);


	return redker;
}









const std::string GMRESGenerator::NormalizeResidualKernel = R"(
__kernel
__attribute__((reqd_work_group_size(256,1,1)))
void NormalizeResidual (const uint size,
                __global double* pR,
				const double rho0)
{
    const uint index = get_global_id(0);

    if (index >= size) return;

    pR[index] /= rho0;
}
)";

const std::string GMRESGenerator::addVectors = R"(
__kernel
__attribute__((reqd_work_group_size(256,1,1)))
void sumVectorKernel (const uint size,
                __global double* resVec,
				__global double* sumVec1,
				__global double* sumVec2)
{
    const uint index = get_global_id(0);

    if (index >= size) return;

    resVec[index] /= sumVec1[index]+sumVec2[index];
}
)";


const std::string GMRESGenerator::addIniGuess = R"(
__kernel
__attribute__((reqd_work_group_size((REPLACE_WORKGROUP_SIZE),1,1)))
void addInitialGuess (const uint size,
                __global double* iniGuess,
				__global double* result)
{
    const uint index = get_global_id(0);

    if (index >= size) return;

    result[index] += iniGuess[index];
}
)";

const std::string GMRESGenerator::NormalizeResidualKernelWithNorm2Base = R"(
__kernel
__attribute__((reqd_work_group_size(256,1,1)))
void NormalizeResidualWithNorm2Base ( 
__global double* input,
const double rho0,
__global double *__restrict__ output,
__local double *sdata)
{
	const uint tid = get_local_id(0);
	const uint bid = get_group_id(0);
	const uint gid = get_global_id(0);
	const uint localSize = 256;
	const uint stride = gid * 2;

	sdata[tid] = input[stride]*input[stride] + input[stride+1]*input[stride+1];
	input[stride] /= rho0;
	input[stride+1] /= rho0;
	barrier(CLK_LOCAL_MEM_FENCE);
	for (uint s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0) { output[bid] = sdata[0]; }
}
)";

const std::string GMRESGenerator::Norm2BaseKernel = R"(
__kernel
__attribute__((reqd_work_group_size(256,1,1)))
void baseNorm2 ( 
__global double* input,
__global double *__restrict__ output,
__local double *sdata)
{
	const uint tid = get_local_id(0);
	const uint bid = get_group_id(0);
	const uint gid = get_global_id(0);
	const uint localSize = 256;
	const uint stride = gid * 2;

	sdata[tid] = input[stride]*input[stride] + input[stride+1]*input[stride+1];
	barrier(CLK_LOCAL_MEM_FENCE);
	for (uint s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0) { output[bid] = sdata[0]; }
}
)";


const std::string GMRESGenerator::Norm2IntermediateKernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void (REPLACE_NAME_MODIFIER) (
__global double* input,
__global double *__restrict__ output,
__local double *sdata)
{
	const uint tid = get_local_id(0);
	const uint gid = get_global_id(0);
	const uint bid = get_group_id(0);
	const uint localSize = WG_SIZE;
	const uint stride = gid * 2;
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
		output[bid] = sdata[0];
	}
}
#undef WG_SIZE
)";


const std::string GMRESGenerator::Norm2FinalKernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void FinalNorm2 (
__global double* input,
__global double *__restrict__ output,
__local double *sdata)
{
	const uint tid = get_local_id(0);
	const uint gid = get_global_id(0);
	const uint localSize = WG_SIZE;
	const uint stride = gid * 2;
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
		output[0] = sdata[0];
	}
}
#undef WG_SIZE
)";


const std::string GMRESGenerator::ReduceRVkKernel = R"(
__kernel
__attribute__((reqd_work_group_size(64,1,1)))
void ReduceRVk (
__global double* input,
__global double *__restrict__ output,
__local double *sdata)
{
	// group size = 64x1
	// global size = 64xkrylov_dims
	//  
	const uint tid = get_local_id(0);
	
	const uint bid = get_group_id(0);
	
	const uint stride = bid * 128 + 2*tid;
	sdata[tid] = input[stride] + input[stride+1];
	barrier(CLK_LOCAL_MEM_FENCE);
	for (unsigned int s = 128 >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		output[bid] = sdata[0];
	}
}
#undef WG_SIZE
)";
