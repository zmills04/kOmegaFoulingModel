#include "BiCGStabGenerator.h"
#include "BiCGStabSolver.h"
#include "ReduceGenerator.h"


const cl_uint BiCGStabGenerator::WG_BITS = 24;
const cl_uint BiCGStabGenerator::ROW_BITS = 32;
const cl_uint BiCGStabGenerator::BLKSIZE = 1024;
const cl_uint BiCGStabGenerator::BLOCK_MULTIPLIER = 3;
const cl_uint BiCGStabGenerator::ROWS_FOR_VECTOR = 1;
const cl_uint BiCGStabGenerator::csrmv_group_size = WORKGROUPSIZE_CSRMV;
const cl_uint BiCGStabGenerator::reduce_group_size = WORKGROUPSIZE_RED;
const cl_uint BiCGStabGenerator::axpby_group_size = WORKGROUPSIZE_AXPBY;
BiCGStabGenerator* BiCGStabGenerator::b_instance = nullptr;


// Reduce kernel generators are more or less duplicated, but i didnt
// want to waste time trying to get them to work together.
void BiCGStabGenerator::ini(const int xsize_, const int xsizefull_, const int ysize_)
{
	// all domain sized arrays will be initialized with x size rounded to 
	// nearest lbm blocksize (usually 128). For reduce kernels, the fact that
	// the x dimension is padded is irrelevant, so we can just pass the padded
	// value to the ReduceGenerator initializer. For the BiCGStabSolver, we 
	// will be initializing variable arrays, so we need the actual x and padded
	// x to make sure we keep the padding set to zero. (important for reductions)


	// Just in case ReduceGenerator hasnt been initialized yet
	// we call ini. It will return immediatly if already done
	ReduceGenerator::ReduceInstance()->ini(xsizefull_, ysize_);
	int wgsize_reduce_final = ReduceGenerator::ReduceInstance()->getFinalWGSize();

	// initializes static variables used in BiCGStab Solver 
	BiCGStabSolver::ini(xsize_, xsizefull_, ysize_);


	//create clsparse solver control. Currently used for debugging, and will remain
	// throughout program. Eventaully will be deleted when sourceGenerator 
	// builds its source which is after adaptive csr meta data has finished being created.
	int status;
	cl_command_queue SPqueuet = clCreateCommandQueueWithProperties(*clEnv::instance()->getContext(),
		*clEnv::instance()->getDevice(), clEnv::instance()->getQueueProperties(), &status);
	ERROR_CHECKING(status, "Error Creating TRqueue", ERROR_OCL_INITIALIZATION);
	cl::CommandQueue SPqueue(SPqueuet, true);
	clsparseStatus clsstatus = clsparseSetup();
	ERROR_CHECKING(clsstatus != clsparseSuccess, "Error Initializing clSparse", ERROR_OCL_INITIALIZATION);
	clSparseControl = clsparseCreateControl(SPqueue());
	ERROR_CHECKING(clSparseControl.status, "Failed to create clsparse control", ERROR_OCL_INITIALIZATION);




	std::string indt = "uint";
	int wgsize_reduce_base = reduce_group_size;
	int wgsize_axpby = axpby_group_size;
	
	addDefine(programSrc, "VALUE_TYPE", "double");
	
	// Get csrmv kernel
	std::string csrstring = genCSRMVKernel(csrmv_group_size, ROW_BITS, WG_BITS, BLKSIZE,
		BLOCK_MULTIPLIER, ROWS_FOR_VECTOR, indt);
	programSrc.append(csrstring);


	//get first modified axpby kernel
	std::string mod1name = "shCalc";
	// vector pY has pM subtracted, pZ has pP added 
	std::string mod1op = "pZ[index] += alpha * pP[index];\n\tpY[index] -= alpha * pM[index];\n";
	std::string mod1in = "__global double *__restrict__ pZ,\n\t__global const double *__restrict__ pM,\n\t__global const double *__restrict__ pP";
	std::string modaxpby1str = genModifiedAxpbyKernel(wgsize_axpby, indt, mod1in, mod1op, mod1name, 0);
	programSrc.append(modaxpby1str);

	//get second modified axpby kernel
	std::string mod2name = "PCalc";
	// all scalars are in same buffer, so only need to pass that buffer and three corresponding offsets
	std::string mod2op = "const VALUE_TYPE omega = *(pAlpha + pOmegaOffset);\n\t"\
		"pY[index] = pR[index] + alpha * (pY[index] - omega * pV[index]);\n";
	std::string mod2in = "const SIZE_TYPE pOmegaOffset,\n\t"\
		"__global const double *__restrict__ pR,\n\t"\
		"__global const double *__restrict__ pV";
	std::string modaxpby2str = genModifiedAxpbyKernel(wgsize_axpby, indt, mod2in, mod2op, mod2name, 0);
	programSrc.append(modaxpby2str);


	//get first modified axpby kernel
	std::string mod3name = "xrCalc";
	// vector pY has pM subtracted, pZ has pP added 
	std::string mod3op = "pZ[index] += alpha * pY[index];\n\tpY[index] -= alpha * pM[index];\n";
	std::string mod3in = "__global double *__restrict__ pZ,\n\t__global const double *__restrict__ pM";
	std::string modaxpby3str = genModifiedAxpbyKernel(wgsize_axpby, indt, mod3in, mod3op, mod3name, 0);
	programSrc.append(modaxpby3str);

	// Sqr Reduce Base
	std::string sqrredstr = genBaseRedKernel(wgsize_reduce_base, Sqr, "Sqr");
	programSrc.append(sqrredstr);

	// Norm Reduce Base
	std::string normredstr = genBaseRedKernel(wgsize_reduce_base, Norm, "Norm");
	programSrc.append(normredstr);

	// Dot Reduce Base
	std::string dotredstr = genBaseRedKernel(wgsize_reduce_base, Dot, "Dot");
	programSrc.append(dotredstr);

	// Generic Sum Reduce Base
	std::string sumredstr = genBaseRedKernel(wgsize_reduce_base, Sum, "Sum");
	programSrc.append(sumredstr);

	// Generic Final reduce kernel
	std::string modred1name = "Sum_Final";

	std::string modred1in = "__global const double *__restrict__ input";
	std::string modred1out = "__global double *__restrict__ output";
	std::string modred1loc = "__local double *sdata";
	std::string modred1off = "const SIZE_TYPE offset";

	std::string modred1opin = "sdata[tid] = input[stride] + input[stride+1];\n";
	std::string modred1oploc = "sdata[tid] += sdata[(tid + s)];\n";
	std::string modred1opout = "output[offset] = sdata[0];";

	std::string sumfinal1str = genFinalRedKernel(wgsize_reduce_final, modred1in, modred1out,
		modred1loc, modred1off, modred1opin, modred1oploc, modred1opout, modred1name);
	programSrc.append(sumfinal1str);

	// Final reduce for Rho and calculation of Beta
	// TODO: try to see if using other threads for beta is faster (would need to sync first)
	std::string modred2name = "Sum_Rho_Calc_Beta";

	std::string modred2in = "__global const double *__restrict__ input";
	std::string modred2out = "__global double *__restrict__ output";
	std::string modred2loc = "__local double *sdata";
	std::string modred2off = "const SIZE_TYPE offsetRho,\n\t"\
		"const SIZE_TYPE offsetBeta,\n\t"\
		"const SIZE_TYPE offsetAlpha,\n\t"\
		"const SIZE_TYPE offsetOmega";

	std::string modred2opin = "sdata[tid] = input[stride] + input[stride+1];\n";
	std::string modred2oploc = "sdata[tid] += sdata[(tid + s)];\n";
	std::string modred2opout = "double rhoOld = output[offsetRho];"\
		"output[offsetRho] = sdata[0];"\
		"output[offsetBeta] = (sdata[0] / rhoOld)*(output[offsetAlpha]/output[offsetOmega]);";

	std::string sumfinal2str = genFinalRedKernel(wgsize_reduce_final, modred2in, modred2out, modred2loc, modred2off,
		modred2opin, modred2oploc, modred2opout, modred2name);
	programSrc.append(sumfinal2str);

	// Final reduce for r dot v and calculation of alpha
	std::string modred3name = "Calc_Alpha";

	std::string modred3in = "__global const double *__restrict__ input";
	std::string modred3out = "__global double *__restrict__ output";
	std::string modred3loc = "__local double *sdata";
	std::string modred3off = "const SIZE_TYPE offsetRho,\n\t"\
		"const SIZE_TYPE offsetAlpha";

	std::string modred3opin = "sdata[tid] = input[stride] + input[stride+1];\n";
	std::string modred3oploc = "sdata[tid] += sdata[(tid + s)];\n";
	std::string modred3opout = "output[offsetAlpha] = output[offsetRho] / sdata[0];";


	std::string sumfinal3str = genFinalRedKernel(wgsize_reduce_final, modred3in, modred3out,
		modred3loc, modred3off, modred3opin, modred3oploc, modred3opout, modred3name);
	programSrc.append(sumfinal3str);

	// It was easier to just create specific kernels to reduce s dot t and t dot t
	addDefine(programSrc, "WG_SIZE", WORKGROUPSIZE_RED);
	addDefine(programSrc, "SIZE_TYPE", "uint");
	programSrc.append(ReduceBaseOmega_kernel);

	addDefine(programSrc, "WG_SIZE", WORKGROUPSIZE_RED);
	addDefine(programSrc, "SIZE_TYPE", "uint");
	programSrc.append(ReduceIntermediateOmega_kernel);

	addDefine(programSrc, "WG_SIZE", wgsize_reduce_final);
	addDefine(programSrc, "SIZE_TYPE", "uint");
	programSrc.append(ReduceFinalOmega_kernel);


	// Add elementwise operations
	programSrc.append(genElementWiseOperation(wgsize_axpby, indt, "(A + B)", "Add"));
	programSrc.append(genElementWiseOperation(wgsize_axpby, indt, "(A - B)", "Subtract"));
	programSrc.append(genElementWiseOperation(wgsize_axpby, indt, "(A * B)", "Multiply"));
	programSrc.append(genElementWiseOperation(wgsize_axpby, indt, "(A / B)", "Divide"));


	// All kernels are generated above, so no need to have a separate call
	// to build the source.
	buildSource();
}


std::string BiCGStabGenerator::genAxpbyKernel(int wgsize, std::string indType, std::string oper,
	std::string repname_, int offset)
{
	std::string axpyker = "";
	addDefine(axpyker, "WG_SIZE", wgsize);
	addDefine(axpyker, "SIZE_TYPE", indType);
	axpyker.append("#define " + oper + "\n");
	axpyker.append(Axpby_kernel);
	std::string offsetstr = "";
	if (offset != 0)
	{
		offsetstr.append(" + ");
		offsetstr.append(std::to_string(offset));
	}
	findAndReplace(axpyker, "(REPLACE_NAME_MODIFIER)", repname_);
	findAndReplace(axpyker, "(REPLACE_OFFSET)", offsetstr);
	return axpyker;
}


std::string BiCGStabGenerator::genAxpyKernel(int wgsize, std::string indType, std::string oper,
	std::string repname_, int offset)
{
	std::string axpyker = "";
	addDefine(axpyker, "WG_SIZE", wgsize);
	addDefine(axpyker, "SIZE_TYPE", indType);

	axpyker.append("#define " + oper + "\n");
	axpyker.append(Axpy_kernel);
	std::string offsetstr = "";
	if (offset != 0)
	{
		offsetstr.append(" + ");
		offsetstr.append(std::to_string(offset));
	}
	findAndReplace(axpyker, "(REPLACE_NAME_MODIFIER)", repname_);
	findAndReplace(axpyker, "(REPLACE_OFFSET)", offsetstr);
	return axpyker;
}


std::string BiCGStabGenerator::genBaseRedKernel(int wgsize, BaseRedType redtype_, std::string name_)
{
	std::string redker = "";
	addDefine(redker, "WG_SIZE", wgsize);
	addDefine(redker, "SIZE_TYPE", "uint");
	redker.append(ReduceBase_kernel);


	std::string inputstr = "";
	std::string operstr = "";
	if (redtype_ == Sum)
	{
		inputstr.append("__global double *__restrict__ input");
		operstr.append("input[stride] + input[stride+1]");
	}
	else if (redtype_ == Dot)
	{
		inputstr.append("__global double *__restrict__ inputA,\n"\
			"__global double *__restrict__ inputB");
		operstr.append("inputA[stride] * inputB[stride] + inputA[stride+1] * inputB[stride+1]");
	}
	else if (redtype_ == Abs || redtype_ == Norm)
	{
		inputstr.append("__global double *__restrict__ input");
		operstr.append("fabs(input[stride]) + fabs(input[stride+1])");
	}
	else if (redtype_ == Sqr)
	{
		inputstr.append("__global double *__restrict__ input");
		operstr.append("(input[stride]*input[stride]) + (input[stride+1]*input[stride+1])");
	}
	else
	{
		printf("error: tried to create a BaseRed kernel with incorrect type\n");
		exit(-10442);
	}
	findAndReplace(redker, "(REPLACE_INPUTS)", inputstr);
	findAndReplace(redker, "(REPLACE_OPERATION)", operstr);
	findAndReplace(redker, "(REPLACE_REDUCE_TYPE)", name_);

	return redker;
}

std::string BiCGStabGenerator::genCSRMVKernel(int wgsize, int rowbits, int wgbits, int blocksize,
	int blockmult, int rows4vec, std::string indt)
{
	std::string csrmvker = "";
	addDefine(csrmvker, "WG_SIZE", wgsize);
	addDefine(csrmvker, "ROWBITS", rowbits);
	addDefine(csrmvker, "WGBITS", wgbits);
	addDefine(csrmvker, "BLOCKSIZE", blocksize);
	addDefine(csrmvker, "BLOCK_MULTIPLIER", blockmult);
	addDefine(csrmvker, "ROWS_FOR_VECTOR", rows4vec);
	addDefine(csrmvker, "INDEX_TYPE", indt);
	std::ifstream in_csrmv("source\\Kernels\\csrmvKernels.cl");
	std::string result_csrmv((std::istreambuf_iterator<char>(in_csrmv)), std::istreambuf_iterator<char>());
	csrmvker.append(result_csrmv);
	return csrmvker;
}

std::string BiCGStabGenerator::genFinalRedKernel(int wgsize, std::string repIn_, std::string repOut_, std::string repLocal_,
	std::string repOff_, std::string repOpIn_, std::string repOpLocal_, std::string repOpOut_, std::string name_)
{
	std::string redker = "";
	addDefine(redker, "WG_SIZE", wgsize);
	addDefine(redker, "SIZE_TYPE", "uint");
	redker.append(ReduceFinal_kernel);
	findAndReplace(redker, "(REPLACE_REDUCE_TYPE)", name_);
	findAndReplace(redker, "(REPLACE_INPUTS)", repIn_);
	findAndReplace(redker, "(REPLACE_OUTPUTS)", repOut_);
	findAndReplace(redker, "(REPLACE_LOCALS)", repLocal_);
	findAndReplace(redker, "(REPLACE_OFFSETS)", repOff_);
	findAndReplace(redker, "(REPLACE_OPERATIONS_INPUT)", repOpIn_);
	findAndReplace(redker, "(REPLACE_OPERATIONS_LOCAL)", repOpLocal_);
	findAndReplace(redker, "(REPLACE_OPERATIONS_OUTPUT)", repOpOut_);

	return redker;
}

std::string BiCGStabGenerator::genModifiedAxpbyKernel(int wgsize, std::string indType,
	std::string repin_, std::string repop_, std::string repname_, int offset)
{
	std::string axpyker = "";
	addDefine(axpyker, "WG_SIZE", wgsize);
	addDefine(axpyker, "SIZE_TYPE", indType);
	axpyker.append(Modified_Axpby_kernel);
	std::string offsetstr = "";
	if (offset != 0)
	{
		offsetstr.append(" + ");
		offsetstr.append(std::to_string(offset));
	}


	findAndReplace(axpyker, "(REPLACE_NAME_MODIFIER)", repname_);
	findAndReplace(axpyker, "(REPLACE_INPUTS)", repin_);
	findAndReplace(axpyker, "(REPLACE_OPERATION)", repop_);
	//findAndReplace(axpyker, "(REPLACE_OFFSET)", offsetstr);
	return axpyker;
}


std::string BiCGStabGenerator::genElementWiseOperation(int wgsize, std::string indType,
	std::string opdefine_, std::string repname_)
{
	std::string elwiseker = "";
	addDefine(elwiseker, "WG_SIZE", wgsize);
	addDefine(elwiseker, "SIZE_TYPE", indType);
	addDefine(elwiseker, "OPERATION(A,B)", opdefine_);



	elwiseker.append(ElementWiseOp_kernel);


	findAndReplace(elwiseker, "(REPLACE_NAME)", repname_);
	return elwiseker;
}




const std::string BiCGStabGenerator::ElementWiseOp_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ElementWise_(REPLACE_NAME) (const SIZE_TYPE size,
                __global VALUE_TYPE* pR,
                __global const VALUE_TYPE* pX,
                __global const VALUE_TYPE* pY)
{
    const SIZE_TYPE index = get_global_id(0);

    if (index >= size) return;

    pR[index] = OPERATION(pX[index], pY[index]);

}

#undef OPERATION
#undef SIZE_TYPE
#undef WG_SIZE
)";


// if using offsets (REPLACE_OFFSET) should be replaced by "+ offsetval"
// otherwise replace it with ""
// pAlphaOffset is kept to allow for alpha to use a vector rather than a single element
// to aid copies. (may need to remove const is this same buffer is not read_only)

// OPERATION, WG_SIZE and SIZE_TYPE need to be defined using #define statement at top
// of kernel;
// by undefining them at the bottom, the kernels for the entire solver can be combined when generating
// the program.
const std::string BiCGStabGenerator::Axpy_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void axpy_(REPLACE_NAME_MODIFIER) (const SIZE_TYPE size,
          __global double *__restrict__ pY,
          __global const double *__restrict__ pAlpha,
         const SIZE_TYPE pAlphaOffset,
          __global const double *__restrict__ pX,
          __global const double *__restrict__ pZ)
{
    const SIZE_TYPE index = get_global_id(0);
    if (index >= size) return;
    const double alpha = *(pAlpha + pAlphaOffset);
    pY[index (REPLACE_OFFSET) ] = OPERATION( ( pZ[index (REPLACE_OFFSET) ] ), ( alpha * pX[index + (REPLACE_OFFSET) ] ) );
}
#undef SIZE_TYPE
#undef WG_SIZE
#undef OPERATION

)";

// Has an offset for Beta as well in this kernel
const std::string BiCGStabGenerator::Axpby_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void axpby_(REPLACE_NAME_MODIFIER) (const SIZE_TYPE size,
          __global double *__restrict__ pY,
          __global const double *__restrict__ pAlpha,
          const SIZE_TYPE pAlphaOffset,
          __global const double *__restrict__ pX,
          __global const double *__restrict__ pBeta,
		  const SIZE_TYPE pBetaOffset,
          __global const double *__restrict__ pZ)
{
    const SIZE_TYPE index = get_global_id(0);
    if (index >= size) return;
    const VALUE_TYPE alpha = *(pAlpha + pAlphaOffset);
    const VALUE_TYPE beta = *(pBeta + pBetaOffset);
    pY[index (REPLACE_OFFSET) ] = OPERATION((alpha * pX[index + (REPLACE_OFFSET) ]), (beta * pZ[index + (REPLACE_OFFSET) ]));
}
#undef SIZE_TYPE
#undef WG_SIZE
#undef OPERATION
)";

// Modified axpby kernel allows for customizations beyond basic Axpby kernel
// Set kernel parameters are size, output buffer, Scalar array (pAlpha) and pAlphaOffset
// remaining kernel parameters and operation (both LHS and RHS are needed in replacement string
// name also has a modifier, so multiple axpby kernels can be compiled in same program source
// may can either add offsets for input/output arrays, or use subbuffers shifted equivalent distance of offset
// sub buffers might save a slight amount of time by reducing operations (iteger ops, so minimal)
// sub-buffers will not work for scalar arrays unless a significant amount of padding was used around each scalar
const std::string BiCGStabGenerator::Modified_Axpby_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void axpby_(REPLACE_NAME_MODIFIER) (const SIZE_TYPE size,
          __global double *__restrict__ pY,
          __global const double *__restrict__ pAlpha,
          const SIZE_TYPE pAlphaOffset,
          (REPLACE_INPUTS) )
{
    const SIZE_TYPE index = get_global_id(0);
    if (index >= size) return;
    const VALUE_TYPE alpha = *(pAlpha + pAlphaOffset);
    (REPLACE_OPERATION) 
}
#undef SIZE_TYPE
#undef WG_SIZE
)";




// Serves as the base kernel for generating various reduces including summation, 
// norms (summation of abs), dot products, etc. 
// (REPLACE_INPUTS) will need to be replaced by __global double ... (1 for reduce, 2 for dot)
// (REPLACE_OPERATION) needs to be replaced by the operation. 
// Make sure that the names in REPLACE_INPUTS matches those in REPLACE_OPERATION
// (REPLACE_OPERATION should include both input[stride] and input[stride+1]
// i.e. for reduce summation REPLACE_OPERATION will be input[stride] + input[stride + 1]
// (REPLACE REDUCE_TYPE) replace with type of reduce. i.e. SUM, or ABS, etc.
// Can use numbering if doing multiple variations of same type of kernel
const std::string BiCGStabGenerator::ReduceBase_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernelBase_(REPLACE_REDUCE_TYPE) ( 
(REPLACE_INPUTS),
__global double *__restrict__ output,
__local double *sdata)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE bid = get_group_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;

	sdata[tid] = (REPLACE_OPERATION);

	barrier(CLK_LOCAL_MEM_FENCE);
	for (SIZE_TYPE s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdata[tid] += sdata[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0) { output[bid] = sdata[0]; }
}

#undef SIZE_TYPE
#undef WG_SIZE
)";



// Similar to reduce base kernel, but number of inputs, outputs (and therefore locals as well)
// can be set. 
// Because of this, there are multiple REPLACE_OPERATION sections which need to include 
// full expressions rather than just the RHS since there can be multiple ones.
// These REPLACE_OPERATION include REPLACE_OPERATION_INPUT (similar to REPLACE_OPERATION in base kernel, 
// REPLACE_OPERATION_LOCAL and REPLACE_OPERATION_OUTPUT
const std::string BiCGStabGenerator::ReduceFinal_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernel_(REPLACE_REDUCE_TYPE) (
(REPLACE_INPUTS),
(REPLACE_OUTPUTS),
(REPLACE_LOCALS),
(REPLACE_OFFSETS))
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;
	(REPLACE_OPERATIONS_INPUT);
	barrier(CLK_LOCAL_MEM_FENCE);
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			(REPLACE_OPERATIONS_LOCAL);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		(REPLACE_OPERATIONS_OUTPUT);
	}
}
#undef SIZE_TYPE
#undef WG_SIZE
)";


const std::string BiCGStabGenerator::ReduceFinalOmega_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernel_Calc_Omega (
__global const double *__restrict__ inputSdT,
__global const double *__restrict__ inputTdT,
__global double *__restrict__ scalarArray,
__local double *sdataSdT,
__local double *sdataTdT,
const SIZE_TYPE offsetOmega)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;
	sdataSdT[tid] = inputSdT[stride] + inputSdT[stride+1];
	sdataTdT[tid] = inputTdT[stride] + inputTdT[stride+1];
	barrier(CLK_LOCAL_MEM_FENCE);
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdataSdT[tid] += sdataSdT[(tid + s)];
			sdataTdT[tid] += sdataTdT[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		scalarArray[offsetOmega] = sdataSdT[0]/sdataTdT[0];
	}
}
#undef SIZE_TYPE
#undef WG_SIZE
)";

const std::string BiCGStabGenerator::ReduceIntermediateOmega_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernel_Intermediate_Omega (
__global const double *__restrict__ inputSdT,
__global const double *__restrict__ inputTdT,
__global double *__restrict__ outputSdT,
__global double *__restrict__ outputTdT,
__local double *sdataSdT,
__local double *sdataTdT)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE bid = get_group_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;
	sdataSdT[tid] = inputSdT[stride] + inputSdT[stride+1];
	sdataTdT[tid] = inputTdT[stride] + inputTdT[stride+1];
	barrier(CLK_LOCAL_MEM_FENCE);
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdataSdT[tid] += sdataSdT[(tid + s)];
			sdataTdT[tid] += sdataTdT[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		outputSdT[bid] = sdataSdT[0];
		outputTdT[bid] = sdataTdT[0];
	}
}
#undef SIZE_TYPE
#undef WG_SIZE
)";


const std::string BiCGStabGenerator::ReduceBaseOmega_kernel = R"(
__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernel_Base_Omega (
__global const double *__restrict__ inputS,
__global const double *__restrict__ inputT,
__global double *__restrict__ outputSdT,
__global double *__restrict__ outputTdT,
__local double *sdataSdT,
__local double *sdataTdT)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE bid = get_group_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;
	sdataSdT[tid] = (inputS[stride]*inputT[stride]) + (inputS[stride+1]*inputT[stride+1]);
	sdataTdT[tid] = (inputT[stride]*inputT[stride]) + (inputT[stride+1]*inputT[stride+1]);
	barrier(CLK_LOCAL_MEM_FENCE);
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			sdataSdT[tid] += sdataSdT[(tid + s)];
			sdataTdT[tid] += sdataTdT[(tid + s)];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		outputSdT[bid] = sdataSdT[0];
		outputTdT[bid] = sdataTdT[0];
	}
}
#undef SIZE_TYPE
#undef WG_SIZE
)";




