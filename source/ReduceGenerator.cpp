#include "ReduceGenerator.h"

bool ReduceGenerator::hostMemFlag = false;
bool ReduceGenerator::deviceMemFlag = false;
int ReduceGenerator::domainRedSize = 0;
int ReduceGenerator::domainRedSteps = 0;
int ReduceGenerator::domainFinalRedSize = 0;
cl_mem* ReduceGenerator::redBuf1 = nullptr;
cl_mem* ReduceGenerator::redBuf2 = nullptr;
int* ReduceGenerator::domainGlobalSizes = nullptr;
ReduceGenerator*  ReduceGenerator::r_instance = nullptr;

// calculates number of reduce steps, sets num_red_steps, and returns
// padded array size
int ReduceGenerator::getNumSteps(int arrsize_, int &num_red_steps)
{
	// get size array needs to be for reductions (padded size)
	int fsize = getPaddedSize(arrsize_);

	// find out how many reduce steps it will take;
	double global_size_cur = (double)(fsize / 2);
	double Num_blocks_cur = global_size_cur / (double)WORKGROUPSIZE_RED;

	num_red_steps = 1;
	while (Num_blocks_cur > 1)
	{
		global_size_cur = Num_blocks_cur / 2.;
		Num_blocks_cur = global_size_cur / (double)WORKGROUPSIZE_RED;
		num_red_steps++;
	}
	return fsize;
}


int ReduceGenerator::getPaddedSize(const int arrsize_)
{
	// get size array needs to be for reductions (padded size)
	int fsize = 2;
	while (fsize < arrsize_)
		fsize *= 2;

	return fsize;
}


// TODO: check if this is correct (use a fully filled array, since a padded array
// might give a false correct result)

// not sure why i used doubles when calculating this, but it works, so why change it
int ReduceGenerator::getReduceSizes(const int arrsize_, int &num_red_steps, int **domfsize_,
	int &intermedSize1_, int& intermedSize2_)
{
	int fsize = getNumSteps(arrsize_, num_red_steps);

	// Initialize array storing global size for each step
	int *globSizes = new int[num_red_steps];


	// repeat previous step, this time filling newly initialized array.
	double global_size_cur = (double)(fsize / 2);
	double Num_blocks_cur = global_size_cur / (double)WORKGROUPSIZE_RED;
	intermedSize1_ = (int)Num_blocks_cur;


	for (int i = 0; i < num_red_steps-1; i++)
	{
		globSizes[i] = (int)global_size_cur;
		global_size_cur = Num_blocks_cur / 2.;
		Num_blocks_cur = global_size_cur / (double)WORKGROUPSIZE_RED;
		if (i == 0)
			intermedSize2_ = (int)Num_blocks_cur;
	}
	if (intermedSize2_ == 0)
		intermedSize2_ = intermedSize1_;
	globSizes[num_red_steps - 1] = (int)global_size_cur;
	*domfsize_ = globSizes;
	return fsize;
}

void ReduceGenerator::ini(int xsize_, int ysize_)
{
	if (deviceMemFlag && hostMemFlag)
		return;

	int intermediateSize1 = 0, intermediateSize2 = 0;
	domainRedSize = iniSizesAndMemory(xsize_*ysize_, domainRedSteps,
		&domainGlobalSizes, &redBuf1, intermediateSize1, intermediateSize2);
	iniDeviceMemory(domainRedSteps, &redBuf2, intermediateSize1, intermediateSize2);


	deviceMemFlag = true;
	hostMemFlag = true;
}

int ReduceGenerator::iniSizesAndMemory(const int arrSize_, int &num_red_steps, 
	int **globSizes_, cl_mem** red_buf_, int &intermedSize1_, int& intermedSize2_)
{

	int fsize = getReduceSizes(arrSize_, num_red_steps, globSizes_,
		intermedSize1_, intermedSize2_);
	iniDeviceMemory(num_red_steps, red_buf_, intermedSize1_, intermedSize2_);
	return fsize;
}

void ReduceGenerator::iniDeviceMemory(const int num_red_steps, cl_mem** red_buf_,
	int intermediateSize1, int intermediateSize2)
{
	cl_mem* redBufTemp = new cl_mem[2];
	int status;
	double zer = 0.;
	redBufTemp[0] = clCreateBuffer(*context, CL_MEM_READ_WRITE,
		sizeof(double) * intermediateSize1, NULL, &status);
	CHECK_ERROR_WITH_EXIT(status, 0, "Error creating intermediate buffer 1");
		
	status = clEnqueueFillBuffer(*ioQue, redBufTemp[0], &zer,
		sizeof(double), 0, sizeof(double)*intermediateSize1, 0, NULL, NULL);
	CHECK_ERROR_WITH_EXIT(status, 0, "Error filling intermediate buffer 1");

	redBufTemp[1] = clCreateBuffer(*context, CL_MEM_READ_WRITE,
		sizeof(double) * intermediateSize2, NULL, &status);
	CHECK_ERROR_WITH_EXIT(status, 0, "Error creating intermediate buffer number 2");

	status = clEnqueueFillBuffer(*ioQue, redBufTemp[1], &zer,
		sizeof(double), 0, sizeof(double)*intermediateSize2, 0, NULL, NULL);
	CHECK_ERROR_WITH_EXIT(status, 0, "Error filling intermediate buffer 2");

	
	*red_buf_ = redBufTemp;
}

std::string ReduceGenerator::getGenericName(reduceType redtype_)
{
	std::string namestr;
	if (redtype_ == Sum) {
		namestr = "Sum";
	}
	else if (redtype_ & MMAS) {
		std::string OPER1 = (redtype_ & MaxAS) ? "fmax" : "fmin";
		std::string OPER2 = (redtype_ & AbsMM) ? "fabs" : "";
		namestr = OPER2 + "_" + OPER1;
	}
	else if (redtype_ == Dot) {
		namestr = "Dot";
	}
	else if (redtype_ & Norm1) {
		namestr = "Norm";
	}
	else if (redtype_ == Sqr) {
		namestr = "Sqr";
	}
	else if (redtype_ == SumNType) {
		namestr = "SumNType";
	}
	else {
		namestr = "";
	}
	return namestr;
}

std::string ReduceGenerator::getFTYPE(varType vtype_)
{
	switch (vtype_)
	{
	case DOUBLE_T:
	{
		return "double";
	}
	case BOOL_T:
	{
		return "bool";
	}
	case CHAR_T:
	{
		return "char";
	}
	case INT_T:
	{
		return "int";
	}
	case UINT_T:
	{
		return "uint";
	}
	}

	return "";
}

size_t ReduceGenerator::getVarTypeSize(varType vtype_)
{
	switch (vtype_)
	{
	case DOUBLE_T:
	{
		return sizeof(double);
	}
	case BOOL_T:
	{
		return sizeof(bool);
	}
	case CHAR_T:
	{
		return sizeof(char);
	}
	case INT_T:
	{
		return sizeof(int);
	}
	case UINT_T:
	{
		return sizeof(unsigned int);
	}
	}

	return 0;
}


ReduceGenerator::hashResult ReduceGenerator::addGenericReduce(reduceType redtype_, std::string &kername, varType vtype_)
{
	std::string redker;
	addDefine(redker, "WG_SIZE", WORKGROUPSIZE_RED);
	addDefine(redker, "FTYPE", getFTYPE(vtype_));
	redker.append(ReduceBase_kernel);

	std::string inputstr = genericInputStr;
	std::string operstr;
	std::string localoperstr = genericForOp;
	std::string namestr;
	if (redtype_ == Sum) {
		operstr = "input[stride] + input[stride+1];";
	}
	else if (redtype_ & MMAS) {
		std::string OPER1 = (redtype_ & MaxAS) ? "fmax" : "fmin";
		std::string OPER2 = (redtype_ & AbsMM) ? "fabs" : "";
		operstr = OPER1 + "(" + OPER2 + "(input[stride]), " + OPER2 + "(input[stride + 1]));";
		localoperstr = "sdata[tid] = " + OPER1 + "(sdata[tid], sdata[(tid + s)]);";
	}
	else if (redtype_ == Dot) {
		inputstr = "__global FTYPE *__restrict__ inputA,\n"\
			"__global FTYPE *__restrict__ inputB,";
		operstr = "inputA[stride] * inputB[stride] + inputA[stride+1] * inputB[stride+1];";
	}
	else if (redtype_ & Norm1) {
		operstr = "fabs(input[stride]) + fabs(input[stride+1]);";
	}
	else if (redtype_ == Sqr) {
		operstr = "(input[stride]*input[stride]) + (input[stride+1]*input[stride+1]);";
	}
	else if (redtype_ == SumNType)
	{
		// if not using in_kernel_ibb, the input is a short and output is an int,
		// so need to specifically set input as a short instread of an ftype.
#ifndef IN_KERNEL_IBB
		inputstr = "__global short *__restrict__ input,";
#endif
		operstr = "((input[stride] & M_FLUID_NODE) ? 1 : 0) + ((input[stride+1] & M_FLUID_NODE) ? 1 : 0);";
	}
	else {
		Gen_Error_Msg(-1593, "tried to create a BaseRed kernel with incorrect type");
	}


	namestr = getGenericName(redtype_) + "_Size" + std::to_string(WORKGROUPSIZE_RED) + "_" + getFTYPE(vtype_);
	findAndReplace(redker, "(REPLACE_INPUTS)", inputstr);
	findAndReplace(redker, "(REPLACE_OPERATION)", operstr);
	findAndReplace(redker, "(REPLACE_REDUCE_TYPE)", namestr);
	findAndReplace(redker, "(REPLACE_OPERATIONS_LOCAL)", localoperstr);

	kername = "ReduceKernelBase_" + namestr;

	return checkHash(redker, kername);
}

ReduceGenerator::hashResult ReduceGenerator::addGenericFinalRedKernel(int wgsize, reduceType redtype_, std::string &kername, varType vtype_)
{
	std::string redker;
	addDefine(redker, "WG_SIZE", wgsize);
	addDefine(redker, "FTYPE", getFTYPE(vtype_));
	redker.append(ReduceFinal_kernel);

	std::string nameAppend = getGenericName(redtype_) + "_Size" + std::to_string(wgsize) + "_" + getFTYPE(vtype_);

	findAndReplace(redker, "(REPLACE_REDUCE_TYPE)", nameAppend);
	findAndReplace(redker, "(REPLACE_INPUTS)", genericInputStr);
	findAndReplace(redker, "(REPLACE_OUTPUTS)", genericOutputStr);
	findAndReplace(redker, "(REPLACE_LOCALS)", genericLocalStr);
	findAndReplace(redker, "(REPLACE_OFFSETS)", genericOffsetStr);
	findAndReplace(redker, "(REPLACE_OPERATIONS_OUTPUT)", "output[Offset] = sdata[0];");

	if (redtype_ & MMAS)
	{
		std::string OPER1 = (redtype_ & MaxAS) ? "fmax" : "fmin";
		findAndReplace(redker, "(REPLACE_OPERATIONS_INPUT)", "sdata[tid] = " + OPER1 + "(" + "input[stride], " + "input[stride + 1]);");
		findAndReplace(redker, "(REPLACE_OPERATIONS_LOCAL)", "sdata[tid] = " + OPER1 + "(sdata[tid], sdata[(tid + s)]);");
	}
	else
	{
		findAndReplace(redker, "(REPLACE_OPERATIONS_INPUT)", "sdata[tid] = input[stride] + input[stride+1];");
		findAndReplace(redker, "(REPLACE_OPERATIONS_LOCAL)", "sdata[tid] += sdata[(tid + s)];");
	}

	kername = "ReduceKernelFinal_" + nameAppend;
	return checkHash(redker, kername);
}

// Could probably be cleaner to template this function instead of passing varType, but eventually need to test what the
// template is to get set define FTYPE and size of local memory. 
void ReduceGenerator::genGenericDomainReduceSet(varType vType, reduceType redtype_, cl_mem **inbuf_, cl_mem **outbuf_,
	RedKernelList *kerptr_, int usemem_)
{
	if (hostMemFlag == false)
	{
		Gen_Error_Msg(-56545, "Domain size must be initialized in"
			"reduce generator before generating source");
	}

	cl_mem **memarr;

	if (usemem_ == 0)
		memarr = &redBuf1;
	else
		memarr = &redBuf2;


	genGenericReduceSet(vType, domainRedSize, redtype_, inbuf_, outbuf_, kerptr_, domainRedSteps, memarr, true);
}


int ReduceGenerator::genGenericReduceSet(varType vType, int &arrSize_, reduceType redtype_, cl_mem **inbuf_, cl_mem **outbuf_,
	RedKernelList *kerptr_, int numker, cl_mem **usemem_temp, bool domainFl)
{
	std::string kername_;
	int *globSizes = nullptr;
	cl_mem* usemem_;
	if (domainFl)
	{// reducing an array thats the same size as domain
		globSizes = domainGlobalSizes;
		usemem_ = *usemem_temp;
	}
	else
	{// reducing an array of any length
	
		// Placeholders to pass into iniSizesAndMemory
		int intermediateSize1 = 0, intermediateSize2 = 0;

		// function will calculate the paddedsize, number of steps, 
		// workgroup sizes, and initialize buffers. The memory in globSizes
		// needs to be freed at the end of this function (only if initialized here,
		// size we dont want to free domainGlobalSizes if domainFl == true).
		// buffers and array holding them will need to be freed by Reducer class.
		int paddedSize = iniSizesAndMemory(arrSize_, numker,
			&globSizes, &usemem_, intermediateSize1, intermediateSize2);
		arrSize_ = paddedSize; 
	}

	int finalwgsize = globSizes[numker - 1];
	kerptr_->ini(numker);

	// creating and adding reduce kernel if hasnt happened yet.
	addGenericReduce(redtype_, kername_, vType);

	// First kernel may have 1 or 2 input kernels
	int numbufs = (redtype_ == Dot) ? 3 : 2;

	cl_mem** bufs;

	bufs = new cl_mem*[numbufs];

	int curbufind = 0;
	bufs[curbufind++] = inbuf_[0];
	if (redtype_ == Dot)
		bufs[curbufind++] = inbuf_[1];

	int intermedBufOut = 0; //will become intermediateBufIn after being used as BufOut, then ^= 1 to switch

	bufs[curbufind++] = &usemem_[intermedBufOut];
		
	kerptr_->operator()(0).setSizes(globSizes[0], WORKGROUPSIZE_RED);
	buildInfo* binfo = new buildInfo(&kerptr_->operator()(0), kername_, vType, numbufs, bufs, NULL, 1, WORKGROUPSIZE_RED);
	kerToBuild.push_back(*binfo);
	
	// All reductions except for Min and Max become simple summations
	// after first step. AbsMin and AbsMax become just Min and Max
	if (redtype_ == AbsMax) { redtype_ = Max; }
	if (redtype_ == AbsMin) { redtype_ = Min; }
	if (!(redtype_ == Max || redtype_ == Min))
	{
		redtype_ = Sum;
	}

	for (int i = 1; i < numker - 1; i++)
	{
		addGenericReduce(redtype_, kername_, vType);
		curbufind = 0;
		
		// intermediate buffer in is intermediate buffer out from last step
		bufs[curbufind++] = &usemem_[intermedBufOut];
		// Switch intermiate buffer to other one to become new intermediate buffer out
		intermedBufOut ^= 1; 
		bufs[curbufind++] = &usemem_[intermedBufOut];
		
		kerptr_->operator()(i).setSizes(globSizes[i], WORKGROUPSIZE_RED);
		binfo = new buildInfo(&kerptr_->operator()(i), kername_, vType, 2, bufs, NULL, 1, WORKGROUPSIZE_RED);
		kerToBuild.push_back(*binfo);
	}

	addGenericFinalRedKernel(finalwgsize, redtype_, kername_, vType);
	curbufind = 0;
	
	// intermediate Buffer input in the last kernel is the output from previous step
	bufs[curbufind++] = &usemem_[intermedBufOut];
	bufs[curbufind++] = *outbuf_;
	kerptr_->operator()(numker - 1).setSizes(globSizes[numker - 1], finalwgsize);
	
	binfo = new buildInfo(&kerptr_->operator()(numker - 1), kername_, vType, 2, bufs, NULL, 1, finalwgsize);
	kerToBuild.push_back(*binfo);
	
	kerptr_->setValueOutputIndex(3);
	kerptr_->setValueOutputArrayIndex(1);
	delete[] bufs;

	// for non-domain sized arrays, memory deallocation will be completed
	// by the Reducer class;
	if (domainFl == false)
	{
		delete[] globSizes;
		*usemem_temp = usemem_;
	}

	return finalwgsize;
}

void ReduceGenerator::callIniKernels()
{
	std::list<buildInfo>::iterator it;
	for (it = kerToBuild.begin(); it != kerToBuild.end(); ++it)
	{
		int status;
		status = it->kerAdd->createKernel(&program, it->kerName);
		CHECK_ERROR_WITH_EXIT(status, 0, "Error creating kernel " + n.kerName);
		
		for (int i = 0; i < it->numBuffer; i++)
		{
			status = it->kerAdd->setArgument(it->bufInds[i], it->bufAdd[i]);
			CHECK_ERROR_WITH_EXIT(status, 0, "Error setting index " + std::to_string(n.bufInds[i]) +
				" to kernel " + n.kerName);
		}
		for (int i = 1; i <= it->numLocalBuf; i++)
		{
			status = it->kerAdd->setLocalMem(it->bufInds[it->numBuffer - 1] + i,
				getVarTypeSize(it->vType)*it->localBufSize);
			CHECK_ERROR_WITH_EXIT(status, 0, "Error setting local buffer at index " +
				std::to_string(n.bufInds[i]) + " to kernel " + n.kerName);
		}
	}
}

const std::string ReduceGenerator::genericInputStr = "__global const FTYPE *__restrict__ input,";
const std::string ReduceGenerator::genericOutputStr = "__global FTYPE *__restrict__ output,";
const std::string ReduceGenerator::genericForOp = "sdata[tid] += sdata[(tid + s)];";
const std::string ReduceGenerator::genericLocalStr = "__local FTYPE *sdata,";
const std::string ReduceGenerator::genericOffsetStr = "const int Offset";


// Serves as the base kernel for generating various reduces including summation, 
// norms (summation of abs), dot products, etc. 
// (REPLACE_INPUTS) will need to be replaced by __global double ... (1 for reduce, 2 for dot)
// (REPLACE_OPERATION) needs to be replaced by the operation. 
// Make sure that the names in REPLACE_INPUTS matches those in REPLACE_OPERATION
// (REPLACE_OPERATION shoul include both input[stride] and input[stride+1]
// i.e. for reduce summation REPLACE_OPERATION will be input[stride] + input[stride + 1]
// (REPLACE REDUCE_TYPE) replace with type of reduce. i.e. SUM, or ABS, etc.
// Can use numbering if doing multiple variations of same type of kernel
const std::string ReduceGenerator::ReduceBase_kernel = R"(
#define SIZE_TYPE	uint

__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernelBase_(REPLACE_REDUCE_TYPE) ( 
(REPLACE_INPUTS)
__global FTYPE *__restrict__ output,
__local FTYPE *sdata)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE bid = get_group_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;

	sdata[tid] = (REPLACE_OPERATION)

	barrier(CLK_LOCAL_MEM_FENCE);
	for (SIZE_TYPE s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			//sdata[tid] += sdata[(tid + s)];
			(REPLACE_OPERATIONS_LOCAL)
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0) { output[bid] = sdata[0]; }
}

#undef SIZE_TYPE
#undef WG_SIZE
#undef FTYPE
)";




// Similar to reduce base kernel, but number of inputs, outputs (and therefore locals as well)
// can be set. 
// Because of this, there are multiple REPLACE_OPERATION sections which need to include 
// full expressions rather than just the RHS since there can be multiple ones.
// These REPLACE_OPERATION include REPLACE_OPERATION_INPUT (similar to REPLACE_OPERATION in base kernel, 
// REPLACE_OPERATION_LOCAL and REPLACE_OPERATION_OUTPUT
const std::string ReduceGenerator::ReduceFinal_kernel = R"(
#define SIZE_TYPE	uint

__kernel
__attribute__((reqd_work_group_size(WG_SIZE,1,1)))
void ReduceKernelFinal_(REPLACE_REDUCE_TYPE) (
(REPLACE_INPUTS)
(REPLACE_OUTPUTS)
(REPLACE_LOCALS)
(REPLACE_OFFSETS)
)
{
	const SIZE_TYPE tid = get_local_id(0);
	const SIZE_TYPE gid = get_global_id(0);
	const SIZE_TYPE localSize = get_local_size(0);
	const SIZE_TYPE stride = gid * 2;

	(REPLACE_OPERATIONS_INPUT)

	barrier(CLK_LOCAL_MEM_FENCE);
	for (unsigned int s = localSize >> 1; s > 0; s >>= 1)
	{
		if (tid < s)
		{
			(REPLACE_OPERATIONS_LOCAL)
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (tid == 0)
	{
		(REPLACE_OPERATIONS_OUTPUT)
	}
}
#undef SIZE_TYPE
#undef WG_SIZE
#undef FTYPE
)";


