// Sparse Matrix class for clSPARSE library
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_REDUCEGENERATOR_H__INCLUDED_)
#define AFX_REDUCEGENERATOR_H__INCLUDED_

#pragma once

#include "SourceGenerator.h"

#ifndef CHECK_KERNEL_STATUS
#define STATUS_CHECK_RUN(codeval, msg) do {} while(0)
#else
#define STATUS_CHECK_RUN(codeval, msg)	if(codeval) { Gen_Error_Msg(codeval, msg); }
#endif

#define STATUS_CHECK_SETUP(codeval, msg)	if(codeval) { Gen_Error_Msg(codeval, msg); }



// TODO: set intermediate buffers to static and create two sets. That way, there
// will only be two sets regardless of the number of reduce class instances that
// exist. Could also add a third set of pointers and create a wrapper which generates
// its own set of buffers if a size other than the regular domain size is needed.
// Also, may want to use just a single class with both base and final kernels, and 
// a method to easily set kernel arguments. Then, a lightweight wrapper can be created
// for reducing the macroscopic kernels, and the BiCGStab kernels can contain methods for
// doing the more intense operations with complicated final reduce steps. Also, the 
// output buffer should not be stored in the reduce class (wrapper class would store it for 
// macroscopic variables, and BiCGStab would store it for FD methods.
// 
// Have a reduce kernel generator class which generates the reduce kernels and stores them in a single
// cl_program. During initialization of code, all arrays needing to be reduced will 
// call generator class to create kernels. Then cl_program is compiled with only kernels needed 
// generator class will track kernels to create once program is compiled, and fill in set arguments for
// local mem, as well as other memories by using pointers to cl_mem stored elsewhere

/// generator class will contain a function that takes addresses of kernels, and buffers, and kernel names
/// these addresses are stored in arrays of pointers. once program is compiled, these arrays of pointers
// are used to initialize the kernels without subsequent calls from classes that will be using the kernels.
// arrays of pointers can then be freed. THe generator class can be instantiated as a static object in reduce class
// so that there is only one instance.



struct buildInfo;






class ReduceGenerator : public sourceGenerator
{
	static ReduceGenerator *r_instance;
	ReduceGenerator(const std::string type_, const clEnv::programType ptype_)
	{
		type = type_;
		pType = ptype_;
	}

	~ReduceGenerator()
	{
		if (deviceMemFlag)
		{
			for (int i = 0; i < domainRedSteps - 1; i++)
			{
				clReleaseMemObject(redBuf1[i]);
				clReleaseMemObject(redBuf1[i]);
			}
		}
		if (hostMemFlag)
		{
			delete[] redBuf1;
			delete[] redBuf2;
		}
	}
public:
	friend class sourceGenerator;
	friend class BiCGStabGenerator;

	static ReduceGenerator *ReduceInstance()
	{
		if (!r_instance)
			r_instance = new ReduceGenerator("ReduceGenerator", clEnv::ReduceType);
		return r_instance;
	}
	
	// char_t is for legacy code (vls.M is an array of chars)
	// NULL_T is for dummy functions (when using template specialization) in ReducerBase 
	enum varType {DOUBLE_T, BOOL_T, CHAR_T, INT_T, UINT_T, NULL_T};

	enum reduceType {
		Sum = 0x1, Max = 0x2, Min = 0x4, Abs = 0x8, 
		Sqr = 0x10, Norm = 0x20, Dot = 0x40, AbsMin = 0x80, AbsMax = 0x100, SumNType = 0x200,
		AbsMM = AbsMin | AbsMax, MinAS = AbsMin | Min, MaxAS = AbsMax | Max,
		SignMM = Min | Max, MMAS = AbsMin | AbsMax | Min | Max, Norm1 = Abs | Norm };

	hashResult addGenericReduce(reduceType redtype_, std::string &kername, varType vtype_);

	static std::string getFTYPE(varType vtype_);
	static size_t getVarTypeSize(varType vtype_);
	static std::string getGenericName(reduceType redtype_);

	void genGenericDomainReduceSet(varType vType, reduceType redtype_, cl_mem **inbuf_, cl_mem **outbuf_,
		RedKernelList *kerptr_, int usemem_ = 0);

	hashResult addGenericFinalRedKernel(int wgsize, reduceType redtype_, std::string &kername, varType vtype_);

	// Returns the size of the work group size of the final reduce kernel
	int genGenericReduceSet(varType vType, int &arrSize_, reduceType redtype_, cl_mem **inbuf_, cl_mem **outbuf_,
		RedKernelList *kerptr_, int numker, cl_mem **usemem_, bool domainFl = false);

	int getNumSteps(int arrsize_, int &num_red_steps);
	
	int getPaddedSize(const int arrsize_);
	
	void callIniKernels();

	void ini(int xsize_, int ysize_);

	int getFinalWGSize() { return domainGlobalSizes[domainRedSteps - 1]; }

	static int domainRedSize;
	static int domainRedSteps;
	static int* domainGlobalSizes;
	static bool hostMemFlag;
	static bool deviceMemFlag;
	static int domainFinalRedSize;


	// These allocate memory, which the caller must free when finished with
	int getReduceSizes(const int arrsize_, int &num_red_steps, int **domfsize_,
		int &intermedSize1_, int& intermedSize2_);
	
	int iniSizesAndMemory(const int arrSize_, int &num_red_steps, int **globSizes_,
		cl_mem** red_buf_, int &intermedSize1_, int& intermedSize2_);

	void iniDeviceMemory(const int num_red_steps, cl_mem** red_buf_, 
		int intermedSize1_, int intermedSize2_);

protected:


	ReduceGenerator(ReduceGenerator const&){};             // copy constructor is private
	ReduceGenerator& operator=(ReduceGenerator const&){};  // assignment operator is private
	std::list<buildInfo> kerToBuild;
	const static std::string ReduceBase_kernel;
	const static std::string ReduceFinal_kernel;
	const static std::string genericInputStr;
	const static std::string genericOutputStr;
	const static std::string wgSizeDefStr;
	const static std::string genericForOp;
	const static std::string genericLocalStr;
	const static std::string genericOffsetStr;
	static cl_mem* redBuf1;
	static cl_mem* redBuf2; //size 2*(domainRedSteps-1)
};


//// TODO: Use only two intermediate reduce buffers (alternate between them. No need for one per step)
//// TODO: delete all unnecessary objects after opencl program is initialized and built
//// TODO: Adjust size of kernel workgroups to eliminate very small final kernel sizes
struct buildInfo
{
	thinKerWrapper *kerAdd;
	std::string kerName;
	cl_mem **bufAdd;
	int numBuffer;
	int *bufInds;
	int numLocalBuf;
	int localBufSize;
	ReduceGenerator::varType vType;

	buildInfo() : numBuffer(0)
	{
		bufInds = nullptr;
		kerAdd = nullptr;
		kerName = "";
		bufAdd = nullptr;
		numLocalBuf = 0;
		localBufSize = 0;
		vType = ReduceGenerator::NULL_T;
	}

	buildInfo(thinKerWrapper *keradd_, std::string &kername_,
		ReduceGenerator::varType vtype_, int numbuffer_,
		cl_mem *buffers_[], int bufferind_[] = NULL, int numlocalbuf_ = 0,
		int localbufsize_ = 0)
		: kerAdd(keradd_), numBuffer(numbuffer_), kerName(kername_),
		numLocalBuf(numlocalbuf_), localBufSize(localbufsize_), vType(vtype_)
	{
		bufInds = new int[numBuffer];
		bufAdd = new cl_mem*[numBuffer];
		for (int i = 0; i < numBuffer; i++)
		{
			if (bufferind_ == NULL)
				bufInds[i] = i;
			bufAdd[i] = buffers_[i];
		}
	}

	~buildInfo()
	{
		delete[] bufAdd;
		delete[] bufInds;
	}

};
//
//template <class Container>
//class BaseReduceKernel
//{
//public:
//	typedef typename Container::type	T;
//
//	//Memory storing reduce results
//	cl_kernel *ker;		// pointers to reduce kernels
//	cl_mem *Buffers;	// Buffers storing intermediate values
//
//	size_t *local_size, *global_size; // pointers to work sizes for each kernel
//	int num_red_steps;
//
//	cl_command_queue *queue;
//	std::string Name; // For printing errors
//	std::string kerName;
//	int FullSize;
//
//	Container *ArrToRed;
//
//	BaseReduceKernel()
//	{
//		queue = NULL;
//		local_size = NULL;
//		global_size = NULL;
//		Buffers = NULL;
//		ker = NULL;
//		FullSize = 0;
//		num_red_steps = 0;
//	}
//
//	~BaseReduceKernel()
//	{
//		free_memory();
//	}
//
//
//	void Gen_Error_Msg(int Errcode, std::string mesg)
//	{
//		FILE *stream;
//		stream = fopen("error.txt", "w+");
//		fprintf(stream, "Error in BaseReduceKernel class: Code: %d, Kernel Name: %s, Array Name: %s\nMessage: %s\n", Errcode, kerName.c_str(), Name.c_str(), mesg.c_str());
//		fclose(stream);
//		printf("Error in BaseReduceKernel class: Code: %d, Kernel Name: %s, Array Name: %s\nMessage: %s\n", Errcode, kerName.c_str(), Name.c_str(), mesg.c_str());
//		exit(Errcode);
//	}
//
//	void setName(std::string redtype)
//	{
//		if (redtype.compare("Sum") == 0)
//		{
//			kerName = "reduce_generic";
//		}
//		else if (redtype.compare("Max") == 0)
//		{
//			kerName = "reduce_generic_max";
//		}
//		else if (redtype.compare("Min") == 0)
//		{
//			kerName = "reduce_generic_min";
//		}
//		else if (redtype.compare("AbsMax") == 0)
//		{
//			kerName = "reduce_generic_absmax";
//		}
//	}
//
//	void free_memory()
//	{
//		FREE(local_size);
//		FREE(global_size);
//		if (Buffers != NULL)
//		{
//			for (int i = 0; i < num_red_steps - 1; i++)
//			{
//				if (Buffers[i] != NULL)
//					clReleaseMemObject(Buffers[i]);
//			}
//			FREE(Buffers);
//		}
//		if (ker != NULL)
//		{
//			for (int i = 0; i < num_red_steps - 1; i++)
//			{
//				if (ker[i] != NULL)
//					clReleaseKernel(ker[i]);
//			}
//			FREE(ker);
//		}
//	}
//
//
//	void ini(cl_command_queue *queue_, int groupSize, Container *arrtored_, cl_program *program, std::string redtype)
//	{
//		ArrToRed = arrtored_;
//		setName(redtype);
//		queue = queue_;
//
//		Name = ArrToRed->getName();
//		Name.append("_Red");
//
//		int size = ArrToRed->getBufferFullSize();
//
//		int status;
//
//		FullSize = 2;
//		while (FullSize < size)
//			FullSize *= 2;
//
//		STATUS_CHECK_SETUP((FullSize != size), "Size of Array must be a power of 2");
//
//
//		double global_size_cur = (double)(FullSize / 2);
//		double Num_blocks_cur = global_size_cur / (double)groupSize;
//
//		num_red_steps = 1;
//		while (Num_blocks_cur > 1)
//		{
//			global_size_cur = Num_blocks_cur / 2.;
//			Num_blocks_cur = global_size_cur / (double)groupSize;
//			num_red_steps++;
//		}
//
//		global_size = (size_t*)calloc(num_red_steps, sizeof(size_t));
//		local_size = (size_t*)calloc(num_red_steps, sizeof(size_t));
//		Buffers = (cl_mem*)malloc((num_red_steps - 1) * sizeof(cl_mem));
//		ker = (cl_kernel*)malloc((num_red_steps - 1)*sizeof(cl_kernel));
//		size_t globsize_cur = FullSize / 2;
//		size_t numblock_cur = globsize_cur / groupSize;
//
//		for (int ii = 0; ii < num_red_steps - 1; ii++)
//		{
//			/////////////Set Kernel Work Sizes//////////////			
//			global_size[ii] = globsize_cur;
//			local_size[ii] = groupSize;
//
//			/////////////Create Kernel//////////////			
//			ker[ii] = clCreateKernel(*program, kerName.c_str(), &status);
//			STATUS_CHECK_SETUP(status, "Error Creating Kernel");
//
//			/////////////Create Intermediate Buffer//////////////	
//			Buffers[ii] = clCreateBuffer(*ArrToRed->getContext(), CL_MEM_READ_WRITE,
//				sizeof(T)* numblock_cur, NULL, &status);
//			STATUS_CHECK_SETUP(status, "Error Creating Temporary Reduce Buffers");
//
//			////////////////Set Kernel Arguments///////////////////			
//			//First kernel is given array to reduce, rest get intermediate buffers
//			if (ii > 0)
//			{
//				status = clSetKernelArg(ker[ii], 0, sizeof(cl_mem), (void *)&Buffers[ii - 1]);
//				STATUS_CHECK_SETUP(status, "Error setting first argument");
//			}
//			status = clSetKernelArg(ker[ii], 1, sizeof(cl_mem), (void *)&Buffers[ii]);
//			STATUS_CHECK_SETUP(status, "Error setting intermediate reduce buffer argument");
//
//			status = clSetKernelArg(ker[ii], 2, local_size[ii] * sizeof(T), NULL);
//			STATUS_CHECK_SETUP(status, "Error setting local buffer");
//
//			/////////////////////// Update WorkSizes for Next Iteration///////////////
//			globsize_cur = numblock_cur / 2;
//			numblock_cur = globsize_cur / groupSize;
//		}
//
//		global_size[num_red_steps - 1] = globsize_cur;
//		local_size[num_red_steps - 1] = globsize_cur;
//	}
//
//	void call_kernels(int num_wait, cl_event *wait, cl_command_queue *que_)
//	{
//		int status = clSetKernelArg(ker[0], 0, sizeof(cl_mem), ArrToRed->get_add_Macro());
//		STATUS_CHECK_RUN(status, "Error Calling Reduce Kernel Number 0");
//
//		status = clEnqueueNDRangeKernel(*que_, ker[0], 1, 0, &global_size[0], &local_size[0], num_wait, wait, NULL);
//		STATUS_CHECK_RUN(status, "Error Calling Reduce Kernel Number 0");
//
//
//		for (int i = 1; i < num_red_steps - 1; i++)
//		{
//			status = clEnqueueNDRangeKernel(*que_, ker[i], 1, 0, &global_size[i], &local_size[i], 0, NULL, NULL);
//			STATUS_CHECK_RUN(status, "Error Calling Base Reduce Kernel");
//		}
//	}
//};













//
//
//
//
//template <typename Container>
//class GenericRedKernelSingle
//{
//public:
//	typedef typename Container::type T;
//
//	Array1D<T> ReduceOneVal; //stores results of reduction
//
//	BaseReduceKernel<Container> BaseRed;
//
//	cl_command_queue *queue;
//
//	std::string Name;
//	std::string kerName, kerNameBase;
//	std::string className;
//
//	cl_kernel ker;
//
//
//
//	GenericRedKernelSingle() : ReduceOneVal(1)
//	{
//		ReduceOneVal.fill((T)0);
//	}
//
//	~GenericRedKernelSingle()
//	{
//		clReleaseKernel(ker);
//	}
//
//	virtual void setNames(std::string redtype)
//	{
//		if (redtype.compare("Sum") == 0)
//		{
//			kerName = "reduce_generic_2_updated";
//		}
//		else if (redtype.compare("Max") == 0)
//		{
//			kerName = "reduce_generic_max_2_updated";
//		}
//		else if (redtype.compare("Min") == 0)
//		{
//			kerName = "reduce_generic_min_2_updated";
//		}
//		else if (redtype.compare("AbsMax") == 0)
//		{
//			kerName = "reduce_generic_absmax_2_updated";
//		}
//		else
//		{
//			printf("%s not a valid reduce type. Choices are Sum, Max, Min, and AbsMax\n", redtype.c_str());
//			exit(-144194);
//		}
//		className = "GenericRedKernelSingle";
//	}
//
//	void Gen_Error_Msg(int Errcode, std::string mesg)
//	{
//		FILE *stream;
//		stream = fopen("error.txt", "w+");
//		fprintf(stream, "Error in %s class: Code: %d, Kernel Name: %s, Array Name: %s\nMessage: %s\n", className.c_str(), Errcode, kerName.c_str(), Name.c_str(), mesg.c_str());
//		fclose(stream);
//		printf("Error in GenericRedKernel class: Code: %d, Kernel Name: %s, Array Name: %s\nMessage: %s\n", Errcode, kerName.c_str(), Name.c_str(), mesg.c_str());
//		exit(Errcode);
//	}
//
//
//	void ini(Container* arrtored_, std::string redtype, int groupSize, cl_program *program, cl_command_queue *queue_ = NULL)
//	{
//		setNames(redtype);
//		if (*queue_ == NULL)
//			queue = arrtored_->getQueue();
//		else
//			queue = queue_;
//
//		Name = arrtored_->getName();
//		Name.append("_Red");
//
//		ReduceOneVal.allocate_buffer_w_copy(*arrtored_->getContext(), *queue);
//
//
//		BaseRed.ini(queue, groupSize, arrtored_, program, redtype);
//
//		int status;
//		ker = clCreateKernel(*program, kerName.c_str(), &status);
//		STATUS_CHECK_SETUP(status, "Error Creating Final Generic Reduction Kernel");
//
//		set_kernel_arguments();
//	}
//
//	virtual void set_kernel_arguments()
//	{
//		int status = clSetKernelArg(ker, 0, sizeof(cl_mem), (void *)&BaseRed.Buffers[BaseRed.num_red_steps - 2]);
//		STATUS_CHECK_SETUP(status, "Error setting first argument");
//
//		status = clSetKernelArg(ker, 1, sizeof(cl_mem), ReduceOneVal.get_buf_add());
//		STATUS_CHECK_SETUP(status, "Error setting output buffer argument");
//
//		status = clSetKernelArg(ker, 2, BaseRed.local_size[BaseRed.num_red_steps - 1] * sizeof(T), NULL);
//		STATUS_CHECK_SETUP(status, "Error setting local buffer");
//	}
//
//
//	virtual T& operator()(const int i)
//	{// Warning, this is a convenience function to get values after they have already been read from device memory
//		return ReduceOneVal(i);
//	}
//
//	virtual bool reduce(int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
//	{
//		Gen_Error_Msg(-1000, "Trying to call reduce() with single output reduce kernel. Use reduceSingle() instead");
//		return false;
//	}
//
//	virtual bool reduce(cl_command_queue *que_, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
//	{
//		Gen_Error_Msg(-1000, "Trying to call reduce() with single output reduce kernel. Use reduceSingle() instead");
//		return false;
//	}
//
//	virtual T reduceSingle(int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
//	{
//		return reduceSingle(queue, num_wait, wait, evt);
//	}
//
//	virtual T reduceSingle(cl_command_queue *que_, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
//	{
//		BaseRed.call_kernels(num_wait, wait, que_);
//
//		int status = clEnqueueNDRangeKernel(*que_, ker, 1, 0, &BaseRed.global_size[BaseRed.num_red_steps - 1],
//			&BaseRed.local_size[BaseRed.num_red_steps - 1], 0, NULL, evt);
//		STATUS_CHECK_RUN(status, "Error calling final reduce kernel save index");
//
//		ReduceOneVal.read_from_buffer_size(*que_, CL_TRUE, 1);
//		return ReduceOneVal(0);
//	}
//
//};
//
//
//
//
//
//template <int size_, typename Container>
//class GenericRedKernel : public GenericRedKernelSingle<Container>
//{
//public:
//	typedef typename Container::type T;
//
//	Array1D<T> ReduceVals; //stores results of reduction
//	int Len = size_;
//	int cur_ind;
//	int next_ind;
//
//	GenericRedKernel() : ReduceVals(size_), cur_ind(0), next_ind(0)
//	{
//		ReduceVals.fill((T)0);
//	}
//
//	~GenericRedKernel()
//	{
//	}
//
//	void setNames(std::string redtype)
//	{
//		if (redtype.compare("Sum") == 0)
//		{
//			kerName = "reduce_generic_2_updated_multi";
//		}
//		else if (redtype.compare("Max") == 0)
//		{
//			kerName = "reduce_generic_max_2_updated_multi";
//		}
//		else if (redtype.compare("Min") == 0)
//		{
//			kerName = "reduce_generic_min_2_updated_multi";
//		}
//		else if (redtype.compare("AbsMax") == 0)
//		{
//			kerName = "reduce_generic_absmax_2_updated_multi";
//		}
//		else
//		{
//			printf("%s not a valid reduce type. Choices are Sum, Max, Min, and AbsMax\n", redtype.c_str());
//			exit(-144194);
//		}
//		className = "GenericRedKernel";
//	}
//
//	void set_kernel_arguments()
//	{
//		int status = clSetKernelArg(ker, 0, sizeof(cl_mem), (void *)&BaseRed.Buffers[BaseRed.num_red_steps - 2]);
//		STATUS_CHECK_SETUP(status, "Error setting first argument");
//
//		status = clSetKernelArg(ker, 1, sizeof(cl_mem), ReduceVals.get_buf_add());
//		STATUS_CHECK_SETUP(status, "Error setting output buffer argument");
//
//		status = clSetKernelArg(ker, 2, BaseRed.local_size[BaseRed.num_red_steps - 1] * sizeof(T), NULL);
//		STATUS_CHECK_SETUP(status, "Error setting local buffer");
//
//		status = clSetKernelArg(ker, 3, sizeof(int), &cur_ind);
//		STATUS_CHECK_SETUP(status, "Error save index");
//	}
//
//
//	T& operator()(const int i)
//	{// Warning, this is a convenience function to get values after they have already been read from device memory
//		return ReduceVals(i);
//	}
//
//	void reset_counter()
//	{
//		cur_ind = 0;
//		next_ind = 0;
//	}
//
//	bool reduce(cl_command_queue *que_, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
//	{
//		cur_ind = next_ind;
//		int status = clSetKernelArg(ker, 3, sizeof(int), &cur_ind);
//		STATUS_CHECK_RUN(status, "Error setting save index");
//
//		BaseRed.call_kernels(num_wait, wait, que_);
//
//		status = clEnqueueNDRangeKernel(*que_, ker, 1, 0, &BaseRed.global_size[BaseRed.num_red_steps - 1],
//			&BaseRed.local_size[BaseRed.num_red_steps - 1], 0, NULL, evt);
//
//		STATUS_CHECK_RUN(status, "Error calling final reduce kernel");
//
//		next_ind++;
//
//		if (next_ind == Len)
//		{
//			ReduceVals.read_from_buffer(*queue, CL_TRUE);
//			next_ind = 0;
//			return true;
//		}
//		return false;
//	}
//
//	bool reduce(int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
//	{
//		return reduce(queue, num_wait, wait, evt);
//	}
//
//	T reduceSingle(int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
//	{
//		return reduceSingle(queue, num_wait, wait, evt);
//	}
//
//	T reduceSingle(cl_command_queue *que_, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
//	{ // this will not increment the counter so that value will be overwritten next time it is called
//
//		int zerval = 0;
//		//// Set kernel to write to ReduceOneVal at position 0///////////////
//		int status = clSetKernelArg(ker, 1, sizeof(cl_mem), ReduceOneVal.get_buf_add());
//		STATUS_CHECK_SETUP(status, "Error setting output buffer argument");
//		status = clSetKernelArg(ker, 3, sizeof(int), &zerval);
//		STATUS_CHECK_SETUP(status, "Error save index");
//
//		/////////// Call Reduce Functions ///////////////
//		BaseRed.call_kernels(num_wait, wait, que_);
//		status = clEnqueueNDRangeKernel(*que_, ker, 1, 0, &BaseRed.global_size[BaseRed.num_red_steps - 1],
//			&BaseRed.local_size[BaseRed.num_red_steps - 1], 0, NULL, evt);
//		STATUS_CHECK_RUN(status, "Error calling final reduce kernel save index");
//
//		/////////// Reset Kernel to write to ReduceVals at position cur_ind/////////////
//		status = clSetKernelArg(ker, 3, sizeof(int), &cur_ind);
//		STATUS_CHECK_RUN(status, "Error setting save index");
//		status = clSetKernelArg(ker, 1, sizeof(cl_mem), ReduceVals.get_buf_add());
//		STATUS_CHECK_SETUP(status, "Error setting output buffer argument");
//
//		///////////////// Return Reduced Value /////////////////////////
//		ReduceOneVal.read_from_buffer_size(*que_, CL_TRUE, 1);
//		return ReduceOneVal(0);
//	}
//
//	T getLastRead()
//	{
//		return ReduceVals(cur_ind);
//	}
//
//	void TransferUpToLast()
//	{
//		ReduceVals.read_from_buffer_size(*que, CL_TRUE, cur_ind);
//	}
//
//	T readAndGetLast()
//	{
//		readUpToLast();
//		return getLastRead();
//	}
//
//	bool savetxt()
//	{
//		return ReduceVals.savetxt(Name.c_str());
//	}
//
//	bool save_txt_from_device()
//	{
//		ReduceVals.read_from_buffer(*que, CL_TRUE);
//		return ReduceVals.save_txt_w_skip(Name.c_str(), skip_val);
//	}
//};






#endif