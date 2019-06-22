// Kernel and DualKernel Class, which help create kernels, set arguments
// and call clenqueuendrangekernel
// (c) Zachary Grant Mills, 2016 
//////////////////////////////////////////////////////////////////////


#if !defined(AFX_KERNELS_H__INCLUDED_)
#define AFX_KERNELS_H__INCLUDED_


#pragma once


#include "StdAfx.h"
#include "HelperFuncs.h"



#ifndef KERNEL_ERROR_CHECKING
#define CHECK_KERNEL_ERROR(test_, msg_, err_)		do {} while(0)
#else
#define CHECK_KERNEL_ERROR(test_, msg_, err_)   if((test_)) { Kernel::Gen_Error_Msg(err_, msg_); }
#endif

#define ALWAYS_CHECK_KERNEL_ERROR(test_, msg_, err_)   if((test_)) { Kernel::Gen_Error_Msg(err_, msg_); }



class Kernel
{
protected:
	cl_kernel ker;
	int dim;
	int optionInd; // index of an argument that changes frequently (used in convenience function setOption(...))

	size_t local_size[3], global_size[3], global_offset[3];
	cl_command_queue *queue;
	std::string Name;
	std::string kerName;


public:
	enum Dimension { X = 0, Y = 1, Z = 2 };

	Kernel()
	{
		optionInd = -1;
		queue = NULL;
		dim = 1;
		memset(local_size, 1, 3*sizeof(size_t));
		memset(global_size, 1, 3*sizeof(size_t));
		memset(global_offset, 0, 3*sizeof(size_t));
	}

	~Kernel()
	{
		free_memory();
	}

	size_t getLocalSize(Dimension dim_ = X)
	{
		return local_size[dim_];
	}

	size_t getGlobalSize(Dimension dim_ = X)
	{
		return global_size[dim_];
	}

	cl_command_queue* getDefaultQueue()
	{
		return queue;
	}

	void free_memory()
	{
		if (ker != NULL)
			clReleaseKernel(ker);
		ker = NULL;
	}

	void setOptionInd(int option_)
	{
		optionInd = option_;
	}

	int getOptionInd() { return optionInd; }

	void Gen_Error_Msg(int Errcode, std::string mesg)
	{
#ifdef LOG_ERROR_IN_FILE
		std::fstream stream;
		stream.open("error.txt", std::ios_base::out);
		stream << "Error in Kernel: " + Name +"\nCode: " + std::to_string(Errcode) +"\nMessage: " + mesg + "\n";
		stream.close();
#endif
		printf("Error in Kernel: %s, Code: %d\nMessage: %s\n", Name.c_str(), Errcode, mesg.c_str());
		exit(Errcode);
	}

	// kerName is name of kernel, Name is name for debugging (can be more descriptive, or to distinguish between 
	// instances calling the same opencl kernel)
	void create_kernel(cl_program program, cl_command_queue *queuet, const std::string kname, const std::string name_ = "")
	{
		kerName = kname;
		if (name_.length() == 0)
			Name = kname;
		else
			Name = name_;
		queue = queuet;

		int status;

		ker = clCreateKernel(program, kerName.c_str(), &status);

		ALWAYS_CHECK_KERNEL_ERROR(status, "Error creating kernel", status);
	}

	void setKernelInfo(cl_command_queue *queuet, const std::string kname, const std::string name_ = "")
	{
		kerName = kname;
		if (name_.length() == 0)
			Name = kname;
		else
			Name = name_;
		queue = queuet;
	}

	void create_kernel(cl_program &program)
	{
		int status;
		ker = clCreateKernel(program, kerName.c_str(), &status);
		ALWAYS_CHECK_KERNEL_ERROR(status, "Error creating kernel " + kerName, status);
	}

	size_t roundGlobalSize(const unsigned int gsize_, Dimension dim_ = X)
	{
#ifdef OPENCL_VERSION_1_2
		return (size_t)ceil((double)gsize_ / (double)local_size[dim_]) * local_size[dim_];
#else
		return (size_t)gsize_;
#endif
	}

	void set_size(globalSize_t gsizex, localSize_t lsizex, globalSize_t gsizey = 1,
		localSize_t lsizey = 1, globalSize_t gsizez = 1, localSize_t lsizez = 1)
	{
		if (gsizez == 1 && gsizey == 1)
			dim = 1;
		else if (gsizez == 1)
			dim = 2;
		else
			dim = 3;

		local_size[0] = lsizex;
		local_size[1] = lsizey;
		local_size[2] = lsizez;
		global_size[0] = roundGlobalSize(gsizex, X);
		global_size[1] = roundGlobalSize(gsizey, Y);
		global_size[2] = roundGlobalSize(gsizez, Z);
	}


	void set_size_with_offset(globalSize_t gsizex, localSize_t lsizex, size_t gwox, 
		globalSize_t gsizey = 1, localSize_t lsizey = 1, size_t gwoy = 0,
		globalSize_t gsizez = 1, localSize_t lsizez = 1, size_t gwoz = 0)
	{
#ifdef OPENCL_VERSION_1_2
		Gen_Error_Msg(ERROR_KERNEL_INITIALIZATION, "Cannot use kernel offsets with OpenCL 1.2");
#endif
		global_offset[0] = gwox;
		global_offset[1] = gwoy;
		global_offset[2] = gwoz;

		if (gsizez == 1 && gsizey == 1)
			dim = 1;
		else if (gsizez == 1)
			dim = 2;
		else
			dim = 3;

		local_size[0] = lsizex;
		local_size[1] = lsizey;
		local_size[2] = lsizez;
		global_size[0] = roundGlobalSize(gsizex, X);
		global_size[1] = roundGlobalSize(gsizey, Y);
		global_size[2] = roundGlobalSize(gsizez, Z);
	}


	void set_local_size(localSize_t lsizex, localSize_t lsizey = 1, localSize_t lsizez = 1)
	{
		if (lsizey == 1 && lsizez == 1)
			dim = 1;
		else if (lsizez == 1)
			dim = 2;
		else
			dim = 3;
		local_size[0] = lsizex;
	}

	void set_global_size(globalSize_t gsizex, globalSize_t gsizey = 1, globalSize_t gsizez = 1)
	{
		if (gsizez == 1 && gsizey == 1)
			dim = 1;
		else if (gsizez == 1)
			dim = 2;
		else
			dim = 3;

		global_size[0] = roundGlobalSize(gsizex, X);
		global_size[1] = roundGlobalSize(gsizey, Y);
		global_size[2] = roundGlobalSize(gsizez, Z);
	}

	void reset_global_size(globalSize_t gsizex, globalSize_t gsizey, globalSize_t gsizez)
	{
		global_size[0] = roundGlobalSize(gsizex, X);
		global_size[1] = roundGlobalSize(gsizey, Y);
		global_size[2] = roundGlobalSize(gsizez, Z);
	}

	void reset_global_size(globalSize_t gsizex, globalSize_t gsizey)
	{
		global_size[0] = roundGlobalSize(gsizex, X);
		global_size[1] = roundGlobalSize(gsizey, Y);
	}

	void reset_global_size(globalSize_t gsizex)
	{
		global_size[0] = roundGlobalSize(gsizex, X);
	}

	void reset_global_size(globalSize_t gsizex, Dimension dim_)
	{
		global_size[dim_] = roundGlobalSize(gsizex, dim_);
	}

	void reset_local_size(localSize_t lsizex, localSize_t lsizey, localSize_t lsizez)
	{
		local_size[0] = lsizex;
		local_size[1] = lsizey;
		local_size[2] = lsizez;
	}
	
	void reset_local_size(localSize_t lsizex, localSize_t lsizey)
	{
		local_size[0] = lsizex;
		local_size[1] = lsizey;
	}

	void reset_local_size(localSize_t lsizex)
	{
		local_size[0] = lsizex;
	}

	void reset_local_size(localSize_t lsize_, Dimension dim_)
	{
		local_size[dim_] = lsize_;
	}

	template <typename T>
	void set_argument(int ind, T *kerArg)
	{
		int status = clSetKernelArg(ker, ind, sizeof(T), (void*)kerArg);
		CHECK_KERNEL_ERROR(status, "Error setting kernel argument number " +\
			std::to_string(ind), status);
	}

	template <typename T>
	void setOption(T *val_)
	{
		CHECK_KERNEL_ERROR(optionInd == -1, "Must set optionInd before using setOptionCallKernel", ERROR_SETTING_KERNEL_ARG);
		set_argument<T>(optionInd, val_);
	}

	// TODO: see if this is necessary for 2.0 kernels since
	//		the actual (not rounded) size can be set to the exact
	//		number of kernel instances to execute.

	// This is for kernels with option indicies corresponding to
	// number of kernels to call
	void setOptionGlobalCallKernel(int val_, cl_command_queue* que = NULL)
	{
		setOption(&val_);
		set_global_call_kernel(val_, que);

	}

	void set_local_memory(int ind, size_t size_)
	{
		int status = clSetKernelArg(ker, ind, size_, NULL);
		CHECK_KERNEL_ERROR(status, "Error setting local mem size for kernel" +\
			" argument number " + std::to_string(ind), status);
	}

	void set_global_call_kernel(int gsizex, cl_command_queue *que = NULL)
	{
		CHECK_KERNEL_ERROR(gsizex < 1, "Global size must be >= 1", ERROR_CALLING_KERNEL);
		set_global_size(gsizex);
		call_kernel(que);
	}

	void call_kernel(cl_command_queue *que = nullptr, int num_list = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		if (que == nullptr)
			que = queue;
		
		int status = clEnqueueNDRangeKernel(*que, ker, dim, global_offset, global_size, local_size, num_list, wait, evt);
		CHECK_KERNEL_ERROR(status, "Error queueing kernel", status);
	}

	template <typename T>
	void setOptionCallKernel(T *val_, cl_command_queue *que = NULL, int num_list = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		CHECK_KERNEL_ERROR(optionInd == -1, "Must set optionInd before using setOptionCallKernel", ERROR_SETTING_KERNEL_ARG);
		set_argument<T>(optionInd, val);
		call_kernel(que, num_list, wait, evt);
	}
};

class DualKernel
{
protected:
	Kernel kerA, kerB;
	int alter;
	
public:
	enum kernelID { kernelA, kernelB, bothKernels };
	Kernel* getKernelAdd(kernelID kerid_)
	{
		if (kerid_ == kernelA)
			return &kerA;
		return &kerB;
	}


	DualKernel()
	{
		alter = 0;
	}

	~DualKernel()
	{
	}

	int getAlter()
	{
		return alter;
	}

	kernelID getAlterID()
	{
		if (alter == 0)
			return kernelA;
		return kernelB;
	}

	int setAlter(int alter_)
	{
		alter = alter_;
	}

	void setOptionInd(int option_)
	{
		kerA.setOptionInd(option_);
		kerB.setOptionInd(option_);
	}


	void setOptionInd(kernelID kID_, int option_)
	{
		if (kID_ == kernelA)
			kerA.setOptionInd(option_);
		else if (kID_ == kernelB)
			kerB.setOptionInd(option_);
		else
		{
			kerA.setOptionInd(option_);
			kerB.setOptionInd(option_);
		}
	}


	void create_kernel(cl_program &program, cl_command_queue *queuet, 
		const std::string knameA, const std::string knameB, std::string name_ = "")
	{
		alter = 0;
		
		kerA.create_kernel(program, queuet, knameA, name_);
		kerB.create_kernel(program, queuet, knameB, name_);
	}

	void setKernelInfo(cl_command_queue *queuet, const std::string knameA, 
		const std::string knameB, std::string name_ = "")
	{
		alter = 0;
		kerA.setKernelInfo(queuet, knameA, name_);
		kerB.setKernelInfo(queuet, knameB, name_);
	}

	void create_kernel(cl_program &program)
	{
		kerA.create_kernel(program);
		kerB.create_kernel(program);
	}

	void set_size(globalSize_t gsizex, localSize_t lsizex, globalSize_t gsizey = 1,
		localSize_t lsizey = 1, globalSize_t gsizez = 1, localSize_t lsizez = 1)
	{
		kerA.set_size(gsizex, lsizex, gsizey, lsizey, gsizez, lsizez);
		kerB.set_size(gsizex, lsizex, gsizey, lsizey, gsizez, lsizez);
	}

	void set_size(kernelID kID_, globalSize_t gsizex, localSize_t lsizex, globalSize_t gsizey = 1,
		localSize_t lsizey = 1, globalSize_t gsizez = 1, localSize_t lsizez = 1)
	{
		if (kID_ == kernelA)
			kerA.set_size(gsizex, lsizex, gsizey, lsizey, gsizez, lsizez);
		else if (kID_ == kernelB)
			kerB.set_size(gsizex, lsizex, gsizey, lsizey, gsizez, lsizez);
		else
		{
			kerA.set_size(gsizex, lsizex, gsizey, lsizey, gsizez, lsizez);
			kerB.set_size(gsizex, lsizex, gsizey, lsizey, gsizez, lsizez);
		}
	}
	
	void set_size_with_offset(globalSize_t gsizex, localSize_t lsizex, size_t gwox,
		globalSize_t gsizey = 1, localSize_t lsizey = 1, size_t gwoy = 0,
		globalSize_t gsizez = 1, localSize_t lsizez = 1, size_t gwoz = 0)
	{
		kerA.set_size_with_offset(gsizex, lsizex, gwox, gsizey, lsizey, gwoy, gsizez, lsizez, gwoz);
		kerB.set_size_with_offset(gsizex, lsizex, gwox, gsizey, lsizey, gwoy, gsizez, lsizez, gwoz);
	}

	void set_size_with_offset(kernelID kID_, globalSize_t gsizex, localSize_t lsizex, size_t gwox,
		globalSize_t gsizey = 1, localSize_t lsizey = 1, size_t gwoy = 0,
		globalSize_t gsizez = 1, localSize_t lsizez = 1, size_t gwoz = 0)
	{
		if (kID_ == kernelA)
			kerA.set_size_with_offset(gsizex, lsizex, gwox, gsizey, lsizey, gwoy, gsizez, lsizez, gwoz);
		else if (kID_ == kernelB)
			kerB.set_size_with_offset(gsizex, lsizex, gwox, gsizey, lsizey, gwoy, gsizez, lsizez, gwoz);
		else
		{
			kerA.set_size_with_offset(gsizex, lsizex, gwox, gsizey, lsizey, gwoy, gsizez, lsizez, gwoz);
			kerB.set_size_with_offset(gsizex, lsizex, gwox, gsizey, lsizey, gwoy, gsizez, lsizez, gwoz);
		}
	}


	void set_local_size(localSize_t lsizex, localSize_t lsizey = 1, localSize_t lsizez = 1)
	{
		kerA.set_local_size(lsizex, lsizey, lsizez);
		kerB.set_local_size(lsizex, lsizey, lsizez);
	}

	void set_local_size(kernelID kID_, localSize_t lsizex, localSize_t lsizey = 1, localSize_t lsizez = 1)
	{
		if (kID_ == kernelA)
			kerA.set_local_size(lsizex, lsizey, lsizez);
		else if (kID_ == kernelB)
			kerB.set_local_size(lsizex, lsizey, lsizez);
		else
		{
			kerA.set_local_size(lsizex, lsizey, lsizez);
			kerB.set_local_size(lsizex, lsizey, lsizez);
		}
	}

	void set_global_size(globalSize_t lsizex, globalSize_t lsizey = 1, globalSize_t lsizez = 1)
	{
		kerA.set_global_size(lsizex, lsizey, lsizez);
		kerB.set_global_size(lsizex, lsizey, lsizez);
	}

	void set_global_size(kernelID kID_, globalSize_t lsizex, globalSize_t lsizey = 1, globalSize_t lsizez = 1)
	{
		if (kID_ == kernelA)
			kerA.set_global_size(lsizex, lsizey, lsizez);
		else if (kID_ == kernelB)
			kerB.set_global_size(lsizex, lsizey, lsizez);
		else
		{
			kerA.set_global_size(lsizex, lsizey, lsizez);
			kerB.set_global_size(lsizex, lsizey, lsizez);
		}
	}

	void reset_global_size(globalSize_t lsizex, globalSize_t lsizey, globalSize_t lsizez)
	{
		kerA.reset_global_size(lsizex, lsizey, lsizez);
		kerB.reset_global_size(lsizex, lsizey, lsizez);
	}

	void reset_global_size(kernelID kID_, globalSize_t lsizex, globalSize_t lsizey, globalSize_t lsizez)
	{
		if (kID_ == kernelA)
			kerA.reset_global_size(lsizex, lsizey, lsizez);
		else if (kID_ == kernelB)
			kerB.reset_global_size(lsizex, lsizey, lsizez);
		else
		{
			kerA.reset_global_size(lsizex, lsizey, lsizez);
			kerB.reset_global_size(lsizex, lsizey, lsizez);
		}
	}


	void reset_global_size(globalSize_t lsizex, globalSize_t lsizey)
	{
		kerA.reset_global_size(lsizex, lsizey);
		kerB.reset_global_size(lsizex, lsizey);
	}

	void reset_global_size(kernelID kID_, globalSize_t lsizex, globalSize_t lsizey)
	{
		if (kID_ == kernelA)
			kerA.reset_global_size(lsizex, lsizey);
		else if (kID_ == kernelB)
			kerB.reset_global_size(lsizex, lsizey);
		else
		{
			kerA.reset_global_size(lsizex, lsizey);
			kerB.reset_global_size(lsizex, lsizey);
		}
	}

	void reset_global_size(globalSize_t lsizex)
	{
		kerA.reset_global_size(lsizex);
		kerB.reset_global_size(lsizex);
	}

	void reset_global_size(kernelID kID_, globalSize_t lsizex)
	{
		if (kID_ == kernelA)
			kerA.reset_global_size(lsizex);
		else if (kID_ == kernelB)
			kerB.reset_global_size(lsizex);
		else
		{
			kerA.reset_global_size(lsizex);
			kerB.reset_global_size(lsizex);
		}
	}

	void reset_global_size(globalSize_t lsizex, Kernel::Dimension dim_)
	{
		kerA.reset_global_size(lsizex, dim_);
		kerB.reset_global_size(lsizex, dim_);
	}

	void reset_global_size(kernelID kID_, globalSize_t lsizex, Kernel::Dimension dim_)
	{
		if (kID_ == kernelA)
			kerA.reset_global_size(lsizex, dim_);
		else if (kID_ == kernelB)
			kerB.reset_global_size(lsizex, dim_);
		else
		{
			kerA.reset_global_size(lsizex, dim_);
			kerB.reset_global_size(lsizex, dim_);
		}
	}

	template <typename T>
	void set_argument(int ind, T *kerArg)
	{
		kerA.set_argument<T>(ind, kerArg);
		kerB.set_argument<T>(ind, kerArg);
	}

	template <typename T>
	void set_argument(kernelID kID_, int ind, T *kerArg)
	{
		if (kID_ == kernelA)
			kerA.set_argument<T>(ind, kerArg);
		else if (kID_ == kernelB)
			kerB.set_argument<T>(ind, kerArg);
		else
		{
			kerA.set_argument<T>(ind, kerArg);
			kerB.set_argument<T>(ind, kerArg);
		}
	}

	template <typename T>
	void setOption(T *val_)
	{
		kerA.setOption<T>(val_);
		kerB.setOption<T>(val_);
	}

	template <typename T>
	void setOption(kernelID kID_, T *val_)
	{
		if (kID_ == kernelA)
			kerA.setOption<T>(val_);
		else if (kID_ == kernelB)
			kerB.setOption<T>(val_);
		else
		{
			kerA.setOption<T>(val_);
			kerB.setOption<T>(val_);
		}
	}

	void set_local_memory(int ind, size_t size_)
	{
		kerA.set_local_memory(ind, size_);
		kerB.set_local_memory(ind, size_);
	}

	void set_argument(kernelID kID_, int ind, size_t size_)
	{
		if (kID_ == kernelA)
			kerA.set_local_memory(ind, size_);
		else if (kID_ == kernelB)
			kerB.set_local_memory(ind, size_);
		else
		{
			kerA.set_local_memory(ind, size_);
			kerB.set_local_memory(ind, size_);
		}
	}

	void set_global_call_kernel(int gsizex, cl_command_queue *que = NULL)
	{
		if (alter == 0)
			kerA.set_global_call_kernel(gsizex, que);
		else
			kerB.set_global_call_kernel(gsizex, que);

		alter ^= 1;
	}

	void set_global_call_kernel(kernelID kID_, int gsizex, 
		cl_command_queue *que = NULL)
	{
		if (kID_ == kernelA)
			kerA.set_global_call_kernel(gsizex, que);
		else if (kID_ == kernelB)
			kerB.set_global_call_kernel(gsizex, que);
		else
		{
			kerA.set_global_call_kernel(gsizex, que);
			kerB.set_global_call_kernel(gsizex, que);
		}
	}

	void call_kernel(cl_command_queue *que = NULL, int num_list = 0, 
		cl_event *wait = NULL, cl_event *evt = NULL)
	{
		if (alter == 0)
			kerA.call_kernel(que, num_list, wait, evt);
		else
			kerB.call_kernel(que, num_list, wait, evt);

		alter ^= 1;
	}

	void call_kernel(kernelID kID_, cl_command_queue *que = NULL, 
		int num_list = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		if (kID_ == kernelA)
			kerA.call_kernel(que, num_list, wait, evt);
		else if (kID_ == kernelB)
			kerB.call_kernel(que, num_list, wait, evt);
		else
		{
			kerA.call_kernel(que, num_list, wait, evt);
			kerB.call_kernel(que, num_list, wait, evt);
		}
	}

	template <typename T>
	void setOptionCallKernel(T *val_, cl_command_queue *que = NULL, 
		int num_list = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		if (alter == 0)
			kerA.setOptionCallKernel<T>(val_, que, num_list, wait, evt);
		else
			kerB.setOptionCallKernel<T>(val_, que, num_list, wait, evt);

		alter ^= 1;
	}

	template <typename T>
	void setOptionCallKernel(kernelID kID_, T *val_, cl_command_queue *que = NULL,
		int num_list = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		if (kID_ == kernelA)
			kerA.setOptionCallKernel<T>(val_, que, num_list, wait, evt);
		else if (kID_ == kernelB)
			kerB.setOptionCallKernel<T>(val_, que, num_list, wait, evt);
		else
		{
			kerA.setOptionCallKernel<T>(val_, que, num_list, wait, evt);
			kerB.setOptionCallKernel<T>(val_, que, num_list, wait, evt);
		}
	}

	template <typename T>
	void set_separate_arguments(int ind, T *kerArgA, T *kerArgB)
	{
		kerA.set_argument<T>(ind, kerArgA);
		kerB.set_argument<T>(ind, kerArgB);
	}


	template <typename T>
	void set_argument_current_kernel(int ind, T *kerArg)
	{
		if (!alter)
			kerA.set_argument<T>(ind, kerArg);
		else
			kerB.set_argument<T>(ind, kerArg);
	}
};

#endif
