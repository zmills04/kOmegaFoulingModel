// Reducer Class which contains all necessary parts to reduce an
// array.
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_REDUCER_H__INCLUDED_)
#define AFX_REDUCER_H__INCLUDED_

#pragma once

#include "StdAfx.h"
#include "ReduceGenerator.h"

template <typename T>
class ReducerBase
{
protected:
	ReduceGenerator::reduceType redType;
	RedKernelList redKer;
	Array1D<T> redOneVal;
	TimeData<T> redVals;


	int arrSize;
	int cur_ind, next_ind;
	cl_mem *intermedBuf;
	bool intermedBufFl;
	ReduceGenerator::varType vType;


	void iniDomain(cl_mem *memloc_)
	{
		cl_mem *Bufout[1] = { redVals.get_buf_add() };
		ReduceGenerator::ReduceInstance()->genGenericDomainReduceSet(vType, redType, &memloc_,
			Bufout, &redKer, redMemSet);

		redMemSet ^= 1;
	}

	void iniDomain(cl_mem *mem1loc_, cl_mem *mem2loc_)
	{
		cl_mem* Bufin[2] = { mem1loc_, mem2loc_ };
		cl_mem *Bufout[1] = { redVals.get_buf_add() };
		ReduceGenerator::ReduceInstance()->genGenericDomainReduceSet(vType, redType, Bufin, Bufout, &redKer, redMemSet);

		redMemSet ^= 1;
	}

	void iniArrays(bool restartFlag)
	{
		redVals.iniFile(restartFlag);
		redVals.zeros(arrSize, arrSize, 1, 1);
		redVals.createTimeArray();
		redVals.allocate_buffer_w_copy();

		redOneVal.zeros(1);
		redOneVal.allocate_buffer_w_copy();
	}

public:
	std::string Name;
	static int redMemSet;
	ReduceGenerator::varType getVarType();

	ReducerBase() : cur_ind(0), next_ind(0)
	{
		intermedBufFl = false;
		vType = getVarType();
	}

	~ReducerBase()
	{
		if (intermedBufFl)
		{
			for (unsigned int i = 0; i < redKer.getSize() - 1; i++)
			{
				clReleaseMemObject(intermedBuf[i]);
			}
			delete[] intermedBuf;
		}
	}

	// Outname will be used for debugging and will also be the output filename
	void ini(ArrayBase<T> &memloc_, bool restartFlag, std::string outname_ = "")
	{
		if (outname_.length() == 0)
			outname_ = memloc_.getName() + "Red" + ReduceGenerator::getGenericName(redType);
		Name = outname_;
		iniArrays(restartFlag);
		// TODO: make it to where I can pass in the string directly without needing to initialize it
		// ahead of time.


		ERROR_CHECKING((redType == ReduceGenerator::Dot), "Trying to initialize " + Name + \
			" as a dot prod reduce but provided only one buffer", ERROR_IN_REDUCE_KERNEL);

		ERROR_CHECKING((redKer.isInitialized()), "Trying to initialize " + Name + \
			" after its already been initialized", ERROR_IN_REDUCE_KERNEL);

		redVals.setName(Name);
		redOneVal.setName(Name + "Single");
		
		// need to make sure buffer is correctly padded,
		// and reallocate if necessary
		int fsize = ReduceGenerator::ReduceInstance()->getPaddedSize(memloc_.getBufferFullSize());
		if (fsize != memloc_.getBufferFullSize())
		{
			memloc_.padBuffer();
		}

		if (ReduceGenerator::ReduceInstance()->domainRedSize == memloc_.getBufferFullSize())
		{
			iniDomain(memloc_.get_buf_add());
			return;
		}

		int numker = 0;

		cl_mem* Bufin[1] = { memloc_.get_buf_add() };
		cl_mem *Bufout[1] = { redVals.get_buf_add() };

		ReduceGenerator::ReduceInstance()->genGenericReduceSet(vType, fsize, redType, Bufin,
			Bufout, &redKer, numker, &intermedBuf);
		intermedBufFl = true;
	}

	void ini(ArrayBase<T> &mem1loc_, ArrayBase<T> &mem2loc_, bool restartFlag, std::string outname_ = "")
	{
		if (outname_.length() == 0)
			outname_ = mem1loc_.getName() + "Dot" + mem2loc_.getName();
		Name = outname_;
		iniArrays(restartFlag);
		ERROR_CHECKING((redType != ReduceGenerator::Dot), "Provided two buffers to non dot prod reduce kernel " + \
			Name + ". Maybe want to create a non-generic reduce kernel set?", ERROR_IN_REDUCE_KERNEL);

		ERROR_CHECKING((redKer.isInitialized()), "Trying to initialize " + Name + \
			" after its already been initialized", ERROR_IN_REDUCE_KERNEL);

		redVals.setName(Name);
		redOneVal.setName(Name + "Single");
		//if (clEnv::instance()->getRestartFlag() == false)
		//	redVals.append2file(-2); // initializes file to be appended to if not a restart

		// need to make sure buffer is correctly padded,
		// and reallocate if necessary
		int fsize1 = ReduceGenerator::ReduceInstance()->getPaddedSize(mem1loc_.getBufferFullSize());
		int fsize2 = ReduceGenerator::ReduceInstance()->getPaddedSize(mem2loc_.getBufferFullSize());

		if ((fsize1 != mem1loc_.getBufferFullSize()) ||
			(fsize2 != mem2loc_.getBufferFullSize()))
		{
			ERROR_CHECKING((fsize1 != fsize2), "trying to do dot prod reduce on "\
				"diff sized arrays using " + Name, ERROR_IN_REDUCE_KERNEL);

			if (mem1loc_.getFullSize() != mem2loc_.getFullSize())
			{
				printf("Warning: trying to do dot prod reduce on arrays with diff "\
					"sizes, but equal padded sizes. Make sure this is not an error\n");
			}

			mem1loc_.padBuffer(fsize1);
			mem2loc_.padBuffer(fsize1);
		}

		if (ReduceGenerator::ReduceInstance()->domainRedSize == mem1loc_.getBufferFullSize())
		{
			iniDomain(mem1loc_.get_buf_add(), mem2loc_.get_buf_add());
			return;
		}

		int numker = 0;
		cl_mem *Bufout[1] = { redVals.get_buf_add() };
		cl_mem* Bufin[2] = { mem1loc_.get_buf_add(), mem2loc_.get_buf_add() };
		ReduceGenerator::ReduceInstance()->genGenericReduceSet(vType, fsize1, redType, Bufin,
			Bufout, &redKer, numker, &intermedBuf);
		intermedBufFl = true;
	}


	T& operator()(const int i)
	{// Warning, this is a convenience function to get values after
		// they have already been read from device memory
		return redVals(i);
	}

	void reset_counter()
	{
		cur_ind = 0;
		next_ind = 0;
	}

	bool reduce(int time_, cl_command_queue *que_ = NULL, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		if (que_ == NULL)
		{
			que_ = clEnv::instance()->getIOqueue();
		}

		cur_ind = next_ind;
		int status = redKer.setOutputIndex(cur_ind);
		ERROR_CHECKING(status, "Error returned setting output index of "\
			+ Name, status);
		for (thinKerWrapper *rk = redKer.begin(); rk != redKer.end(); ++rk)
		{
			int status = rk->operator()(que_, num_wait, wait, evt);
			ERROR_CHECKING(status, "Error returned running"\
				" reduce kernels of " + Name, status);
		}
		bool retval = redVals.setTimeAndIncrement(p.Time, que_);
		if (retval)
			next_ind = 0;
		
	}


	void appendToFileFromDevice(cl_command_queue* que_ = NULL, cl_bool block_flag = CL_TRUE, 
		int num_wait = 0, cl_event* wait = NULL, cl_event* evt = NULL)
	{
		if (cur_ind == 0)
			return;
		
		if (que_ == NULL)
		{
			que_ = clEnv::instance()->getIOqueue();
		}
		
		redVals.appendData_from_device(que_, block_flag, num_wait, wait, evt);
		next_ind = 0;
	}


	T reduceSingle(cl_command_queue *que_ = NULL, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{ // this will not increment the counter so that value will be overwritten next time it is called
		if (que_ == NULL)
		{
			que_ = clEnv::instance()->getIOqueue();
		}

		//// Set kernel to write to ReduceOneVal at position 0///////////////
		int status = redKer.setOutputIndex(0);
		ERROR_CHECKING(status, "Error returned setting output index of "\
			+ Name, status);
		status = redKer.setOutputArray(redOneVal.get_buf_add());
		ERROR_CHECKING(status, "Error returned changing output of "\
			+ Name + " to redOneVal buffer", status);

		/////////// Call Reduce Functions ///////////////
		for (thinKerWrapper *rk = redKer.begin(); rk != redKer.end(); ++rk)
		{
			int status = rk->operator()(que_, num_wait, wait, evt);
			ERROR_CHECKING(status, "Error returned running"\
				" reduce kernels of " + Name, status);
		}

		/////////// Reset Kernel to write to ReduceVals at position cur_ind/////////////
		redKer.setOutputIndex(cur_ind);
		ERROR_CHECKING(status, "Error returned setting output index of "\
			+ Name, status);
		redKer.setOutputArray(redVals.get_buf_add());
		ERROR_CHECKING(status, "Error returned changing output of "\
			+ Name + " to redVals buffer", status);

		redOneVal.read_from_buffer_size(1, que_);
		return redOneVal(0);
	}

	T getLastRead()
	{
		return redVals(cur_ind);
	}

	void transferUpToLast(cl_command_queue *que_ = NULL)
	{
		redVals.read_from_buffer_size(cur_ind, que_);
	}

	T readAndGetLast(cl_command_queue *que_ = NULL)
	{
		transferUpToLast(que_);
		return getLastRead();
	}
};

template <typename T> int ReducerBase<T>::redMemSet = 0;


// This is to allow for passing reducer without instantiating multiple functions
// for each type of reducer (argument will be ReducerBase, which will accept all
// versions of Reducer)
// For domain variables, which use the same intermediate buffers, only variable types which have sizes
// less than or equal to doubles can be used (opencl buffers are type agnostic (addressing based on the 
// type declared in the kernel, so a buffer for 100 doubles can store 100 ints with only have the buffer being used)
template <typename T, ReduceGenerator::reduceType redT_ = ReduceGenerator::Sum, int ArrSize_ = REDUCE_RESULTS_SIZE>
class Reducer : public ReducerBase<T>
{
public:
	Reducer()
	{
		redType = redT_;
		arrSize = ArrSize_;
	}
};


#endif