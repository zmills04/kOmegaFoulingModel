// Array Classes which hold both Host and Device memory and provide
// convenience functions for transfering memory between host and device
// (c) Zachary Grant Mills, 2016 
//////////////////////////////////////////////////////////////////////

// TODO: add pinned memory array
// TODO: check if necessary to set name in all derived array classes

#if !defined(AFX_ARRAY_H__INCLUDED_)
#define AFX_ARRAY_H__INCLUDED_

#pragma once

#include "StdAfx.h"
#include "oclEnvironment.h"
#include "logger.h"
#include "StreamFuncs.h"


#ifdef _DEBUG
#define ARRAY_DEBUG
#endif //_DEBUG


#ifndef BIN_SAVE_DIR
#define BIN_SAVE_DIR	""
#endif //BIN_SAVE_DIR


#ifndef ARRAY_ERROR_CHECKING
#define CHECK_ARRAY_ERROR(test_, msg_, err_)		do {} while(0)
#else
#define CHECK_ARRAY_ERROR(test_, msg_, err_)   if((test_)) { ArrayBase::Gen_Error_Msg(err_, msg_); }
#endif


// Simple array with pointer based iterator
// would probably be better to use something
// from std library like std::list or vector,
// but was curious about how to implement 
// a simple iterator. 
template <class T>
class iterList
{
protected:
	T* list_;
	unsigned int size;
public:
	iterList()
	{
		list_ = nullptr;
		size = 0;
	}
	~iterList()
	{
		delete[] list_;
	}
	void ini(const unsigned int size_)
	{
		size = size_;
		if (list_ != nullptr)
			delete[] list_;
		
		list_ = new T[size];
	}
	T* begin() { return &list_[0]; }
	T* end() { return &list_[size]; }

	T& operator()(int i_)
	{
		return list_[i_];
	}
	T operator[](int i_)
	{
		return list_[i_];
	}
	bool isInitialized()
	{
		if (size > 0)
			return true;
		return false;
	}
	unsigned int getSize() { return size; }
};

// small kernel class for reduce kernels, do not use for main kernels
// as this class has limited functionality
class thinKerWrapper
{
private:
	cl_kernel ker;
	bool kerCreated = false;
	size_t gsize, lsize;
public:
	~thinKerWrapper()
	{
		if (kerCreated)
			clReleaseKernel(ker);

	}
	void setSizes(int gsize_, int lsize_)
	{
		gsize = (size_t)gsize_;
		lsize = (size_t)lsize_;
	}

	cl_kernel& getKernel() { return ker; }

	int operator()(const cl_command_queue *que_, const int num_wait = 0,
		cl_event* wait = NULL, cl_event* evt = NULL)
	{
		return clEnqueueNDRangeKernel(*que_, ker, 1, NULL, &gsize, &lsize, num_wait, wait, evt);
	}

	int createKernel(cl_program *program, std::string name_)
	{
		kerCreated = true;
		int status;
		ker = clCreateKernel(*program, name_.c_str(), &status);
		return status;
	}

	template <typename T>
	int setArgument(const int ind_, T *val_)
	{
		return clSetKernelArg(ker, ind_, sizeof(T), (void*)val_);
	}

	template <typename T>
	int setArgument(const int ind_, const T val_)
	{
		return clSetKernelArg(ker, ind_, sizeof(T), (void*)&val_);
	}

	int setLocalMem(const int ind_, size_t size_loc)
	{
		return clSetKernelArg(ker, ind_, size_loc, NULL);
	}
};

// main reason for implementation was for storing reduce kernels
class RedKernelList : public iterList<thinKerWrapper>
{
private:
	int output_index;
	int output_array_index;
public:
	RedKernelList() : output_index(-1), output_array_index(-1) {}
	void setValueOutputIndex(const int outi_) { output_index = outi_; }
	void setValueOutputArrayIndex(const int outi_) { output_array_index = outi_; }
	// set argument will try error, so no need to check for one here
	int setOutputIndex(const int outval_)
	{
		return list_[size - 1].setArgument(output_index, &outval_);
	}
	int setOutputArray(cl_mem* buf_)
	{
		return list_[size - 1].setArgument(output_array_index, buf_);
	}


	int getOutputIndex() { return output_index; }
};

// Forward declarations 
template <typename T>
class Array1D;

template <typename T>
class Array2D;

template <typename T>
class Array3D;


class staticBaseVar
{
public:
	static cl_context* context;
	static std::string curDumpSaveDir;



	static void setArrayContext(cl_context *context_)
	{
		context = context_;// clEnv::instance()->getContext();
	}


	static void setCurDumpDir(std::string curDump_)
	{
		curDumpSaveDir = curDump_;
	}
};


template <class T> 
class ArrayBase : public staticBaseVar
{
public:
	cl_command_queue *ioQue;

	typedef T		TYPE;
protected:
	T *Array;
	cl_mem Buffer;
	int SizeX, SizeY, SizeZ;
	int FullSizeX, FullSizeY, FullSizeZ;
	int FullSize, BufferFullSize;
	int Host_Alloc_Flag, Device_Alloc_Flag;
	int clFlags;

	std::string Name = "";
	

	void Gen_Error_Msg(int Errcode, std::string mesg)
	{
#ifdef LOG_ERROR_IN_FILE
		std::fstream stream;
		fileopen(stream, std::string("error.txt"), FileOut);
		stream << "Error in Array: " + Name + "\nCode: " + std::to_string(Errcode) + "\nMessage: " + mesg + "\n";
		stream.close();
#endif
		printf("Error in Array: %s, Code: %d\nMessage: %s\n", Name.c_str(), Errcode, mesg.c_str());
		exit(Errcode);
	}
	
	void ini(std::string name_)
	{
		SizeX = 0;
		FullSizeX = 0;
		SizeY = 0;
		FullSizeY = 0;
		SizeZ = 0;
		FullSizeZ = 0;
		FullSize = 0;
		Host_Alloc_Flag = 0;
		Device_Alloc_Flag = 0;
		Array = nullptr;
		Buffer = nullptr;
		clFlags = 0;
		ioQue = clEnv::instance()->getIOqueue();
//		context = clEnv::instance()->getContext();
		Name = name_;
	}

	void ini()
	{
		SizeX = 0;
		FullSizeX = 0;
		SizeY = 0;
		FullSizeY = 0;
		SizeZ = 0;
		FullSizeZ = 0;
		FullSize = 0;
		Host_Alloc_Flag = 0;
		Device_Alloc_Flag = 0;
		Array = nullptr;
		Buffer = nullptr;
		clFlags = 0;
		ioQue = clEnv::instance()->getIOqueue();
//		context = clEnv::instance()->getContext();
	}

	ArrayBase(const ArrayBase &other)
	{
		SizeX = other.SizeX;
		FullSizeX = other.FullSizeX;
		SizeY = other.SizeX;
		FullSizeY = other.FullSizeY;
		SizeZ = other.SizeX;
		FullSizeZ = other.FullSizeZ;
		FullSize = other.FullSize;
		Host_Alloc_Flag = other.Host_Alloc_Flag;
		Device_Alloc_Flag = other.Device_Alloc_Flag;
		clFlags = other.clFlags;
		Name = other.Name;
		BufferFullSize = other.BufferFullSize;
		ioQue = clEnv::instance()->getIOqueue();
		//context = clEnv::instance()->getContext();
		Array = nullptr;
		Buffer = nullptr;
		if (Host_Alloc_Flag)
		{
			AllocHost();
			memcpy(Array, other.get_array(), FullSize*sizeof(T));
		}
		if (Device_Alloc_Flag)
		{
			allocate_buffer_size(BufferFullSize);
			copy_to_buffer();
		}
	}

	//ArrayBase<T>& operator=(const ArrayBase<T> &other)
	//{
	//	SizeX = other.SizeX;
	//	FullSizeX = other.FullSizeX;
	//	SizeY = other.SizeX;
	//	FullSizeY = other.FullSizeY;
	//	SizeZ = other.SizeX;
	//	FullSizeZ = other.FullSizeZ;
	//	FullSize = other.FullSize;
	//	Host_Alloc_Flag = other.Host_Alloc_Flag;
	//	Device_Alloc_Flag = other.Device_Alloc_Flag;
	//	clFlags = other.clFlags;
	//	Name = other.Name;
	//	BufferFullSize = other.BufferFullSize;
	//	ioQue = clEnv::instance()->getIOqueue();
	//	context = clEnv::instance()->getContext();
	//	Array = nullptr;
	//	Buffer = nullptr;
	//	if (Host_Alloc_Flag)
	//	{
	//		AllocHost();
	//		memcpy(Array, other.get_array(), FullSize*sizeof(T));
	//	}
	//	if (Device_Alloc_Flag)
	//	{
	//		allocate_buffer_size(BufferFullSize);
	//		copy_to_buffer();
	//	}
	//}


	void AllocHost()
	{
		if (Array != nullptr)
			delete[] Array;
		Host_Alloc_Flag = 1;
		Array = new T[FullSize];
	}

	void AllocHostZeros()
	{
		AllocHost();
		zeros();
	}

	void AllocDevice(cl_mem_flags clflags, int buffersize = -1)
	{
		if (buffersize == -1)
			buffersize = FullSize;

		if (Device_Alloc_Flag) { clReleaseMemObject(Buffer); }
		clFlags = clflags;
		int status;
		Buffer = clCreateBuffer(*context, clFlags, sizeof(T)* buffersize, NULL, &status);
		CHECK_ARRAY_ERROR(status, "Error Allocating Buffer", status);

		Device_Alloc_Flag = 1;
	}
	

	void Copy_HtoD(cl_command_queue *queue = NULL, cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		CHECK_ARRAY_ERROR((Device_Alloc_Flag == 0 || Host_Alloc_Flag == 0),
			"Copying unallocated memory to device", ERROR_BUFFER_ALLOCATION);

		if (queue == NULL) { queue = ioQue;  }
		
		int status = clEnqueueWriteBuffer(*queue, Buffer, block_flag, 0, FullSize*sizeof(T), Array, num_wait, wait, evt);
		
		CHECK_ARRAY_ERROR(status, "Copying to device buffer", status);
	}

	void Copy_DtoH(cl_command_queue *queue = NULL, cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		CHECK_ARRAY_ERROR((Device_Alloc_Flag == 0 || Host_Alloc_Flag == 0),
			"Copying unallocated memory to device", ERROR_BUFFER_ALLOCATION);

		if (queue == NULL) { queue = ioQue; }

		int status = clEnqueueReadBuffer(*queue, Buffer, block_flag, 0, FullSize*sizeof(T), Array, num_wait, wait, evt);

		CHECK_ARRAY_ERROR(status, "Copying from buffer to device", status);
	}

	void Copy_HtoD_Size(int num_el, cl_command_queue *queue = NULL, cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		CHECK_ARRAY_ERROR((Device_Alloc_Flag == 0 || Host_Alloc_Flag == 0),
			"Copying unallocated memory to device", ERROR_BUFFER_ALLOCATION);

		if (queue == NULL) { queue = ioQue; }

		CHECK_ARRAY_ERROR((num_el > BufferFullSize),
			"Trying to copy memory region larger than buffer length", ERROR_BUFFER_ALLOCATION);

		if (num_el > FullSize)
		{
			int status = clEnqueueWriteBuffer(*queue, Buffer, block_flag, 0, FullSize*sizeof(T), Array, num_wait, wait, NULL);
			CHECK_ARRAY_ERROR(status, "Copying to device buffer", status);
			T zer{ 0 };
			FillBufferFunc((void*)&zer, queue, sizeof(T), num_el - FullSize, 0, NULL, evt);
		}
		int status = clEnqueueWriteBuffer(*queue, Buffer, block_flag, 0, FullSize*sizeof(T), Array, num_wait, wait, NULL);
		CHECK_ARRAY_ERROR(status, "Copying to device buffer", status);
	}

	void Copy_DtoH_Size(int num_el, cl_command_queue *queue = NULL, cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		CHECK_ARRAY_ERROR((Device_Alloc_Flag == 0 || Host_Alloc_Flag == 0),
			"Copying unallocated memory to device", ERROR_BUFFER_ALLOCATION);

		if (queue == NULL) { queue = ioQue; }

		int status = clEnqueueReadBuffer(*queue, Buffer, block_flag, 0, num_el*sizeof(T), Array, num_wait, wait, evt);

		CHECK_ARRAY_ERROR(status, "Copying from buffer to device", status);
	}

	void enqueue_copy_buffer_to_buffer(cl_mem src_buf, cl_mem dest_buf, int copy_length = -1,
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event* wait = NULL, cl_event *blk_evt = NULL)
	{
		if (copy_length = -1) { copy_length = FullSize; }
		if (queue == NULL) { queue = ioQue; }

		int status = clEnqueueCopyBuffer(*queue, src_buf, dest_buf, 0, 0, copy_length * sizeof(T), num_wait, wait, blk_evt);

		CHECK_ARRAY_ERROR(status, "Copying between buffers", status);
	}


	void FillBufferFunc(void* fillval, cl_command_queue *queue, size_t elsize = -1, int size = -1, int offset = 0,
		int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		CHECK_ARRAY_ERROR((Device_Alloc_Flag == 0),
			"Filling Uninitialized buffer\n", ERROR_MEMORY_COPY);

		if (queue == NULL) { queue = ioQue; }
		if (elsize == -1) { elsize = sizeof(T); }
		if (size == -1) { size = FullSize; }

		int	status = clEnqueueFillBuffer(*queue, Buffer, fillval, elsize, offset*elsize, size*elsize, num_wait, wait, evt);
		
		CHECK_ARRAY_ERROR(status, "Error filling device buffer", status);
	}

public:
	typedef T type;

	ArrayBase()
	{
		if (Name.length() == 0)
			Name = "Default";
	}

	virtual ~ArrayBase()
	{
		FreeHost();
		FreeDevice();
	}

	virtual ArrayBase& operator=(ArrayBase &other)
	{
		SizeX = other.getSizeX();
		FullSizeX = other.getFullSizeX();
		SizeY = other.getSizeY();
		FullSizeY = other.getFullSizeY();
		SizeZ = other.getSizeZ();
		FullSizeZ = other.getFullSizeZ();
		FullSize = other.getFullSize();
		Host_Alloc_Flag = other.getHostAllocFlag();
		Device_Alloc_Flag = other.getDeviceAllocFlag();
		clFlags = other.getclFlags();
		Name = other.getName();
		BufferFullSize = other.getBufferFullSize();
		ioQue = clEnv::instance()->getIOqueue();
		//context = clEnv::instance()->getContext();
		Array = nullptr;
		Buffer = nullptr;
		if (Host_Alloc_Flag)
		{
			AllocHost();
			memcpy(Array, other.get_array(), FullSize*sizeof(T));
		}
		if (Device_Alloc_Flag)
		{
			allocate_buffer_size(BufferFullSize);
			copy_to_buffer();
		}
		return *this;
	}

	virtual ArrayBase& operator=(ArrayBase *other)
	{
		SizeX = other->getSizeX();
		FullSizeX = other->getFullSizeX();
		SizeY = other->getSizeY();
		FullSizeY = other->getFullSizeY();
		SizeZ = other->getSizeZ();
		FullSizeZ = other->getFullSizeZ();
		FullSize = other->getFullSize();
		Host_Alloc_Flag = other->getHostAllocFlag();
		Device_Alloc_Flag = other->getDeviceAllocFlag();
		clFlags = other->getclFlags();
		Name = other->getName();
		BufferFullSize = other->getBufferFullSize();
		ioQue = clEnv::instance()->getIOqueue();
		//context = clEnv::instance()->getContext();
		Array = nullptr;
		Buffer = nullptr;
		if (Host_Alloc_Flag)
		{
			AllocHost();
			memcpy(Array, other->get_array(), FullSize*sizeof(T));
		}
		if (Device_Alloc_Flag)
		{
			allocate_buffer_size(BufferFullSize);
			copy_to_buffer();
		}
		return *this;
	}

	
	int getHostAllocFlag() { return Host_Alloc_Flag; }
	int getDeviceAllocFlag() { return Device_Alloc_Flag; }
	int getclFlags() { return clFlags; }

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
///////////                                               ///////////
///////////    Initialization and allocation function     ///////////
///////////                                               ///////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

	void setName(std::string name_)
	{
		Name = name_;
	}
	
	void allocate_buffer(cl_mem_flags clflags = CL_MEM_READ_WRITE)
	{
		BufferFullSize = FullSize;
		AllocDevice(clflags, FullSize);
	}

	void allocate_buffer_size(cl_mem_flags clflags, int buffersize)
	{
		BufferFullSize = buffersize;
		AllocDevice(clflags, buffersize);

		// Will set padding to zeros if the buffer is larger than the host array.
		if(BufferFullSize > FullSize)
			Copy_HtoD_Size(BufferFullSize);
	}

	void resetSizesDebug(int newsize)
	{
		FullSize = newsize;
		SizeX = newsize;
		BufferFullSize = newsize;
	}

	void allocate_buffer_size(int buffersize)
	{
		allocate_buffer_size(CL_MEM_READ_WRITE, buffersize);
	}

	void allocate_buffer_w_copy(cl_mem_flags clflags = CL_MEM_READ_WRITE, cl_bool block_flag = CL_TRUE)
	{
		BufferFullSize = FullSize;
		AllocDevice(clflags);
		Copy_HtoD(NULL, block_flag);
	}

	T reduce()
	{
		T ret_ = (T)0; 
		for (int i = 0; i < SizeX; i++)
		{
			for (int j = 0; j < SizeY; j++)
			{
				for (int k = 0; k < SizeZ; k++)
				{
					ret_ += Array[i + j*FullSizeX + k*FullSizeY*FullSizeX];
				}
			}
		}
		return ret_;
	}

	void padBuffer(cl_mem_flags clflags = CL_MEM_READ_WRITE)
	{
		int fsize = 2;
		while (fsize < FullSize)
			fsize *= 2;

		Copy_DtoH(NULL, CL_TRUE);
		allocate_buffer_size(clflags, fsize);
	}

	void zeros()
	{
		if (FullSize == 0)
			return;
		T* p = (T*)calloc(FullSize, sizeof(T));
		memcpy(Array, p, FullSize*sizeof(T));
		FREE(p);
	}

	//Free Host Memory
	void FreeHost()
	{
		if (Array != nullptr)
			delete[] Array;
		Array = nullptr;
		Host_Alloc_Flag = 0;
	}

	void FreeDevice(void)
	{
		if (Device_Alloc_Flag == 1)
		{
			int status = clReleaseMemObject(Buffer);
			CHECK_ARRAY_ERROR(status, "Error freeing buffer", status);
		}
		Device_Alloc_Flag = 0;
	}

	void delete_array()
	{
		ArrayBase<T>::FreeHost();
		ArrayBase<T>::FreeDevice();
	}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
///////////                                               ///////////
///////////  Convenience Functions for testing/debugging  ///////////
///////////                                               ///////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
	
	bool Test_Negatives()
	{
		for (int i = 0; i < FullSize; i++)
		{
			if (Array[i] < 0)
				return true;
		}
		return false;
	}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
///////////                                               ///////////
/////////// Functions for accessing private class members ///////////
///////////                                               ///////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

	std::string getName()
	{
		return Name;
	}

	void* get_add_Macro()
	{
		return get_buf_add();
	}

	T* get_array()
	{
		return Array;
	}

	cl_mem get_buffer()
	{
		return Buffer;
	}

	cl_mem* get_buf_add()
	{
		return &Buffer;
	}

	int getSizeX(void)
	{
		return SizeX;
	}

	int getSizeY(void)
	{
		return SizeY;
	}

	int getSizeZ(void)
	{
		return SizeZ;
	}

	int getFullSize(void)
	{
		return FullSize;
	}

	int getFullSizeX(void)
	{
		return FullSizeX;
	}

	int getFullSizeY(void)
	{
		return FullSizeY;
	}

	int getFullSizeZ(void)
	{
		return FullSizeZ;
	}

	int getSize(int n)
	{
		if (n == 1)
			return SizeX;
		else if (n == 2)
			return SizeY;
		else if (n == 3)
			return SizeZ;

		return NULL;
	}

	int getBufferFullSize()
	{
		return BufferFullSize;
	}

	int getFullSize(int n)
	{
		if (n == 1)
			return FullSizeX;
		else if (n == 2)
			return FullSizeY;
		else if (n == 3)
			return FullSizeZ;

		return NULL;
	}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////                                               /////////////////
/////////////////  Functions for writing between mem locations  /////////////////
/////////////////                                               /////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

	//void copyHostMem(ArrayBase* other)
	//{
	//	memcpy(Array, other.get_array(), FullSize * sizeof(T));
	//}

/////////////////////////////////////////////////////////////////////////////////
///////////////////    Filling memory with a single value    ////////////////////
/////////////////////////////////////////////////////////////////////////////////

	void FillBuffer(T fillval, cl_command_queue *queue = NULL, int num_wait = 0,
		cl_event *wait = NULL, cl_event *evt = NULL)
	{
		FillBufferFunc((void*)&fillval, queue, sizeof(T), FullSize, 0, num_wait, wait, evt);
	}
	
	void FillBuffer(T fillval, int fillsize, cl_command_queue *queue = NULL, int num_wait = 0,
		cl_event *wait = NULL, cl_event *evt = NULL)
	{
		FillBufferFunc((void*)&fillval, queue, sizeof(T), fillsize, 0, num_wait, wait, evt);
	}

	// Fills with a generic type (i.e. if filling Array1Dv2d with cl_double)
	// fillsize needs to be provided here since we are most likely filling with different
	// sized types
	void FillBuffer(void* fillval, size_t typesize, int fillsize, int offset = 0, cl_command_queue *queue = NULL, int num_wait = 0,
		cl_event *wait = NULL, cl_event *evt = NULL)
	{
		FillBufferFunc(fillval, queue, typesize, fillsize, offset, num_wait, wait, evt);
	}
	
	// dont use frequently as this sets entire padded array to zero and then fills it
	void fill(const T val)
	{
		zeros();
		for (int i = 0; i < SizeX; i++)
		{
			for (int j = 0; j < SizeY; j++)
			{
				for (int k = 0; k < SizeZ; k++)
				{
					Array[i + FullSizeX*(j + k*FullSizeZ)] = val;
				}
			}
		}
	}
	
	// does not change values in padded sections
	// TODO: test to see if this is faster than fillAll (it probably isnt)
	void fastfill(const T val)
	{
		for (int i = 0; i < SizeX; i++)
		{
			for (int j = 0; j < SizeY; j++)
			{
				for (int k = 0; k < SizeZ; k++)
				{
					Array[i + FullSizeX*(j + k*FullSizeZ)] = val;
				}
			}
		}
	}

	// sets all elements in array (including padded regions) to val
	void fillAll(const T val)
	{
		T *p = Array;
		for (int i = 0; i < FullSize; i++) 
			*(p++) = val;
	}



	void fillByNodeType(const T val, Array2D<NTYPE_TYPE>& Map, NTYPE_TYPE nodeVal = M_FLUID_NODE)
	{
		for (int i = 0; i < SizeX; i++)
		{
			for (int j = 0; j < SizeY; j++)
			{
				if ((Map(i, j) & nodeVal) == 0)
					continue;

				for (int k = 0; k < SizeZ; k++)
				{
					Array[i + FullSizeX * (j + k * FullSizeZ)] = val;
				}
			}
		}
	}

/////////////////////////////////////////////////////////////////////////////////
///////////////////   Writing Array member to buffer member  ////////////////////
/////////////////////////////////////////////////////////////////////////////////

	void copy_to_buffer(cl_command_queue *queue = NULL, cl_bool block_flag = CL_TRUE,
		int num_wait = 0, cl_event *wait = NULL, cl_event *evt = NULL)
	{
		Copy_HtoD(queue, block_flag, num_wait, wait, evt);
	}

	void copy_to_buffer_size(int size_temp, cl_command_queue *queue = NULL, 
		cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event *wait = NULL,
		cl_event *evt = NULL)
	{
		Copy_HtoD_Size(size_temp, queue, block_flag, num_wait, wait, evt);
	}

/////////////////////////////////////////////////////////////////////////////////
////////////////     Writing member buffer to member array    ///////////////////
/////////////////////////////////////////////////////////////////////////////////
	void read_from_buffer(cl_command_queue *queue = NULL,
		cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event *wait = NULL,
		cl_event *evt = NULL)
	{
		if (Host_Alloc_Flag == 0)
		{
			AllocHost();
		}

		Copy_DtoH(queue, block_flag, num_wait, wait, evt);
	}

	void read_from_buffer_size(int numel, cl_command_queue *queue = NULL, 
		cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event *wait = NULL,
		cl_event *evt = NULL)
	{
		if (Host_Alloc_Flag == 0)
		{
			AllocHost();
		}

		Copy_DtoH_Size(numel, queue, block_flag, num_wait, wait, evt);
	}


/////////////////////////////////////////////////////////////////////////////////
/////////////   Writing/reading non member array and/or non member buffer ///////
/////////////////////////////////////////////////////////////////////////////////
	void read_from_buffer_to_array(cl_mem &buf, cl_bool block_flag = CL_TRUE, int size_ = -1,
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event *wait = NULL,
		cl_event *evt = NULL)
	{
		if (queue == NULL)
			queue = ioQue;
		if (size_ == -1)
			size_ = FullSize;

		int status = clEnqueueReadBuffer(*queue, buf, block_flag, 0, size_*sizeof(T), Array,
			num_wait, wait, evt);

		CHECK_ARRAY_ERROR(status, "Error reading non-member buffer into array", status);
	}

	void write_array_to_buffer(cl_mem &buf, cl_bool block_flag = CL_TRUE, int size_ = -1,
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event *wait = NULL,
		cl_event *evt = NULL)
	{
		if (queue == NULL)
			queue = ioQue;
		if (size_ == -1)
			size_ = FullSize;

		int status = clEnqueueWriteBuffer(*queue, buf, block_flag, 0,
			size_*sizeof(T), Array, num_wait, wait, evt);
		CHECK_ARRAY_ERROR(status, "Error copying array to non member buffer", status);
	}



	void write_array_to_buffer(T *Array_, cl_bool block_flag = CL_TRUE, int size_ = -1,
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event *wait = NULL, 
		cl_event *evt = NULL)
	{
		if (queue == NULL)
			queue = ioQue;
		if (size_ == -1)
			size_ = FullSize;

		int status = clEnqueueWriteBuffer(*queue, Buffer, block_flag, 0, 
			size_*sizeof(T), Array_, num_wait, wait, evt);
		CHECK_ARRAY_ERROR(status, "Error copying non-member array to device", status);
	}

	void read_from_buffer_to_array(T *Array_, cl_bool block_flag = CL_TRUE, int size_ = -1,
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event *wait = NULL,
		cl_event *evt = NULL)
	{
		if (queue == NULL)
			queue = ioQue;
		if (size_ == -1)
			size_ = FullSize;

		int status = clEnqueueReadBuffer(*queue, Buffer, block_flag, 0, size_*sizeof(T), Array_,
			num_wait, wait, evt);

		CHECK_ARRAY_ERROR(status, "Error reading buffer into non-member array", status);
	}

/////////////////////////////////////////////////////////////////////////////////
///////////////// Copying between buffer member and non member///////////////////
/////////////////////////////////////////////////////////////////////////////////
	void enqueue_copy_to_buffer_blocking(cl_mem src_buf, int copy_length = -1, 
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event* wait = NULL)
	{
		cl_event copyEvent;
		enqueue_copy_buffer_to_buffer(src_buf, Buffer, copy_length, queue,
			num_wait, wait, &copyEvent);
		clWaitForEvents(1, &copyEvent);
	}

	void enqueue_copy_from_buffer_blocking(cl_mem dest_buf, int copy_length = -1, 
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event* wait = NULL)
	{
		cl_event copyEvent;
		enqueue_copy_buffer_to_buffer(Buffer, dest_buf, copy_length, queue, num_wait, wait, &copyEvent)
			clWaitForEvents(1, &copyEvent);
		clWaitForEvents(1, &copyEvent);
	}
	

	void enqueue_copy_to_buffer(cl_mem src_buf, int copy_length = -1, 
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event* wait = NULL,
		cl_event *blk_evt = NULL)
	{
		enqueue_copy_buffer_to_buffer(src_buf, Buffer, copy_length, queue, num_wait, wait, blk_evt);
	}

	void enqueue_copy_from_buffer(cl_mem dest_buf, int copy_length = -1, 
		cl_command_queue *queue = NULL, int num_wait = 0, cl_event* wait = NULL,
		cl_event *blk_evt = NULL)
	{
		enqueue_copy_buffer_to_buffer(Buffer, dest_buf, copy_length, queue, num_wait, wait, blk_evt);
	}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////                                               /////////////////
/////////////////    Functions reading from/writing to file     /////////////////
/////////////////                                               /////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// TODO: Make sure that any function appending FileName does not receive FileName
// passed by address
/////////////////////////////////////////////////////////////////////////////////
///////// Generic Save/Load functions which read/write text and binaries ////////
/////////////////////////////////////////////////////////////////////////////////
	bool save(std::string FileName = "")
	{
		if (FileName.length() == 0)
			FileName = Name;
		return savetxt(FileName) && ArrayBase<T>::savebin(FileName);
	}

	bool save_from_device(std::string FileName = "")
	{
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}

		read_from_buffer();
		bool ret = save(FileName);
		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}


	//SFINAE to avoid trying to read structure arrays from text files
	//https://www.fluentcpp.com/2018/05/15/make-sfinae-pretty-1-what-value-sfinae-brings-to-code/
	template<typename T_ = T>
	bool load(std::string FileName = "", typename std::enable_if< !std::is_class<T_>::value, std::nullptr_t >::type = nullptr)
	{
		return !((ArrayBase<T>::loadbin(FileName) == false) && (loadtxt(FileName) == false));
	}

	template<typename T_ = T>
	bool load(std::string FileName = "", typename std::enable_if< std::is_class<T_>::value, std::nullptr_t >::type = nullptr)
	{
		return !((loadbin(FileName) == false));
	}

	// Can not load structs using plain text (with generic programming and ease atleast, so this will cause the function to
	// always return false for arrays of structs (will still be able to load binaries for structures)



/////////////////////////////////////////////////////////////////////////////////
///////// Save text functions, which call save2file (or variation on it) ////////
/////////  (save2file is imlemented in derived classes i.e pure virtual) ////////
/////////////////////////////////////////////////////////////////////////////////

	bool savetxt(std::string FileName = "")
	{
		return save2file(FileName);
	}

	bool savetxt_full(std::string FileName = "")
	{
		return save2file_full(FileName);
	}

	bool save_txt_from_device(std::string FileName = "", cl_command_queue *queue = NULL)
	{ // savetxt will convert to Name if FileName = "", so no need to convert it here
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}

		read_from_buffer(queue);
		bool ret = savetxt(FileName);

		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}

	bool save_txt_from_device_full(std::string FileName = "", cl_command_queue *queue = NULL)
	{// savetxt will convert to Name if FileName = "", so no need to convert it here
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}

		read_from_buffer(queue);
		bool ret = savetxt_full(FileName);

		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}

	virtual bool save2file(std::string FileName = "") = 0;
	
	virtual bool save2file_full(std::string FileName = "")
	{
		return save2file(FileName);
	}


/////////////////////////////////////////////////////////////////////////////////
///////// Save Binary functions, (not virtual as it works for all arrays) ///////
/////////////////////////////////////////////////////////////////////////////////
	bool savebin(std::string FileName = "")
	{//appends .bin to end of file name and directory to beginning
#ifndef CREATE_BIN_FILES
		return true;
#endif //CREATE_BIN_FILES
		if (FileName.length() == 0)
			FileName = Name;
		std::string Buf = BIN_SAVE_DIR;
		Buf.append(FileName + ".bin");

		std::fstream stream;
		if (fileopen(stream, Buf, BinaryOut) == false)
			return false;
		stream.write((char*)Array, FullSize*sizeof(T));
		stream.close();
		return true;
	};

	bool savebin_ignore_flag(std::string FileName = "")
	{//appends .bin to end of file name and directory to beginning
		if (FileName.length() == 0)
			FileName = Name;
		std::string Buf = BIN_SAVE_DIR;
		Buf.append(FileName + ".bin");

		std::fstream stream;
		if (fileopen(stream, Buf, BinaryOut) == false)
			return false;
		stream.write((char*)Array, FullSize*sizeof(T));
		stream.close();
		return true;
	}

	bool save_bin_from_device(std::string FileName = "", cl_command_queue *queue = NULL)
	{ //calls savebin, which will set FileName correctly if defaults to ""
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}

		read_from_buffer(queue);
		bool ret = savebin(FileName);
		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}

	bool save_bin_from_device_ignore_flag(std::string FileName = "", cl_command_queue *queue = NULL)
	{ //calls savebin, which will set FileName correctly if defaults to ""
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}
		if (FileName.length() == 0)
			FileName = Name;
		read_from_buffer(queue);
		bool ret = savebin_ignore_flag(FileName);
		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}

// TODO: Make sure that saving only FullSize is sufficient for all arrays
//		  and that it wont affect arrays using padding
/////////////////////////////////////////////////////////////////////////////////
///////// Save Binary functions, (not virtual as it works for all arrays) ///////
/////////////////////////////////////////////////////////////////////////////////
	bool loadbin(std::string FileName = "")
	{//appends .bin to end of file name
		if (FileName.length() == 0)
			FileName = Name;
		std::string Buf = BIN_SAVE_DIR;
		Buf.append(FileName + ".bin");

		std::fstream stream;
		if (fileopen(stream, Buf, BinaryIn) == false)
			return false;
		stream.read((char*)Array, FullSize*sizeof(T));
		stream.close();
		return true;
	};

	virtual bool loadtxt(std::string FileName = "") = 0;

// TODO: test all binary read/write/count functions
/////////////////////////////////////////////////////////////////////////////////
/////////             Count functions to determine array size             ///////
/////////////////////////////////////////////////////////////////////////////////
	int count(std::string FileName)
	{
		int counter = countbin(FileName);
		if (counter > 0) return counter;
		return counttxt(FileName);
	}

	int countbin(std::string FileName)
	{
		if (FileName.length() == 0)
			FileName = Name;
		std::string Buf = BIN_SAVE_DIR;
		Buf.append(FileName + ".bin");


		int length = 0;
		std::fstream stream;
		if (fileopen(stream, Buf, BinaryIn) == true)
		{
			// get length of file:
			stream.seekg(0, stream.end);
			length = stream.tellg();
			length /= sizeof(T);
			stream.close();
		}

		return length;
	}

	virtual int counttxt(std::string FileName = "")
	{ 
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");

		int counter = 0;
		double Buf;
				
		std::fstream stream;
		if (fileopen(stream, FileName, FileIn) == false)
			return false;

		while (stream >> Buf)
			counter++;

		stream.close();
		return counter;
	};

/////////////////////////////////////////////////////////////////////////////////
/////////      Functions for renaming/cleaning/copying output files       ///////
/////////////////////////////////////////////////////////////////////////////////

	

	
	void RenameTxtFile(std::string SourceName = "")
	{
		if (SourceName.length() == 0)
			SourceName = Name + ".txt";
		std::string DestFileName = curDumpSaveDir + SLASH + SourceName;
		RenameFile(SourceName, DestFileName);
	}

	void RenameBinFile(std::string SourceName = "")
	{
		if (SourceName.length() == 0)
			SourceName = Name + ".bin";
		std::string DestFileName = curDumpSaveDir + SLASH + SourceName;
		RenameFile(SourceName, DestFileName);
	}


	void CleanFile(std::string SourceName = "")
	{
		int maxAttempt = 100;
		if (SourceName.length() == 0)
			SourceName = Name;
		std::fstream stream;
		if(fileopen(stream, SourceName, FileOut, maxAttempt))
		{
			stream.close();
		}
		else
		{
			LOGMESSAGE("could not clean file " + SourceName + ", so this file will retain old information");
		}
	}


protected:
	void delay(double sec)
	{
		time_t StartTime, CurrTime;
		time(&StartTime);
		do
		{
			time(&CurrTime);
		} while (difftime(CurrTime, StartTime) < sec);
	};

	bool fileopen(std::fstream &stream, std::string &FileName, FileMode Mode)
	{
		if (Mode == BinaryIn || Mode == FileIn)
			return fileopen(stream, FileName, Mode, 1);
		else
			return fileopen(stream, FileName, Mode, 100);
	};

	bool fileopen(std::fstream &stream, std::string &FileName, FileMode mode_, const int att)
	{
		for (int i = 0; i < att; i++)
		{
			switch (mode_)
			{
			case FileIn:
			{
				stream.open(FileName, std::ios_base::in);
				break;
			}
			case FileOut:
			{
				stream.open(FileName, std::ios_base::out);
				break;
			}
			case FileAppend:
			{
				stream.open(FileName, std::ios_base::out | std::ios_base::app);
				break;
			}
			case BinaryIn:
			{
				stream.open(FileName, std::ios_base::in | std::ios_base::binary);
				break;
			}
			case BinaryOut:
			{
				stream.open(FileName, std::ios_base::out | std::ios_base::binary);
				break;
			}
			}

			if (stream.is_open())
				return true;
			delay(0.1);
		}
		return false;
	};
};




template <class T> class Array1D : public ArrayBase <T>
{
private:
	void set_excess_to_zero()
	{
		if (SizeX == FullSize)
			return;

		size_t extra = (FullSize - SizeX);
		memset(&Array[SizeX], 0, extra*sizeof(T));
	}

protected:
	void reallocate_host(const int n1)
	{
		T *p = new T[FullSize];
		memcpy(p, ArrayBase<T>::Array, FullSize*sizeof(T));
		int oldsize = ArrayBase<T>::FullSize;

		oldsize = MIN(oldsize, n1);

		delete[] ArrayBase::Array;
		ArrayBase<T>::Array = new T[n1];
		ArrayBase<T>::Host_Alloc_Flag = 1;
		memcpy(ArrayBase<T>::Array, p, oldsize*sizeof(T));
		delete[] p;
	}
	
	//memory reallocation on device
	void reallocate_device(const int n1)
	{
		T *p = new T[FullSize];
		read_from_buffer_to_array(p, CL_TRUE, FullSize);
		FreeDevice();

		int oldsize = ArrayBase<T>::FullSize;
		oldsize = MIN(oldsize, n1);

		allocate_buffer_size(n1);
		write_array_to_buffer(p, CL_TRUE, n1);
		delete[] p;
	}

	void testIndex(const int n1)
	{
		using std::to_string;
		CHECK_ARRAY_ERROR((0 > n1 || n1 >= SizeX),
			"Error in 1D array! Index n1=" + to_string(n1) +\
			" out of range [0..." + to_string(SizeX) + "]\n",
			ERROR_MEMORY_OUT_OF_BOUNDS);
		CHECK_ARRAY_ERROR((n1 >= FullSize),
			"Error in 1D array!", ERROR_MEMORY_OUT_OF_BOUNDS);
	}

public:

	Array1D(const Array1D<T>& other) : ArrayBase<T>(other)
	{}

	Array1D(std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);
	}

	Array1D(const int n1, std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);

		SizeX = n1;
		FullSize = n1;
		FullSizeX = n1;

		SizeY = 1;
		SizeZ = 1;
		FullSizeY = 1;
		FullSizeZ = 1;

		ArrayBase<T>::AllocHost();
	}

	Array1D(const int n1, const int fs, std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);

		SizeX = n1;
		FullSize = fs;
		FullSizeX = fs;

		SizeY = 1;
		SizeZ = 1;
		FullSizeY = 1;
		FullSizeZ = 1;

		ArrayBase<T>::AllocHost();
	}

	Array1D& operator=(ArrayBase<T> &other)
	{
		SizeX = other.getSizeX();
		FullSizeX = other.getFullSizeX();
		SizeY = other.getSizeY();
		FullSizeY = other.getFullSizeY();
		SizeZ = other.getSizeZ();
		FullSizeZ = other.getFullSizeZ();
		FullSize = other.getFullSize();
		Host_Alloc_Flag = other.getHostAllocFlag();
		Device_Alloc_Flag = other.getDeviceAllocFlag();
		clFlags = other.getclFlags();
		Name = other.getName();
		BufferFullSize = other.getBufferFullSize();
		ioQue = clEnv::instance()->getIOqueue();
		//context = clEnv::instance()->getContext();
		Array = nullptr;
		Buffer = nullptr;
		if (Host_Alloc_Flag)
		{
			AllocHost();
			memcpy(Array, other.get_array(), FullSize*sizeof(T));
		}
		if (Device_Alloc_Flag)
		{
			allocate_buffer_size(BufferFullSize);
			copy_to_buffer();
		}
		return *this;
	}

	Array1D& operator=(ArrayBase *other)
	{
		SizeX = other->getSizeX();
		FullSizeX = other->getFullSizeX();
		SizeY = other->getSizeY();
		FullSizeY = other->getFullSizeY();
		SizeZ = other->getSizeZ();
		FullSizeZ = other->getFullSizeZ();
		FullSize = other->getFullSize();
		Host_Alloc_Flag = other->getHostAllocFlag();
		Device_Alloc_Flag = other->getDeviceAllocFlag();
		clFlags = other->getclFlags();
		Name = other->getName();
		BufferFullSize = other->getBufferFullSize();
		ioQue = clEnv::instance()->getIOqueue();
		//context = clEnv::instance()->getContext();
		Array = nullptr;
		Buffer = nullptr;
		if (Host_Alloc_Flag)
		{
			AllocHost();
			memcpy(Array, other->get_array(), FullSize*sizeof(T));
		}
		if (Device_Alloc_Flag)
		{
			allocate_buffer_size(BufferFullSize);
			copy_to_buffer();
		}
		return *this;
	}

	Array1D& operator=(Array1D *other)
	{
		ArrayBase<T>* arrtemp = static_cast<ArrayBase<T>*>(other);
		this->operator=(arrtemp);
		return *this;
	}
	
	Array1D& operator=(Array2D<T> *other)
	{
		ArrayBase<T>* arrtemp = static_cast<ArrayBase<T>*>(other);
		this->operator=(arrtemp);

		// need to fold y dimension into x dimension 
		// (FullSize and BufferFullSize are unchanged)
		this->SizeX = this->FullSizeX*this->SizeY;
		this->FullSizeX *= this->FullSizeY;
		this->SizeY = 1;
		this->FullSizeY = 1;
		return *this;
	}



	T& operator[](const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif 
		return *(ArrayBase<T>::Array + n1);
	}

	T& operator()(const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		return *(ArrayBase<T>::Array + n1);
	}

	T getEl(const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		return *(ArrayBase<T>::Array + n1);
	}

	void setEl(const T val, const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		*(ArrayBase<T>::Array + n1) = val;
	}

	//memory allocation
	void allocate(const int n1)
	{
		ArrayBase<T>::FreeHost();

		SizeX = n1;
		FullSize = n1;
		FullSizeX = n1;

		SizeY = 1;
		SizeZ = 1;
		FullSizeY = 1;
		FullSizeZ = 1;

		ArrayBase<T>::AllocHost();
	}

	void iniFromStdVec(std::vector<T> &stvec_)
	{
		allocate(stvec_.size());
		for (int i = 0; i < SizeX; i++)
			setEl(stvec_[i], i);
	}

	void exportToStdVec(std::vector<T> &stvec_)
	{
		if (stvec_.size() != FullSize)
			stvec_.resize(FullSize);
		std::memcpy(stvec_.data(), this->Array, FullSize*sizeof(T));
	}



	void allocate(const int n1, const int fs)
	{
		ArrayBase<T>::FreeHost();

		SizeX = n1;
		FullSize = fs;
		FullSizeX = fs;

		SizeY = 1;
		SizeZ = 1;
		FullSizeY = 1;
		FullSizeZ = 1;

		ArrayBase<T>::AllocHost();
	}

	void zeros(const int n1)
	{
		allocate(n1);
		ArrayBase<T>::zeros();
	}

	void zeros(const int n1, const int fs1)
	{
		allocate(n1, fs1);
		ArrayBase<T>::zeros();
	}

	void fillas2D_nopad(T filval, const int nx_, const int fx_, const int ny_)
	{
		for (int i = 0; i < nx_; i++)
		{
			for (int j = 0; j < ny_; j++)
			{
				setEl(filval, i + j*fx_);
			}
		}
	}


	//	maintains length of array at fsize, but changes 
	//	SizeX and sets excess to 0. (called afer reallocate)
	void reset_sizes(int actsize, int fsize)
	{
		CHECK_ARRAY_ERROR((actsize > fsize || fsize != FullSize),
			"Error in reset_sizes of Array1D", ERROR_MEMORY_OUT_OF_BOUNDS);

		FullSizeX = fsize;
		SizeX = actsize;
		FullSize = fsize;

		set_excess_to_zero();
	}

	void reset_sizeX()
	{
		SizeX = FullSize;
	}


	void reallocate(const int n1)
	{
		if (Host_Alloc_Flag)
		{
			reallocate_host(n1);
		}

		if (Device_Alloc_Flag)
			reallocate_device(n1);

		SizeX = n1;
		FullSize = n1;
		FullSizeX = n1;
	}

	void reallocate_host_only(const int n1)
	{
		if (n1 == FullSize)
			return;
	
		if (Host_Alloc_Flag)
			reallocate_host(n1);

		SizeX = n1;
		FullSize = n1;
		FullSizeX = n1;
	}

	void reallocate_device_only(const int n1)
	{
		if (Device_Alloc_Flag)
			reallocate_device(n1);

		SizeX = n1;
		FullSize = n1;
		FullSizeX = n1;
	}

	
	bool savetxt_as_2D(const int Xsize_, const int XsizeFull_, const int Ysize_, 
		std::string FileName = "")
	{
		return save2file_as_2D(Xsize_, XsizeFull_, Ysize_, FileName);
	}
	
	bool save_txt_from_device_as_2D(const int Xsize_, const int XsizeFull_, const int Ysize_, 
		std::string FileName = "", cl_command_queue *queue = NULL)
	{
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}
		read_from_buffer(queue);
		bool ret = savetxt_as_2D(Xsize_, XsizeFull_, Ysize_, FileName);
		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}
	
	bool savetxt_w_skip(int skip_val, std::string FileName = "")
	{
		return save2file_w_skip(Buf, skip_val);
	}


	bool save2file(std::string FileName = "")
	{
		int i;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");

		std::fstream stream;
		if (fileopen(stream, FileName, FileOut) == false)
			return false;

		for (i = 0; i < ArrayBase<T>::SizeX; i++)
		{
			stream << getEl(i) << "\n";
		}
		stream.close();
		return true;
	}


	virtual bool append2file_from_device(int toind_ = -1, std::string FileName = "",
		cl_command_queue *queue = NULL)
	{
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}

		read_from_buffer(queue);
		bool ret = append2file(toind_, FileName);

		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}

	// Kind of hacky, but passing -2 will overwrite existing file
	// with an empty one.
	virtual bool append2file(int toind_ = -1, std::string FileName = "")
	{
		if (toind_ == -1)
			toind_ = SizeX;
		else if (toind_ == -2)
		{
			if (FileName.length() == 0)
				FileName = Name;
			FileName.append(".txt");

			std::fstream stream;
			if (fileopen(stream, FileName, FileOut) == false)
				return false;
			stream.close();

			return true;
		}


		int i;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");

		std::fstream stream;
		if (fileopen(stream, FileName, FileAppend) == false)
			return false;

		for (i = 0; i < toind_; i++)
		{
			stream << getEl(i) << "\n";
		}
		stream.close();
		return true;
	}


	bool save2file_full(std::string FileName = "")
	{
		int i;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");
		std::fstream stream;
		if (fileopen(stream, FileName, FileOut) == false)
			return false;

		for (i = 0; i < ArrayBase<T>::SizeX; i++)
		{
			stream << getEl(i) << "\n";
		}
		stream.close();
		return true;
	}


	bool save2file_as_2D(const int Xsize_, const int XsizeFull_, const int Ysize_, std::string FileName = "")
	{
		int i, j;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");
		std::fstream stream;
		if (fileopen(stream, FileName, FileOut) == false)
			return false;

		for (i = 0; i < Xsize_; i++)
		{
			for (j = 0; j < Ysize_; j++)
			{
				stream << getEl(i + XsizeFull_*j) << "\t";
			}
			stream << "\n";
		}
		stream.close();
		return true;
	}

	// Could have been a function that I used for debugging, but is no
	// longer needed
	// TODO: check if this is still used and remove if it isnt.
	bool save2file_w_skip(int skip, std::string FileName = "")
	{
		int i;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");
		std::fstream stream;
		if (fileopen(stream, FileName, FileOut) == false)
			return false;

		for (i = 0; i < ArrayBase<T>::SizeX; i++)
		{
			stream << getEl(i*skip) << "\n";
		}
		stream.close();
		return true;
	}


	virtual bool loadtxt(std::string FileName = "") override
	{
		if (FileName.length() == 0)
			FileName = Name;
		
		FileName.append(".txt");
		int i;
		T Buf;

		std::fstream stream;
		if (fileopen(stream, FileName, FileIn) == false)
			return false;
		for (i = 0; i < ArrayBase<T>::SizeX; i++)
		{	
			stream >> Buf;
			setEl(Buf, i);
		}
		stream.close();
		return true;
	}


	virtual bool loadtxtnew(std::string FileName = "")
	{
		zeros(counttxt(FileName));
		return loadtxt(FileName);
	}

	bool loadbinnew(std::string FileName = "")
	{
		zeros(countbin(FileName));
		return loadbin(FileName);
	}

	bool loadnew(std::string FileName = "")
	{
		int counter = countbin(FileName);
		if (counter > 0)
		{
			zeros(counter);
			return loadbin(FileName);
		}
		return loadtxtnew(FileName);
	}



};


typedef  Array1D <double> Array1Dd;
typedef  Array1D <int> Array1Di;
typedef  Array1D <signed char> Array1Dc;
typedef  Array1D <float> Array1Df;
typedef  Array1D <unsigned int> Array1Du;
typedef  Array1D <cl_short> Array1Ds;

//typedef  Array1D <par> Array1DP;
//typedef  Array1D <bLinks> Array1DBL;
//typedef  Array1D <Pparam> Array1DPP;
//typedef  Array1D <rampI> Array1DRI;
//typedef  Array1D <foulI> Array1DFI;


// Variable size array (use append to add element to array)
// will increase array size by set multiple (resizeMult) of wgSize
// when element is appended to end of current array size (append function
// will return true if array is appended to inform that buffer needs to be
// reallocated as this is not done when host array is reallocated to avoid
// multiple reallocations of buffer).

template <typename T>
class DynamicArray1D : public Array1D<T>
{
protected:
	int wgSize;
	int curEl;
	int resizeMult;

public:
	DynamicArray1D(int inisize, std::string name_ = "default", const int wgsize_ = 128,
		const int resizem_ = RESIZE_MULTIPLIER) : 
		Array1D<T>(name_), wgSize(wgsize_), curEl(0), resizeMult(resizem_)
	{
		if (inisize == 0)
			inisize = wgSize;

		getGlobalSizeMacro(inisize, wgSize);

		Array1D<T>::zeros(inisize);
	}

	DynamicArray1D(std::string name_ = "default", const int wgsize_ = 128,
		const int resizem_ = RESIZE_MULTIPLIER) :
		Array1D<T>(name_), wgSize(wgsize_), curEl(0), resizeMult(resizem_)
	{
	}

	void allocateDynamic(const int inisize)
	{
		getGlobalSizeMacro(inisize, wgSize);
		Array1D<T>::zeros(inisize);
	}

	bool append(const T val_)
	{
		Array1D<T>::operator()(curEl) = val_;

		curEl++;
		if (curEl < FullSize)
			return false;

		int newSize = FullSize + resizeMult*wgSize;
		reallocate_host(newSize);

		SizeX = newSize;
		FullSize = newSize;
		FullSizeX = newSize;
		return true;
	}

	void zeros()
	{
		curEl = 0;
		ArrayBase<T>::zeros();
	}

	void reset()
	{
		curEl = 0;
	}

	void copy_dynamic_to_buffer(cl_command_queue* queue = NULL,
		cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event* wait = NULL,
		cl_event* evt = NULL)
	{
		ArrayBase<T>::copy_to_buffer_size(curEl, queue, block_flag, num_wait, wait, evt);
	}

	// memory reallocation on device, does not save existing data
	// on device to copy over after reallocation
	void reallocate_device_dynamic()
	{
		FreeDevice();
		allocate_buffer_size(FullSize);
	}

	// memory reallocation on device, does not save existing data
	// on device to copy over after reallocation
	void resize_device_dynamic()
	{
		FreeDevice();
		allocate_buffer_size(BufferFullSize + resizeMult * wgSize);
	}

	// memory reallocation on device, which copies current host
	// contents to newly reallocated buffer rather than copying 
	// values held previously in buffer.
	void reallocate_device_dynamic_w_copy()
	{
		FreeDevice();
		allocate_buffer_size(FullSize);
		ArrayBase<T>::copy_to_buffer_size(curEl)
	}

	int curSize() { return curEl; }
	int* curSizeAdd() { return &curEl; }
	int curFullSize() { return FullSize; }
	int* curFullSizeAdd() { return &FullSize; }
	int* curBufferSizeAdd() { return &BufferFullSize; }

};

typedef  DynamicArray1D <double> DynArray1Dd;
typedef  DynamicArray1D <int> DynArray1Di;
typedef  DynamicArray1D <signed char> DynArray1Dc;
typedef  DynamicArray1D <float> DynArray1Df;
typedef  DynamicArray1D <unsigned int> DynArray1Du;
typedef  DynamicArray1D <cl_short> DynArray1Ds;


template <typename T, typename Tbase, int VecSize_> 
class Array1Dv : public Array1D <T>
{
private:
	void set_excess_to_zero()
	{
		if (SizeX == FullSize)
			return;
		T zero_ = { (Tbase)0 };
		for (int i = SizeX; i < FullSize; i++)
		{
			Array[i] = zero_;
		}
	}
public:
	Array1Dv(std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);
	}
	//bool loadtxtnew(std::string FileName = "") override
	//{
	//	zeros(counttxt(FileName) / 2);
	//	return loadtxt(FileName);
	//};

};

template <typename T>
void operator >> (const YAML::Node& node, Array1D<T>& arr) 
{
	for (int i = 0; i < arr.getSizeX(); i++)
		node[i] >> arr(i);
}

template <typename T>
YAML::Emitter& operator << (YAML::Emitter& out, Array1D<T>& arr)
{
	out << YAML::Flow;
	out << YAML::BeginSeq;
	for (int i = 0; i < arr.getSizeX(); i++)
		out << arr(i);
	out << YAML::EndSeq;
	return out;
}

typedef  Array1Dv <cl_short2, cl_short, 2> Array1Dv2s;
typedef  Array1Dv <cl_short3, cl_short, 3> Array1Dv3s;
typedef  Array1Dv <cl_short4, cl_short, 4> Array1Dv4s;
typedef  Array1Dv <cl_short8, cl_short, 8> Array1Dv8s;
typedef  Array1Dv <cl_short16, cl_short, 16> Array1Dv16s;

typedef  Array1Dv <cl_ushort2, cl_ushort, 2> Array1Dv2us;
typedef  Array1Dv <cl_ushort3, cl_ushort, 3> Array1Dv3us;
typedef  Array1Dv <cl_ushort4, cl_ushort, 4> Array1Dv4us;
typedef  Array1Dv <cl_ushort8, cl_ushort, 8> Array1Dv8us;
typedef  Array1Dv <cl_ushort16, cl_ushort, 16> Array1Dv16us;


typedef  Array1Dv <cl_double2, double, 2> Array1Dv2d;
typedef  Array1Dv <cl_double3, double, 3> Array1Dv3d;
typedef  Array1Dv <cl_double4, double, 4> Array1Dv4d;
typedef  Array1Dv <cl_double8, double, 8> Array1Dv8d;
typedef  Array1Dv <cl_double16, double, 16> Array1Dv16d;

typedef  Array1Dv <cl_float2, float, 2> Array1Dv2f;
typedef  Array1Dv <cl_float3, float, 3> Array1Dv3f;
typedef  Array1Dv <cl_float4, float, 4> Array1Dv4f;
typedef  Array1Dv <cl_float8, float, 8> Array1Dv8f;
typedef  Array1Dv <cl_float16, float, 16> Array1Dv16f;

typedef  Array1Dv <cl_int2, int, 2> Array1Dv2i;
typedef  Array1Dv <cl_int3, int, 3> Array1Dv3i;
typedef  Array1Dv <cl_int4, int, 4> Array1Dv4i;
typedef  Array1Dv <cl_int8, int, 8> Array1Dv8i;
typedef  Array1Dv <cl_int16, int, 16> Array1Dv16i;

typedef  Array1Dv <cl_uint2, unsigned int, 2> Array1Dv2u;
typedef  Array1Dv <cl_uint3, unsigned int, 3> Array1Dv3u;
typedef  Array1Dv <cl_uint4, unsigned int, 4> Array1Dv4u;
typedef  Array1Dv <cl_uint8, unsigned int, 8> Array1Dv8u;
typedef  Array1Dv <cl_uint16, unsigned int, 16> Array1Dv16u;

template <class T> class Array2D : public ArrayBase<T>
{
private:
	void testIndex(const int n1, const int n2)
	{
		CHECK_ARRAY_ERROR((0 > n1 || n1 >= SizeX),
			("Tried to access element " + std::to_string(n1) + 
			" in X dim of array of size [" + std::to_string(SizeX) +
			", " + std::to_string(SizeY) + "]") , ERROR_MEMORY_OUT_OF_BOUNDS);

		CHECK_ARRAY_ERROR((0 > n2 || n2 >= SizeY),
			"Tried to access element " + std::to_string(n2) +
			" in Y dim of array of size [" + std::to_string(SizeX) +
			", " + std::to_string(SizeY) + "]", ERROR_MEMORY_OUT_OF_BOUNDS);

		CHECK_ARRAY_ERROR((FullSizeY * n1 + n2 < 0 || FullSizeY * n1 + n2 >= FullSize),
			"Error in 2D array!", ERROR_MEMORY_OUT_OF_BOUNDS);
	}

protected:
	void reallocate_host(const int n1)
	{
		T *p = new T[FullSize];
		memcpy(p, ArrayBase<T>::Array, FullSize*sizeof(T));
		int oldsize = ArrayBase<T>::FullSize;

		oldsize = MIN(oldsize, n1*FullSizeY);

		delete[] ArrayBase::Array;
		ArrayBase<T>::Array = new T[n1*FullSizeY];
		ArrayBase<T>::Host_Alloc_Flag = 1;
		memcpy(ArrayBase<T>::Array, p, oldsize*sizeof(T));
		delete[] p;
	}




	//memory reallocation on device
	void reallocate_device(const int n1)
	{
		T *p = new T[FullSize];
		read_from_buffer_to_array(p, CL_TRUE, FullSize);
		FreeDevice();

		int oldsize = ArrayBase<T>::FullSize;
		oldsize = MIN(oldsize, n1*FullSizeY);

		allocate_buffer_size(n1*FullSizeY);
		write_array_to_buffer(p, CL_TRUE, n1*FullSizeY);
		delete[] p;
	}

public:

	Array2D()
	{
		ArrayBase<T>::ini();
	}

	Array2D(std::string name_)
	{
		ArrayBase<T>::ini(name_);
	}

	Array2D(const int n1, const int n2, std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);

		SizeX = n1;
		FullSizeX = n1;

		SizeY = n2;
		FullSizeY = n2;

		SizeZ = 1;
		FullSizeZ = 1;

		FullSize = FullSizeX*FullSizeY*FullSizeZ;

		ArrayBase<T>::AllocHost();
	}

	Array2D(const int n1, const int fs1, const int n2, const int fs2, std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);

		SizeX = n1;
		FullSizeX = fs1;

		SizeY = n2;
		FullSizeY = fs2;

		SizeZ = 1;
		FullSizeZ = 1;

		FullSize = FullSizeX*FullSizeY*FullSizeZ;

		ArrayBase<T>::AllocHost();
	}

	virtual ~Array2D()
	{
		ArrayBase<T>::FreeHost();
		ArrayBase<T>::FreeDevice();
	}

	T& operator()(const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif 
		return *(ArrayBase<T>::Array + n1 + FullSizeX * n2);
	}

	T& operator[](const int n1)
	{
		return *(ArrayBase<T>::Array + n1);
	}

	T& operator()(cl_int2 n)
	{
		return operator()(n.x, n.y);
	}

	T getEl(const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif
		return *(ArrayBase<T>::Array + n1 + FullSizeX * n2);
	}

	void setEl(const T val, const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif
		*(ArrayBase<T>::Array + n1 + FullSizeX * n2) = val;
	}

	void reallocate_xdim(const int n1)
	{
		if (Host_Alloc_Flag)
		{
			reallocate_host(n1);
		}

		if (Device_Alloc_Flag)
			reallocate_device(n1);

		SizeX = n1;
		FullSizeX = n1;

		FullSize = FullSizeX*FullSizeY*FullSizeZ;
	}

	//memory allocation
	void allocate(const int n1, const int n2)
	{
		ArrayBase<T>::ini();

		SizeX = n1;
		FullSizeX = n1;

		SizeY = n2;
		FullSizeY = n2;

		SizeZ = 1;
		FullSizeZ = 1;

		FullSize = FullSizeX*FullSizeY*FullSizeZ;

		ArrayBase<T>::AllocHost();
	}

	void allocate(cl_int2 n)
	{
		allocate(n.x, n.y);
	}

	void allocate(const int n1, const int fs1, const int n2, const int fs2)
	{
		ArrayBase<T>::ini();

		SizeX = n1;
		FullSizeX = fs1;

		SizeY = n2;
		FullSizeY = fs2;

		SizeZ = 1;
		FullSizeZ = 1;

		FullSize = FullSizeX*FullSizeY*FullSizeZ;

		ArrayBase<T>::AllocHost();
	}

	void zeros(const int n1, const int n2)
	{
		allocate(n1, n2);
		ArrayBase<T>::zeros();
	}

	void zeros(cl_int2 n)
	{
		allocate(n);
		ArrayBase<T>::zeros();
	}

	void zeros(const int n1, const int fs1, const int n2, const int fs2)
	{
		allocate(n1, fs1, n2, fs2);
		ArrayBase<T>::zeros();
	}

	void fill_noBuffer(const T val)
	{
		for (int i = 0; i < SizeX; i++)
		{
			for (int j = 0; j < SizeY; j++)
			{
				setEl(val, i, j);
			}
		}
	}

	void reallocate_host(const int n1, const int n2)
	{
		T *p = new T[FullSize];
		memcpy(p, ArrayBase<T>::Array, FullSize*sizeof(T));
		int oldsize = ArrayBase<T>::FullSize;

		delete[] ArrayBase::Array;
		ArrayBase<T>::Array = new T[n1*n2];
		ArrayBase<T>::Host_Alloc_Flag = 1;

		int Xlim = MIN(SizeX, n1);
		int Ylim = MIN(SizeY, n2);

		for (int i = 0; i < Xlim; i++)
		{
			for (int j = 0; j < Ylim; j++)
			{
				Array[i*n2 + j] = p[i*FullSizeX + j];
			}
		}

		delete[] p; 

		SizeX = n1;
		SizeY = n2;
		FullSizeX = n1;
		FullSizeY = n2;
		FullSize = n1*n2;
	};

	virtual bool save2file(std::string FileName = "")
	{
		int i,j;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");
		std::fstream stream;
		if (fileopen(stream, FileName, FileOut) == false)
			return false;

		for (i = 0; i < ArrayBase<T>::SizeX; i++)
		{
			for (j = 0; j < ArrayBase<T>::SizeY; j++)
			{
				stream << getEl(i, j) << "\t";
			}
			stream << "\n";
		}
		stream.close();
		return true;
	}

	virtual bool savetxt_as_1D(std::string FileName = "")
	{
		int i, j;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");
		std::fstream stream;
		if (fileopen(stream, FileName, FileOut) == false)
			return false;

		for (i = 0; i < ArrayBase<T>::FullSizeX; i++)
		{
			for (j = 0; j < ArrayBase<T>::SizeY; j++)
			{
				stream << std::setw(16) << Array[i + j*FullSizeX] << "\n";
			}
			stream << "\n";
		}
		stream.close();
		return true;
	}

	bool save_txt_from_device_as_1D(std::string FileName = "", cl_command_queue *queue = NULL)
	{ // savetxt will convert to Name if FileName = "", so no need to convert it here
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}

		read_from_buffer(queue);
		bool ret = savetxt_as_1D(FileName);

		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}
	
	virtual bool save2file_full(std::string FileName = "")
	{
		int i,j;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");
		std::fstream stream;
		if (fileopen(stream, FileName, FileOut) == false)
			return false;

		for (i = 0; i < ArrayBase<T>::SizeX; i++)
		{
			for (j = 0; j < ArrayBase<T>::SizeY; j++)
			{
				stream << getEl(i, j) << "\t";
			}
			stream << "\n";
		}
		stream.close();
		return true;
	}

	bool loadtxt(std::string FileName = "")
	{
		int i, j;
		T Buf;
		if (FileName.length() == 0)
			FileName = Name;

		FileName.append(".txt");

		std::fstream stream;
		if (fileopen(stream, FileName, FileIn) == false)
			return false;

		for (i = 0; i < ArrayBase<T>::SizeX; i++)
		{
			for (j = 0; j < ArrayBase<T>::SizeY; j++)
			{
				stream >> Buf;
				setEl(Buf, i, j);
			}
		}
		return true;
	}


	virtual bool loadtxtnew(const int s2, std::string FileName = "")
	{
		zeros(counttxt(FileName) / s2, s2);
		return loadtxt(FileName);
	}

	bool loadbinnew(const int s2, std::string FileName = "")
	{
		zeros(countbin(FileName) / s2, s2);
		return loadbin(FileName);
	}

	bool loadnew(const int s2, std::string FileName = "")
	{
		int counter = ArrayBase<T>::countbin(FileName);
		if (counter > 0)
		{
			zeros(counter / s2, s2);
			return ArrayBase<T>::loadbin(FileName);
		}
		return loadtxtnew(s2, FileName);
	}
};

typedef  Array2D <double> Array2Dd;
typedef  Array2D <int> Array2Di;
typedef  Array2D <signed char> Array2Dc;
typedef  Array2D <float> Array2Df;
typedef  Array2D <unsigned int> Array2Du;

// This is used for saving time data sets (e.g. shear stress,
// averaged values, Nu, etc), so it is specifically written
// for these arrays.
template <class T> class TimeData : public Array2D<T>
{
private:
	bool fileCreated;
	int* timeArr;
	int curInd;

	// If data is processed before being written, 
	// saveDataFunc needs to be set to point to a
	// processing and output function. This function
	// needs to take the file stream and index as 
	// arguments, and calculate and write processed
	// values to the file.
	void (*saveDataFunc)(TimeData<T>*, std::fstream&, int);



public:

	TimeData()
	{
		// If this constructor is called, the name must be set separately before
		// the file can be initialized or an error will be thrown.
		ArrayBase<T>::ini("defaultName");
		fileCreated = false;
		saveDataFunc = &TimeData::genericSaveFunc;

	}

	TimeData(std::string name_)
	{
		ArrayBase<T>::ini(name_);
		fileCreated = false;
		saveDataFunc = &TimeData::genericSaveFunc;
	}

	virtual ~TimeData()
	{
		delete[] timeArr;
		timeArr = nullptr;
		
		ArrayBase<T>::FreeHost();
		ArrayBase<T>::FreeDevice();
	}

	void iniFile(bool restartRun_)
	{
		bool testName = Name.compare("defaultName") == 0;

		CHECK_ARRAY_ERROR(testName, "TimeData arrays must be constructed with a\n"\
			"string containing a name, or have the name set separately before\n"\
			"iniFile can be called", ERROR_IN_TIMEDATA_ARRAY);


		std::string FileName = Name + ".txt";
		

		std::string logMess;

		if (restartRun_)
		{
			logMess = "Not initializing new file for " + Name + " array, because "\
				"run is being restarted. Warning: This will result in error if "\
				"file doesn't currently exist.";
		}
		else
		{
			std::fstream stream;

			logMess = "Initializing file for " + Name + " array. Warning: "\
				"If restarting run, old data has been deleted";

			if (fileopen(stream, FileName, FileOut) == false)
			{
				CHECK_ARRAY_ERROR(true, "Error initializing file for "
					"saving TimeData Array", ERROR_IN_TIMEDATA_ARRAY);
			}

			stream.close();
		}


		LOGMESSAGE(logMess);
		fileCreated = true;
	}

	static void genericSaveFunc(TimeData<T>* td_, std::fstream& stream_, int j_)
	{
		for (int i = 0; i < td_->SizeX; i++)
		{
			stream_ << td_->getEl(i, j_) << "\t";
		}
	}

	void setFunction(void (*funPtr_)(TimeData<T>*, std::fstream&, int))
	{
		saveDataFunc = funPtr_;
	}

	void createTimeArray()
	{
		CHECK_ARRAY_ERROR(ArrayBase::SizeY == 0, "Array size must be set "
			"before calling createTimeArray", ERROR_IN_TIMEDATA_ARRAY);
		
		timeArr = new int[ArrayBase::SizeY];
		curInd = 0;
	}

	int getCurInd() { return curInd; }
	int* getCurIndAdd() { return &curInd; }
	// Should never be called, but including just in case
	void resetCurInd() { curInd = 0; }


	// make sure that the curInd is set as a kernel argument before 
	// calling this function.
	bool setTimeAndIncrement(int time_, cl_command_queue* queue = NULL,
		cl_bool block_flag = CL_TRUE, int num_wait = 0, cl_event* wait = NULL,
		cl_event* evt = NULL)
	{
		// Write current time to timeArr
		timeArr[curInd++] = time_;
		
		// Check if array is full and save if so
		if (curInd == ArrayBase<T>::SizeY)
		{
			appendData_from_device(queue, block_flag, num_wait, wait, evt);
			return true;
		}
		return false;
	}

	bool appendData_from_device(cl_command_queue* queue = NULL, 
		cl_bool block_flag = CL_TRUE, int num_wait = 0, 
		cl_event* wait = NULL, cl_event* evt = NULL)
	{
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}

		ArrayBase<T>::read_from_buffer_size(SizeX*curInd, queue, block_flag, num_wait, wait, evt);
		bool ret = appendData();

		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}

	bool appendData()
	{
		int j;

		CHECK_ARRAY_ERROR(!fileCreated, "Trying to append to unitialized "
			"TimeData save file", ERROR_IN_TIMEDATA_ARRAY);

		std::string FileName = Name + ".txt";
		std::fstream stream;

		// Handle file being unable to open
		if (fileopen(stream, FileName, FileAppend) == false)
		{
			std::string logMess = "Unable to open output file for " + Name;

			if (curInd < SizeY)
			{
				// If called at a save step when array has not been completely filled,
				// the data will not need to be overwritten yet, so a message is logged,
				// and the function returns false to indicate that any counters should
				// not be reset.
				logMess += ". Data will only be lost if simulation "\
					"crashes before next save attempt.";
			}
			else
			{
				// If called when the array is full, curInd is decremented so that
				// on next call to setTimeAndIncrement another attempt to open the file
				// will occur. This will result in the last timestep of saved data being
				// overwritten, but is better than throwing out all data currently in the
				// array.
				logMess += ". Some Data has been lost.";
				curInd--;
			}
			LOGMESSAGE(logMess);
			return false;
		}

		for (j = 0; j < curInd; j++)
		{
			stream << timeArr[j] << "\t";
			(*saveDataFunc)(this, stream, j);
			stream << "\n";
		}

		curInd = 0;
		stream.close();
		return true;
	}
};


// Keeping just in case I do need to specialize anything
template <typename T, typename Tbase, int VecSize_>
class Array2Dv : public Array2D <T>
{
public:
	Array2Dv(std::string name_)
	{
		ArrayBase<T>::ini(name_);
	}
	Array2Dv()
	{
	}

};


typedef  Array2Dv <cl_double2, double, 2> Array2Dv2d;
typedef  Array2Dv <cl_double3, double, 3> Array2Dv3d;
typedef  Array2Dv <cl_double4, double, 4> Array2Dv4d;
typedef  Array2Dv <cl_double8, double, 8> Array2Dv8d;
typedef  Array2Dv <cl_double16, double, 16> Array2Dv16d;

typedef  Array2Dv <cl_float2, float, 2> Array2Dv2f;
typedef  Array2Dv <cl_float3, float, 3> Array2Dv3f;
typedef  Array2Dv <cl_float4, float, 4> Array2Dv4f;
typedef  Array2Dv <cl_float8, float, 8> Array2Dv8f;
typedef  Array2Dv <cl_float16, float, 16> Array2Dv16f;

typedef  Array2Dv <cl_int2, int, 2> Array2Dv2i;
typedef  Array2Dv <cl_int3, int, 3> Array2Dv3i;
typedef  Array2Dv <cl_int4, int, 4> Array2Dv4i;
typedef  Array2Dv <cl_int8, int, 8> Array2Dv8i;
typedef  Array2Dv <cl_int16, int, 16> Array2Dv16i;

typedef  Array2Dv <cl_uint2, unsigned int, 2> Array2Dv2u;
typedef  Array2Dv <cl_uint3, unsigned int, 3> Array2Dv3u;
typedef  Array2Dv <cl_uint4, unsigned int, 4> Array2Dv4u;
typedef  Array2Dv <cl_uint8, unsigned int, 8> Array2Dv8u;
typedef  Array2Dv <cl_uint16, unsigned int, 16> Array2Dv16u;




template <class T> class Array3D : public ArrayBase<T>
{
private:
	void testIndex(const int n1, const int n2, const int n3)
	{
#ifdef WARN_OUTSIDE_MAIN_BOUNDS
		if (n1 >= ArrayBase<T>::SizeX && n1 < ArrayBase<T>::FullSizeX)
		{
			printf("Accessing element %d, which is > SizeX (%d), but within FullSizeX bounds (%d)\n", n1, ArrayBase<T>::SizeX, ArrayBase<T>::FullSizeX);
		}
#endif
		CHECK_ARRAY_ERROR((0 > n1 || n1 >= ArrayBase<T>::FullSizeX),
			("Tried to access element n1 = " + std::to_string(n1) +
				" in X dim of array of size [" + std::to_string(SizeX) +
				", " + std::to_string(SizeY) + ", " + std::to_string(SizeZ) +
				"]"), ERROR_MEMORY_OUT_OF_BOUNDS);

#ifdef WARN_OUTSIDE_MAIN_BOUNDS
		if (n2 >= ArrayBase<T>::SizeY && n2 < ArrayBase<T>::FullSizeY)
		{
			printf("Accessing element %d, which is > SizeY (%d), but within FullSizeY bounds (%d)\n", n2, ArrayBase<T>::SizeY, ArrayBase<T>::FullSizeY);
		}
#endif
		CHECK_ARRAY_ERROR((0 > n2 || n2 >= ArrayBase<T>::FullSizeY),
			("Tried to access element n2 = " + std::to_string(n2) +
				" in X dim of array of size [" + std::to_string(SizeX) +
				", " + std::to_string(SizeY) + ", " + std::to_string(SizeZ) +
				"]"), ERROR_MEMORY_OUT_OF_BOUNDS);

#ifdef WARN_OUTSIDE_MAIN_BOUNDS
		if (n3 >= ArrayBase<T>::SizeZ && n3 < ArrayBase<T>::FullSizeZ)
		{
			printf("Accessing element %d, which is > SizeZ (%d), but within FullSizeZ bounds (%d)\n", n3, ArrayBase<T>::SizeZ, ArrayBase<T>::FullSizeZ);
		}
#endif
		CHECK_ARRAY_ERROR((0 > n3 || n3 >= ArrayBase<T>::FullSizeZ),
			("Tried to access element n3 = " + std::to_string(n3) +
				" in X dim of array of size [" + std::to_string(SizeX) +
				", " + std::to_string(SizeY) + ", " + std::to_string(SizeZ) +
				"]"), ERROR_MEMORY_OUT_OF_BOUNDS);
	}

public:


	Array3D(std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);
	}

	Array3D(const int n1, const int n2, const int n3, std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);

		ArrayBase<T>::SizeX = n1;
		ArrayBase<T>::FullSizeX = n1;

		ArrayBase<T>::SizeY = n2;
		ArrayBase<T>::FullSizeY = n2;

		ArrayBase<T>::SizeZ = n3;
		ArrayBase<T>::FullSizeZ = n3;

		ArrayBase<T>::FullSize = FullSizeX*FullSizeY*FullSizeZ;

		ArrayBase<T>::AllocHost();
	}

	Array3D(const int n1, const int fs1, const int n2, const int fs2,
		const int n3, const int fs3, std::string name_ = "default")
	{
		ArrayBase<T>::ini(name_);

		ArrayBase<T>::SizeX = n1;
		ArrayBase<T>::FullSizeX = fs1;

		ArrayBase<T>::SizeY = n2;
		ArrayBase<T>::FullSizeY = fs2;

		ArrayBase<T>::SizeZ = n3;
		ArrayBase<T>::FullSizeZ = fs3;

		ArrayBase<T>::FullSize = FullSizeX*FullSizeY*FullSizeZ;

		ArrayBase<T>::AllocHost();
	}

	virtual ~Array3D()
	{
		ArrayBase<T>::FreeHost();
		ArrayBase<T>::FreeDevice();
	}

	T& operator[](const int n1)
	{
		return *(ArrayBase<T>::Array + n1);
	}

	T& operator()(const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2, n3);
#endif 
		return *(ArrayBase<T>::Array + n1 + ArrayBase<T>::FullSizeX * n2 + ArrayBase<T>::FullSizeX *ArrayBase<T>::FullSizeY * n3);
	}

	T& operator()(cl_int3 n)
	{
		return operator()(n.x, n.y, n.z);
	}

	T& operator()(cl_int2 n1, int n2)
	{
		return operator()(n1.x, n1.y, n2);
	}

	T getEl(const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2, n3);
#endif 
		return *(ArrayBase<T>::Array + n1 + ArrayBase<T>::FullSizeX * n2 +
			ArrayBase<T>::FullSizeX *ArrayBase<T>::FullSizeY * n3);
	}

	T getEl(const int n1, const int n2, const int n3, bool testfl)
	{
#ifdef ARRAY_DEBUG
		if (testfl)
			testIndex(n1, n2, n3);
#endif 
		return *(ArrayBase<T>::Array + n1 + ArrayBase<T>::FullSizeX * n2 +
			ArrayBase<T>::FullSizeX *ArrayBase<T>::FullSizeY * n3);
	}

	void setEl(const T val, const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2, n3);
#endif
		*(ArrayBase<T>::Array + n1 + ArrayBase<T>::FullSizeX * n2 + 
			ArrayBase<T>::FullSizeX *ArrayBase<T>::FullSizeY * n3) = val;
	}

	//memory allocation
	void allocate(const int n1, const int n2, const int n3)
	{
		ArrayBase<T>::ini();

		ArrayBase<T>::SizeX = n1;
		ArrayBase<T>::FullSizeX = n1;

		ArrayBase<T>::SizeY = n2;
		FullSizeY = n2;

		ArrayBase<T>::SizeZ = n3;
		ArrayBase<T>::FullSizeZ = n3;

		ArrayBase<T>::FullSize = FullSizeX*FullSizeY*FullSizeZ;

		ArrayBase<T>::AllocHost();
	}

	void allocate(cl_int3 n)
	{
		allocate(n.x, n.y, n.z)
	}

	void allocate(const int n1, const int fs1, const int n2, const int fs2, const int n3, const int fs3)
	{
		ArrayBase<T>::ini();

		ArrayBase<T>::SizeX = n1;
		ArrayBase<T>::FullSizeX = fs1;

		ArrayBase<T>::SizeY = n2;
		ArrayBase<T>::FullSizeY = fs2;

		ArrayBase<T>::SizeZ = n3;
		ArrayBase<T>::FullSizeZ = fs3;

		ArrayBase<T>::FullSize = FullSizeX*FullSizeY*FullSizeZ;

		ArrayBase<T>::AllocHost();
	}

	void zeros(const int n1, const int n2, const int n3)
	{
		allocate(n1, n2, n3);
		ArrayBase<T>::zeros();
	}

	void zeros(cl_int3 n)
	{
		allocate(n);
		ArrayBase<T>::zeros();
	}

	void zeros(const int n1, const int fs1, const int n2, const int fs2, const int n3, const int fs3)
	{
		allocate(n1, fs1, n2, fs2, n3, fs3);
		ArrayBase<T>::zeros();
	}

	bool save_txt_from_device_as_multi2D(std::string FileName = "", bool FullFlag = false)
	{
		int delflag = 0;
		if (ArrayBase<T>::Host_Alloc_Flag == 0)
		{
			delflag = 1;
		}
		read_from_buffer();
		bool ret = savetxt_as_multi2D(FileName, FullFlag);
		if (delflag)
			ArrayBase<T>::FreeHost();
		return ret;
	}

	virtual bool save2file(std::string FileName = "")
	{
		int i, j, k;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");
		std::fstream stream;
		if (fileopen(stream, FileName, FileOut) == false)
			return false;

		for (i = 0; i < SizeX; i++)
		{
			for (j = 0; j < SizeY; j++)
			{
				for (k = 0; k < SizeZ; k++)
				{
					stream << getEl(i, j, k) << "\t";
				}
				stream << "\n";
			}
			stream << "\n";
		}
		stream.close();
		return true;
	}

	virtual bool savetxt_as_multi2D(std::string FileName = "", bool FullFlag = false)
	{
		if (FileName.length() == 0)
			FileName = Name;
		
		int i, j, k;
		int Xbounds = (FullFlag) ? FullSizeX : SizeX;
		int Ybounds = (FullFlag) ? FullSizeY : SizeY;
		int Zbounds = (FullFlag) ? FullSizeZ : SizeZ;

		for (k = 0; k < Zbounds; k++)
		{
			std::string FileNameOut = FileName + "_" + std::to_string(k) + ".txt";
			std::fstream stream;
			if (fileopen(stream, FileNameOut, FileOut) == false)
				return false;


			for (i = 0; i < Xbounds; i++)
			{
				for (j = 0; j < Ybounds; j++)
				{
					stream << getEl(i, j, k) << "\t";
				}
				stream << "\n";
			}
			stream.close();
		}

		return true;
	}
	
	virtual bool loadtxt(std::string FileName = "")
	{
		int i, j, k;
		T Buf;
		if (FileName.length() == 0)
			FileName = Name;
		FileName.append(".txt");

		std::fstream stream;
		if (fileopen(stream, FileName, FileIn) == false)
			return false;

		for (i = 0; i < ArrayBase<T>::SizeX; i++)
		{
			for (j = 0; j < ArrayBase<T>::SizeY; j++)
			{
				for (k = 0; k < ArrayBase<T>::SizeZ; k++)
				{
					stream >> Buf; 
					setEl(Buf, i, j, k);
				}
			}
		}
		stream.close();
		return true;
	}

	void reallocate_host(const int n1, const int n2, int n3)
	{
		T * p = new T[FullSize];
		memcpy(p, ArrayBase<T>::Array, FullSize*sizeof(T));


		delete [] Array;
		Array = new T[n1*n2*n3];
		ArrayBase<T>::Host_Alloc_Flag = 1;

		int Xlim = MIN(SizeX, n1);
		int Ylim = MIN(SizeY, n2);
		int Zlim = MIN(SizeZ, n3);

		for (int i = 0; i < Xlim; i++)
		{
			for (int j = 0; j < Ylim; j++)
			{
				for (int k = 0; k < Zlim; k++)
				{
					Array[i*n2*n3 + j*n3 + k] = p[i*FullSizeY*FullSizeZ + j*FullSizeZ + k];
				}
			}
		}

		SizeX = n1;
		SizeY = n2;
		SizeZ = n3;
		FullSizeX = n1;
		FullSizeY = n2;
		FullSizeY = n3;
		FullSize = n1*n2*n3;

		delete [] p;
	};

	void fill(const T val_)
	{
		std::fill(Array, Array + FullSize, val_);
	}
};

typedef  Array3D <double> Array3Dd;
typedef  Array3D <int> Array3Di;
typedef  Array3D <signed char> Array3Dc;
typedef  Array3D <float> Array3Df;
typedef  Array3D <unsigned int> Array3Du;





class Array1DGL : public Array1D <float>
{
public:
	GLuint GLpos, GLcol;
	cl_float3 ColorVec;
	cl_mem ColorBuf;
	int ColorVBO_flag;
	int Point_flag;
	int Device_Col_Alloc_Flag;
	float *ColArray;
	int FullSize_col;
	int Host_Color_Alloc_Flag;
	Array1DGL(std::string name_ = "")
	{
		if (name_.length() == 0)
			name_ = "default";
		Name = name_;
		GLpos = -1;
		GLcol = -1;
		ColorVec = { { -1.f, -1.f, -1.f } };
		ColorBuf = nullptr;
		ColorVBO_flag = -1;
		Point_flag = -1;
		Device_Col_Alloc_Flag = 0;
		Host_Color_Alloc_Flag = 0;
		ArrayBase<float>::ini();
	}

	~Array1DGL()
	{
		Array1DGL::delete_array();
	}

	void FreeDeviceColor()
	{
		if (Device_Col_Alloc_Flag)
		{
			clReleaseMemObject(ColorBuf);
			Device_Col_Alloc_Flag = 0;
			ColorVBO_flag = -1;
			ColorVec = { { -1.f, -1.f, -1.f } };
		}
	}

	void FreeHostColor()
	{
		delete[] ColArray;
		Host_Color_Alloc_Flag = 0;
	}

	void FreeGLColor()
	{
		if (Device_Col_Alloc_Flag)
			FreeDeviceColor();

		if (GLcol != -1)
		{
			glDeleteBuffers(1, &GLcol);
			GLcol = -1;
		}

	}

	void FreeGLArray()
	{
		if (GLpos != -1)
		{
			glDeleteBuffers(1, &GLpos);
			GLpos = -1;
		}
		Point_flag = -1;
	}

	void delete_array()
	{
		FreeGLColor();
		FreeHostColor();
		ArrayBase<float>::FreeDevice();
		FreeGLArray();
		ArrayBase<float>::FreeHost();
	}

	void CreateVBO_position(bool vbo_type, int clflags = CL_MEM_READ_WRITE)
	{
		Point_flag = vbo_type;
		CHECK_ARRAY_ERROR(Host_Alloc_Flag == 0, 
			"trying to create VBO with empty array", ERROR_IN_OPENGL_ARRAY);
		GLenum target = GL_ARRAY_BUFFER;
		GLenum usage = GL_DYNAMIC_DRAW;
		int dataSize = sizeof(cl_float)*FullSize;
		GLuint pid = 0;
		int err;
		glGenBuffers(1, &pid);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error generating glbuffer", err);
		
		glBindBuffer(target, pid);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error binding glbuffer", err);
		
		glBufferData(target, dataSize, (void*)Array, usage);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error writing to glbuffer", err);
		
		int bsize = 0;
		glGetBufferParameteriv(target, GL_BUFFER_SIZE, &bsize);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error calling glGetBufferParameteriv", err);
		
		if (dataSize != bsize)
		{
			glDeleteBuffers(1, &pid);
			CHECK_ARRAY_ERROR(true, "[createVBO()] Data size is mismatched with input array",
				ERROR_IN_OPENGL_ARRAY);
			pid = 0;
		}
		glBindBuffer(target, 0);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error error binding buffer", err);
		
		GLpos = pid;
		glFinish();
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "glFinish returned error", err);

		create_pos_buffer_from_gl(*staticBaseVar::context, clflags);
	}

	void AllocateColor(int csize)
	{
		ColArray = new float[csize];
		FullSize_col = csize;
		Host_Color_Alloc_Flag = 1;
	}

	void CreateVBO_color(int clflags = CL_MEM_READ_WRITE)
	{
		ColorVBO_flag = VAR_COLOR_VBO;
		CHECK_ARRAY_ERROR(Host_Color_Alloc_Flag == 0, 
			"trying to create color VBO with empty array", ERROR_IN_OPENGL_ARRAY);
	
		int err = 0;
		GLenum target = GL_ARRAY_BUFFER;
		GLenum usage = GL_DYNAMIC_DRAW;
		int dataSize = sizeof(cl_float)*FullSize_col;
		GLuint cid = 0;
		
		glGenBuffers(1, &cid);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error generating VBO color buffer", err);

		glBindBuffer(target, cid);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error binding VBO color buffer", err);

		glBufferData(target, dataSize, (void*)ColArray, usage);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error copying to VBO color buffer", err);


		int bcsize = 0;
		glGetBufferParameteriv(target, GL_BUFFER_SIZE, &bcsize);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error calling glGetBufferParameteriv", err);

		if (dataSize != bcsize)
		{
			glDeleteBuffers(1, &cid);
			CHECK_ARRAY_ERROR(true, "[createVBO()] Data size is mismatched with input array",
				ERROR_IN_OPENGL_ARRAY);
			cid = 0;
		}
		glBindBuffer(target, 0);
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error binding buffer", err);

		GLcol = cid;
		glFinish();
		err = glGetError();
		CHECK_ARRAY_ERROR(err, "Error call glFinish()", err);
		
		create_col_buffer_from_gl(*context, clflags);
	}

	void release_gl_objects(cl_command_queue queue)
	{
		int status = clEnqueueReleaseGLObjects(queue, 1, &Buffer, 0, 0, 0);
		CHECK_ARRAY_ERROR(status, "Error releasing GL objects", status);
		if (VAR_COLOR_VBO)
		{
			status = clEnqueueReleaseGLObjects(queue, 1, &ColorBuf, 0, 0, 0);
			CHECK_ARRAY_ERROR(status, "Error releasing GL Color objects", status);
		}
	}

	void acquire_gl_objects(cl_command_queue queue)
	{
		int status = clEnqueueAcquireGLObjects(queue, 1, &Buffer, 0, 0, 0);
		CHECK_ARRAY_ERROR(status, "Error acquiring GL objects", status);
		if (VAR_COLOR_VBO)
		{
			status += clEnqueueAcquireGLObjects(queue, 1, &ColorBuf, 0, 0, 0);
			CHECK_ARRAY_ERROR(status, "Error releasing GL color objects", status);
		}
	}

	void CreateVBO_color(cl_float3 color_temp)
	{
		ColorVBO_flag = CONST_COLOR_VBO;
		ColorVec = color_temp;
	}


	void create_pos_buffer_from_gl(cl_context context_, int clflags)
	{
		Device_Alloc_Flag = 1;
		clFlags = clflags;
		int status;
		Buffer = clCreateFromGLBuffer(context_, clFlags, GLpos, &status);
		CHECK_ARRAY_ERROR(status, "Error creating pos buffer from GL Buffer", status);
	}

	void create_col_buffer_from_gl(cl_context context_, int clflags)
	{
		CHECK_ARRAY_ERROR(GLcol == -1, "Color VBO must be created before clBuffer",
			ERROR_IN_OPENGL_ARRAY);
		int status;
		ColorBuf = clCreateFromGLBuffer(context_, clflags, GLcol, &status);
		CHECK_ARRAY_ERROR(status, "Error creating color buffer from GL Buffer", status);
		Device_Col_Alloc_Flag = 1;
	}

	void Render_VBO(double size)
	{
		if (ColorVBO_flag == VAR_COLOR_VBO)
		{
			if (Point_flag == POINT_VBO)
				Render_Points_Variable(size);
			else
				Render_Lines_Variable(size);
		}
		else
		{
			if (Point_flag == POINT_VBO)
				Render_Points_Const(size);
			else
				Render_Lines_Const(size);
		}
	}

	cl_mem* get_col_buf_add()
	{
		return &ColorBuf;
	}

	void Render_Lines_Const(double size)
	{
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(size);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

		glBindBuffer(GL_ARRAY_BUFFER, GLpos);
		glVertexPointer(2, GL_FLOAT, 0, 0);

		glEnableClientState(GL_VERTEX_ARRAY);
		glColor3f(ColorVec.x, ColorVec.y, ColorVec.z);


		glDisableClientState(GL_NORMAL_ARRAY);

		glDrawArrays(GL_LINE_STRIP, 0, FullSize / 2);

		glDisableClientState(GL_VERTEX_ARRAY);

		glDisableClientState(GL_LINE_SMOOTH);
	}


	void Render_Individual_Lines_Const(double size)
	{
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(size);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

		glBindBuffer(GL_ARRAY_BUFFER, GLpos);
		glVertexPointer(2, GL_FLOAT, 0, 0);

		glEnableClientState(GL_VERTEX_ARRAY);
		glColor3f(ColorVec.x, ColorVec.y, ColorVec.z);


		glDisableClientState(GL_NORMAL_ARRAY);

		glDrawArrays(GL_LINES, 0, FullSize / 2);

		glDisableClientState(GL_VERTEX_ARRAY);

		glDisableClientState(GL_LINE_SMOOTH);
	}

	void set_color_array(cl_float3 cv)
	{
		CHECK_ARRAY_ERROR(Host_Color_Alloc_Flag == 0,
			"Trying to write to unallocated color array", ERROR_IN_OPENGL_ARRAY);
		
		for (int i = 0; i < FullSize_col / 3; i++)
		{
			ColArray[i * 3] = cv.x;
			ColArray[i * 3 + 1] = cv.y;
			ColArray[i * 3 + 2] = cv.z;
		}
	}

	void set_color_vector(cl_float3 cv)
	{
		ColorVec = { { cv.x, cv.y, cv.z } };
	}


	void Render_Lines_Variable(double size)
	{
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(size);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

		glBindBuffer(GL_ARRAY_BUFFER, GLcol);
		glColorPointer(3, GL_FLOAT, 0, 0);

		glBindBuffer(GL_ARRAY_BUFFER, GLpos);
		glVertexPointer(2, GL_FLOAT, 0, 0);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		glDisableClientState(GL_NORMAL_ARRAY);
		glDrawArrays(GL_LINE_STRIP, 0, FullSize / 2);

		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);

		glDisable(GL_LINE_SMOOTH);
	}

	void Render_Points_Const(double size)
	{
		glEnable(GL_POINT_SMOOTH);
		glPointSize(size);

		glBindBuffer(GL_ARRAY_BUFFER, GLpos);
		glVertexPointer(2, GL_FLOAT, 0, 0);

		glEnableClientState(GL_VERTEX_ARRAY);
		glColor3f(ColorVec.x, ColorVec.y, ColorVec.z);

		glDisableClientState(GL_NORMAL_ARRAY);

		glDrawArrays(GL_POINTS, 0, FullSize / 2);

		glDisableClientState(GL_VERTEX_ARRAY);

		glDisableClientState(GL_POINT_SMOOTH);
	}

	void Render_Points_Variable(double size)
	{
		glEnable(GL_POINT_SMOOTH);
		glPointSize(size);

		glBindBuffer(GL_ARRAY_BUFFER, GLcol);
		glColorPointer(3, GL_FLOAT, 0, 0);

		glBindBuffer(GL_ARRAY_BUFFER, GLpos);
		glVertexPointer(2, GL_FLOAT, 0, 0);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		glDisableClientState(GL_NORMAL_ARRAY);

		glDrawArrays(GL_POINTS, 0, FullSize / 2);

		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);

		glDisableClientState(GL_POINT_SMOOTH);
	}


};

#endif // !defined(AFX_ARRAY_H__INCLUDED_)
