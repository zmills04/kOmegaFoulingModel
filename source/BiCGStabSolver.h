// Sparse Matrix class for clSPARSE library
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BICGSTABSOLVER_H__INCLUDED_)
#define AFX_BICGSTABSOLVER_H__INCLUDED_
#pragma once

#include "StdAfx.h"
#include "SparseMatrix.h"
#include "BiCGStabGenerator.h"
#include "Kernels.h"
#include "Array.h"

class BiCGStabDebugger;
// Easier to just use clsparse's functions to calculate meta info for matrix
// and extract it from there. The matrix_meta class is not accessable through
// the include files, so it has been copied and pasted here.
struct matrix_meta
{
	matrix_meta() : rowBlockSize(0), offRowBlocks(0)
	{
	}

	void clear()
	{
		offRowBlocks = rowBlockSize = 0;
		rowBlocks = ::cl::Buffer();
	}

	::cl::Buffer rowBlocks;  /*!< Meta-data used for csr-adaptive algorithm; can be NULL */

	clsparseIdx_t rowBlockSize;  /*!< Size of array used by the rowBlocks handle */

	clsparseIdx_t offRowBlocks;
};




// TODO: use subbuffers to create offsets in order to avoid lower section
//		of domain (and reduce worksize to avoid upper section). This may provide
//		an additional speed up since each of these kernel may be called multiple
//		times per timestep. 

class BiCGStabSolver
{
public:

	cl_command_queue *calcQue; // can be different for each sparse solver
	
	std::string Name;

	static int fullSize; // size of vectors including padding(>= XsizeFull*nY)
	static int colSize; // actual size of vectors w/out padding, and size(IA) - 1
	int NNZ; // number of nonzeros sinze of A and JA
	static int Xsize, Ysize, XsizeFull; // size of domain
	static int numRedKer, *globalWorkSizeRed;
	int maxIters;
	double relTol, absTol, initialResidual;
	// variables used for adaptive CSRMV kernels. Calculated by clSparse
	// Keeping this unique to each instance just in case there may be slight
	// differences based on sparse layout.
	cl_uint rowBlockSize, offRowBlocks;
	cl_mem rowBlocks;

	int globalSizeCSRMV, globalSizeAXPBY;


	// IA and JA arrays, pointer used to allow for same CSR_Inds to be shared across multiple
	// instances, but right now, its not currently working, so use a different one for each
	CSR_Inds *Inds;

	// A matrix
	Array1Dd Amat;

	Array1Dd bVec; // also used as storage for r and s vectors during iterations

	// Actual values of macroscopic variable being solved for (pointer to array in containing class)
	Array2Dd *xVec; // also used as storage for h vector during iterations

	// Stores scalars 
	// { norm_b, norm_r, norm_s, rho, omega, alpha, beta }
	Array1Dd scalarBuf;


	// Sharing buffers between each solver instance to save on memory. This requires
	// all instances to share the same matrix/vector sizes and does not allow for
	// solving operations to overlap.
	static cl_mem r0Hat, pVec, vVec, tVec;

	// Using two buffers for each intermediate set (alternates between the two, instead of one per reduce step)
	static cl_mem *reduceBufSet1, *reduceBufSet2;
	static bool staticVarInitialized;

	// Using a different kernel for each step will reduce number of steps during calculation,
	// but increase total number of kernels. Not sure if this may cause issues with GPU (wasting
	// memory or possibly reaching limit on number of kernels). 
	// Kernels for sparse solvers. axpbyKernelSingle refers to one used to calculate P vector
	// axpbyKernelDouble is used to calculate s and h and x and r

	//Step 1
	RedKernelList bNorm;

	//Step 2
	thinKerWrapper rCSRMV;

	// Step 3: copy r_0 to rhat_0 and p1

	// Step 4
	RedKernelList rhoSqr;

	// Step f1/f2: 
	RedKernelList rhoDot;

	// Step f3
	thinKerWrapper pAXPBY;

	// Step f4
	thinKerWrapper vCSRMV;

	// Step f5
	RedKernelList alphaDot;

	// Step f6
	thinKerWrapper shAXPBY;

	// Step f7: uses bNorm from step 1

	// Step f8
	thinKerWrapper tCSRMV;

	// Step f9
	RedKernelList omegaDotSqr; // Two sets of reduction base kernels and a single final reduce
	// calculating omega. length of list is num_red_ker*2 - 1
	// Step f10
	thinKerWrapper xrAXPBY;


	// ElementWise Operations
	thinKerWrapper addKer, subKer, multKer, divKer;




	// Step f11: uses bNorm from step 1

	enum ScalarIndex { indNormB = 0, indNormR = 1, indNormS = 2, indRho = 3, indOmega = 4, indAlpha = 5, indBeta = 6 };
	enum { C, E, W, N, S };

	BiCGStabSolver() : scalarBuf(7, "scalarBuf")
	{
		rowBlocks = nullptr;
	}

	~BiCGStabSolver()
	{
		FREE_OCL_BUFFER(rowBlocks);
		if (staticVarInitialized)
		{
			FREE_OCL_BUFFER(r0Hat);
			FREE_OCL_BUFFER(pVec);
			FREE_OCL_BUFFER(vVec);
			FREE_OCL_BUFFER(tVec);
			FREE_OCL_BUFFER_ARRAY(reduceBufSet1, 2);
			FREE_OCL_BUFFER_ARRAY(reduceBufSet2, 2);
			delete[] globalWorkSizeRed;
			staticVarInitialized = false;
		}
	}

	void setMaxIters(int iters_) { maxIters = iters_; }
	void setAbsTol(int tol_) { absTol = tol_; }
	void setRelTol(int tol_) { relTol = tol_; }

	void iniAXPBYKernels();

	// y = alpha*A*x + beta * y;
	void createCSRMVKernel(thinKerWrapper &ker, double alphaval, double betaval, cl_mem *xvec_, cl_mem* yvec_);

	void iniCSRMVKernels();

	void createElementwiseKernels();

	// initializes kernels and sets arguments for everthing but input buffer(s) in first kernel,
	// and output, local memory and offset values in last kernel.
	// startInd is the index of the intermediate buffer in the first kernel (1 for all but dot reduce)
	void createReduceListBases(RedKernelList &redlist, const int startInd, std::string base1ker,
		std::string base2ker, std::string finalker, cl_mem* mem_use);

	void createBNorm();

	void createRhoSqr();

	void createRhoDot();

	void createAlphaDot();

	void createOmegaDotSqr();

	static void iniBuffer(cl_mem &buf_, const int size_, const std::string name_);

	// initializes all static variables
	static void ini(int xsize_, int xsizefull_, int ysize_);

	// creates CSR matrix and clsparsecontrol object, generates meta information, copies
	// it over and deletes clsparse objects since they are no longer needed.
	void getCSRMeta();

	//size_t globalSizeCSRMV, globalSizeAXPBY;
	void CreateSolver(Array2Dd *macro_, CSR_Inds *inds_, cl_command_queue *calcque_, int maxiters_, double reltol_, double abstol_);

	void runReduce(RedKernelList &redlist_, cl_command_queue* que_, int num_wait = 0,
		cl_event *wait = NULL, cl_event* evt = NULL);

	bool reduceAndCheckConvergence(cl_command_queue *que_, bool setInitialRes = false, int num_wait = 0,
		cl_event *wait = NULL, cl_event* evt = NULL);

	void copy_buffers(cl_mem *src_buf, cl_mem *dest_buf);

	void solve();

	// for debugging
	void printScalars();

	// Functions for transfering between host and device and outputting data

	//void setInitialValue(double inival, bool fullArrFlag = false);

	//void setInitialValueRows(double inival, std::vector<int> &rowi);
	//void setInitialValueCols(double inival, std::vector<int> &coli);

	double& BiCGStabSolver::operator()(int i, int j)
	{
		return xVec->operator()(i, j);
	}

	int getBufferFullSize();

	std::string getName();
	///////// Load and CheckPoint Methods ///////////



	//////// Save Methods /////////////////////
	bool saveAxbCSR();

	bool saveAxb_w_indicies();

	bool saveAxbCSR_from_device();
	bool saveAxb_w_indicies_from_device();

	cl_mem* get_add_A();
	cl_mem* get_add_IndArr();

	Array1Dd* getBVec() { return &bVec; }
	Array2Dd* getMacroArray() { return xVec; }

	bool savetxt(std::string outname = "");

	bool savetxt_from_device(std::string outname = "");
	bool save_bvec(std::string outname = "");

	bool save_bvec_from_device(std::string outname = "");

	bool saveCheckPoint(std::string outname = "");

	cl_mem* get_add_b() { return bVec.get_buf_add(); }

	cl_mem* get_add_Macro() { return xVec->get_buf_add(); }

	void copy_to_device(const int blFlag = true);

	void copy_to_host(const int blFlag = true);

	void copy_inds_to_host(const int blFlag = true);

	void copy_inds_to_device(const int blFlag = true);

	void copy_to_host_all(const int blFlag = true);

	void copy_to_device_all(const int blFlag = true);

	int rows() { return Inds->rows; }
	int cols() { return Inds->cols; }
	int nnz() { return Inds->nnz(); }

	double A(const int i, const int j, int dir = C);

	bool testInd(const int i, const int j, const int dir = C);

	bool save_w_indicies(std::string Name, bool fromDevFlag = false);
	bool saveCSR(std::string outname, bool fromDevFlag = false);
	bool saveCSR_row_col_val(std::string outname, bool fromDevFlag = false);

	bool saveAxb_w_indicies_from_device_as_bin();
	bool save_w_indicies_as_bin(std::string Name, bool fromDevFlag);
	void saveTempArray(cl_mem *mem_, std::string savename_)
	{
		Array1Dd tempArr(fullSize, savename_);
		tempArr.read_from_buffer_to_array(*mem_);
		tempArr.savebin();
	}

};



class BiCGStabDebugger
{
public:
	clsparseCsrMatrix clsCSR;
	Array1Dd clA, clB, clX, clP, clS, clT, clH;
	clsparseScalar clOne, clZero, clNegOne;
	cldenseVector clXVec, clBVec, clPVec, clSVec, clTVec, clHVec, clRVec;

	enum ArrayTag { B, X, P, S, T, H, R };

	~BiCGStabDebugger()
	{ // class destructors handle freeing remaining memory
		clReleaseMemObject(clOne.value);
		clReleaseMemObject(clZero.value);
		clReleaseMemObject(clNegOne.value);
	}

	BiCGStabDebugger(Array1Dd &A_, Array1Dd &B_, Array2Dd &X_, CSR_Inds *inds_) :
		clA(A_.getBufferFullSize(), "clA"),
		clB(B_.getBufferFullSize(), "clB"), clX(X_.getBufferFullSize(), "clX"),
		clP(B_.getBufferFullSize(), "clP"), clS(X_.getBufferFullSize(), "clS"),
		clT(B_.getBufferFullSize(), "clT"), clH(X_.getBufferFullSize(), "clH")
	{
		clsparseInitCsrMatrix(&clsCSR);

		clsCSR.num_nonzeros = inds_->nnz();
		clsCSR.num_rows = inds_->rows;
		clsCSR.num_cols = inds_->cols;

		clA.allocate_buffer_w_copy();
		clA.FillBuffer(0.);
		clsCSR.values = clA.get_buffer();
		clsCSR.col_indices = inds_->JA.get_buffer();
		clsCSR.row_pointer = inds_->IA.get_buffer();
		clsparseCsrMetaCreate(&clsCSR, BiCGStabGenerator::BiCGStabInstance()->clSparseControl.control);
		clA.write_array_to_buffer(A_.get_array(), CL_TRUE, inds_->nnz());

		int fullSize = B_.getBufferFullSize();

		clX.allocate_buffer_w_copy();
		clX.FillBuffer(0.);
		clXVec.num_values = fullSize;
		clXVec.values = clX.get_buffer();
		clX.write_array_to_buffer(X_.get_array(), CL_TRUE, X_.getFullSize());



		clB.allocate_buffer_w_copy();
		clB.FillBuffer(0.);
		clBVec.num_values = fullSize;
		clBVec.values = clB.get_buffer();
		clB.write_array_to_buffer(B_.get_array(), CL_TRUE, X_.getFullSize());



		clP.allocate_buffer_w_copy();
		clPVec.num_values = fullSize;
		clPVec.values = clP.get_buffer();

		clS.allocate_buffer_w_copy();
		clSVec.num_values = fullSize;
		clSVec.values = clS.get_buffer();

		clT.allocate_buffer_w_copy();
		clTVec.num_values = fullSize;
		clTVec.values = clT.get_buffer();

		clH.allocate_buffer_w_copy();
		clHVec.num_values = fullSize;
		clTVec.values = clH.get_buffer();


		clsparseInitScalar(&clOne);
		clOne.value = clCreateBuffer(CLCONTEXT, CL_MEM_READ_ONLY, sizeof(double), nullptr, nullptr);
		clsparseInitScalar(&clZero);
		clZero.value = clCreateBuffer(CLCONTEXT, CL_MEM_READ_ONLY, sizeof(double), nullptr, nullptr);
		clsparseInitScalar(&clNegOne);
		clNegOne.value = clCreateBuffer(CLCONTEXT, CL_MEM_READ_ONLY, sizeof(double), nullptr, nullptr);


		double *halpha = (double*)clEnqueueMapBuffer(IOQUEUE, clOne.value, CL_TRUE, CL_MAP_WRITE, 0, sizeof(double), 0, nullptr, nullptr, nullptr);
		*halpha = 1.;
		clEnqueueUnmapMemObject(IOQUEUE, clOne.value, halpha, 0, nullptr, nullptr);

		halpha = (double*)clEnqueueMapBuffer(IOQUEUE, clZero.value, CL_TRUE, CL_MAP_WRITE, 0, sizeof(double), 0, nullptr, nullptr, nullptr);
		*halpha = 0.;
		clEnqueueUnmapMemObject(IOQUEUE, clZero.value, halpha, 0, nullptr, nullptr);

		halpha = (double*)clEnqueueMapBuffer(IOQUEUE, clNegOne.value, CL_TRUE, CL_MAP_WRITE, 0, sizeof(double), 0, nullptr, nullptr, nullptr);
		*halpha = -1.;
		clEnqueueUnmapMemObject(IOQUEUE, clNegOne.value, halpha, 0, nullptr, nullptr);

	}

	clsparseCsrMatrix* getMat() { return &clsCSR; }

	bool saveMat() { return clA.save_bin_from_device(); }

	cldenseVector* getVec(ArrayTag at_)
	{
		switch (at_)
		{
		case X:
		{
			return &clXVec;
		}
		case B:
		{
			return &clBVec;
		}
		case P:
		{
			return &clPVec;
		}
		case S:
		{
			return &clSVec;
		}
		case T:
		{
			return &clTVec;
		}
		case H:
		{
			return &clHVec;
		}
		case R:
		{
			return &clRVec;
		}
		}
		return nullptr;
	}

	bool saveVec(ArrayTag at_)
	{
		switch (at_)
		{
		case X:
		{
			return clX.save_bin_from_device();
		}
		case B:
		{
			return clB.save_bin_from_device();
		}
		case P:
		{
			return clP.save_bin_from_device();
		}
		case S:
		{
			return clS.save_bin_from_device();
		}
		case T:
		{
			return clT.save_bin_from_device();
		}
		case H:
		{
			return clH.save_bin_from_device();
		}
		case R:
		{
			return clB.save_bin_from_device("clR");
		}
		}
		return false;
	}

	clsparseControl getControl()
	{
		return BiCGStabGenerator::BiCGStabInstance()->clSparseControl.control;
	}
};






#endif // AFX_BICGSTAB_H__INCLUDED_