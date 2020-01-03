// Sparse Matrix class for clSPARSE library
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GMRESSOLVER_H__INCLUDED_)
#define AFX_GMRESSOLVER_H__INCLUDED_
#pragma once

#include "StdAfx.h"
#include "SparseMatrix.h"
#include "BiCGStabGenerator.h"
#include "GMRESGenerator.h"
#include "Kernels.h"
#include "Array.h"



// TODO: use subbuffers to create offsets in order to avoid lower section
//		of domain (and reduce worksize to avoid upper section). This may provide
//		an additional speed up since each of these kernel may be called multiple
//		times per timestep. 

class GMRESSolver
{
public:
	enum { C, E, W, N, S };
	   
	cl_command_queue* calcQue; // can be different for each sparse solver

	std::string Name;

	static cl_uint fullSize; // size of vectors including padding(>= XsizeFull*nY)
	static cl_uint colSize; // actual size of vectors w/out padding, and size(IA) - 1
	cl_uint NNZ; // number of nonzeros sinze of A and JA

	static cl_uint Xsize, Ysize, XsizeFull; // size of domain
	static int numRedKer, *globalWorkSizeRed, *localWorkSizeRed;

	int maxIters;
	cl_uint maxRestarts;
	double relTol, absTol, initialResidual;

	// variables used for adaptive CSRMV kernels. Calculated by clSparse
	// Keeping this unique to each instance just in case there may be slight
	// differences based on sparse layout.
	cl_uint rowBlockSize, offRowBlocks;
	Array1Du rowBlocks;
	cl_uint numRowBlocks;

	int globalSizeMatrixVector, globalSizeVectorVector;

#ifdef _DEBUG
	Array2Dd bVec_copy, xVec_copy;
#endif

	// IA and JA arrays, pointer used to allow for same CSR_Inds to be shared across multiple
	// instances, but right now, its not currently working, so use a different one for each
	CSR_Inds* Inds;


	// A matrix
	Array1Dd Amat;

	Array1Dd bVec; // also used as storage for r and s vectors during iterations

	// Actual values of macroscopic variable being solved for (pointer to array in containing class)
	Array2Dd* xVec; // also used as storage for h vector during iterations

	// used for when resetting time after solution blows up
	Array2Dd* xVecPrev;
	bool resetTimeFlag;
	std::function<void(void)> resetTimeFunc;

	Array1Dd Rmat, values_xi_k, rho_0, update_coefficients, norm2TempVals, r_dot_vk_buffer;
	//Array1Dd values_eta_k;

	// Sharing buffers between each solver instance to save on memory. This requires
	// all instances to share the same matrix/vector sizes and does not allow for
	// solving operations to overlap.
	static cl_mem device_inner_prod_buffer, resVec, krylov_basis, vi_in_vk_buffer, iniXVal;


	// Using two buffers for each intermediate set (alternates between the two, instead of one per reduce step)
	static cl_mem* reduceBufSet1;
	static bool staticVarInitialized;
	
	static cl_uint krylovDimSize, krylovMatSize, bufSizePerVector, numBufChunks;
	// Using a different kernel for each step will reduce number of steps during calculation,
	// but increase total number of kernels. Not sure if this may cause issues with GPU (wasting
	// memory or possibly reaching limit on number of kernels). 
	// Kernels for sparse solvers. axpbyKernelSingle refers to one used to calculate P vector
	// axpbyKernelDouble is used to calculate s and h and x and r

	// Norm2 of residual for calculation of rho_0
	RedKernelList resNorm2, resNorm2WithNormalize;

	// calculation of RHS at beginning and residual after each restart
	thinKerWrapper updateRHS, updateResidual;

	// scales residual
	thinKerWrapper scaleRes;

	// pipelined gmres product for k = 0 and k > 0
	thinKerWrapper gmresProdK0, gmresProdK;

	// pipelined gmres normalize vk
	thinKerWrapper normalizeVk;

	thinKerWrapper gramSchmidtStage1, gramSchmidtStage2;

	// reduction of r dor vk values
	thinKerWrapper redRDotVk;

	// final step per iteration
	thinKerWrapper updateResult;

	// vclReduce Implementation
	thinKerWrapper vclReduce;

	thinKerWrapper addIniGuess;

	double vclReduceMethod();

#ifdef _DEBUG
	Array1Dd DebugOut, rho0Monitor;
	thinKerWrapper sumKer, iniMonitorResid, curMonitorResid;
	RedKernelList redMonitor;
	double iniResMonitor;
	static cl_mem xMonitor, residMonitor;
#endif

	double avgRes, norm_rhs, prevErr;
	int iterCount;

	void resetAmatrixArguments();
	enum debugSaveType { innerProd, residualVec, krylovBasis, rDotVK, viInVK};

	void debugSaveBuffer(debugSaveType dbgType_, int kval = 0)
	{
		int size_ = 0;
		std::string name_;
		switch (dbgType_)
		{
		case innerProd:
		{
			size_ = numBufChunks * bufSizePerVector;
			DebugOut.read_from_buffer_to_array(device_inner_prod_buffer, TRUE, size_, calcQue);
			name_ = "vlbm_InnerProdVec_" + std::to_string(kval);
			DebugOut.save2fileSize(size_, name_);
			break;
		}
		case residualVec:
		{
			size_ = fullSize;
			DebugOut.read_from_buffer_to_array(resVec, TRUE, size_, calcQue);
			name_ = "vlbm_residual";
			DebugOut.save2fileSize(size_, name_);
			break;
		}
		case krylovBasis:
		{
			size_ = colSize * (kval+1);
			DebugOut.read_from_buffer_to_array(krylov_basis, TRUE, size_, calcQue);
			name_ = "vlbm_krylov_basis";
			DebugOut.save2fileSize2D(fullSize, kval+1, name_);
			break;
		}
		case rDotVK:
		{
			size_ = (kval+1) * bufSizePerVector;
			r_dot_vk_buffer.read_from_buffer();
			name_ = "vlbm_RdotVk";
			r_dot_vk_buffer.save2fileSize2D(bufSizePerVector, kval+1, name_);
			break;
		}
		case viInVK:
		{
			size_ = (kval+1) * bufSizePerVector;
			DebugOut.read_from_buffer_to_array(vi_in_vk_buffer, TRUE, size_, calcQue);
			name_ = "vlbm_viInVk";
			DebugOut.save2fileSize2D(bufSizePerVector, kval+1, name_);
			break;
		}
		}
	}

	void iniMonitorResidual()
	{
		bVec_copy.enqueue_copy_to_buffer_blocking(bVec.get_buffer());
		xVec_copy.enqueue_copy_to_buffer_blocking(xVec->get_buffer());
		iniMonitorResid(calcQue);
		runReduce(redMonitor, calcQue);
		rho0Monitor.read_from_buffer(calcQue);
		iniResMonitor = rho0Monitor(0);
	}

	void debugPrintResidual(double prevResid_)
	{
		sumKer(calcQue);
		curMonitorResid(calcQue);
		runReduce(redMonitor, calcQue);
		rho0Monitor.read_from_buffer(calcQue);
		std::cout << "Residual estimate vs. true residual: " << prevResid_ <<
			" vs. " << rho0Monitor(0) / iniResMonitor << std::endl;

	}

	GMRESSolver() 
	{
	}

	~GMRESSolver()
	{
		if (staticVarInitialized)
		{
			FREE_OCL_BUFFER(device_inner_prod_buffer);
			FREE_OCL_BUFFER(resVec);
			FREE_OCL_BUFFER(krylov_basis);
			FREE_OCL_BUFFER(vi_in_vk_buffer)
			FREE_OCL_BUFFER_ARRAY(reduceBufSet1, 2);
#ifdef _DEBUG
			FREE_OCL_BUFFER(xMonitor);
			FREE_OCL_BUFFER(residMonitor);
#endif
			delete[] globalWorkSizeRed;
			delete[] localWorkSizeRed;
			staticVarInitialized = false;
		}
	}

	void callScaleRes()
	{
		scaleRes.setArgument(2, &rho_0[0]);
		scaleRes(calcQue);
	}

	void callReduceVI(cl_uint kval)
	{
		if (kval == 0)
			return;
		FINISH_QUEUES;
		r_dot_vk_buffer.read_from_buffer_size(kval * bufSizePerVector);
		for (int i = 0; i < (int)kval; ++i)
		{
			values_xi_k[i] = 0.0;
			for (int j = 0; j < (int)bufSizePerVector; ++j)
				values_xi_k[i] += r_dot_vk_buffer[i * (int)bufSizePerVector + j];
		}


		//redRDotVk.setSizes(64 * kval, 64);
		//redRDotVk(calcQue);
		//values_xi_k.read_from_buffer(calcQue);
	}

	void setMaxIters(int iters_) { maxIters = iters_; }
	void setAbsTol(int tol_) { absTol = tol_; }
	void setRelTol(int tol_) { relTol = tol_; }
	void getRowBlockInfo();
	// Fills diagonal elements of rows associated with solid boundary
// nodes with 1 (otherwise solver will encounter div by 0)
	void fillSolidBoundaryNodes();

	// y = alpha*A*x + beta * y;
	void createKernels();

	static void iniBuffer(cl_mem& buf_, const int size_, std::string name_);

	// initializes all static variables
	static void ini(int xsize_, int xsizefull_, int ysize_);

	//size_t globalSizeCSRMV, globalSizeAXPBY;
	void CreateSolver(Array2Dd* macro_, CSR_Inds* inds_, cl_command_queue* calcque_, int maxiters_, double reltol_, double abstol_);

	//size_t globalSizeCSRMV, globalSizeAXPBY;
	void CreateSolverWithVarTime(Array2Dd* macro_, Array2Dd* macroPrev_, std::function<void(void)>& resetTimeFunc_,
		CSR_Inds* inds_, cl_command_queue* calcque_, int maxiters_, double reltol_, double abstol_);


	void runReduce(RedKernelList& redlist_, cl_command_queue* que_, int num_wait = 0,
		cl_event* wait = NULL, cl_event* evt = NULL);

	bool reduceAndCheckConvergence(cl_command_queue* que_, bool setInitialRes = false, int num_wait = 0,
		cl_event* wait = NULL, cl_event* evt = NULL);

	bool checkConvergence();
	
	void fillArraysWithInitialValues();

	void copy_buffers(cl_mem* src_buf, cl_mem* dest_buf);

	void copyToPrevSolution();

	void copyFromPrevSolution();

	void solve();

	// Functions for transfering between host and device and outputting data

	//void setInitialValue(double inival, bool fullArrFlag = false);

	//void setInitialValueRows(double inival, std::vector<int> &rowi);
	//void setInitialValueCols(double inival, std::vector<int> &coli);

	double& GMRESSolver::operator()(int i, int j)
	{
		return xVec->operator()(i, j);
	}

	void calcNorm2Residual()
	{
		runReduce(resNorm2, calcQue);
		rho_0.read_from_buffer(calcQue);
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
	bool saveCSR_As_MM(std::string outname, bool fromDevFlag = false);
	bool saveVec_As_MM(double* arr_, std::string outname, bool fromDevFlag = false);

	bool saveAxb_w_indicies_from_device_as_bin();
	bool save_w_indicies_as_bin(std::string Name, bool fromDevFlag);
	void saveTempArray(cl_mem* mem_, std::string savename_)
	{
		Array1Dd tempArr(fullSize, savename_);
		tempArr.read_from_buffer_to_array(*mem_);
		tempArr.savebin();
	}


	void debugCheckForNans();
	void debugTestResidual(double residual_);
};


#endif