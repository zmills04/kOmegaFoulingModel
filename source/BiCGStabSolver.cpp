#include "BiCGStabSolver.h"
//#include "clVariablesLS.h"
#include "clProblem.h"
#include "SourceGenerator.h"
#include "ReduceGenerator.h"
//#include "BiCGStabGenerator.h"



int BiCGStabSolver::fullSize = 0;
int BiCGStabSolver::colSize = 0; 
int BiCGStabSolver::Xsize = 0;
int BiCGStabSolver::Ysize = 0;
int BiCGStabSolver::XsizeFull = 0; 
int BiCGStabSolver::numRedKer = 0;
int* BiCGStabSolver::globalWorkSizeRed = nullptr;
cl_mem BiCGStabSolver::r0Hat = nullptr;
cl_mem BiCGStabSolver::pVec = nullptr;
cl_mem BiCGStabSolver::vVec = nullptr;
cl_mem BiCGStabSolver::tVec = nullptr;

cl_mem* BiCGStabSolver::reduceBufSet1 = nullptr;
cl_mem* BiCGStabSolver::reduceBufSet2 = nullptr;
bool BiCGStabSolver::staticVarInitialized = false;

void BiCGStabSolver::iniAXPBYKernels()
{
	// r -= alpha*v (r becomes s) -> pY -= alpha*pM
	// x += alpha*p (x becomes h) -> pZ += alpha*pP
	int alphaInd = indAlpha;
	shAXPBY.setSizes(globalSizeAXPBY, WORKGROUPSIZE_AXPBY);
	shAXPBY.createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "axpby_shCalc");
	shAXPBY.setArgument(0, &colSize);
	shAXPBY.setArgument(1, bVec.get_buf_add());
	shAXPBY.setArgument(2, scalarBuf.get_buf_add());
	shAXPBY.setArgument(3, &alphaInd);
	shAXPBY.setArgument(4, xVec->get_buf_add());
	shAXPBY.setArgument(5, &vVec);
	shAXPBY.setArgument(6, &pVec);


	// h += omega*s (h becomes x) -> pZ += alpha*pY
	// x -= omega*s (r becomes s) -> pY -= alpha*pM
	int omegaInd = indOmega;
	xrAXPBY.setSizes(globalSizeAXPBY, WORKGROUPSIZE_AXPBY);
	xrAXPBY.createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "axpby_xrCalc");
	xrAXPBY.setArgument(0, &colSize);
	xrAXPBY.setArgument(1, bVec.get_buf_add());
	xrAXPBY.setArgument(2, scalarBuf.get_buf_add());
	xrAXPBY.setArgument(3, &omegaInd);
	xrAXPBY.setArgument(4, xVec->get_buf_add());
	xrAXPBY.setArgument(5, &tVec);

	// In kernel, variable alpha is beta
	int betaInd = indBeta;
	pAXPBY.setSizes(globalSizeAXPBY, WORKGROUPSIZE_AXPBY);
	pAXPBY.createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "axpby_PCalc");
	pAXPBY.setArgument(0, &colSize);
	pAXPBY.setArgument(1, &pVec);
	pAXPBY.setArgument(2, scalarBuf.get_buf_add());
	pAXPBY.setArgument(3, &betaInd);
	pAXPBY.setArgument(4, &omegaInd);
	pAXPBY.setArgument(5, bVec.get_buf_add());
	pAXPBY.setArgument(6, &vVec);
}


void BiCGStabSolver::createElementwiseKernels()
{
	subKer.setSizes(globalSizeAXPBY, WORKGROUPSIZE_AXPBY);
	subKer.createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "ElementWise_Subtract");
	subKer.setArgument(0, colSize);

	addKer.setSizes(globalSizeAXPBY, WORKGROUPSIZE_AXPBY);
	addKer.createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "ElementWise_Add");
	addKer.setArgument(0, colSize);

	multKer.setSizes(globalSizeAXPBY, WORKGROUPSIZE_AXPBY);
	multKer.createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "ElementWise_Multiply");
	multKer.setArgument(0, colSize);

	divKer.setSizes(globalSizeAXPBY, WORKGROUPSIZE_AXPBY);
	divKer.createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "ElementWise_Divide");
	divKer.setArgument(0, colSize);

}



// y = alpha*A*x + beta * y;
void BiCGStabSolver::createCSRMVKernel(thinKerWrapper &ker, double alphaval, double betaval, cl_mem *xvec_, cl_mem* yvec_)
{
	ker.setSizes(globalSizeCSRMV, WORKGROUPSIZE_CSRMV);
	ker.createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "csrmv_adaptive");
	ker.setArgument(0, Amat.get_buf_add());
	ker.setArgument(1, Inds->JA.get_buf_add());
	ker.setArgument(2, Inds->IA.get_buf_add());
	ker.setArgument(3, xvec_);
	ker.setArgument(4, yvec_);
	ker.setArgument(5, &rowBlocks);
	ker.setArgument<double>(6, alphaval);
	ker.setArgument<double>(7, betaval);
}

void BiCGStabSolver::iniCSRMVKernels()
{
	// b = b - A*x_0   -> b become r0
	createCSRMVKernel(rCSRMV, -1., 1., xVec->get_buf_add(), bVec.get_buf_add());

	// v = A*p
	createCSRMVKernel(vCSRMV, 1., 0., &pVec, &vVec);

	// t = A*s
	createCSRMVKernel(tCSRMV, 1., 0., bVec.get_buf_add(), &tVec);
}

// initializes kernels and sets arguments for everthing but input buffer(s) in first kernel,
// and output, local memory and offset values in last kernel.
// startInd is the index of the intermediate buffer in the first kernel (1 for all but dot reduce)
void BiCGStabSolver::createReduceListBases(RedKernelList &redlist, const int startInd, std::string base1ker,
	std::string base2ker, std::string finalker, cl_mem* mem_use)
{
	int iBufOut = 0;
	redlist.ini(numRedKer);
	redlist(0).setSizes(globalWorkSizeRed[0], WORKGROUPSIZE_RED);
	redlist(0).createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), base1ker);
	redlist(0).setArgument(startInd, &mem_use[0]);
	redlist(0).setLocalMem(startInd + 1, WORKGROUPSIZE_RED*sizeof(double));

	for (int i = 1; i < numRedKer - 1; i++)
	{
		redlist(i).setSizes(globalWorkSizeRed[i], WORKGROUPSIZE_RED);
		redlist(i).createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), base2ker);
		redlist(i).setArgument(0, &mem_use[iBufOut]);
		iBufOut ^= 1;
		redlist(i).setArgument(1, &mem_use[iBufOut]);
		redlist(i).setLocalMem(2, WORKGROUPSIZE_RED*sizeof(double));
	}
	redlist(numRedKer - 1).setSizes(globalWorkSizeRed[numRedKer - 1], globalWorkSizeRed[numRedKer - 1]);
	redlist(numRedKer - 1).createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), finalker);
	redlist(numRedKer - 1).setArgument(0, &mem_use[iBufOut]);
}

void BiCGStabSolver::createBNorm()
{
	createReduceListBases(bNorm, 1, "ReduceKernelBase_Norm",
		"ReduceKernelBase_Sum", "ReduceKernel_Sum_Final", reduceBufSet1);

	bNorm(0).setArgument(0, bVec.get_buf_add());
	bNorm(numRedKer - 1).setArgument(1, scalarBuf.get_buf_add());
	bNorm(numRedKer - 1).setLocalMem(2, globalWorkSizeRed[numRedKer - 1] * sizeof(double));
	bNorm.setValueOutputIndex(3);
}

void BiCGStabSolver::createRhoSqr()
{
	createReduceListBases(rhoSqr, 1, "ReduceKernelBase_Sqr",
		"ReduceKernelBase_Sum", "ReduceKernel_Sum_Final", reduceBufSet2);
	rhoSqr(0).setArgument(0, bVec.get_buf_add());
	rhoSqr(numRedKer - 1).setArgument(1, scalarBuf.get_buf_add());
	rhoSqr(numRedKer - 1).setLocalMem(2, globalWorkSizeRed[numRedKer - 1] * sizeof(double));
	int scalarInd = indRho;
	rhoSqr(numRedKer - 1).setArgument(3, &scalarInd);
}

void BiCGStabSolver::createRhoDot()
{
	createReduceListBases(rhoDot, 2, "ReduceKernelBase_Dot",
		"ReduceKernelBase_Sum", "ReduceKernel_Sum_Rho_Calc_Beta", reduceBufSet2);

	rhoDot(0).setArgument(0, bVec.get_buf_add());
	rhoDot(0).setArgument(1, &r0Hat);

	rhoDot(numRedKer - 1).setArgument(1, scalarBuf.get_buf_add());
	rhoDot(numRedKer - 1).setLocalMem(2, globalWorkSizeRed[numRedKer - 1] * sizeof(double));

	int scalarInd = indRho;
	rhoDot(numRedKer - 1).setArgument(3, &scalarInd);
	scalarInd = indBeta;
	rhoDot(numRedKer - 1).setArgument(4, &scalarInd);
	scalarInd = indAlpha;
	rhoDot(numRedKer - 1).setArgument(5, &scalarInd);
	scalarInd = indOmega;
	rhoDot(numRedKer - 1).setArgument(6, &scalarInd);
}

void BiCGStabSolver::createAlphaDot()
{
	createReduceListBases(alphaDot, 2, "ReduceKernelBase_Dot",
		"ReduceKernelBase_Sum", "ReduceKernel_Calc_Alpha", reduceBufSet1);

	alphaDot(0).setArgument(0, &vVec);
	alphaDot(0).setArgument(1, &r0Hat);

	alphaDot(numRedKer - 1).setArgument(1, scalarBuf.get_buf_add());
	alphaDot(numRedKer - 1).setLocalMem(2, globalWorkSizeRed[numRedKer - 1] * sizeof(double));

	int scalarInd = indRho;
	alphaDot(numRedKer - 1).setArgument(3, &scalarInd);
	scalarInd = indAlpha;
	alphaDot(numRedKer - 1).setArgument(4, &scalarInd);
}

void BiCGStabSolver::createOmegaDotSqr()
{
	int iBufOut = 0;
	omegaDotSqr.ini(numRedKer);
	omegaDotSqr(0).setSizes(globalWorkSizeRed[0], WORKGROUPSIZE_RED);
	omegaDotSqr(0).createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "ReduceKernel_Base_Omega");
	omegaDotSqr(0).setArgument(0, bVec.get_buf_add());
	omegaDotSqr(0).setArgument(1, &tVec);
	omegaDotSqr(0).setArgument(2, &reduceBufSet1[iBufOut]);
	omegaDotSqr(0).setArgument(3, &reduceBufSet2[iBufOut]);
	omegaDotSqr(0).setLocalMem(4, WORKGROUPSIZE_RED*sizeof(double));
	omegaDotSqr(0).setLocalMem(5, WORKGROUPSIZE_RED*sizeof(double));

	for (int i = 0; i < numRedKer - 1; i++)
	{
		omegaDotSqr(i).setSizes(globalWorkSizeRed[i], WORKGROUPSIZE_RED);
		omegaDotSqr(i).createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "ReduceKernel_Intermediate_Omega");
		omegaDotSqr(i).setArgument(0, &reduceBufSet1[iBufOut]);
		omegaDotSqr(i).setArgument(1, &reduceBufSet2[iBufOut]);
		iBufOut ^= 1;
		omegaDotSqr(i).setArgument(2, &reduceBufSet1[iBufOut]);
		omegaDotSqr(i).setArgument(3, &reduceBufSet2[iBufOut]);
		omegaDotSqr(i).setLocalMem(4, WORKGROUPSIZE_RED*sizeof(double));
		omegaDotSqr(i).setLocalMem(5, WORKGROUPSIZE_RED*sizeof(double));
	}

	omegaDotSqr(numRedKer - 1).setSizes(globalWorkSizeRed[numRedKer - 1], globalWorkSizeRed[numRedKer - 1]);
	omegaDotSqr(numRedKer - 1).createKernel(BiCGStabGenerator::BiCGStabInstance()->getProgram(), "ReduceKernel_Calc_Omega");
	omegaDotSqr(numRedKer - 1).setArgument(0, &reduceBufSet1[iBufOut]);
	omegaDotSqr(numRedKer - 1).setArgument(1, &reduceBufSet2[iBufOut]);
	omegaDotSqr(numRedKer - 1).setArgument(2, scalarBuf.get_buf_add());
	omegaDotSqr(numRedKer - 1).setLocalMem(3, globalWorkSizeRed[numRedKer - 1] * sizeof(double));
	omegaDotSqr(numRedKer - 1).setLocalMem(4, globalWorkSizeRed[numRedKer - 1] * sizeof(double));
	int scalarInd = indOmega;
	omegaDotSqr(numRedKer - 1).setArgument(5, &scalarInd);
}

void BiCGStabSolver::iniBuffer(cl_mem &buf_, const int size_, std::string name_)
{
	int status;
	double zer = 0.;
	buf_ = clCreateBuffer(*clEnv::instance()->getContext(), CL_MEM_READ_WRITE, sizeof(double) * size_, NULL, &status);
	ERROR_CHECKING_OCL(status, "Error creating " + name_ + " in BiCGStabSolver Class", ERROR_CREATING_BICGSTAB_SOLVER);

	status = clEnqueueFillBuffer(*clEnv::instance()->getIOqueue(), buf_, &zer,
		sizeof(double), 0, sizeof(double)*size_, 0, NULL, NULL);
	ERROR_CHECKING_OCL(status, "Error filling " + name_ + " in BiCGStabSolver Class", ERROR_CREATING_BICGSTAB_SOLVER);
}

// initializes all static variables
void BiCGStabSolver::ini(int xsize_, int xsizefull_, int ysize_)
{
	if (staticVarInitialized)
		return;
	Xsize = xsize_;
	Ysize = ysize_;
	XsizeFull = xsizefull_;
	colSize = XsizeFull*Ysize;
	int iSize1 = 0, iSize2 = 0;

	fullSize = ReduceGenerator::ReduceInstance()->iniSizesAndMemory(colSize, numRedKer,
		&globalWorkSizeRed, &reduceBufSet1, iSize1, iSize2);

	ReduceGenerator::ReduceInstance()->iniDeviceMemory(numRedKer, &reduceBufSet2, iSize1, iSize2);
	iniBuffer(r0Hat, fullSize, "r0Hat");
	iniBuffer(pVec, fullSize, "pVec");
	iniBuffer(vVec, fullSize, "vVec");
	iniBuffer(tVec, fullSize, "tVec");

	staticVarInitialized = true;
}

// creates CSR matrix, generates meta information, copies
// it over and deletes clsparse objects since they are no longer needed.
// clSparse control object created in BiCGStabGenerator and destroyed
// when sourceGenerator compiles its program (final step in initialization)
void BiCGStabSolver::getCSRMeta()
{
	int status;
	// generate meta info for adaptive csrmv using clSparse
	clsparseCsrMatrix csrMat;
	clsparseInitCsrMatrix(&csrMat);

	csrMat.num_nonzeros = Inds->nnz();
	csrMat.num_rows = Inds->rows;
	csrMat.num_cols = Inds->cols;

	csrMat.values = Amat.get_buffer();
	csrMat.col_indices = Inds->JA.get_buffer();
	csrMat.row_pointer = Inds->IA.get_buffer();

	clsparseCsrMetaCreate(&csrMat, BiCGStabGenerator::BiCGStabInstance()->clSparseControl.control);
	matrix_meta* meta_ptr = static_cast< matrix_meta* >(csrMat.meta);
	rowBlockSize = meta_ptr->rowBlockSize;
	offRowBlocks = meta_ptr->offRowBlocks;

	// Creating own buffer to store rowBlocks array in device memory, since the clsparse object is temporary
	rowBlocks = clCreateBuffer(*clEnv::instance()->getContext(), CL_MEM_READ_WRITE, rowBlockSize*sizeof(cl_ulong), NULL, &status);
	ERROR_CHECKING_OCL(status, "error creating rowBlocks buffer in BiCGStab", ERROR_CREATING_BICGSTAB_SOLVER);

	// copying data into new rowBlocks buffer. Now clsparse objects are unneeded.
	status = clEnqueueCopyBuffer(*clEnv::instance()->getIOqueue(), meta_ptr->rowBlocks(), rowBlocks, 0, 0, rowBlockSize*sizeof(cl_ulong), 0, NULL, NULL);
	ERROR_CHECKING_OCL(status, "error copying rowBlocks into BiCGStab class buffer object", ERROR_CREATING_BICGSTAB_SOLVER);

	// just to ensure that meta owned buffer is freed
	meta_ptr->clear();
}

void BiCGStabSolver::CreateSolverWithVarTime(Array2Dd* macro_, Array2Dd* macroPrev_, std::function<void(void)>& resetTimeFunc_,
	CSR_Inds* inds_, cl_command_queue* calcque_, int maxiters_, double reltol_, double abstol_)
{
	resetTimeFunc = resetTimeFunc_;
	xVecPrev = macroPrev_;
	resetTimeFlag = true;
	CreateSolver(macro_, inds_, calcque_, maxiters_, reltol_, abstol_);
}

void BiCGStabSolver::setArraysWithoutCreation(Array2Dd* macro_, CSR_Inds* inds_, cl_command_queue* calcque_,
	int maxiters_, double reltol_, double abstol_)
{
	calcQue = calcque_;
	Name = macro_->getName() + "Solver";
	Inds = inds_;
	xVec = macro_;

	Xsize = Inds->Xsize;
	Ysize = Inds->Ysize;
	XsizeFull = Inds->XsizeFull;
	colSize = XsizeFull * Ysize;

	NNZ = Inds->nnz();
	maxIters = maxiters_;
	relTol = reltol_;
	absTol = abstol_;


	ini(Inds->Xsize, Inds->XsizeFull, Inds->Ysize);

	NNZ = Inds->nnz();
	maxIters = maxiters_;
	relTol = reltol_;
	absTol = abstol_;


	std::string Aname = Name + "_Amat";
	std::string bname = Name + "_bVec";
	std::string xname = Name + "_xVec";
	std::string scalarname = Name + "_Scalar";

	// Only want to call these once when a CSR_Inds instance is shared
	// between BiCGStabSolver instances
	std::string ianame = Name + "_IA";
	std::string janame = Name + "_JA";
	std::string raname = Name + "_RA";

	Amat.zeros(Inds->nnz());
	Amat.setName(Aname);
	fillSolidBoundaryNodes();

	Inds->IA.setName(ianame);
	Inds->JA.setName(janame);
	Inds->RA.setName(raname);

	Inds->allocate_buffers(*clEnv::instance()->getContext());
	Inds->copy_RA_to_device(clEnv::instance()->getIOqueue());
	Inds->copy_to_device(clEnv::instance()->getIOqueue());
	Inds->iniFlag = true;
	
	
	bVec.zeros(fullSize);
	bVec.setName(bname);
	bVec.allocate_buffer_w_copy();
	Amat.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
	xVec->allocate_buffer_size(p.FullSize);
	xVec->copy_to_buffer();


	vclA.setWithExistingBuffers(Inds->IA.get_buffer(), Inds->JA.get_buffer(),
		Amat.get_buffer(), Inds->rows, Inds->cols, Inds->nnz(), viennacl::ocl::current_context());

	vclX.setWithExistingBuffer(xVec->get_buffer(), Inds->rows, 0, 0, viennacl::ocl::current_context());
	vclB.setWithExistingBuffer(bVec.get_buffer(), Inds->rows, 0, 0, viennacl::ocl::current_context());
}


//size_t globalSizeCSRMV, globalSizeAXPBY;
void BiCGStabSolver::CreateSolver(Array2Dd *macro_, CSR_Inds *inds_, cl_command_queue *calcque_, int maxiters_, double reltol_, double abstol_)
{
	calcQue = calcque_;
	Name = macro_->getName() + "Solver";
	Inds = inds_;
	xVec = macro_;

	// Just incase the static variables have not been
	// initialized.
	ini(Inds->Xsize, Inds->XsizeFull, Inds->Ysize);
	
	NNZ = Inds->nnz();
	maxIters = maxiters_;
	relTol = reltol_;
	absTol = abstol_;


	std::string Aname = Name + "_Amat";
	std::string bname = Name + "_bVec";
	std::string xname = Name + "_xVec";
	std::string scalarname = Name + "_Scalar";

#ifdef _DEBUG
	std::string bVecName_copy = Name + "_bVec_copy";
	std::string xVecName_copy = Name + "_xVec_copy";
	bVec_copy.zeros(Xsize, XsizeFull, Ysize, Ysize);
	xVec_copy.zeros(Xsize, XsizeFull, Ysize, Ysize);
	bVec_copy.setName(bVecName_copy);
	xVec_copy.setName(xVecName_copy);
	bVec_copy.allocate_buffer_w_copy();
	xVec_copy.allocate_buffer_w_copy();
#endif


	// Only want to call these once when a CSR_Inds instance is shared
	// between BiCGStabSolver instances
	if (!(Inds->iniFlag))
	{
		std::string ianame = Name + "_IA";
		std::string janame = Name + "_JA";
		std::string raname = Name + "_RA";

		Inds->IA.setName(ianame);
		Inds->JA.setName(janame);
		Inds->RA.setName(raname);

		Inds->allocate_buffers(*clEnv::instance()->getContext());
		Inds->copy_RA_to_device(clEnv::instance()->getIOqueue());
		Inds->copy_to_device(clEnv::instance()->getIOqueue());
		Inds->iniFlag = true;
	}
	Amat.zeros(Inds->nnz());
	Amat.setName(Aname);
	fillSolidBoundaryNodes();
	
	Amat.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
	scalarBuf.setName(scalarname);
	scalarBuf.fill(1.);
	scalarBuf.allocate_buffer_w_copy();
	
	bVec.zeros(fullSize);
	bVec.setName(bname);
	bVec.allocate_buffer_w_copy();

	// Will be reducing this possibly, so allocating correct size
	xVec->allocate_buffer_size(p.FullSize);
	xVec->copy_to_buffer();

	if (resetTimeFlag)
	{
		xVecPrev->allocate_buffer_size(p.FullSize);
		xVecPrev->FreeHost();
	}

	// AXPBY kernel sizes
	int blocksNum = (colSize + WORKGROUPSIZE_AXPBY - 1) / WORKGROUPSIZE_AXPBY;
	globalSizeAXPBY = size_t(blocksNum * WORKGROUPSIZE_AXPBY);

	getCSRMeta();

	// csrmv kernel sizes
	globalSizeCSRMV = size_t(((rowBlockSize / 2) - 1) * WORKGROUPSIZE_CSRMV);

	createOmegaDotSqr();
	createAlphaDot();
	createRhoDot();
	createRhoSqr();
	createBNorm();
	iniCSRMVKernels();
	iniAXPBYKernels();

}


void BiCGStabSolver::fillSolidBoundaryNodes()
{
	for (int i = 0; i < Inds->Xsize; i++)
	{
		for (int j = 0; j < Inds->Ysize; j++)
		{
			if (Inds->testSolidBoundaryNode(i, j))
			{
				Amat(Inds->getInd(i, j, CSR_Inds::C)) = 1.;
			}
		}
	}
}

void BiCGStabSolver::runReduce(RedKernelList &redlist_, cl_command_queue* que_, int num_wait,
	cl_event *wait, cl_event* evt)
{
	for (thinKerWrapper *rk = redlist_.begin(); rk != redlist_.end(); ++rk)
	{
		int status = rk->operator()(que_, num_wait, wait, evt);
		clEnv::instance()->finishQueues();
		ERROR_CHECKING(status, ("Error returned running reduce kernels of " + Name), status);
	}
}

bool BiCGStabSolver::reduceAndCheckConvergence(cl_command_queue *que_, bool setInitialRes, int num_wait,
	cl_event *wait, cl_event* evt)
{
	runReduce(bNorm, que_, num_wait, wait, evt);
	scalarBuf.read_from_buffer_size(2);
	double residual = scalarBuf(1) / scalarBuf(0);
	//avgRes = (avgRes * (double)(iterCount - 1) + residual) / ((double)iterCount);

	if (setInitialRes)
	{
		initialResidual = residual;
	}


	if (resetTimeFlag)
	{
		if (isnan(residual) || (iterCount > 1 && residual > 10. * avgRes))
		{
			resetTimeFunc();
			return true;
		}
	}
#ifdef _DEBUG
	debugTestResidual(residual);
#endif

#ifdef PRINT_BICGSTAB_RESIDUALS
		printf("Iteration %d Residual = %g\n", iterCount, residual);
#endif

	if (residual <= relTol || residual <= absTol * initialResidual)
	{
		return true;
	}
	avgRes = residual;
	return false;
}


void BiCGStabSolver::copy_buffers(cl_mem *src_buf, cl_mem *dest_buf)
{
	int status = clEnqueueCopyBuffer(*clEnv::instance()->getIOqueue(),
		*src_buf, *dest_buf, 0, 0, colSize*sizeof(double), 0, NULL, NULL);
	ERROR_CHECKING(status, ("Error copying kernels in BiCGStabSolver " + Name), status);
}

void BiCGStabSolver::printScalars()
{
	using std::cout;
	static bool namesPrinted = false;
	if (!namesPrinted)
	{
		cout << "\nScalar Array: [NormB, NormR, NormS, Rho, Omega, Alpha, Beta]\n";
	}

	scalarBuf.read_from_buffer();
	cout << "Scalar Array: [" << scalarBuf(0) << ", " << scalarBuf(1) << ", ";
	cout << scalarBuf(2) << ", " << scalarBuf(3) << ", " << scalarBuf(4);
	cout << ", " << scalarBuf(5) << ", " << scalarBuf(6) << "]\n";


}

void BiCGStabSolver::copyToPrevSolution()
{
	if (resetTimeFlag)
	{
		xVecPrev->enqueue_copy_to_buffer_blocking(xVec->get_buffer());
	}
}

void BiCGStabSolver::copyFromPrevSolution()
{
	if (resetTimeFlag)
	{
		xVec->enqueue_copy_to_buffer_blocking(xVecPrev->get_buffer());
	}
}

void BiCGStabSolver::debugTestResidual(double residual_)
{
	if (isnan(residual_))
	{
		bVec_copy.save_txt_from_device();
		xVec_copy.save_txt_from_device();
		vlb.saveDebug();
	}
}

void BiCGStabSolver::debugCheckForNans()
{
	bVec_copy.enqueue_copy_to_buffer_blocking(bVec.get_buffer());
	xVec_copy.enqueue_copy_to_buffer_blocking(xVec->get_buffer());
	if (p.Time > 40000)
	{
		if (bVec.checkForNans() || xVec->checkForNans() || Amat.checkForNans())
		{
			int iout, jout, kout;
			if (bVec.checkForNans(iout, jout, kout, false))
			{
				std::cout << "nan in bVec at (" << iout << ", " <<
					jout << ", " << kout << ") at time" << p.Time << std::endl;
			}
			if (Amat.checkForNans(iout, jout, kout, false))
			{
				std::cout << "nan in Amat at (" << iout << ", " <<
					jout << ", " << kout << ") at time" << p.Time << std::endl;
			}
			if (xVec->checkForNans(iout, jout, kout, false))
			{
				std::cout << "nan in xVec at (" << iout << ", " <<
					jout << ", " << kout << ") at time" << p.Time << std::endl;
			}
			vlb.saveDebug();
		}
	}
}


void BiCGStabSolver::solve()
{
#ifdef _DEBUG
	//debugCheckForNans();
#endif

	avgRes = 0.;
	iterCount = 0;

	scalarBuf.FillBuffer(1.);
	bNorm.setOutputIndex(indNormB);

	// Step 1
	runReduce(bNorm, calcQue);
	clFlush(*calcQue);
	
	//using save Queue to read to avoid reading before finishing
	// blocking will ensure checking after reading
	scalarBuf.read_from_buffer_size(1);

	double h_norm_b = scalarBuf(0);

#ifdef _DEBUG
	debugTestResidual(h_norm_b);
#endif

	if (h_norm_b <= 1.0e-16)
	{
		copy_buffers(bVec.get_buf_add(), xVec->get_buf_add());
		return;
	}

	bNorm.setOutputIndex(indNormR);

	// Step 2
	rCSRMV(calcQue);
	clFinish(*calcQue);

	// Step 3
	copy_buffers(bVec.get_buf_add(), &r0Hat);
	copy_buffers(bVec.get_buf_add(), &pVec);

	// Step 4
	runReduce(rhoSqr, calcQue);


	bool firstIter = true;
	while (iterCount < maxIters)
	{
		if (iterCount > 0)
		{
			// Step f1, f2
			runReduce(rhoDot, calcQue);


			// Step f3
			pAXPBY(calcQue);
		}

		// Step f4
		vCSRMV(calcQue);

		// Step f5
		runReduce(alphaDot, calcQue);

		// Step f6
		shAXPBY(calcQue);

		// Step f7

		if (reduceAndCheckConvergence(calcQue, firstIter))
		{
			//converged, so h is correct result and already stored in x
			return;
		}

		firstIter = false;

		// Step f8
		tCSRMV(calcQue);

		// Step f9
		runReduce(alphaDot, calcQue);

		// Step f10
		xrAXPBY(calcQue);

		if (reduceAndCheckConvergence(calcQue, firstIter))
		{
			//converged, so x is correct result
			return;
		}

		iterCount++;
	}
}

//void BiCGStabSolver::setInitialValue(double inival, bool fullArrFlag)
//{
//	if (fullArrFlag)
//	{
//		xVec->FillBuffer(inival);
//		return;
//	}
//
//	for (int i = 0; i < Xsize; i++)
//	{
//		for (int j = 0; j < Ysize; j++)
//		{
//			if (vls.nType(i, j) & Inds->fluidFlag)
//			{
//				xVec->operator()(i, j) = inival;
//			}
//		}
//	}
//
//	xVec->copy_to_buffer();
//}

//void BiCGStabSolver::setInitialValueRows(double inival, std::vector<int> &rowi)
//{
//	xVec->read_from_buffer(NULL, CL_TRUE);
//	for (int i = 0; i < Xsize; i++)
//	{
//		for (int jj = 0; jj < rowi.size(); jj++)
//		{
//			int j = rowi[jj];
//			if (vls.nType(i, j) & Inds->fluidFlag)
//			{
//				xVec->operator()(i, j) = inival;
//			}
//		}
//	}
//	xVec->copy_to_buffer();
//}
//
//void BiCGStabSolver::setInitialValueCols(double inival, std::vector<int> &coli)
//{
//	xVec->read_from_buffer(NULL, CL_TRUE);
//	for (int ii = 0; ii < coli.size(); ii++)
//	{
//		int i = coli[ii];
//		for (int j = 0; j < Ysize; j++)
//		{
//			if (vls.nType(i, j) & Inds->fluidFlag)
//			{
//				xVec->operator()(i, j) = inival;
//			}
//		}
//	}
//
//	xVec->copy_to_buffer(NULL, CL_TRUE);
//}



int BiCGStabSolver::getBufferFullSize()
{
	return fullSize;
}

std::string BiCGStabSolver::getName()
{
	return Name;
}

///////// Load and CheckPoint Methods ///////////



//////// Save Methods /////////////////////
bool BiCGStabSolver::saveAxbCSR()
{
	bool ret = savetxt();
	ret &= save_bvec();
	ret &= saveCSR(Name);
	return ret;
}

bool BiCGStabSolver::saveAxb_w_indicies()
{
	bool ret = savetxt();
	ret &= save_bvec();
	ret &= save_w_indicies(Name.append("_A"), false);
	return ret;
}

bool BiCGStabSolver::saveAxbCSR_from_device()
{
	bool ret = savetxt_from_device();
	ret &= save_bvec_from_device();
	ret &= saveCSR(Name, true);
	return ret;
}

bool BiCGStabSolver::saveAxb_w_indicies_from_device()
{
	bool ret = savetxt_from_device();
	ret &= save_bvec_from_device();
	std::string nameOut = Name;
	nameOut.append("_A");
	ret &= save_w_indicies(nameOut, true);
	return ret;
}

bool BiCGStabSolver::saveAxb_w_indicies_from_device_as_bin()
{
	bool ret = xVec->save_bin_from_device(Name);
	ret &= bVec.save_bin_from_device();
	std::string nameOut = Name;
	nameOut.append("_A");
	ret &= save_w_indicies_as_bin(nameOut, true);
	return ret;
}

cl_mem* BiCGStabSolver::get_add_A()
{
	return Amat.get_buf_add();
}

cl_mem* BiCGStabSolver::get_add_IndArr()
{
	return Inds->IndArray.get_buf_add();
}

bool BiCGStabSolver::savetxt(std::string outname)
{
	if (outname.length() == 0)
		outname = Name;
	return xVec->savetxt(outname);
}

bool BiCGStabSolver::savetxt_from_device(std::string outname)
{
	if (outname.length() == 0)
		outname = Name;
	return xVec->save_txt_from_device(outname);
}

bool BiCGStabSolver::save_bvec(std::string outname)
{
	if (outname.length() == 0)
		outname = Name + "_bvec";
	return bVec.savetxt_as_2D(Xsize, XsizeFull, Ysize, outname);
}

bool BiCGStabSolver::save_bvec_from_device(std::string outname)
{
	if (outname.length() == 0)
		outname = Name + "_bvec";
	return bVec.save_txt_from_device_as_2D(Xsize, XsizeFull, Ysize, outname);
}

// Amat (both values and indicies), bVec and other info can be 
// recreated, so no need to waste space and time saving it.
bool BiCGStabSolver::saveCheckPoint(std::string outname)
{
	return xVec->save_bin_from_device(outname);
}

void BiCGStabSolver::copy_to_device(const int blFlag)
{
	Amat.copy_to_buffer(NULL, blFlag);
}

void BiCGStabSolver::copy_to_host(const int blFlag)
{
	Amat.read_from_buffer(NULL, blFlag);
}

void BiCGStabSolver::copy_inds_to_host(const int blFlag)
{
	Inds->copy_to_host(NULL, blFlag);
}

void BiCGStabSolver::copy_inds_to_device(const int blFlag)
{
	Inds->copy_to_device(NULL, blFlag);
}

void BiCGStabSolver::copy_to_host_all(const int blFlag)
{
	copy_to_host(blFlag);
	copy_inds_to_host(blFlag);
}

void BiCGStabSolver::copy_to_device_all(const int blFlag)
{
	copy_to_device(blFlag);
	copy_inds_to_device(blFlag);
}

double BiCGStabSolver::A(const int i, const int j, int dir)
{
	int ind = Inds->getInd(i, j, dir);
	if (ind > -1)
		return Amat(ind);
	else
		return 0;
}

bool BiCGStabSolver::testInd(const int i, const int j, const int dir)
{
	if (Inds->getInd(i, j, dir) == -1)
		return false;
	return true;
}


bool BiCGStabSolver::save_w_indicies(std::string Name, bool fromDevFlag)
{
	if (fromDevFlag)
	{
		copy_to_host_all(true);
	}
	Array2Dd Outarray(nnz(), 3);

	for (int i = 0; i < nnz(); i++)
	{
		Outarray(i, 0) = (double)(Inds->RowIndex(i));
		Outarray(i, 1) = (double)(Inds->ColIndex(i));
		Outarray(i, 2) = Amat(i);
	}

	return Outarray.savetxt(Name);
}

bool BiCGStabSolver::save_w_indicies_as_bin(std::string Name, bool fromDevFlag)
{
	if (fromDevFlag)
	{
		copy_to_host_all(true);
	}
	Array2Dd Outarray(nnz(), 3);

	for (int i = 0; i < nnz(); i++)
	{
		Outarray(i, 0) = (double)(Inds->RowIndex(i));
		Outarray(i, 1) = (double)(Inds->ColIndex(i));
		Outarray(i, 2) = Amat(i);
	}

	return Outarray.savebin(Name);
}




bool BiCGStabSolver::saveCSR(std::string outname, bool fromDevFlag )
{
	bool ret;
	if (fromDevFlag)
	{
		ret = Inds->saveIA(outname);
		ret &= Inds->saveJA(outname);
		outname.append("_A");
		ret &= Amat.save_txt_from_device(outname);
	}
	else
	{
		ret = Inds->saveIA(outname);
		ret &= Inds->saveJA(outname);
		outname.append("_A");
		ret &= Amat.savetxt(outname);
	}
	return ret;
}

bool BiCGStabSolver::saveCSR_row_col_val(std::string outname, bool fromDevFlag)
{
	bool ret;
	if (fromDevFlag)
	{
		ret &= Inds->saveJA(outname);
		ret &= Inds->saveRA(outname);
		outname.append("_A");
		ret = Amat.save_txt_from_device(outname);

	}
	else
	{
		ret = Inds->saveJA(outname);
		ret &= Inds->saveRA(outname);
		outname.append("_A");
		ret = Amat.savetxt(outname);
	}
	return ret;
}