#include "GMRESSolver.h"
//#include "clVariablesLS.h"
#include "clProblem.h"
#include "SourceGenerator.h"
#include "ReduceGenerator.h"
//#include "BiCGStabGenerator.h"

cl_uint GMRESSolver::bufSizePerVector = 128;
cl_uint GMRESSolver::numBufChunks = 3;
cl_uint GMRESSolver::krylovMatSize = 0;
cl_uint GMRESSolver::krylovDimSize = GMRES_KRYLOV_DIM_SIZE;
cl_uint GMRESSolver::fullSize = 0;
cl_uint GMRESSolver::colSize = 0;
cl_uint GMRESSolver::Xsize = 0;
cl_uint GMRESSolver::Ysize = 0;
cl_uint GMRESSolver::XsizeFull = 0;
int GMRESSolver::numRedKer = 0;
int* GMRESSolver::globalWorkSizeRed = nullptr;
int* GMRESSolver::localWorkSizeRed = nullptr;
cl_mem GMRESSolver::device_inner_prod_buffer = nullptr;
cl_mem GMRESSolver::resVec = nullptr;
cl_mem GMRESSolver::krylov_basis = nullptr;
//cl_mem GMRESSolver::r_dot_vk_buffer = nullptr;
cl_mem GMRESSolver::vi_in_vk_buffer = nullptr;
cl_mem GMRESSolver::iniXVal = nullptr;
cl_mem* GMRESSolver::reduceBufSet1 = nullptr;

#ifdef _DEBUG
cl_mem GMRESSolver::xMonitor = nullptr;
cl_mem GMRESSolver::residMonitor = nullptr;
#endif

bool GMRESSolver::staticVarInitialized = false;


class vecWrapper
{
public:
	std::vector<double> hostVec;
	viennacl::vector<double> devVec;
	unsigned int size;
	std::string name;

	vecWrapper(unsigned int size_, std::string name_) : size(size_), name(name_), hostVec(size_), devVec(size_)
	{
		//size = size_;
		//name = name_;
	}

	vecWrapper(unsigned int size_, std::string name_, double* arr) : size(size_), name(name_), hostVec(size_), devVec(size_)
	{
//		size = size_;
//		name = name_;
		copyDataTo(arr, size_);
	}

	vecWrapper(unsigned int size_, std::string name_, double fillVal) : size(size_), name(name_), hostVec(size_), devVec(size_)
	{
//		size = size_;
//		name = name_;
		fill(fillVal);
	}

	vecWrapper(viennacl::vector<double> vec_, std::string name_) : name(name_), devVec(vec_), hostVec(vec_.size())
	{
		size = (unsigned int)vec_.size();
//		name = name_;
		readFromDevice();
	}

	void fill(double val_)
	{
		std::fill(hostVec.begin(), hostVec.end(), val_);
		writeToDevice();
	}

	void zeros()
	{
		fill(0.);
	}

	void readFromDevice()
	{
		viennacl::fast_copy(devVec.begin(), devVec.end(), hostVec.begin());
	}

	void writeToDevice()
	{
		viennacl::fast_copy(hostVec.begin(), hostVec.end(), devVec.begin());
	}

	void copyDataTo(double* arr, unsigned int len)
	{
		if (len == 0)
			len = size;
		else if (len > size)
			printf("trying to copy data from larger array %s array\n", name.c_str());

		for (unsigned int i = 0; i < len; i++)
		{
			hostVec[i] = arr[i];
		}

		writeToDevice();
	}


	//void compareToArray(double* arr, const unsigned int start = 0, unsigned int stop = -1, unsigned int skip = 1)
	//{
	//	stop = (stop == -1) ? size : stop;
	//	double sum0 = 0, sum1 = 0, diffSum = 0;
	//	readFromDevice();
	//	if (stop > size)
	//	{
	//		printf("range to compare larger size of %s array\n", name.c_str());
	//		return;
	//	}
	//	for (unsigned int i = start; i < stop; i += skip)
	//	{
	//		sum0 += fabs(hostVec[i]);
	//		sum1 += fabs(arr[i]);
	//		diffSum += fabs(hostVec[i] - arr[i]);
	//	}
	//	printf("Comparison for %s from %d to %d: vcl = %g, lbm = %g, diff = %g\n", name.c_str(), start, stop, sum0, sum1, diffSum);

	//}

	void compareToArrayNoCopy(double* arr, const unsigned int start = 0, unsigned int stop = -1, unsigned int skip = 1)
	{
		stop = (stop == -1) ? size : stop;
		double sum0 = 0, sum1 = 0, diffSum = 0;

		if (stop > size)
		{
			printf("range to compare larger size of %s array\n", name.c_str());
			return;
		}
		for (unsigned int i = start; i < stop; i += skip)
		{
			sum0 += fabs(hostVec[i]);
			sum1 += fabs(arr[i]);
			diffSum += fabs(hostVec[i] - arr[i]);
		}
		printf("Comparison for %s from %d to %d: vcl = %g, lbm = %g, diff = %g\n", name.c_str(), start, stop, sum0, sum1, diffSum);

	}

	void compareToArrayNoCopy(cl_mem& arrMem, const unsigned int start = 0, unsigned int stop = -1, unsigned int skip = 1)
	{
		stop = (stop == -1) ? size : stop;
		Array1Dd arr("compArr");
		arr.zeros(stop);

		FINISH_QUEUES;
		arr.read_from_buffer_to_array(arrMem, TRUE, stop, IOQUEUE_REF);

		compareToArrayNoCopy(arr.get_array(), start, stop, skip);
	}

	template <typename T>
	void compareToArray(T arr, const unsigned int start = 0, unsigned int stop = -1, unsigned int skip = 1)
	{
		readFromDevice();
		compareToArrayNoCopy(arr, start, stop, skip);
	}

	//void compareToArray(double* arr, const unsigned int start = 0, unsigned int stop = -1, unsigned int skip = 1)
	//{
	//	readFromDevice();
	//	compareToArrayNoCopy(arr, start, stop, skip);
	//}


	viennacl::vector<double>& operator()()
	{
		return devVec;
	}

	double& operator()(const int n1)
	{
		return *(hostVec.begin() + n1);
	}

	double& operator[](const int n1)
	{
		return *(hostVec.begin() + n1);
	}

	void printHostVec(unsigned int colSize=1)
	{
		std::fstream stream;
		std::string nameout = name + ".txt";
		stream.open(nameout, std::ios_base::out);
		unsigned int rowSize;
		if (!stream.is_open())
			std::cout << "error opening " + nameout << "\n";
		
		rowSize = size / colSize;
		

		for (unsigned int i = 0; i < colSize; i++)
		{
			for (unsigned int j = 0; j < rowSize; j++)
			{
				stream << hostVec[j * colSize + i] << "\t";
			}
			stream << "\n";
		}
		stream.close();
	}

	void printDeviceVec(unsigned int colSize=1)
	{
		readFromDevice();
		printHostVec(colSize);
	}

};


void createVCLMatrix(viennacl::compressed_matrix<double> &A, int colSize_, int nnz_, unsigned int *IA, unsigned int *JA, double *AA)
{
	A.set(IA, JA, AA, colSize_, colSize_, nnz_);
}


void GMRESSolver::solve()
{
	// get rhs updated with initial guess
	updateRHS(calcQue);

	fillArraysWithInitialValues();
	
	unsigned int buffer_size_per_vector = 128;
	unsigned int num_buffer_chunks = 3;

	norm_rhs = rho_0(0);


	double rho = 1.0;
	iterCount = 0;

	for (unsigned int restart_count = 0; restart_count < maxRestarts; ++restart_count)
	{
		if (restart_count > 0)
		{
			// compute new residual without introducing a temporary for A*x:
			updateResidual(calcQue);			
			rho_0(0) = vclReduceMethod();
		}


		callScaleRes();


		rho = 1.0;

		if (checkConvergence())
			break;
	
		cl_uint k;
		for (k = 0; k < krylovDimSize; ++k)
		{
			if (k == 0)
			{
				gmresProdK0(calcQue);
				cl_uint zer = 0;

				normalizeVk.setArgument(1, &zer);
				normalizeVk.setArgument(4, &zer);
				normalizeVk.setArgument(8, &zer);
				normalizeVk(calcQue);
			}
			else
			{
				///// GMRES Product Kernel
				// My Implementation
				cl_uint vkm1Start = (k - 1) * colSize;
				cl_uint vkStart = vkm1Start + colSize;
				gmresProdK.setArgument(6, &vkm1Start);
				gmresProdK.setArgument(8, &vkStart);
				gmresProdK(calcQue);

				gramSchmidtStage1.setArgument(3, &k);
				gramSchmidtStage1(calcQue);


				gramSchmidtStage2.setArgument(3, &k);
				gramSchmidtStage2(calcQue);

				cl_uint offsetInR = k * (krylovDimSize + 1);
				cl_uint chunkOffset = k * bufSizePerVector;
				normalizeVk.setArgument(1, &vkStart);
				normalizeVk.setArgument(4, &offsetInR);
				normalizeVk.setArgument(8, &chunkOffset);
				normalizeVk(calcQue);
			}
		}

		callReduceVI(k);
		Rmat.read_from_buffer(calcQue);

		for (cl_uint i = 0; i < k; ++i)
		{
			if (std::fabs(Rmat[i + i * k]) < (relTol * Rmat[0]))
			{
				k = i;
				break;
			}
		}
		cl_uint full_krylov_dim = k;
		
		// Compute error estimator:
		for (cl_uint i = 0; i < k; ++i)
		{
			iterCount += 1; //increase iteration counter

			// check for accumulation of round-off errors for poorly conditioned systems
			if (values_xi_k[i] >= rho || values_xi_k[i] <= -rho)
			{
				k = i;
				break;  // restrict Krylov space at this point. No gain from using additional basis vectors, since orthogonality is lost.
			}

			// update error estimator
			rho *= std::sin(std::acos(values_xi_k[i] / rho));
		}


		// Solve minimization problem:
		for (int i2 = static_cast<int>(k) - 1; i2 > -1; --i2)
		{
			for (cl_uint j = i2 + 1; j < k; ++j)
				values_xi_k[i2] -= Rmat[i2 + j * full_krylov_dim] * values_xi_k[j];

			values_xi_k[i2] /= Rmat[i2 + i2 * full_krylov_dim];
		}

		for (cl_uint i = 0; i < k; ++i)
			values_xi_k[i] *= rho_0[0];


		values_xi_k.copy_to_buffer(calcQue);

		updateResult.setArgument(6, &k);
		updateResult(calcQue);

		FINISH_QUEUES;

		prevErr = fabs(rho * rho_0[0] / norm_rhs);

	}
	addIniGuess(calcQue);

}



//void GMRESSolver::solve()
//{
//	typedef double ScalarType;
//	viennacl::compressed_matrix<double> vclA;
//	Amat.read_from_buffer();
//	createVCLMatrix(vclA, colSize, Inds->nnz(), Inds->IA.get_array(), Inds->JA.get_array(), Amat.get_array());
//	//Amat.FreeHost();
//	//Amat.FreeDevice();
//	//Amat.setPreAllocatedBuffer(vclA.handle().opencl_handle().get());
//
//	//Inds->JA.FreeHost();
//	//Inds->JA.FreeDevice();
//	//Inds->JA.setPreAllocatedBuffer(vclA.handle2().opencl_handle().get());
//
//	//Inds->JA.FreeHost();
//	//Inds->JA.FreeDevice();
//	//Inds->JA.setPreAllocatedBuffer(vclA.handle1().opencl_handle().get());
//
//	//resetAmatrixArguments();
//	bVec.read_from_buffer();
//	xVec->read_from_buffer();
//	vecWrapper vclB(colSize, "vclB", bVec.get_array()), vclX(colSize, "vclX", xVec->get_array());
//
//	vecWrapper vclRHS(colSize, "vclRHS");
//	vclRHS.devVec = viennacl::linalg::prod(vclA, vclX.devVec);
//	vclRHS.devVec = vclB.devVec - vclRHS.devVec;
//
//	// get rhs updated with initial guess
//	updateRHS(calcQue);
//
//	fillArraysWithInitialValues();
//	vclRHS.compareToArray(resVec, 0, colSize);
//
//	vecWrapper vclResidual(vclRHS.devVec, "vclResidual");
//	vecWrapper vclResult(vclRHS.size, "vclResult", 0.0);
//
//	vecWrapper vclKrylovBasis(colSize * krylovDimSize, "vclKrylovBasis");
//	vecWrapper vclR(krylovDimSize * krylovDimSize, "vclR");
//
//	unsigned int buffer_size_per_vector = 128;
//	unsigned int num_buffer_chunks = 3;
//
//	vecWrapper vclInnerProd(num_buffer_chunks * buffer_size_per_vector, "vclInnerProd", 0.0);
//	vecWrapper vclRDotVk(buffer_size_per_vector * krylovDimSize, "vclRDotVk", 0.0); // holds result of first reduction stage for <r, v_k> on device
//	vecWrapper vclViInVk(buffer_size_per_vector * krylovDimSize, "vclViInVk", 0.0); 
//	vecWrapper vclValuesXiK(krylovDimSize, "vclValuesXiK", 0.0); // holds values \xi_k = <r, v_k> on device
//	vecWrapper vclEtaK(krylovDimSize, "vclEtaK", 0.0);
//	vecWrapper vclUpdateCoeff(krylovDimSize, "vclUpdateCoeff", 0.0);
//
//
//
//	ScalarType vclNormRHS = viennacl::linalg::norm_2(vclResidual.devVec);
//	ScalarType vclRho0 = vclNormRHS;
//	ScalarType vclRho = ScalarType(1);
//
//	
//	//double norm_rhs2 = vclReduceMethod();
//	norm_rhs = rho_0(0);
//
//	
//	double rho = 1.0;
//	iterCount = 0;
//	int vclEndFlag = 0;
//
//	for (unsigned int restart_count = 0; restart_count < maxRestarts; ++restart_count)
//	{
//		if (restart_count > 0)
//		{
//			// compute new residual without introducing a temporary for A*x:
//			updateResidual(calcQue);
//
//			vclResidual.devVec = viennacl::linalg::prod(vclA, vclResult.devVec);
//			vclResidual() = vclRHS() - vclResidual();
//
//			vclRho0 = viennacl::linalg::norm_2(vclResidual());
//
//
//			rho_0(0) = vclReduceMethod();
//			std::cout << "restarting GMRES with rho0 = " << vclRho0 << " (vcl) and " << rho_0(0) << " (lbm)\n";
//		}
//
//		//rho_0[0] = vclNormRHS;
//		callScaleRes();
//		vclResidual() /= vclRho0;
//
//		vclRho = ScalarType(1.0);
//		rho = 1.0;
//
//		if (vclRho0 / vclNormRHS < relTol || vclRho0 < absTol)
//		{
//			std::cout << "vcl would break here\n";
//			vclX.compareToArray(iniXVal);
//			FINISH_QUEUES;
//			xVec->read_from_buffer();
//			vclResult.compareToArray(xVec->get_array());
//			if(vclEndFlag == 0)
//				vclX.devVec = vclResult.devVec + vclX.devVec;
//			vclEndFlag++;
//		}
//
//		if (checkConvergence())
//		{
//			break;
//		}
//		cl_uint k;
//		for (k = 0; k < krylovDimSize; ++k)
//		{
//			if (k == 0)
//			{
//				gmresProdK0(calcQue);
//
//				viennacl::vector_range<viennacl::vector<ScalarType> > v0(vclKrylovBasis(), viennacl::range(0, vclRHS.size));
//				viennacl::linalg::pipelined_gmres_prod(vclA, vclResidual(), v0, vclInnerProd());
//
//				//vclKrylovBasis.compareToArray(krylov_basis, 0, (k + 1) * colSize);
//				//vclInnerProd.compareToArray(device_inner_prod_buffer);
//
//				cl_uint zer = 0;
//
//				normalizeVk.setArgument(1, &zer);
//				normalizeVk.setArgument(4, &zer);
//				normalizeVk.setArgument(8, &zer);
//				normalizeVk(calcQue);
//
//				// Normalize v_1 and compute first reduction stage for <r, v_0> in device_r_dot_vk_buffer:
//				viennacl::linalg::pipelined_gmres_normalize_vk(v0, vclResidual(),
//					vclR(), k * krylovDimSize + k,
//					vclInnerProd(), vclRDotVk(),
//					buffer_size_per_vector, k * buffer_size_per_vector);
//				//vclR.compareToArray(Rmat.get_array(), 0, 1+ krylovDimSize *(k+1), krylovDimSize);
//				//vclInnerProd.compareToArray(device_inner_prod_buffer);
//				//vclKrylovBasis.compareToArray(krylov_basis, 0, (k+1)*colSize);
//				//vclRDotVk.compareToArray(r_dot_vk_buffer, 0, (k + 1) * buffer_size_per_vector);
//
//			}
//			else
//			{
//				///// GMRES Product Kernel
//				// My Implementation
//				cl_uint vkm1Start = (k - 1) * colSize;
//				cl_uint vkStart = vkm1Start + colSize;
//				gmresProdK.setArgument(6, &vkm1Start);
//				gmresProdK.setArgument(8, &vkStart);
//				gmresProdK(calcQue);
//
//				// VCL Implementation
//				viennacl::vector_range<viennacl::vector<ScalarType> > vk(vclKrylovBasis(), viennacl::range(k * vclRHS.size, k * vclRHS.size + vclRHS.size));
//				viennacl::vector_range<viennacl::vector<ScalarType> > vk_minus_1(vclKrylovBasis(), viennacl::range(k * vclRHS.size - vclRHS.size, k*vclRHS.size));
//				viennacl::linalg::pipelined_gmres_prod(vclA, vk_minus_1, vk, vclInnerProd());
//				
//				//vclKrylovBasis.compareToArray(krylov_basis, 0, (k + 1)* colSize);
//				//vclInnerProd.compareToArray(device_inner_prod_buffer);
//
//				gramSchmidtStage1.setArgument(3, &k);
//				gramSchmidtStage1(calcQue);
//
//				viennacl::linalg::pipelined_gmres_gram_schmidt_stage1(vclKrylovBasis(), vclRHS.size, vclRHS.size, k, vclViInVk(), buffer_size_per_vector);
//				//vclViInVk.compareToArray(vi_in_vk_buffer);
//
//				gramSchmidtStage2.setArgument(3, &k);
//				gramSchmidtStage2(calcQue);
//
//				viennacl::linalg::pipelined_gmres_gram_schmidt_stage2(vclKrylovBasis(), vclRHS.size, vclRHS.size, k,
//					vclViInVk(), vclR(), krylovDimSize, vclInnerProd(), buffer_size_per_vector);
//
//				//vclR.compareToArray(Rmat.get_array(), 0, 1 + krylovDimSize * (k + 1), krylovDimSize);
//				//vclInnerProd.compareToArray(device_inner_prod_buffer);
//				//vclKrylovBasis.compareToArray(krylov_basis, 0, (k + 1)* colSize);
//				//vclRDotVk.compareToArray(r_dot_vk_buffer, 0, (k + 1)* buffer_size_per_vector);
//				//vclViInVk.compareToArray(vi_in_vk_buffer);
//
//				cl_uint offsetInR = k * (krylovDimSize + 1);
//				cl_uint chunkOffset = k * bufSizePerVector;
//				normalizeVk.setArgument(1, &vkStart);
//				normalizeVk.setArgument(4, &offsetInR);
//				normalizeVk.setArgument(8, &chunkOffset);
//				normalizeVk(calcQue);
//
//
//				viennacl::linalg::pipelined_gmres_normalize_vk(vk, vclResidual(),
//					vclR(), k* krylovDimSize + k,
//					vclInnerProd(), vclRDotVk(),
//					buffer_size_per_vector, k* buffer_size_per_vector);
//				
//				//vclR.compareToArray(Rmat.get_array(), 0, 1 + krylovDimSize * (k + 1), krylovDimSize);
//				//vclInnerProd.compareToArray(device_inner_prod_buffer);
//				//vclKrylovBasis.compareToArray(krylov_basis, 0, (k + 1)* colSize);
//				//vclRDotVk.compareToArray(r_dot_vk_buffer, 0, (k + 1)* buffer_size_per_vector);
//
//			}
//
//			//vclR.printDeviceVec(krylovDimSize);
//			//Rmat.save_txt_from_device_as_2D(krylovDimSize, krylovDimSize, krylovDimSize, "lbmR");
//			//#ifdef _DEBUG
//			//			debugSaveBuffer(innerProd, k);
//			//#endif
//		}
//
//		unsigned int vclk = k;
//		
//		callReduceVI(k);
//		Rmat.read_from_buffer(calcQue);
//
//		vclRDotVk.readFromDevice();
//		for (std::size_t i = 0; i < vclk; ++i)
//		{
//			vclValuesXiK[(unsigned int)i] = ScalarType(0);
//			for (std::size_t j = 0; j < buffer_size_per_vector; ++j)
//				vclValuesXiK[(unsigned int)i] += vclRDotVk[(unsigned int)(i * buffer_size_per_vector + j)];
//		}
//
//		
//
//		vclR.readFromDevice();
//		unsigned int vcl_full_krylov_dim = vclk; //needed for proper access to R
//		for (unsigned int i = 0; i < vclk; ++i)
//		{
//			if (std::fabs(vclR[i + i * vclk]) < relTol * vclR[0])
//			{
//				vclk = i;
//				break;
//			}
//		}
//
//
//		cl_uint vclfull_krylov_dim = vclk;
//
//		for (cl_uint i = 0; i < k; ++i)
//		{
//			if (std::fabs(Rmat[i + i * k]) < (relTol * Rmat[0]))
//			{
//				k = i;
//				break;
//			}
//		}
//		cl_uint full_krylov_dim = k;
//		// Compute error estimator:
//		for (cl_uint i = 0; i < k; ++i)
//		{
//			iterCount += 1; //increase iteration counter
//
//			// check for accumulation of round-off errors for poorly conditioned systems
//			if (values_xi_k[i] >= rho || values_xi_k[i] <= -rho)
//			{
//				k = i;
//				break;  // restrict Krylov space at this point. No gain from using additional basis vectors, since orthogonality is lost.
//			}
//
//			// update error estimator
//			rho *= std::sin(std::acos(values_xi_k[i] / rho));
//		}
//
//		// Compute error estimator:
//		for (cl_uint i = 0; i < vclk; ++i)
//		{
//
//			// check for accumulation of round-off errors for poorly conditioned systems
//			if (vclValuesXiK[i] >= vclRho || vclValuesXiK[i] <= -vclRho)
//			{
//				vclk = i;
//				break;  // restrict Krylov space at this point. No gain from using additional basis vectors, since orthogonality is lost.
//			}
//
//			// update error estimator
//			vclRho *= std::sin(std::acos(vclValuesXiK[i] / vclRho));
//		}
//
//		// Solve minimization problem:
//		//host_values_eta_k_buffer = host_values_xi_k;
//
//		for (int i2 = static_cast<int>(k) - 1; i2 > -1; --i2)
//		{
//			for (cl_uint j = i2 + 1; j < k; ++j)
//				values_xi_k[i2] -= Rmat[i2 + j * full_krylov_dim] * values_xi_k[j];
//
//			values_xi_k[i2] /= Rmat[i2 + i2 * full_krylov_dim];
//		}
//
//		for (unsigned int i2 = 0; i2 < krylovDimSize; i2++)
//			vclEtaK[i2] = vclValuesXiK[i2];
//		for (int i2 = static_cast<int>(vclk) - 1; i2 > -1; --i2)
//		{
//			unsigned int i = static_cast<unsigned int>(i2);
//			for (unsigned int j = static_cast<unsigned int>(i) + 1; j < vclk; ++j)
//				vclEtaK[i] -= vclR[i + j * vclfull_krylov_dim] * vclEtaK[j];
//
//			vclEtaK[i] /= vclR[i + i * vclfull_krylov_dim];
//		}
//
//		vclEtaK.compareToArrayNoCopy(values_xi_k.get_array(), 0, vclk);
//
//		for (cl_uint i = 0; i < k; ++i)
//			values_xi_k[i] *= rho_0[0];
//
//		for (cl_uint i = 0; i < vclk; ++i)
//			vclUpdateCoeff[i] = vclRho0 * vclEtaK[i];
//
//		vclUpdateCoeff.writeToDevice();
//		
//		values_xi_k.copy_to_buffer(calcQue);
//		
//		updateResult.setArgument(6, &k);
//		updateResult(calcQue);
//
//		viennacl::linalg::pipelined_gmres_update_result(vclResult(), vclResidual(),
//			vclKrylovBasis(), vclRHS().size(), vclRHS().internal_size(),
//			vclUpdateCoeff(), vclk);
//
//		FINISH_QUEUES;
//
//		prevErr = fabs(rho * rho_0[0] / norm_rhs);
//#ifdef _DEBUG
//		debugPrintResidual(prevErr);
//#endif
//	}
//
//	addIniGuess(calcQue);
//	if (vclEndFlag == 0)
//	{
//		std::cout << "\n\nLBM has converged, vcl has not. Final Comparison below\n";
//	}
//	else if (vclEndFlag == 1)
//	{
//		std::cout << "\n\nFinished at same time as vcl. Final Comparison below\n";
//	}
//	else
//	{
//		std::cout << "\n\nvcl finished " << vclEndFlag - 1 << " restarts ago. Final Comparison below\n";
//	}
//	xVec->read_from_buffer(calcQue);
//	vclX.compareToArray(xVec->get_array());
//}


double GMRESSolver::vclReduceMethod()
{
	vclReduce(calcQue);
	FINISH_QUEUES;
	norm2TempVals.read_from_buffer();
	double retval = 0.;
	for (int i = 0; i < norm2TempVals.getFullSize(); i++)
		retval += norm2TempVals(i);
	return sqrt(retval);
}

//void GMRESSolver::solve()
//{
//#ifdef _DEBUG
//	//saveVec_As_MM(xVec->get_array(), "xVec_MM.mtx", true);
//	//saveVec_As_MM(bVec.get_array(), "bVec_MM.mtx", true);
//	//saveCSR_As_MM("Amatrix.mtx", true);
//	iniMonitorResidual();
//#endif
//
//
//	// get rhs updated with initial guess
//	updateRHS(calcQue);
//
//	fillArraysWithInitialValues();
//	double* vclData;
//
//	vclData = (double*)malloc(sizeof(double) * colSize);
//
//	std::fstream stream;
//
//	stream.open("C:\\Users\\Grant\\Documents\\Research_local\\kOmega\\vcl_RHS.bin", std::ios_base::binary | std::ios_base::in);
//	stream.read((char*)vclData, colSize * sizeof(double));
//	stream.close();
//	Array1Dd diffArr("diffArr");
//	diffArr.zeros(colSize);
//	double cumDiff = 0.;
//
//	DebugOut.read_from_buffer_to_array(resVec, TRUE, colSize, calcQue);
//	double sumVCL = 0., sumLBM = 0.;
//	for (int iii = 0; iii < colSize; iii++)
//	{
//		diffArr(iii) = vclData[iii] - DebugOut(iii);
//		cumDiff += fabs(diffArr(iii));
//		sumVCL += fabs(vclData[iii]);
//		sumLBM += fabs(DebugOut(iii));
//	}
//	diffArr.savetxt();
//
//
//	//#ifdef _DEBUG
//	//	bVec.save_txt_from_device("", calcQue);
//	//#endif
//
//
//	norm_rhs = rho_0(0);
//	double rho = 1.0;
//	iterCount = 0;
//	for (unsigned int restart_count = 0; restart_count < maxRestarts; ++restart_count)
//	{
//		if (restart_count > 0)
//		{
//			// compute new residual without introducing a temporary for A*x:
//			updateResidual(calcQue);
//			calcNorm2Residual();
//		}
//
//		callScaleRes();
//
//		rho = 1.0;
//		if (checkConvergence())
//			break;
//		cl_uint k;
//		for (k = 0; k < krylovDimSize; ++k)
//		{
//			if (k == 0)
//			{
//				gmresProdK0(calcQue);
//				//#ifdef _DEBUG
//				//				debugSaveBuffer(innerProd);
//				//				debugSaveBuffer(krylovBasis, k);
//				//#endif
//
//
//
//				cl_uint zer = 0;
//
//				normalizeVk.setArgument(1, &zer);
//				normalizeVk.setArgument(4, &zer);
//				normalizeVk.setArgument(8, &zer);
//				normalizeVk(calcQue);
//
//				//#ifdef _DEBUG
//				//				debugSaveBuffer(rDotVK,k);
//				//				debugSaveBuffer(krylovBasis,k);
//				//				Rmat.save_txt_from_device("", calcQue);
//				//#endif
//			}
//			else
//			{
//				cl_uint vkm1Start = (k - 1) * colSize;
//				cl_uint vkStart = vkm1Start + colSize;
//				gmresProdK.setArgument(6, &vkm1Start);
//				gmresProdK.setArgument(8, &vkStart);
//				gmresProdK(calcQue);
//
//				//#ifdef _DEBUG
//				//				debugSaveBuffer(innerProd);
//				//				debugSaveBuffer(krylovBasis,k);
//				//#endif
//
//				gramSchmidtStage1.setArgument(3, &k);
//				gramSchmidtStage1(calcQue);
//
//				//#ifdef _DEBUG
//				//				debugSaveBuffer(innerProd);
//				//				debugSaveBuffer(krylovBasis,k);
//				//#endif
//
//				gramSchmidtStage2.setArgument(3, &k);
//				gramSchmidtStage2(calcQue);
//
//				//#ifdef _DEBUG
//				//				debugSaveBuffer(innerProd);
//				//				debugSaveBuffer(krylovBasis,k);
//				//#endif
//
//				cl_uint offsetInR = k * (krylovDimSize + 1);
//				cl_uint chunkOffset = k * bufSizePerVector;
//				normalizeVk.setArgument(1, &vkStart);
//				normalizeVk.setArgument(4, &offsetInR);
//				normalizeVk.setArgument(8, &chunkOffset);
//				normalizeVk(calcQue);
//
//				//#ifdef _DEBUG
//				//				debugSaveBuffer(rDotVK, k);
//				//				debugSaveBuffer(krylovBasis, k);
//				//				Rmat.save_txt_from_device("", calcQue);
//				//#endif
//
//			}
//			//#ifdef _DEBUG
//			//			debugSaveBuffer(innerProd, k);
//			//#endif
//		}
//		callReduceVI(k);
//		Rmat.read_from_buffer(calcQue);
//		//#ifdef _DEBUG
//		//		//debugSaveBuffer(krylovBasis, k-1);
//		//		//Rmat.save2file_as_2D(krylovDimSize, krylovDimSize, krylovDimSize, "vlbm_Rmatrix");
//		//		debugSaveBuffer(rDotVK, k-1);
//		//		//values_xi_k.save2file("vlbm_Xik");
//		//#endif
//
//
//		cl_uint full_krylov_dim = k;
//		for (cl_uint i = 0; i < k; ++i)
//		{
//			if (std::fabs(Rmat[i + i * k]) < (relTol * Rmat[0]))
//			{
//				k = i;
//				break;
//			}
//		}
//
//		// Compute error estimator:
//		for (cl_uint i = 0; i < k; ++i)
//		{
//			iterCount += 1; //increase iteration counter
//
//			// check for accumulation of round-off errors for poorly conditioned systems
//			if (values_xi_k[i] >= rho || values_xi_k[i] <= -rho)
//			{
//				k = i;
//				break;  // restrict Krylov space at this point. No gain from using additional basis vectors, since orthogonality is lost.
//			}
//
//			// update error estimator
//			rho *= std::sin(std::acos(values_xi_k[i] / rho));
//		}
//
//		// Solve minimization problem:
//		//host_values_eta_k_buffer = host_values_xi_k;
//
//		for (int i2 = static_cast<int>(k) - 1; i2 > -1; --i2)
//		{
//			for (cl_uint j = i2 + 1; j < k; ++j)
//				values_xi_k[i2] -= Rmat[i2 + j * full_krylov_dim] * values_xi_k[j];
//
//			values_xi_k[i2] /= Rmat[i2 + i2 * full_krylov_dim];
//		}
//
//		for (cl_uint i = 0; i < k; ++i)
//			values_xi_k[i] *= rho_0[0];
//		values_xi_k.copy_to_buffer(calcQue);
//		values_xi_k.save2file();
//		updateResult.setArgument(6, &k);
//		updateResult(calcQue);
//
//		prevErr = fabs(rho * rho_0[0] / norm_rhs);
//#ifdef _DEBUG
//		debugPrintResidual(prevErr);
//#endif
//	}
//
//}



void GMRESSolver::resetAmatrixArguments()
{
	curMonitorResid.setArgument(0, Inds->IA.get_buf_add());
	curMonitorResid.setArgument(1, Inds->JA.get_buf_add());
	curMonitorResid.setArgument(3, Amat.get_buf_add());

	iniMonitorResid.setArgument(0, Inds->IA.get_buf_add());
	iniMonitorResid.setArgument(1, Inds->JA.get_buf_add());
	iniMonitorResid.setArgument(3, Amat.get_buf_add());

	updateRHS.setArgument(0, Inds->IA.get_buf_add());
	updateRHS.setArgument(1, Inds->JA.get_buf_add());
	updateRHS.setArgument(3, Amat.get_buf_add());

	updateResidual.setArgument(0, Inds->IA.get_buf_add());
	updateResidual.setArgument(1, Inds->JA.get_buf_add());
	updateResidual.setArgument(3, Amat.get_buf_add());

	gmresProdK0.setArgument(0, Inds->IA.get_buf_add());
	gmresProdK0.setArgument(1, Inds->JA.get_buf_add());
	gmresProdK0.setArgument(3, Amat.get_buf_add());

	gmresProdK.setArgument(0, Inds->IA.get_buf_add());
	gmresProdK.setArgument(1, Inds->JA.get_buf_add());
	gmresProdK.setArgument(3, Amat.get_buf_add());
}


void GMRESSolver::createKernels()
{
	norm2TempVals.zeros(128);
	norm2TempVals.allocate_buffer_w_copy();
	vclReduce.setSizes(128*128, 128);
	vclReduce.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "normVCL");
	vclReduce.setArgument(0, &resVec);
	vclReduce.setArgument(1, &colSize);
	vclReduce.setLocalMem(2, 128 * sizeof(double));
	vclReduce.setArgument(3, norm2TempVals.get_buf_add());

	addIniGuess.setSizes(colSize, WORKGROUPSIZE_GMRESPROD);
	addIniGuess.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "addInitialGuess");
	addIniGuess.setArgument(0, &colSize);
	addIniGuess.setArgument(1, &iniXVal);
	addIniGuess.setArgument(2, xVec->get_buf_add());





	cl_uint4 layoutx;
	int localWorkSize, globalWorkSize;
	int iBufOut;
#ifdef _DEBUG
	rho0Monitor.zeros(1);
	rho0Monitor.allocate_buffer_w_copy();
	layoutx = { {0, 1, colSize, fullSize} };
	localWorkSize = WORKGROUPSIZE_VCLSPARSE;
	globalWorkSize = localWorkSize * numRowBlocks;
	curMonitorResid.setSizes(globalWorkSize, localWorkSize);
	curMonitorResid.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "updateResidual");
	curMonitorResid.setArgument(0, Inds->IA.get_buf_add());
	curMonitorResid.setArgument(1, Inds->JA.get_buf_add());
	curMonitorResid.setArgument(2, rowBlocks.get_buf_add());
	curMonitorResid.setArgument(3, Amat.get_buf_add());
	curMonitorResid.setArgument(4, &numRowBlocks);
	curMonitorResid.setArgument(5, &xMonitor);
	curMonitorResid.setArgument(6, &layoutx);
	curMonitorResid.setArgument(7, &residMonitor);
	curMonitorResid.setArgument(8, &layoutx);
	curMonitorResid.setArgument(9, bVec_copy.get_buf_add());

	iniMonitorResid.setSizes(globalWorkSize, localWorkSize);
	iniMonitorResid.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "updateResidual");
	iniMonitorResid.setArgument(0, Inds->IA.get_buf_add());
	iniMonitorResid.setArgument(1, Inds->JA.get_buf_add());
	iniMonitorResid.setArgument(2, rowBlocks.get_buf_add());
	iniMonitorResid.setArgument(3, Amat.get_buf_add());
	iniMonitorResid.setArgument(4, &numRowBlocks);
	iniMonitorResid.setArgument(5, xVec->get_buf_add());
	iniMonitorResid.setArgument(6, &layoutx);
	iniMonitorResid.setArgument(7, &residMonitor);
	iniMonitorResid.setArgument(8, &layoutx);
	iniMonitorResid.setArgument(9, bVec.get_buf_add());

	sumKer.setSizes(fullSize, WORKGROUPSIZE_VCLSPARSE);
	sumKer.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "sumVectorKernel");
	sumKer.setArgument(0, &colSize);
	sumKer.setArgument(1, &xMonitor);
	sumKer.setArgument(2, xVec->get_buf_add());
	sumKer.setArgument(3, xVec_copy.get_buf_add());

	iBufOut = 0;
	redMonitor.ini(numRedKer);
	redMonitor(0).setSizes(globalWorkSizeRed[0], localWorkSizeRed[0]);
	redMonitor(0).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "baseNorm2");
	redMonitor(0).setArgument(0, &residMonitor);
	redMonitor(0).setArgument(1, &reduceBufSet1[iBufOut]);
	redMonitor(0).setLocalMem(2, WORKGROUPSIZE_RED * sizeof(double));

	for (int i = 1; i < numRedKer - 1; i++)
	{
		int lworksizeTmp = localWorkSizeRed[i];
		redMonitor(i).setSizes(globalWorkSizeRed[i], lworksizeTmp);
		if (lworksizeTmp == WORKGROUPSIZE_RED)
			redMonitor(i).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "intermediateNorm2Full");
		else
			redMonitor(i).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "intermediateNorm2");
		redMonitor(i).setArgument(0, &reduceBufSet1[iBufOut]);
		iBufOut ^= 1;
		redMonitor(i).setArgument(1, &reduceBufSet1[iBufOut]);
		redMonitor(i).setLocalMem(2, lworksizeTmp * sizeof(double));
	}
	redMonitor(numRedKer - 1).setSizes(globalWorkSizeRed[numRedKer - 1], localWorkSizeRed[numRedKer - 1]);
	redMonitor(numRedKer - 1).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "FinalNorm2");
	redMonitor(numRedKer - 1).setArgument(0, &reduceBufSet1[iBufOut]);
	redMonitor(numRedKer - 1).setArgument(1, rho0Monitor.get_buf_add());
	redMonitor(numRedKer - 1).setLocalMem(2, localWorkSizeRed[numRedKer - 1] * sizeof(double));

#endif


	layoutx = { {0, 1, colSize, fullSize} };
	localWorkSize = WORKGROUPSIZE_VCLSPARSE;
	globalWorkSize = localWorkSize * numRowBlocks;
	updateRHS.setSizes(globalWorkSize, localWorkSize);
	updateRHS.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "updateRHS");
	updateRHS.setArgument(0, Inds->IA.get_buf_add());
	updateRHS.setArgument(1, Inds->JA.get_buf_add());
	updateRHS.setArgument(2, rowBlocks.get_buf_add());
	updateRHS.setArgument(3, Amat.get_buf_add());
	updateRHS.setArgument(4, &numRowBlocks);
	updateRHS.setArgument(5, xVec->get_buf_add());
	updateRHS.setArgument(6, &layoutx);
	updateRHS.setArgument(7, bVec.get_buf_add());
	updateRHS.setArgument(8, &layoutx);

	updateResidual.setSizes(globalWorkSize, localWorkSize);
	updateResidual.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "updateResidual");
	updateResidual.setArgument(0, Inds->IA.get_buf_add());
	updateResidual.setArgument(1, Inds->JA.get_buf_add());
	updateResidual.setArgument(2, rowBlocks.get_buf_add());
	updateResidual.setArgument(3, Amat.get_buf_add());
	updateResidual.setArgument(4, &numRowBlocks);
	updateResidual.setArgument(5, xVec->get_buf_add());
	updateResidual.setArgument(6, &layoutx);
	updateResidual.setArgument(7, &resVec);
	updateResidual.setArgument(8, &layoutx);
	updateResidual.setArgument(9, bVec.get_buf_add());

	scaleRes.setSizes(fullSize, WORKGROUPSIZE_VCLSPARSE);
	scaleRes.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "NormalizeResidual");
	scaleRes.setArgument(0, &colSize);
	scaleRes.setArgument(1, &resVec);
	
	cl_uint zer = 0;
	gmresProdK0.setSizes(128 * WORKGROUPSIZE_GMRESPROD, WORKGROUPSIZE_GMRESPROD);
	gmresProdK0.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "gmres_csr_prod");
	gmresProdK0.setArgument(0, Inds->IA.get_buf_add());
	gmresProdK0.setArgument(1, Inds->JA.get_buf_add());
	gmresProdK0.setArgument(2, rowBlocks.get_buf_add());
	gmresProdK0.setArgument(3, Amat.get_buf_add());
	gmresProdK0.setArgument(4, &numRowBlocks);
	gmresProdK0.setArgument(5, &resVec);
	gmresProdK0.setArgument(6, &zer);
	gmresProdK0.setArgument(7, &krylov_basis);
	gmresProdK0.setArgument(8, &zer);
	gmresProdK0.setArgument(9, &colSize);
	gmresProdK0.setArgument(10, &device_inner_prod_buffer);
	gmresProdK0.setArgument(11, &bufSizePerVector);
	gmresProdK0.setLocalMem(12, WORKGROUPSIZE_GMRESPROD * sizeof(double));
	gmresProdK0.setLocalMem(13, WORKGROUPSIZE_GMRESPROD * sizeof(double));
	gmresProdK0.setLocalMem(14, 1024 * sizeof(double));
	


	gmresProdK.setSizes(128 * WORKGROUPSIZE_GMRESPROD, WORKGROUPSIZE_GMRESPROD);
	gmresProdK.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "gmres_csr_prod");
	gmresProdK.setArgument(0, Inds->IA.get_buf_add());
	gmresProdK.setArgument(1, Inds->JA.get_buf_add());
	gmresProdK.setArgument(2, rowBlocks.get_buf_add());
	gmresProdK.setArgument(3, Amat.get_buf_add());
	gmresProdK.setArgument(4, &numRowBlocks);
	gmresProdK.setArgument(5, &krylov_basis);
	//gmresProdK.setArgument(6, &zer);//must be set before call
	gmresProdK.setArgument(7, &krylov_basis);
	//gmresProdK.setArgument(8, &zer);// must be set before call
	gmresProdK.setArgument(9, &colSize);
	gmresProdK.setArgument(10, &device_inner_prod_buffer);
	gmresProdK.setArgument(11, &bufSizePerVector);
	gmresProdK.setLocalMem(12, WORKGROUPSIZE_GMRESPROD * sizeof(double));
	gmresProdK.setLocalMem(13, WORKGROUPSIZE_GMRESPROD * sizeof(double));
	gmresProdK.setLocalMem(14, 1024 * sizeof(double));



	normalizeVk.setSizes(128 * WORKGROUPSIZE_GMRESPROD, WORKGROUPSIZE_GMRESPROD);
	cl_uint chunk_size = 128;
	normalizeVk.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "gmres_normalize_vk");
	normalizeVk.setArgument(0, &krylov_basis);
	//normalizeVk.setArgument(1, ); // must be set before call
	normalizeVk.setArgument(2, &resVec);
	normalizeVk.setArgument(3, Rmat.get_buf_add());
	//normalizeVk.setArgument(4, );//offset in Rmatrix
	normalizeVk.setArgument(5, &device_inner_prod_buffer);
	normalizeVk.setArgument(6, &chunk_size);
	normalizeVk.setArgument(7, r_dot_vk_buffer.get_buf_add());
	//normalizeVk.setArgument(8, &zer);// must be set before call
	normalizeVk.setArgument(9, &colSize);
	normalizeVk.setLocalMem(10, WORKGROUPSIZE_GMRESPROD * sizeof(double));
	

	gramSchmidtStage1.setSizes(128 * WORKGROUPSIZE_GMRESPROD, WORKGROUPSIZE_GMRESPROD);
	gramSchmidtStage1.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "gmres_gram_schmidt_1");
	gramSchmidtStage1.setArgument(0, &krylov_basis);
	gramSchmidtStage1.setArgument(1, &colSize);
	gramSchmidtStage1.setArgument(2, &colSize);
	//gramSchmidtStage1.setArgument(3, );//set before call
	gramSchmidtStage1.setArgument(4, &vi_in_vk_buffer);
	gramSchmidtStage1.setArgument(5, &chunk_size);

	gramSchmidtStage2.setSizes(128 * WORKGROUPSIZE_GMRESPROD, WORKGROUPSIZE_GMRESPROD);
	gramSchmidtStage2.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "gmres_gram_schmidt_2");
	gramSchmidtStage2.setArgument(0, &krylov_basis);
	gramSchmidtStage2.setArgument(1, &colSize);
	gramSchmidtStage2.setArgument(2, &colSize);
	//gramSchmidtStage2.setArgument(3, );//set before call
	gramSchmidtStage2.setArgument(4, &vi_in_vk_buffer);
	gramSchmidtStage2.setArgument(5, &chunk_size);
	gramSchmidtStage2.setArgument(6, Rmat.get_buf_add());
	gramSchmidtStage2.setArgument(7, &krylovDimSize);
	gramSchmidtStage2.setArgument(8, &device_inner_prod_buffer);
	gramSchmidtStage2.setLocalMem(9, 7* WORKGROUPSIZE_GMRESPROD*sizeof(double));



	//redRDotVk.setSizes(64* krylovDimSize, 64); // needs to be set before each call
	//redRDotVk.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "ReduceRVk");
	//redRDotVk.setArgument(0, &r_dot_vk_buffer);
	//redRDotVk.setArgument(1, values_xi_k.get_buf_add());
	//redRDotVk.setLocalMem(2, 64*krylovDimSize*sizeof(double));

	updateResult.setSizes(128 * WORKGROUPSIZE_GMRESPROD, WORKGROUPSIZE_GMRESPROD);
	updateResult.createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "gmres_update_result");
	updateResult.setArgument(0, xVec->get_buf_add());
	updateResult.setArgument(1, &resVec);
	updateResult.setArgument(2, &krylov_basis);
	updateResult.setArgument(3, &colSize);
	updateResult.setArgument(4, &colSize);
	updateResult.setArgument(5, values_xi_k.get_buf_add());
	//updateResult.setArgument(6, );// must be set before call



	iBufOut = 0;
	resNorm2WithNormalize.ini(numRedKer);
	resNorm2WithNormalize(0).setSizes(globalWorkSizeRed[0], localWorkSizeRed[0]);
	resNorm2WithNormalize(0).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "NormalizeResidualWithNorm2Base");
	resNorm2WithNormalize(0).setArgument(0, &resVec);
	//resNorm2WithNormalize(0).setArgument(1, &rho_0[0]);//need to set before calling
	resNorm2WithNormalize(0).setArgument(2, &reduceBufSet1[iBufOut]);
	resNorm2WithNormalize(0).setLocalMem(3, WORKGROUPSIZE_RED * sizeof(double));

	for (int i = 1; i < numRedKer - 1; i++)
	{
		int lworksizeTmp = localWorkSizeRed[i];
		resNorm2WithNormalize(i).setSizes(globalWorkSizeRed[i], lworksizeTmp);
		if(lworksizeTmp == WORKGROUPSIZE_RED)
			resNorm2WithNormalize(i).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "intermediateNorm2Full");
		else
			resNorm2WithNormalize(i).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "intermediateNorm2");
		resNorm2WithNormalize(i).setArgument(0, &reduceBufSet1[iBufOut]);
		iBufOut ^= 1;
		resNorm2WithNormalize(i).setArgument(1, &reduceBufSet1[iBufOut]);
		resNorm2WithNormalize(i).setLocalMem(2, lworksizeTmp * sizeof(double));
	}
	resNorm2WithNormalize(numRedKer - 1).setSizes(globalWorkSizeRed[numRedKer - 1], localWorkSizeRed[numRedKer - 1]);
	resNorm2WithNormalize(numRedKer - 1).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "FinalNorm2");
	resNorm2WithNormalize(numRedKer - 1).setArgument(0, &reduceBufSet1[iBufOut]);
	resNorm2WithNormalize(numRedKer - 1).setArgument(1, rho_0.get_buf_add());
	resNorm2WithNormalize(numRedKer - 1).setLocalMem(2, localWorkSizeRed[numRedKer-1] * sizeof(double));


	iBufOut = 0;
	resNorm2.ini(numRedKer);
	resNorm2(0).setSizes(globalWorkSizeRed[0], localWorkSizeRed[0]);
	resNorm2(0).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "baseNorm2");
	resNorm2(0).setArgument(0, &resVec);
	resNorm2(0).setArgument(1, &reduceBufSet1[iBufOut]);
	resNorm2(0).setLocalMem(2, WORKGROUPSIZE_RED * sizeof(double));

	for (int i = 1; i < numRedKer - 1; i++)
	{
		int lworksizeTmp = localWorkSizeRed[i];
		resNorm2(i).setSizes(globalWorkSizeRed[i], lworksizeTmp);
		if (lworksizeTmp == WORKGROUPSIZE_RED)
			resNorm2(i).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "intermediateNorm2Full");
		else
			resNorm2(i).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "intermediateNorm2");
		resNorm2(i).setArgument(0, &reduceBufSet1[iBufOut]);
		iBufOut ^= 1;
		resNorm2(i).setArgument(1, &reduceBufSet1[iBufOut]);
		resNorm2(i).setLocalMem(2, lworksizeTmp * sizeof(double));
	}
	resNorm2(numRedKer - 1).setSizes(globalWorkSizeRed[numRedKer - 1], localWorkSizeRed[numRedKer - 1]);
	resNorm2(numRedKer - 1).createKernel(GMRESGenerator::GMRESInstance()->getProgram(), "FinalNorm2");
	resNorm2(numRedKer - 1).setArgument(0, &reduceBufSet1[iBufOut]);
	resNorm2(numRedKer - 1).setArgument(1, rho_0.get_buf_add());
	resNorm2(numRedKer - 1).setLocalMem(2, localWorkSizeRed[numRedKer - 1] * sizeof(double));
}

void GMRESSolver::iniBuffer(cl_mem& buf_, const int size_, std::string name_)
{
	int status;
	double zer = 0.;
	buf_ = clCreateBuffer(*clEnv::instance()->getContext(), CL_MEM_READ_WRITE, sizeof(double) * size_, NULL, &status);
	ERROR_CHECKING_OCL(status, "Error creating " + name_ + " in BiCGStabGenerator Class", ERROR_CREATING_BICGSTAB_SOLVER);

	status = clEnqueueFillBuffer(*clEnv::instance()->getIOqueue(), buf_, &zer,
		sizeof(double), 0, sizeof(double) * size_, 0, NULL, NULL);
	ERROR_CHECKING_OCL(status, "Error filling " + name_ + " in BiCGStabGenerator Class", ERROR_CREATING_BICGSTAB_SOLVER);
}

// initializes all static variables
void GMRESSolver::ini(int xsize_, int xsizefull_, int ysize_)
{
	if (staticVarInitialized)
		return;
	Xsize = xsize_;
	Ysize = ysize_;
	XsizeFull = xsizefull_;
	colSize = XsizeFull * Ysize;
	int iSize1 = 0, iSize2 = 0;

	fullSize = ReduceGenerator::ReduceInstance()->iniSizesAndMemory(colSize, numRedKer,
		&globalWorkSizeRed, &reduceBufSet1, iSize1, iSize2);

	localWorkSizeRed = new int[numRedKer];
	for (int i = 0; i < numRedKer-2; i++)
		localWorkSizeRed[i] = GMRESGenerator::GMRESInstance()->reduce_group_size;
	localWorkSizeRed[numRedKer - 2] = GMRESGenerator::GMRESInstance()->reduceNextToLastSize;
	localWorkSizeRed[numRedKer - 1] = globalWorkSizeRed[numRedKer-1];

	iniBuffer(device_inner_prod_buffer, numBufChunks * bufSizePerVector, "devInnerProdVec");
	//iniBuffer(r_dot_vk_buffer, krylovDimSize* bufSizePerVector, "RdotVk");
	iniBuffer(vi_in_vk_buffer, krylovDimSize* bufSizePerVector, "viInVk");
	iniBuffer(resVec, fullSize, "resVec");
	iniBuffer(krylov_basis, colSize* krylovDimSize, "krylovBasis");
	iniBuffer(iniXVal, colSize, "iniXVal");

#ifdef _DEBUG
	iniBuffer(xMonitor, fullSize, "xMonitor");
	iniBuffer(residMonitor, fullSize, "residMonitor");
#endif



	staticVarInitialized = true;
}


void GMRESSolver::getRowBlockInfo()
{
	int num_entries_in_current_batch = 0;

	const size_t shared_mem_size = 1024; // number of column indices loaded to shared memory, number of floating point values loaded to shared memory

	int row_block_num_ = 0;
	Array1Du rowBlocks_;
	rowBlocks_.zeros(colSize + 1);
	for (cl_uint i = 0; i < colSize; ++i)
	{
		int entries_in_row = int(Inds->IA[i + 1]) - int(Inds->IA[i]);
		num_entries_in_current_batch += entries_in_row;
		//if (entries_in_row > 0)
		//	std::cout << "i = " << i << ", entries in row = " << entries_in_row << ", num current batch = " << num_entries_in_current_batch << "\n";

		if (num_entries_in_current_batch > shared_mem_size)
		{
			int rows_in_batch = i - rowBlocks_[row_block_num_];
			if (rows_in_batch > 0) // at least one full row is in the batch. Use current row in next batch.
				rowBlocks_[++row_block_num_] = i--;
			else // row is larger than buffer in shared memory
				rowBlocks_[++row_block_num_] = i + 1;
			num_entries_in_current_batch = 0;
		}
	}
	if (num_entries_in_current_batch > 0)
		rowBlocks_[++row_block_num_] = colSize;

	numRowBlocks = row_block_num_;

	if (row_block_num_ > 0) //matrix might be empty...
	{
		rowBlocks.zeros(numRowBlocks + 1);
		for (cl_uint i = 0; i < numRowBlocks + 1; i++)
		{
			rowBlocks[i] = rowBlocks_[i];
		}
		rowBlocks.allocate_buffer_w_copy(CL_MEM_READ_ONLY);
	}	
}

void GMRESSolver::CreateSolverWithVarTime(Array2Dd* macro_, Array2Dd* macroPrev_, std::function<void(void)>& resetTimeFunc_,
	CSR_Inds* inds_, cl_command_queue* calcque_, int maxiters_, double reltol_, double abstol_)
{
	resetTimeFunc = resetTimeFunc_;
	xVecPrev = macroPrev_;
	resetTimeFlag = true;
	CreateSolver(macro_, inds_, calcque_, maxiters_, reltol_, abstol_);
}

//size_t globalSizeCSRMV, globalSizeAXPBY;
void GMRESSolver::CreateSolver(Array2Dd* macro_, CSR_Inds* inds_, cl_command_queue* calcque_, int maxiters_, double reltol_, double abstol_)
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
	DebugOut.zeros(fullSize * krylovDimSize);
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

	std::string rname = Name + "_Rmat";
	std::string xi_k_name = Name + "_xi_k";
	//std::string eta_K_name = Name + "_eta_k";
	std::string updateCoeffname = Name + "_updateCoeffs";
	std::string RhoName = Name + "_Rho0";
	
	Rmat.setName(rname);
	values_xi_k.setName(xi_k_name);
	//values_eta_k.setName(eta_K_name);
	update_coefficients.setName(updateCoeffname);
	rho_0.setName(RhoName);

	Rmat.zeros(krylovDimSize * krylovDimSize);
	values_xi_k.zeros(krylovDimSize);
	//values_eta_k.zeros(krylovDimSize);
	update_coefficients.zeros(krylovDimSize);
	r_dot_vk_buffer.zeros(krylovDimSize * bufSizePerVector);
	r_dot_vk_buffer.allocate_buffer_w_copy();
	rho_0.zeros(1);

	Rmat.allocate_buffer_w_copy();
	values_xi_k.allocate_buffer_w_copy();
	//values_eta_k.allocate_buffer_w_copy();
	//update_coefficients.allocate_buffer_w_copy();
	rho_0.allocate_buffer_w_copy();

	maxRestarts = (cl_uint)(maxiters_ / krylovDimSize);

	getRowBlockInfo();
	createKernels();


}


void GMRESSolver::fillSolidBoundaryNodes()
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


void GMRESSolver::runReduce(RedKernelList & redlist_, cl_command_queue * que_, int num_wait,
	cl_event * wait, cl_event * evt)
{
	for (thinKerWrapper* rk = redlist_.begin(); rk != redlist_.end(); ++rk)
	{
		int status = rk->operator()(que_, num_wait, wait, evt);
		ERROR_CHECKING(status, ("Error returned running reduce kernels of " + Name), status);
	}
	clEnv::instance()->finishQueues();
}

bool GMRESSolver::reduceAndCheckConvergence(cl_command_queue * que_, bool setInitialRes, int num_wait,
	cl_event * wait, cl_event * evt)
{
	return true;
}


void GMRESSolver::copy_buffers(cl_mem * src_buf, cl_mem * dest_buf)
{
	int status = clEnqueueCopyBuffer(*clEnv::instance()->getIOqueue(),
		*src_buf, *dest_buf, 0, 0, colSize * sizeof(double), 0, NULL, NULL);
	ERROR_CHECKING(status, ("Error copying kernels in BiCGStabSolver " + Name), status);
}


void GMRESSolver::copyToPrevSolution()
{
	if (resetTimeFlag)
	{
		xVecPrev->enqueue_copy_to_buffer_blocking(xVec->get_buffer());
	}
}

void GMRESSolver::copyFromPrevSolution()
{
	if (resetTimeFlag)
	{
		xVec->enqueue_copy_to_buffer_blocking(xVecPrev->get_buffer());
	}
}

void GMRESSolver::debugTestResidual(double residual_)
{
	if (isnan(residual_))
	{
		bVec_copy.save_txt_from_device();
		xVec_copy.save_txt_from_device();
		vlb.saveDebug();
	}
}

void GMRESSolver::debugCheckForNans()
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


void GMRESSolver::fillArraysWithInitialValues()
{
	Rmat.FillBuffer(0, calcQue);
	bVec.enqueue_copy_from_buffer(resVec, fullSize, calcQue);
	xVec->enqueue_copy_from_buffer(iniXVal, fullSize, calcQue);
	// Not sure if this is necessary
	xVec->FillBuffer(0, calcQue);
	values_xi_k.FillBuffer(0, calcQue);
	rho_0(0) = vclReduceMethod();
	//runReduce(resNorm2, calcQue);
	//rho_0.read_from_buffer();
	//rho_0(0) = sqrt(rho_0(0));
}




bool GMRESSolver::checkConvergence()
{
	double residual = rho_0(0) / norm_rhs;

	if (resetTimeFlag)
	{
		if (isnan(residual))
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

	if (residual <= relTol || rho_0(0) <= absTol)
	{
		return true;
	}
	return false;
}


int GMRESSolver::getBufferFullSize()
{
	return fullSize;
}

std::string GMRESSolver::getName()
{
	return Name;
}

///////// Load and CheckPoint Methods ///////////



//////// Save Methods /////////////////////
bool GMRESSolver::saveAxbCSR()
{
	bool ret = savetxt();
	ret &= save_bvec();
	ret &= saveCSR(Name);
	return ret;
}

bool GMRESSolver::saveAxb_w_indicies()
{
	bool ret = savetxt();
	ret &= save_bvec();
	ret &= save_w_indicies(Name.append("_A"), false);
	return ret;
}

bool GMRESSolver::saveAxbCSR_from_device()
{
	bool ret = savetxt_from_device();
	ret &= save_bvec_from_device();
	ret &= saveCSR(Name, true);
	return ret;
}

bool GMRESSolver::saveAxb_w_indicies_from_device()
{
	bool ret = savetxt_from_device();
	ret &= save_bvec_from_device();
	std::string nameOut = Name;
	nameOut.append("_A");
	ret &= save_w_indicies(nameOut, true);
	return ret;
}

bool GMRESSolver::saveAxb_w_indicies_from_device_as_bin()
{
	bool ret = xVec->save_bin_from_device(Name);
	ret &= bVec.save_bin_from_device();
	std::string nameOut = Name;
	nameOut.append("_A");
	ret &= save_w_indicies_as_bin(nameOut, true);
	return ret;
}

cl_mem* GMRESSolver::get_add_A()
{
	return Amat.get_buf_add();
}

cl_mem* GMRESSolver::get_add_IndArr()
{
	return Inds->IndArray.get_buf_add();
}

bool GMRESSolver::savetxt(std::string outname)
{
	if (outname.length() == 0)
		outname = Name;
	return xVec->savetxt(outname);
}

bool GMRESSolver::savetxt_from_device(std::string outname)
{
	if (outname.length() == 0)
		outname = Name;
	return xVec->save_txt_from_device(outname);
}

bool GMRESSolver::save_bvec(std::string outname)
{
	if (outname.length() == 0)
		outname = Name + "_bvec";
	return bVec.savetxt_as_2D(Xsize, XsizeFull, Ysize, outname);
}

bool GMRESSolver::save_bvec_from_device(std::string outname)
{
	if (outname.length() == 0)
		outname = Name + "_bvec";
	return bVec.save_txt_from_device_as_2D(Xsize, XsizeFull, Ysize, outname);
}

// Amat (both values and indicies), bVec and other info can be 
// recreated, so no need to waste space and time saving it.
bool GMRESSolver::saveCheckPoint(std::string outname)
{
	return xVec->save_bin_from_device(outname);
}

void GMRESSolver::copy_to_device(const int blFlag)
{
	Amat.copy_to_buffer(NULL, blFlag);
}

void GMRESSolver::copy_to_host(const int blFlag)
{
	Amat.read_from_buffer(NULL, blFlag);
}

void GMRESSolver::copy_inds_to_host(const int blFlag)
{
	Inds->copy_to_host(NULL, blFlag);
}

void GMRESSolver::copy_inds_to_device(const int blFlag)
{
	Inds->copy_to_device(NULL, blFlag);
}

void GMRESSolver::copy_to_host_all(const int blFlag)
{
	copy_to_host(blFlag);
	copy_inds_to_host(blFlag);
}

void GMRESSolver::copy_to_device_all(const int blFlag)
{
	copy_to_device(blFlag);
	copy_inds_to_device(blFlag);
}

double GMRESSolver::A(const int i, const int j, int dir)
{
	int ind = Inds->getInd(i, j, dir);
	if (ind > -1)
		return Amat(ind);
	else
		return 0;
}

bool GMRESSolver::testInd(const int i, const int j, const int dir)
{
	if (Inds->getInd(i, j, dir) == -1)
		return false;
	return true;
}


bool GMRESSolver::save_w_indicies(std::string Name, bool fromDevFlag)
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


bool GMRESSolver::saveVec_As_MM(double* arr_, std::string Name, bool fromDevFlag)
{
	std::fstream stream;
	if (fileopen(stream, Name, FileOut) == false)
		return false;
	for (cl_uint i = 0; i < colSize; i++)
	{
		stream << arr_[i] << "\n";
	}

	stream.close();
	return true;
}

bool GMRESSolver::saveCSR_As_MM(std::string Name, bool fromDevFlag)
{
	if (fromDevFlag)
	{
		copy_to_host_all(true);
	}

	std::fstream stream;
	if (fileopen(stream, Name, FileOut) == false)
		return false;
	stream << "%% MatrixMarket matrix coordinate real general\n";
	stream << colSize << " " << colSize << " " << nnz() << "\n";
	for (int i = 0; i < nnz(); i++)
	{
		stream << (Inds->RowIndex(i) + 1) << " ";
		stream << (Inds->ColIndex(i) + 1) << " ";
		stream << Amat(i) << "\n";
	}
	stream.close();
	return true;
}

bool GMRESSolver::save_w_indicies_as_bin(std::string Name, bool fromDevFlag)
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




bool GMRESSolver::saveCSR(std::string outname, bool fromDevFlag)
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

bool GMRESSolver::saveCSR_row_col_val(std::string outname, bool fromDevFlag)
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
		ret &= Inds->saveJA(outname);
		ret &= Inds->saveRA(outname);
		outname.append("_A");
		ret = Amat.savetxt(outname);
	}
	return ret;
}