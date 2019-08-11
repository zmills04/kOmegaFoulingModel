// Sparse Matrix class for clSPARSE library
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SPARSEMATRIX_H__INCLUDED_)
#define AFX_SPARSEMATRIX_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

//#include <stdio.h>
//#include <memory.h>
//#include <CL\cl.h>
#include "StdAfx.h"
#include "Array.h"
#include <clsparse_export.h>


#ifdef _DEBUG
//#define ARRAY_DEBUG
#endif //_DEBUG


#ifndef BIN_SAVE_DIR
#define BIN_SAVE_DIR	""
#endif //BIN_SAVE_DIR

class CSR_Inds
{
public:
	enum { C, E, W, N, S }; // Enum of indicies in *IndArray (i.e. east neighbor of element (i,j) = neigh_index(i+j*XsizeFull,E)) 

	Array1Di IA; // row information
	Array1Di JA; // column information (i.e. the location in the domain that the element corresponds to (in 1D i.e. i+XsizeFull*j))
	Array1Di RA; // RA(i) gives the full matrix row index of element i in CSR matrix (provides easier access to info, mainly used for saving files)

	NTYPE_TYPE solidFlag, fluidFlag, solidBoundFlag;

	int nnz_, rows, cols; //number of nonzeros, rows in full matrix, cols in full matrix. (rows=cols=vlb.XsizeFull*vlb.nY)
	int XsizeFull, YsizeFull; // size of macro variable arrays
	int Xsize, Ysize; // will most likely not be needed, but could provide for bounds checking

	int IA_memflag, JA_memflag, AR_memflag;
	
// gives quick access to a given element and its neighbors in A array 
// (index of (i,j) in A array is IndArray(i+j*XsizeFull, 0 or C), its E neigh is IndArray(i,j,E))
	Array2Di IndArray;  

// Pointer to vls.nType
#ifdef IN_KERNEL_IBB
	Array2Di* Map;
#else
	Array2D<cl_short>* Map;
#endif


	CSR_Inds(NTYPE_TYPE solidflag_, NTYPE_TYPE fluidflag_,
		NTYPE_TYPE solidboundflag_) : solidFlag(solidflag_), 
	fluidFlag(fluidflag_), solidBoundFlag(solidboundflag_)
	{}

	~CSR_Inds() {}

	void ini(const int nx_, const int fx_, const int ny_, const int fy_, Array2D<NTYPE_TYPE> *map_)
	{
		ERROR_CHECKING((solidFlag == 0x0), "solidFlag in CSR_Inds is set to 0x0, indicating "\
			"that the nType flags were not set using the class constructor",\
			ERROR_INITIALIZING_CSR_CLASS);

		cols = fx_*fy_;
		Xsize = nx_;
		Ysize = ny_;
		XsizeFull = fx_;
		YsizeFull = ny_;
		rows = cols;
		IA.zeros(cols + 1);

		IndArray.zeros(cols, 5);
		Map = map_;

		ERROR_CHECKING((Map == nullptr), "NULL value passed for Map in CSR_Inds", \
			ERROR_INITIALIZING_CSR_CLASS);

		fill_IA();
		fill_JA();

		IA_memflag = CL_MEM_READ_WRITE;
		JA_memflag = CL_MEM_READ_WRITE;
		AR_memflag = CL_MEM_READ_WRITE;
		
	}

	void set_IAmemflag(int memflag_)
	{
		IA_memflag = memflag_;
	}

	void set_ARmemflag(int memflag_)
	{
		AR_memflag = memflag_;
	}

	void set_JAmemflag(int memflag_)
	{
		IA_memflag = memflag_;
	}

	void set_memflags(int memflag_)
	{
		IA_memflag = memflag_;
		JA_memflag = memflag_;
		AR_memflag = memflag_;
	}


	/// Fills Loc with the location in the sparse matrix, 
	void fill_IA()
	{
		//Array2Di countarr;
		//countarr.zeros(XsizeFull, YsizeFull);
		IndArray.fill(-1);
		nnz_ = 0;
		int row_num = 1;
		for (int j = 0; j < YsizeFull; j++)
		{
			for (int i = 0; i < XsizeFull; i++)
			{
				int localcount = testRow(i, j);
				//countarr(i, j) = localcount;
				IA(i + j*XsizeFull + 1) = nnz_;
			}
		}
		//countarr.savetxt("count");
		JA.zeros(nnz_);
		RA.zeros(nnz_);
	}

	//virtual int testRow(const int i, const int j)
	//{
	//	int localcount = 0;
	//	if (i >= Xsize || j >= Ysize)
	//		return 0;

	//	// If node is solid and not a (solid) boundary node, return 0
	//	if ((Map->operator()(i, j) & solidFlag) & !(Map->operator()(i, j) & solidFlag))
	//		return 0;
	//	
	//	int ind = i + j*XsizeFull;
	//	
	//	if (testNeigh(i, j, S) == true)
	//	{
	//		IndArray(ind, S) = nnz_++;
	//		localcount++;
	//	}
	//	if (testNeigh(i, j, W) == true)
	//	{
	//		IndArray(ind, W) = nnz_++;
	//		localcount++;
	//	}
	//	localcount++;
	//	IndArray(ind, C) = nnz_++;
	//	if (testNeigh(i, j, E) == true)
	//	{
	//		IndArray(ind, E) = nnz_++;
	//		localcount++;
	//	}
	//	if (testNeigh(i, j, N) == true)
	//	{
	//		IndArray(ind, N) = nnz_++;
	//		localcount++;
	//	}
	//	return localcount;
	//}

	virtual int getEastIndex(int i)
	{
		return i+1;
	}

	virtual int getWestIndex(int i)
	{
		return i-1;
	}

	virtual int testRow(const int i, const int j)
	{
		int localcount = 0;
		if (i >= Xsize || j >= Ysize)
			return 0;

		// If node is solid and not a (solid) boundary node, return 0
		if ((Map->operator()(i, j) & solidFlag) & !(Map->operator()(i, j) & solidFlag))
			return 0;

		int ind = i + j * XsizeFull;

		if (testNeigh(i, j, S) == true)
		{
			IndArray(ind, S) = nnz_++;
			localcount++;
		}

		// This is done to correctly order the indicies in ascending order when using periodic
		// boundaries (an inherited class)
		int ie = getEastIndex(i);
		int iw = getWestIndex(i);

		if (iw < i && i < ie)
		{
			if (testNeigh(i, j, W) == true)
			{
				IndArray(ind, W) = nnz_++;
				localcount++;
			}

			localcount++;
			IndArray(ind, C) = nnz_++;

			if (testNeigh(i, j, E) == true)
			{
				IndArray(ind, E) = nnz_++;
				localcount++;
			}
		}
		else if (ie < iw && iw < i)
		{
			if (testNeigh(i, j, E) == true)
			{
				IndArray(ind, E) = nnz_++;
				localcount++;
			}

			if (testNeigh(i, j, W) == true)
			{
				IndArray(ind, W) = nnz_++;
				localcount++;
			}

			localcount++;
			IndArray(ind, C) = nnz_++;
		}
		else
		{
			localcount++;
			IndArray(ind, C) = nnz_++;

			if (testNeigh(i, j, E) == true)
			{
				IndArray(ind, E) = nnz_++;
				localcount++;
			}

			if (testNeigh(i, j, W) == true)
			{
				IndArray(ind, W) = nnz_++;
				localcount++;
			}

		}

		if (testNeigh(i, j, N) == true)
		{
			IndArray(ind, N) = nnz_++;
			localcount++;
		}
		return localcount;
	}

	bool testNeigh(int i, int j, const int dir)
	{
		// if the center node is a solid boundary node,
		// nothing but the c coefficient is necessary,
		// so testNeigh should always return false.
		if (Map->operator()(i, j) & solidBoundFlag)
			return false;

		if (dir == E)
			getEastIndex(i);
		else if (dir == W)
			getWestIndex(i);
		else if (dir == N)
			j++;
		else if (dir == S)
			j--;
		else
			printf("error, calling testNeigh with incorrect direction\n");
				
		if (i < 0 || i >= Xsize || j < 0 || j >= Ysize)
			return false;

		// if neighbor is a fluid node or solid boundary, return true
		if (Map->operator()(i, j) & (solidBoundFlag | fluidFlag))
		{
			return true;
		}

		return false;
	}

	virtual void fill_JA()
	{
		int curel = 0;
		for (int j = 0; j < YsizeFull; j++)
		{
			for (int i = 0; i < XsizeFull; i++)
			{
				int ind = i + j*XsizeFull;
				if (IndArray(ind, C) == -1)
					continue;

				if (IndArray(ind, S) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i + (j - 1)*XsizeFull;
				}
				if (IndArray(ind, W) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i - 1 + (j)*XsizeFull;
				}

				RA(curel) = ind;
				JA(curel++) = ind;

				if (IndArray(ind, E) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i + 1 + (j)*XsizeFull;
				}

				if (IndArray(ind, N) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i + (j + 1)*XsizeFull;
				}
			}
		}



		if (curel != nnz_)
			printf("Error in initialization of CSR indicies");

	}

	/////////////////////////  Access and Utility functions ////////////////////
	int nnz()
	{
		return nnz_;
	}

	int getInd(const int i, const int j, int dir = C)
	{
		return IndArray(i+j*XsizeFull, dir);
	}


	void copy_to_device(cl_command_queue *que, int blockflag = CL_TRUE)
	{
		IA.copy_to_buffer(que, blockflag);
		JA.copy_to_buffer(que, blockflag);
		IndArray.copy_to_buffer(que, blockflag);
	}

	void copy_to_host(cl_command_queue *que, int blockflag = CL_TRUE)
	{
		IA.read_from_buffer(que, blockflag);
		JA.read_from_buffer(que, blockflag);
		IndArray.read_from_buffer(que, blockflag);
	}

	void copy_RA_to_device(cl_command_queue *que, int blockflag = CL_TRUE)
	{
		RA.copy_to_buffer(que, blockflag);
	}

	void allocate_buffers(cl_context context)
	{
		IA.allocate_buffer(IA_memflag);
		JA.allocate_buffer(JA_memflag);
		IndArray.allocate_buffer(AR_memflag);
		RA.allocate_buffer(AR_memflag);
	}

	int RowIndex(const int i)
	{
		return RA(i);
	}

	int ColIndex(const int i)
	{
		return JA(i);
	}


	bool saveall(std::string prefix_, cl_command_queue *devque = NULL)
	{
		bool ret = saveIA(prefix_, devque);
		ret &= saveJA(prefix_, devque);
		ret &= saveRA(prefix_, devque);
		ret &= saveIndArray(prefix_, devque);
		return ret;
	}

	bool saveIA(std::string append_, cl_command_queue *devque = NULL)
	{
		append_.append("_IA");
		if (!(devque == NULL))
			return IA.save_txt_from_device(append_, devque);
		else
			return IA.savetxt(append_);
	}

	bool saveRA(std::string append_, cl_command_queue *devque = NULL)
	{
		append_.append("_RA");
		if (!(devque == NULL))
			return RA.save_txt_from_device(append_, devque);
		else
			return RA.savetxt(append_);
	}

	bool saveJA(std::string append_, cl_command_queue *devque = NULL)
	{
		append_.append("_JA");
		if (!(devque == NULL))
			return JA.save_txt_from_device(append_, devque);
		else
			return JA.savetxt(append_);
	}

	bool saveIndArray(std::string append_, cl_command_queue *devque = NULL)
	{
		append_.append("_Inds");
		if (!(devque == NULL))
			return IndArray.save_txt_from_device(append_, devque);
		else
			return IndArray.savetxt(append_);
	}

	bool saveIndArray2DFormat(std::string append_, cl_command_queue *devque = NULL)
	{
		if (devque != NULL)
		{
			IndArray.read_from_buffer(devque);
		}

		Array2Di outarray(Xsize, Ysize);
		for (int i = 0; i < Xsize; i++)
		{
			for (int j = 0; j < Ysize; j++)
			{
				outarray(i, j) = getInd(i, j, C);
			}
		}
		append_.append("_Inds");
		return outarray.savetxt(append_); 

	}
};


class CSR_Inds_Periodic : public CSR_Inds
{
public:


	CSR_Inds_Periodic(NTYPE_TYPE solidflag_, NTYPE_TYPE fluidflag_,
		NTYPE_TYPE solidboundflag_) :
		CSR_Inds(solidflag_, fluidflag_, solidboundflag_)
	{};


	~CSR_Inds_Periodic() {};



	int getEastIndex(int i)
	{
		return (i < Xsize - 1) ? (i + 1) : (0);
	}

	int getWestIndex(int i)
	{
		return (i > 0) ? (i - 1) : (Xsize - 1);
	}

	void fill_JA()
	{
		int curel = 0;
		for (int j = 0; j < YsizeFull; j++)
		{
			for (int i = 0; i < XsizeFull; i++)
			{
				int ind = i + j*XsizeFull;
				if (IndArray(ind, C) == -1)
					continue;


				if (IndArray(ind, S) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i + (j - 1)*XsizeFull;
				}

				if (i == 0)
				{
					RA(curel) = ind;
					JA(curel++) = ind;

					if (IndArray(ind, E) >= 0)
					{
						RA(curel) = ind;
						JA(curel++) = i + 1 + j*XsizeFull;
					}

					if (IndArray(ind, W) >= 0)
					{
						RA(curel) = ind;
						JA(curel++) = Xsize - 1 + j*XsizeFull;;
					}
				}
				else if (i == Xsize-1)
				{
					if (IndArray(ind, E) >= 0)
					{
						RA(curel) = ind;
						JA(curel++) = j*XsizeFull;
					}

					if (IndArray(ind, W) >= 0)
					{
						RA(curel) = ind;
						JA(curel++) = i - 1 + j*XsizeFull;;
					}

					RA(curel) = ind;
					JA(curel++) = ind;
				}
				else
				{
					if (IndArray(ind, W) >= 0)
					{
						RA(curel) = ind;
						JA(curel++) = i - 1 + j*XsizeFull;;
					}

					RA(curel) = ind;
					JA(curel++) = ind;

					if (IndArray(ind, E) >= 0)
					{
						RA(curel) = ind;
						JA(curel++) = i + 1 + j*XsizeFull;
					}
				}

				if (IndArray(ind, N) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i + (j + 1)*XsizeFull;
				}
			}
		}

		if (curel != nnz_)
			printf("Error in initialization of CSR indicies");
	}

};


//
//
//template <NTYPE_TYPE solidFlag = M0_SOLID_NODE, NTYPE_TYPE fluidFlag = M0_FLUID_NODE,
//	NTYPE_TYPE solidBoundFlag = SOLID_BOUNDARY_NODE0, NTYPE_TYPE neswBoundFlag = NESW_BOUNDARY_NODE0>
//class CSR_Inds_Periodic_BL : public CSR_Inds_Periodic<solidFlag, fluidFlag, solidBoundFlag, neswBoundFlag>
//{
//public:
//
//	CSR_Inds_Periodic_BL() {}; // only using default constructor and destructor
//	~CSR_Inds_Periodic_BL() {};
//
//	/// Fills Loc with the location in the sparse matrix, 
//	void fill_IA()
//	{
//		IndArray.fill(-1);
//		nnz_ = 0;
//		int row_num = 1;
//		for (int j = 1; j < YsizeFull-1; j++)
//		{
//			for (int i = 0; i < XsizeFull; i++)
//			{
//				int localcount = testRow(i, j);
//				IA(i + j*XsizeFull + 1) = nnz_;
//			}
//		}
//		JA.zeros(nnz_);
//		RA.zeros(nnz_);
//	}
//
//	int testRow(const int i, const int j)
//	{
//		int localcount = 0;
//		if (i >= Xsize || j >= Ysize)
//			return 0;
//
//		if (Map->operator()(i, j) & solidFlag)
//			return 0;
//		int ind = i + j*XsizeFull;
//		if ((Map->operator()(i, j - 1) & solidFlag) || (Map->operator()(i, j + 1) & solidFlag))
//		{
//			localcount++;
//			IndArray(ind, C) = nnz_++;
//			return localcount;
//		}
//		
//		
//		int ie = (i < Xsize - 1) ? (i + 1) : (i);
//		int iw = (i > 0) ? (i - 1) : (Xsize - 1);
//		int jn = j + 1;
//		int js = j - 1;
//
//		if (js >= 0)
//		{
//			if (Map->operator()(i, js) & fluidFlag)
//			{
//				IndArray(ind, S) = nnz_++;
//				localcount++;
//			}
//		}
//
//		if (iw < i && i < ie)
//		{
//			if (Map->operator()(iw, j) & fluidFlag)
//			{
//				IndArray(ind, W) = nnz_++;
//				localcount++;
//			}
//
//			localcount++;
//			IndArray(ind, C) = nnz_++;
//
//			if (Map->operator()(ie, j) & fluidFlag)
//			{
//				IndArray(ind, E) = nnz_++;
//				localcount++;
//			}
//		}
//		else if (ie < iw && iw < i)
//		{
//			if (Map->operator()(ie, j) & fluidFlag)
//			{
//				IndArray(ind, E) = nnz_++;
//				localcount++;
//			}
//
//			if (Map->operator()(iw, j) & fluidFlag)
//			{
//				IndArray(ind, W) = nnz_++;
//				localcount++;
//			}
//
//			localcount++;
//			IndArray(ind, C) = nnz_++;
//		}
//		else
//		{
//			localcount++;
//			IndArray(ind, C) = nnz_++;
//
//			if (Map->operator()(ie, j) & fluidFlag)
//			{
//				IndArray(ind, E) = nnz_++;
//				localcount++;
//			}
//
//			if (Map->operator()(iw, j) & fluidFlag)
//			{
//				IndArray(ind, W) = nnz_++;
//				localcount++;
//			}
//		}
//
//		if (jn < Ysize - 1)
//		{
//			if (Map->operator()(i, jn) & fluidFlag)
//			{
//				IndArray(ind, N) = nnz_++;
//				localcount++;
//			}
//		}
//
//		return localcount;
//	}
//
//	void fill_JA()
//	{
//		int curel = 0;
//		for (int j = 0; j < YsizeFull; j++)
//		{
//			for (int i = 0; i < XsizeFull; i++)
//			{
//				int ind = i + j*XsizeFull;
//				if (IndArray(ind, C) == -1)
//					continue;
//
//
//				if (IndArray(ind, S) >= 0)
//				{
//					RA(curel) = ind;
//					JA(curel++) = i + (j - 1)*XsizeFull;
//				}
//
//				if (i == 0)
//				{
//					RA(curel) = ind;
//					JA(curel++) = ind;
//
//					if (IndArray(ind, E) >= 0)
//					{
//						RA(curel) = ind;
//						JA(curel++) = i + 1 + j*XsizeFull;
//					}
//
//					if (IndArray(ind, W) >= 0)
//					{
//						RA(curel) = ind;
//						JA(curel++) = Xsize - 1 + j*XsizeFull;;
//					}
//				}
//				else if (i == Xsize - 1)
//				{
//					if (IndArray(ind, E) >= 0)
//					{
//						RA(curel) = ind;
//						JA(curel++) = j*XsizeFull;
//					}
//
//					if (IndArray(ind, W) >= 0)
//					{
//						RA(curel) = ind;
//						JA(curel++) = i - 1 + j*XsizeFull;;
//					}
//
//					RA(curel) = ind;
//					JA(curel++) = ind;
//				}
//				else
//				{
//					if (IndArray(ind, W) >= 0)
//					{
//						RA(curel) = ind;
//						JA(curel++) = i - 1 + j*XsizeFull;;
//					}
//
//					RA(curel) = ind;
//					JA(curel++) = ind;
//
//					if (IndArray(ind, E) >= 0)
//					{
//						RA(curel) = ind;
//						JA(curel++) = i + 1 + j*XsizeFull;
//					}
//				}
//
//				if (IndArray(ind, N) >= 0)
//				{
//					RA(curel) = ind;
//					JA(curel++) = i + (j + 1)*XsizeFull;
//				}
//			}
//		}
//
//		if (curel != nnz_)
//			printf("Error in initialization of CSR indicies");
//	}
//
//};

//
//
//template <typename T>
//class CSR
//{
//public:
//	enum { C, E, W, N, S };
//	
//	// IA and JA arrays, pointer used to allow for same CSR_Inds to be shared across multiple
//	// instances, but right now, its not currently working, so use a different one for each
//	CSR_Inds *Inds;
//	
//	//int Xsize, Ysize; // domain sizes (i.e. rows >= Xsize*Ysize)
//	//int	XsizeFull, YsizeFull; // padded domain sizes (i.e. rows = XsizeFull*YsizeFull)
//	//int rows, cols; // sizes of full matrix (which CSR represents)
//	//int nnz, fullsize; //number of non-zeros and actual size of matrix a with padding
//
//	Array1D<T> Acpu; // A vector of CSR matrix
//
//	cl_context *context; // address of context and default que to avoid need to pass in multiple functions
//	cl_command_queue *que; 
//	
//	clsparseCsrMatrix clsCSR;
//	
//	CSR() {}; 
//	~CSR() 
//	{
//		clsparseCsrMetaDelete(&clsCSR);
//	};
//
//	void ini(CSR_Inds *inds_, cl_command_queue *queue_, cl_context* context_, clsparseCreateResult *clSparseControl_)
//	{
//		Inds = inds_;
//		context = context_;
//		que = queue_;
//
//		clsparseInitCsrMatrix(&clsCSR);
//		
//		clsCSR.num_nonzeros = Inds->nnz();
//		clsCSR.num_rows = Inds->rows;
//		clsCSR.num_cols = Inds->cols;
//
//		Inds->allocate_buffers(*context_);
//		Inds->copy_RA_to_device(*queue_);
//		Inds->copy_to_device(*queue_);
//
//		
//		Acpu.zeros(Inds->nnz());
//		Acpu.allocate_buffer_w_copy( CL_MEM_READ_WRITE);
//
//		clsCSR.values = Acpu.get_buffer();
//		clsCSR.col_indices = Inds->JA.get_buffer();
//		clsCSR.row_pointer = Inds->IA.get_buffer();
//
//		clsparseCsrMetaCreate(&clsCSR, clSparseControl_->control);
//		//Ind buffer allocation must be done separately to avoid multiple allocations
//	}
//
//	void copy_to_device(const int blFlag = true)
//	{
//		Acpu.copy_to_buffer(*que, blFlag);
//	}
//
//	void copy_to_host(const int blFlag = true)
//	{
//		Acpu.read_from_buffer(*que, blFlag);
//	}
//
//	void copy_inds_to_host(const int blFlag = true)
//	{
//		Inds->copy_to_host(*que, blFlag);
//	}
//
//	void copy_inds_to_device(const int blFlag = true)
//	{
//		Inds->copy_to_device(*que, blFlag);
//	}
//
//	void copy_to_host_all(const int blFlag = true)
//	{
//		copy_to_host(blFlag);
//		copy_inds_to_host(blFlag);
//	}
//
//	void copy_to_device_all(const int blFlag = true)
//	{
//		copy_to_device(blFlag);
//		copy_inds_to_device(blFlag);
//	}
//
//	int rows() { return Inds->rows };
//	int cols() { return Inds->cols };
//	int nnz() {	return Inds->nnz(); }
//
//	T& A(const int i, const int j, int dir = C)
//	{
//		int ind = Inds->getInd(i, j, dir);
//		if (ind > -1)
//			return &Aarr(ind);
//		else
//			return NULL;
//	}
//
//	bool testInd(const int i, const int j, const int dir = C)
//	{
//		if (Inds->getInd(i, j, dir) == -1)
//			return false;
//		return true;
//	}
//
//
//	bool save_w_indicies(std::string Name, bool fromDevFlag = false)
//	{
//		if (fromDevFlag)
//		{
//			copy_to_host_all(true);
//		}
//		Array2D<T> Outarray(nnz(), 3);
//
//		for (int i = 0; i < nnz(); i++)
//		{
//			Outarray(i, 0) = (T)(Inds->RowIndex(i));
//			Outarray(i, 1) = (T)(Inds->ColIndex(i));
//			Outarray(i, 2) = Acpu(i);
//		}
//
//		return Outarray.savetxt(Name);
//	}
//
//	bool saveCSR(std::string outname, bool fromDevFlag = false)
//	{
//		bool ret;
//		if (fromDevFlag)
//		{
//			ret &= Inds->saveIA(outname, que);
//			ret &= Inds->saveJA(outname, que);
//			outname.append("_A");
//			ret = Acpu.save_txt_from_device(outname, que);
//
//		}
//		else
//		{
//			ret &= Inds->saveIA(outname);
//			ret &= Inds->saveJA(outname);
//			outname.append("_A");
//			ret = Acpu.savetxt(outname);
//		}
//		return ret;
//	}
//
//
//	bool saveCSR_row_col_val(std::string outname, bool fromDevFlag = false)
//	{
//		bool ret;
//		if (fromDevFlag)
//		{
//			ret &= Inds->saveJA(outname, *que);
//			ret &= Inds->saveRA(outname, *que);
//			outname.append("_A");
//			ret = Acpu.save_txt_from_device(outname, *que);
//
//		}
//		else
//		{
//			ret &= Inds->saveJA(outname);
//			ret &= Inds->saveRA(outname);
//			outname.append("_A");
//			ret = Acpu.savetxt(outname);
//		}
//		return ret;
//	}
//
//};
//
//
//
//template <typename T>
//class FDEqSetBase
//{
//public:
//	typedef T type;
//
//	CSR<T> Amat;
//	cldenseVector xvec, bvec;
//	Array1D<T> Xcpu, Bcpu;
//
//	std::string Name;
//	
//	cl_command_queue *que;
//	cl_context* context;
//	int Xsize, Ysize, XsizeFull, YsizeFull, FullSize, nnz;
//	int FullSize_Red;
//	Array2Dc *Map;
//
//	clsparseCreateResult *clSparseResult;
//
//
//	FDEqSetBase(){};
//	~FDEqSetBase(){};
//
//	virtual void ini(std::string name_, CSR_Inds *inds_, clsparseCreateResult *clSparseControl_, cl_command_queue *queue_, cl_context *context_) = 0;
//
//
//	void setInitialValue(T inival, bool fullArrFlag = false)
//	{
//		if (fullArrFlag)
//		{
//			Xcpu.FillBuffer(inival);
//			return;
//		}
//
//		for (int i = 0; i < Xsize; i++)
//		{
//			for (int j = 0; j < Ysize; j++)
//			{
//				if (Map->operator()(i, j) == 1)
//				{
//					Xcpu(i + XsizeFull*j) = inival;
//				}
//			}
//		}
//				
//		Xcpu.copy_to_buffer();
//	}
//
//
//	void setInitialValueRows(T inival, std::vector<int> &rowi)
//	{
//		for (int i = 0; i < Xsize; i++)
//		{
//			for (int jj = 0; jj < rowi.size(); jj++)
//			{
//				int j = rowi[jj];
//				if (Map->operator()(i, j) == 1)
//				{
//					Xcpu(i + XsizeFull*j) = inival;
//				}
//			}
//		}
//		Xcpu.copy_to_buffer();
//	}
//
//	void setInitialValueCols(T inival, std::vector<int> &coli)
//	{
//		Xcpu.read_from_buffer(*que, CL_TRUE);
//		for (int ii = 0; ii < coli.size(); ii++)
//		{
//			int i = coli[ii];
//			for (int j = 0; j < Ysize; j++)
//			{
//				if (Map->operator()(i, j) == 1)
//				{
//					Xcpu(i + XsizeFull*j) = inival;
//				}
//			}
//		}
//
//		Xcpu.copy_to_buffer(*que, CL_TRUE);
//	}
//
//	void initializeCpuVecs()
//	{
//		Xcpu.zeros(FullSize_Red);
//		Bcpu.zeros(FullSize_Red);
//
//		Xcpu.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
//		Bcpu.allocate_buffer_w_copy(CL_MEM_READ_WRITE);
//	}
//
//	virtual T& operator()(int i, int j) = 0;
//
//	int getBufferFullSize()
//	{
//		return FullSize_Red;
//	}
//
//	virtual void solve() = 0;
//
//
//	cl_context* getContext()
//	{
//		return context;
//	}
//
//	cl_command_queue* getQueue()
//	{
//		return que;
//	}
//
//	std::string getName()
//	{
//		return Name;
//	}
//
/////////// Load and CheckPoint Methods ///////////
//
//
//
////////// Save Methods /////////////////////
//	bool saveAxbCSR()
//	{
//		bool ret = savetxt();
//		ret &= save_bvec();
//		ret &= Amat.saveCSR(Name);
//		return ret;
//	}
//
//	bool saveAxb_w_indicies()
//	{
//		bool ret = savetxt();
//		ret &= save_bvec();
//		ret &= Amat.save_w_indicies(Name.append("_A"), false);
//		return ret;
//	}
//
//	bool saveAxbCSR_from_device()
//	{
//		bool ret = savetxt_from_device();
//		ret &= save_bvec_from_device();
//		ret &= Amat.saveCSR(Name, true);
//		return ret;
//	}
//
//	bool saveAxb_w_indicies_from_device()
//	{
//		bool ret = savetxt_from_device();
//		ret &= save_bvec_from_device();
//		ret &= Amat.save_w_indicies(Name.append("_A"), true);
//		return ret;
//	}
//
//	cl_mem* get_add_A()
//	{
//		return Amat.Acpu.get_buf_add();
//	}
//
//	cl_mem* get_add_IndArr()
//	{
//		return Amat.Inds->IndArray.get_buf_add();
//	}
//
//	virtual bool savetxt() = 0;
//
//	virtual bool savetxt_from_device() = 0;
//
//	virtual bool save_bvec() = 0;
//
//	virtual bool save_bvec_from_device() = 0;
//
//	virtual bool save_bvec_as_vector_from_device() = 0;
//
//	virtual cl_mem* get_add_b() = 0;
//
//	virtual cl_mem* get_add_Macro() = 0;
//
//};
//
//template <typename T>
//class FDEqSetExplicit : public FDEqSetBase<T>
//{
//public:
//	typedef type T;
//
//	int alter;
//	clsparseScalar alpha;
//	clsparseScalar beta;
//
//	FDEqSetExplicit(){};
//	~FDEqSetExplicit(){};
//
//	void ini(std::string name_, CSR_Inds *inds_, clsparseCreateResult *clSparseControl_, cl_command_queue *queue_, cl_context *context_)
//	{
//		context = context_;
//		clSparseResult = clSparseControl_;
//		que = queue_;
//		Name = name_;
//
//		Amat.ini(inds_, queue_, context_, clSparseControl_);
//		Map = inds_->Map;
//		Xsize = inds_->Xsize;
//		Ysize = inds_->Ysize;
//		XsizeFull = inds_->XsizeFull;
//		YsizeFull = inds_->YsizeFull;
//		FullSize = XsizeFull*YsizeFull;
//
//		FullSize_Red = 2;
//		while (FullSize_Red < FullSize)
//			FullSize_Red *= 2;
//
//		alter = 0;
//
//		nnz = Amat.nnz();
//
//		que = queue_;
//
//		initializeCpuVecs();
//		initializeSpMvConstants();
//
//		xvec.num_values = FullSize;
//		xvec.values = Xcpu.get_buffer();
//
//		bvec.num_values = FullSize;
//		bvec.values = Bcpu.get_buffer();
//	}
//
//
//	void initializeSpMvConstants()
//	{
//		clsparseInitScalar(&alpha);
//		alpha.value = clCreateBuffer(*context, CL_MEM_READ_ONLY, sizeof(T), nullptr, nullptr);
//		clsparseInitScalar(&beta);
//		beta.value = clCreateBuffer(*context, CL_MEM_READ_ONLY, sizeof(T), nullptr, nullptr);
//
//		T one = 1.0;
//		T zero = 0.0;
//		// alpha = 1;
//		T *halpha = (T *)clEnqueueMapBuffer(*que, alpha.value, CL_TRUE, CL_MAP_WRITE, 0, sizeof(T), 0, nullptr, nullptr, nullptr);
//		*halpha = one;
//		clEnqueueUnmapMemObject(*que, alpha.value, halpha, 0, nullptr, nullptr);
//		//beta = 1;
//		T* hbeta = (T*)clEnqueueMapBuffer(*que, beta.value, CL_TRUE, CL_MAP_WRITE, 0, sizeof(T), 0, nullptr, nullptr, nullptr);
//		*hbeta = one;
//		clEnqueueUnmapMemObject(*que, beta.value, hbeta, 0, nullptr, nullptr);
//	}
//
//	T& operator()(int i, int j)
//	{
//		if (alter == 0)
//			return Xcpu(i + XsizeFull*j);
//		else
//			return Bcpu(i + XsizeFull*j);
//	}
//
//	void solve()
//	{
//		clsparseStatus status;
//		if (alter == 0)
//			status = clsparseDcsrmv(&alpha, &Amat.clsCSR, &xvec, &beta, &bvec, clSparseResult->control);
//		else
//			status = clsparseDcsrmv(&alpha, &Amat.clsCSR, &bvec, &beta, &xvec, clSparseResult->control);
//
//
//		if (status != clsparseSuccess)
//		{
//			std::cout << "Problem with execution SpMV algorithm."
//				<< " Error: " << status << std::endl;
//		}
//
//		alter ^= 1;
//	}
//
//
//
//	///////// Load and CheckPoint Methods ///////////
//
//
//
//	//////// Save Methods /////////////////////
//	bool savetxt()
//	{
//		if (alter == 0)
//			return Xcpu.savetxt_as_2D(Name, Xsize, XsizeFull, Ysize);
//		else
//			return Bcpu.savetxt_as_2D(Name, Xsize, XsizeFull, Ysize);
//	}
//
//	bool savetxt_from_device()
//	{
//		if (alter == 0)
//			return Xcpu.save_txt_from_device_as_2D(Name, *que, Xsize, XsizeFull, Ysize);
//		else
//			return Bcpu.save_txt_from_device_as_2D(Name, *que, Xsize, XsizeFull, Ysize);
//	}
//
//	bool save_bvec()
//	{
//		std::string outname = Name;
//		outname.append("_bvec");
//		if (alter == 0)
//			return Bcpu.savetxt_as_2D(outname, Xsize, XsizeFull, Ysize);
//		else
//			return Xcpu.savetxt_as_2D(outname, Xsize, XsizeFull, Ysize);
//	}
//
//	bool save_bvec_from_device()
//	{
//		std::string outname = Name;
//		outname.append("_bvec");
//		if (alter == 0)
//			return Bcpu.save_txt_from_device_as_2D(outname, *que, Xsize, XsizeFull, Ysize);
//		else
//			return Xcpu.save_txt_from_device_as_2D(outname, *que, Xsize, XsizeFull, Ysize);
//	}
//
//	bool save_bvec_as_vector_from_device()
//	{
//		std::string outname = Name;
//		outname.append("_bvec");
//		if (alter == 0)
//			return Bcpu.save_txt_from_device(outname, *que);
//		else
//			return Xcpu.save_txt_from_device(outname, *que);
//	}
//
//	cl_mem* get_add_Macro()
//	{
//		if (alter == 0)
//			return Xcpu.get_buf_add();
//		else
//			return Bcpu.get_buf_add();
//	}
//
//	cl_mem* get_add_b()
//	{
//		if (alter == 0)
//			return Bcpu.get_buf_add();
//		else
//			return Xcpu.get_buf_add();
//	}
//};
//
//
//
//
//template <typename T>
//class FDEqSetImplicit : public FDEqSetBase<T>
//{
//public:
//	typedef type T;
//
//	clsparseCreateSolverResult solverResult;
//
//
//	FDEqSetImplicit() {};
//	~FDEqSetImplicit() {};
//	
//
//	void ini(std::string name_, CSR_Inds *inds_, clsparseCreateResult *clSparseControl_, cl_command_queue *queue_, cl_context *context_)
//	{
//		clSparseResult = clSparseControl_;
//		que = queue_;
//		context = context_;
//		Name = name_;
//		Amat.ini(inds_, queue_, context_, clSparseControl_);
//		Map = inds_->Map;
//		Xsize = inds_->Xsize;
//		Ysize = inds_->Ysize;
//		XsizeFull = inds_->XsizeFull;
//		YsizeFull = inds_->YsizeFull;
//		FullSize = XsizeFull*YsizeFull;
//		
//		FullSize_Red = 2;
//		while (FullSize_Red < FullSize)
//			FullSize_Red *= 2;
//
//		nnz = Amat.nnz();
//
//		initializeCpuVecs();
//		
//		xvec.num_values = FullSize;
//		xvec.values = Xcpu.get_buffer();
//
//		bvec.num_values = FullSize;
//		bvec.values = Bcpu.get_buffer();
//	}
//
//
//	void createSolverControl(PRECONDITIONER pretype, double reltol, double abstol, int maxiters)
//	{
//		solverResult = clsparseCreateSolverControl(pretype, maxiters, reltol, abstol);
//		CLSPARSE_V(solverResult.status, "Failed to create clsparse solver control");
//		clsparseSolverPrintMode(solverResult.control, CLSPARSE_PRINT);
//	}
//
//	void solve()
//	{
//		clsparseStatus status = clsparseDcsrbicgStab(&xvec, &Amat.clsCSR, &bvec, solverResult.control, clSparseResult->control);
//		if (status != clsparseSuccess)
//		{
//			std::cout << "Problem with execution SpMV algorithm."
//				<< " Error: " << status << std::endl;
//		}
//		//status = clsparseReleaseSolverControl(solverResult.control);
//
//		//if (status != clsparseSuccess)
//		//{
//		//	std::cout << "Problem with releasing control object."
//		//		<< " Error: " << status << std::endl;
//		//}
//	}
//
//	T& operator()(int i, int j)
//	{
//		return Xcpu(i + XsizeFull*j);
//	}
//	
//	bool savetxt()
//	{
//		return Xcpu.savetxt_as_2D(Xsize, XsizeFull, Ysize, Name);
//	}
//
//	bool savetxt_from_device()
//	{
//		return Xcpu.save_txt_from_device_as_2D(Xsize, XsizeFull, Ysize, Name, que);
//	}
//
//	bool save_bvec()
//	{
//		std::string outname = Name;
//		outname.append("_bvec");
//		return Bcpu.savetxt_as_2D(Xsize, XsizeFull, Ysize, outname);
//	}
//
//	bool save_bvec_from_device()
//	{
//		std::string outname = Name;
//		outname.append("_bvec");
//		return Bcpu.save_txt_from_device_as_2D(Xsize, XsizeFull, Ysize, outname, que);
//	}
//
//	bool save_bvec_as_vector_from_device()
//	{
//		std::string outname = Name;
//		outname.append("_bvec");
//		return Bcpu.save_txt_from_device(outname, que);
//	}
//
//	cl_mem* get_add_b()
//	{
//		return Bcpu.get_buf_add();
//	}
//
//	cl_mem* get_add_Macro()
//	{
//		return Xcpu.get_buf_add();
//	}
//
//};
//
//
//
//
//
#endif