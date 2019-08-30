#include "SparseMatrix.h"


void CSR_Inds::ini(const int nx_, const int fx_, const int ny_, const int fy_, Array2D<NTYPE_TYPE>* map_)
{
	ERROR_CHECKING((solidFlag == 0x0), "solidFlag in CSR_Inds is set to 0x0, indicating "\
		"that the nType flags were not set using the class constructor", \
		ERROR_INITIALIZING_CSR_CLASS);

	cols = fx_ * fy_;
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

void CSR_Inds::set_IAmemflag(int memflag_)
{
	IA_memflag = memflag_;
}

void CSR_Inds::set_ARmemflag(int memflag_)
{
	AR_memflag = memflag_;
}

void CSR_Inds::set_JAmemflag(int memflag_)
{
	IA_memflag = memflag_;
}

void CSR_Inds::set_memflags(int memflag_)
{
	IA_memflag = memflag_;
	JA_memflag = memflag_;
	AR_memflag = memflag_;
}


/// Fills Loc with the location in the sparse matrix, 
void CSR_Inds::fill_IA()
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
			IA(i + j * XsizeFull + 1) = nnz_;
		}
	}
	//countarr.savetxt("count");
	JA.zeros(nnz_);
	RA.zeros(nnz_);
}

int CSR_Inds::testRow(const int i, const int j)
{
	int localcount = 0;
	if (i >= Xsize || j >= Ysize)
		return 0;

	int ind = i + j * XsizeFull;

	// Must handle if node is a solid boundary node, so first test for
	// if node is solid
	if ((Map->operator()(i, j) & solidFlag))
	{
		// if node is a solid boundary node, only C element (diagonal) is
		// a non-zero, and its coefficient is set to 1.
		if (Map->operator()(i, j) & solidBoundFlag)
		{
			IndArray(ind, C) = nnz_++;
			return 1;
		}
		else
		{// if not a solid boundary node, return
			return 0;
		}
	}



	if (testNeigh(i, j, S) == true)
	{
		IndArray(ind, S) = nnz_++;
		localcount++;
	}

	// This is done to correctly order the indicies in ascending order when using periodic
	// boundaries (an inherited class)
	int ie = getEastIndex(i);
	int iw = getWestIndex(i);

	if (i == 0 && j == 3)
		int a = 3;

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

bool CSR_Inds::testNeigh(int i, int j, const int dir)
{
	// if the center node is a solid boundary node,
	// nothing but the c coefficient is necessary,
	// so testNeigh should always return false.
	if (Map->operator()(i, j) & solidBoundFlag)
		return false;

	if (dir == E)
		i = getEastIndex(i);
	else if (dir == W)
		i = getWestIndex(i);
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

void CSR_Inds::fill_JA()
{
	int curel = 0;
	for (int j = 0; j < YsizeFull; j++)
	{
		for (int i = 0; i < XsizeFull; i++)
		{
			int ind = i + j * XsizeFull;
			if (IndArray(ind, C) == -1)
				continue;

			if (IndArray(ind, S) >= 0)
			{
				RA(curel) = ind;
				JA(curel++) = i + (j - 1) * XsizeFull;
			}
			if (IndArray(ind, W) >= 0)
			{
				RA(curel) = ind;
				JA(curel++) = i - 1 + (j)* XsizeFull;
			}

			RA(curel) = ind;
			JA(curel++) = ind;

			if (IndArray(ind, E) >= 0)
			{
				RA(curel) = ind;
				JA(curel++) = i + 1 + (j)* XsizeFull;
			}

			if (IndArray(ind, N) >= 0)
			{
				RA(curel) = ind;
				JA(curel++) = i + (j + 1) * XsizeFull;
			}
		}
	}



	ERROR_CHECKING((curel != nnz_), "Error in initialization of CSR indicies",
		ERROR_INITIALIZING_CSR_CLASS);

}

bool CSR_Inds::testSolidNode(const int i, const int j)
{
	if (Map->operator()(i, j) & solidFlag)
		return true;
	return false;
}
bool CSR_Inds::testSolidBoundaryNode(const int i, const int j)
{
	if (Map->operator()(i, j) & solidBoundFlag)
		return true;
	return false;
}
bool CSR_Inds::testFluidNode(const int i, const int j)
{
	if (Map->operator()(i, j) & fluidFlag)
		return true;
	return false;
}

void CSR_Inds::copy_to_device(cl_command_queue* que, int blockflag)
{
	IA.copy_to_buffer(que, blockflag);
	JA.copy_to_buffer(que, blockflag);
	IndArray.copy_to_buffer(que, blockflag);
}

void CSR_Inds::copy_to_host(cl_command_queue* que, int blockflag)
{
	IA.read_from_buffer(que, blockflag);
	JA.read_from_buffer(que, blockflag);
	IndArray.read_from_buffer(que, blockflag);
}

void CSR_Inds::copy_RA_to_device(cl_command_queue* que, int blockflag)
{
	RA.copy_to_buffer(que, blockflag);
}

void CSR_Inds::allocate_buffers(cl_context context)
{
	IA.allocate_buffer(IA_memflag);
	JA.allocate_buffer(JA_memflag);
	IndArray.allocate_buffer(AR_memflag);
	RA.allocate_buffer(AR_memflag);
}


bool CSR_Inds::saveall(std::string prefix_, cl_command_queue* devque)
{
	bool ret = saveIA(prefix_, devque);
	ret &= saveJA(prefix_, devque);
	ret &= saveRA(prefix_, devque);
	ret &= saveIndArray(prefix_, devque);
	return ret;
}

bool CSR_Inds::saveIA(std::string append_, cl_command_queue* devque)
{
	append_.append("_IA");
	if (!(devque == NULL))
		return IA.save_txt_from_device(append_, devque);
	else
		return IA.savetxt(append_);
}

bool CSR_Inds::saveRA(std::string append_, cl_command_queue * devque)
{
	append_.append("_RA");
	if (!(devque == NULL))
		return RA.save_txt_from_device(append_, devque);
	else
		return RA.savetxt(append_);
}

bool CSR_Inds::saveJA(std::string append_, cl_command_queue * devque)
{
	append_.append("_JA");
	if (!(devque == NULL))
		return JA.save_txt_from_device(append_, devque);
	else
		return JA.savetxt(append_);
}

bool CSR_Inds::saveIndArray(std::string append_, cl_command_queue * devque)
{
	append_.append("_Inds");
	if (!(devque == NULL))
		return IndArray.save_txt_from_device(append_, devque);
	else
		return IndArray.savetxt(append_);
}

bool CSR_Inds::saveIndArray2DFormat(std::string append_, cl_command_queue * devque)
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


void CSR_Inds_Periodic::fill_JA()
{
	int curel = 0;
	for (int j = 0; j < YsizeFull; j++)
	{
		for (int i = 0; i < XsizeFull; i++)
		{
			int ind = i + j * XsizeFull;
			if (IndArray(ind, C) == -1)
				continue;


			if (IndArray(ind, S) >= 0)
			{
				RA(curel) = ind;
				JA(curel++) = i + (j - 1) * XsizeFull;
			}

			if (i == 0)
			{
				RA(curel) = ind;
				JA(curel++) = ind;

				if (IndArray(ind, E) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i + 1 + j * XsizeFull;
				}

				if (IndArray(ind, W) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = Xsize - 1 + j * XsizeFull;;
				}
			}
			else if (i == Xsize - 1)
			{
				if (IndArray(ind, E) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = j * XsizeFull;
				}

				if (IndArray(ind, W) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i - 1 + j * XsizeFull;;
				}

				RA(curel) = ind;
				JA(curel++) = ind;
			}
			else
			{
				if (IndArray(ind, W) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i - 1 + j * XsizeFull;;
				}

				RA(curel) = ind;
				JA(curel++) = ind;

				if (IndArray(ind, E) >= 0)
				{
					RA(curel) = ind;
					JA(curel++) = i + 1 + j * XsizeFull;
				}
			}

			if (IndArray(ind, N) >= 0)
			{
				RA(curel) = ind;
				JA(curel++) = i + (j + 1) * XsizeFull;
			}
		}
	}

	ERROR_CHECKING((curel != nnz_), "Error in initialization of CSR indicies",
		ERROR_INITIALIZING_CSR_CLASS);

}
