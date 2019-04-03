// Sparse Matrix class for clSPARSE library
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BICGSTABGENERATOR_H__INCLUDED_)
#define AFX_BICGSTABGENERATOR_H__INCLUDED_
#pragma once

#include "SourceGenerator.h"




/////////////////////////////////////////////////////////////////////////////
/*
Overview of the BiConjugate Gradient Stabilized method implementation

1) norm(b) - uses standard abs reduce method if norm(b) < machine precision,
no need to calculate x ( = zero vector). norm(b) also used in convergence checking

2) r_0 = b - A*x_0 - uses csrmv kernel (no modifications are made to this kernel,
so only one needed (here alpha = -1, beta = 1). r_0 stored in b array during solve since csrmv calculates
y = alpha*A*x + beta*y and b vector not needed during solution.


3) r_0 copied to rhat_0 and p_1

4) rho_1 calculated from <rhat_0, r_0> = <r_0, r_0> - uses Sqr reduce base kernel

5) rho_0 = alpha = omega_0 = 1 (scalars, which will share an array)

5) for i = 1,2,3,...     (note im1 = i-1)
if( i == 1) Skip f1-f3
f1) rho_i = <rhat_0, r_im1> - uses standard dot product reduce base kernel

f2) beta = (rho_i/rho_im1) * (alpha/omega_im1) - Final reduce step for step f1 will calculate beta as well

f3) p_i = r_im1 + beta*(p_im1 - omega*v_im1) - will use modified axpby kernel
(p_0 = v_0 = zero vector, so p_1 = r_0). This kernel will only be used here, so
name needs to be appended to be unique.

f4) v_i = A*p_i - will use standard csrmv kernel (alpha = 1, beta = 0)

f5) alpha = rho_i / <rhat_0, v_i> - will use standard dot product kernel w/ modified final kernel finalizining <>
diving rho_i by it to calculate alpha

f6) s = r_im1 - alpha*v_i
h = x_im1 + alpha*p_i
Will both be done in single modified axpby where ops are r_im1 -= alpha*v_i and x_im1 += alpha*p_i
with r buffer now storing s_i values and x buffer now storing h values.
Kernel name will again need to be unique.

f7) norm(s) -> copy to host and check convergence - if converged, h is result, which is already stored in x buffer

f8) t = A*s - uses standard csrmv kernel (alpha = 1, beta = 0)

f9) omega_i = <t,s>/<t,t> - Base kernels already made (dot and sqr), modified final kernel will calculate quotient
of resulting reduces and store in omega_i

f10) x_i = h + omega_i*s
r_i = s - omega_i*t - Uses modified axpby kernel, which ensures that xi is calculated before calculating r_i (since
that will modify s and x_i is a f(s)

f11) calc norm(r_i) -> copy to host and check convergence. If converged, x is solution. Uses standard norm reduce.

Note: reduce is performed by sequential kernels reducing array storing values down by a factor equal to the Work group size
each step. For most reduce operations these kernels will consist of a base reduce kernel where the reduce operator is applied
(i.e. the dot product is performed, or abs(), etc) along with the first reduction round. A generic summing reduce kernel is
used to continue reducing the array until it reaches a size less than the workgroup. Then a final reduce kernel is used to
sum the remaining elements of the array, and if necessary a modified version of this final reduce kernel can perform an
operation on the final summed values (i.e. divide it by another scalar, etc). This final reduce can also be modified to
perform the final reduction step on multiple arrays and perform use their final sums to calculate a single value to store
(which is done in step f5). This will work for any reduce operations that only need to be performed in the initial
reduction step. For some, such as max and min reductions, a modified kernel is needed for each reduction step since
max and min operations need to be applied at each step.

Note: any arrays that will be reduced need to be padded with zeros to give a length equal to a power of 2. If not, the
kernel will try to access out of bounds elements.

Kernels needed -
1) csrmv - only one needed for steps 2, f4, f8
2) modified axpby - one for step f3 and 1 for steps f6 and f10
1) norm base reduce for 1, f7 and f11
1) dot base reduce for f1, f5 and f9
1) sqr base reduce of 4, and f9
1) generic sum base reduce for intermediate steps
1) generic final reduce kernel for final norm reduction and step 4
3) modified final reduce for summing rho_i and calculating beta (step f1 and f2), calculation of alpha from dot product
reduction of rhat and v (step f5), and single final kernel to finish reduction of both dot products and calculation
of omega in step f9

scalars will be held in single array with order { norm_b, norm_r, norm_s, rho, omega, alpha, beta }  *note that
padding may be used between the scalars for reading/writing performance

Parameters for csrmv (clsparse's csrmv adaptive) need to be obtained by creating a clsparsecsrmatrix and calling
clsparsecsrmetacreate to generate the parameter values.
*/
///////////////////////////////////////////////////////////////////////////


// Note that very large domains may require the kernel generation functions to use ulong for the index type
// in order to avoid overflow at large element numbers


// Not deriving from reduce generator since there is a lot of overlap. A lot of
// methods are very similar with reduceGenerator ones and OOP could be used
// to cut down on code, but not worth the effort.
class BiCGStabGenerator : public sourceGenerator
{
private:
	static BiCGStabGenerator *b_instance;
	BiCGStabGenerator(const std::string type_, const clEnv::programType ptype_)
	{
		type = type_;
		pType = ptype_;
	}

	~BiCGStabGenerator()
	{
	}

public:

	static BiCGStabGenerator *BiCGStabInstance()
	{
		if (!b_instance)
			b_instance = new BiCGStabGenerator("BiCGStabGenerator", clEnv::BiCGSType);
		return b_instance;
	}



	enum BaseRedType {Sum, Dot, Norm, Abs, Sqr};
	
	static const std::string Axpy_kernel;
	static const std::string Axpby_kernel;
	static const std::string ReduceBase_kernel;
	static const std::string ReduceFinal_kernel;
	static const std::string Modified_Axpby_kernel;
	static const std::string ReduceBaseOmega_kernel;
	static const std::string ReduceIntermediateOmega_kernel;
	static const std::string ReduceFinalOmega_kernel;
	static const std::string ElementWiseOp_kernel;

	// These are the constants used by clsparse when calculating the matrix meta
	// defined in cpp file
	static const cl_uint WG_BITS;
	static const cl_uint ROW_BITS;
	static const cl_uint BLKSIZE;
	static const cl_uint BLOCK_MULTIPLIER;
	static const cl_uint ROWS_FOR_VECTOR;

	static const cl_uint csrmv_group_size;
	static const cl_uint reduce_group_size;
	static const cl_uint axpby_group_size;

	std::string index_t = "uint";
	clsparseCreateResult clSparseControl;

	// Generates the string containing all kernels necessary for solving with BiCGS
	void ini(const int xsize_, const int xsizefull_, const int ysize_);

	std::string genCSRMVKernel(int wgsize, int rowbits, int wgbits, int blocksize,
		int blockmult, int rows4vec, std::string indt);

	std::string genBaseRedKernel(int wgsize, BaseRedType redtype_, std::string name_);

	std::string genFinalRedKernel(int wgsize, std::string repIn_, std::string repOut_, std::string repLocal_,
		std::string repOff_, std::string repOpIn_, std::string repOpLocal_, std::string repOpOut_, std::string name_);

	std::string genAxpyKernel(int wgsize, std::string indType, std::string oper,
		std::string repname_, int offset = 0);

	std::string genAxpbyKernel(int wgsize, std::string indType, std::string oper,
		std::string repname_, int offset = 0);

	std::string genModifiedAxpbyKernel(int wgsize, std::string indType,
		std::string repin_, std::string repop_, std::string repname_, int offset = 0);

	std::string genElementWiseOperation(int wgsize, std::string indType,
		std::string opdefine_, std::string repname_);

	void deleteClSparseObjects()
	{
		clsparseReleaseControl(BiCGStabGenerator::BiCGStabInstance()->clSparseControl.control);
		clsparseTeardown();
	}
};




#endif // AFX_BICGSTABGENERATOR_H__INCLUDED_