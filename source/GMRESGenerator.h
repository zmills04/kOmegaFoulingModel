// Sparse Matrix class for clSPARSE library
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GMRESGENERATOR_H__INCLUDED_)
#define AFX_GMRESGENERATOR_H__INCLUDED_
#pragma once

#include "StdAfx.h"
#include "SourceGenerator.h"
#include "Kernels.h"
#include "Array.h"

/*
passed: A, rhs = 
Variables:

	Matricies/Vectors:
		residual: initally set = rhs (bVec)
		result: initially set to 0's
		device_krylov_basis: size of bVec*krylov dimension
		device_buffer_R: size of krylov dimension^2
		device_inner_prod_buffer: num_buffer_chunks*buffer_size_per_vector
		device_r_dot_vk_buffer: buffer_size_per_vector*kryloc dimension
		device_vi_in_vk_buffer: buffer_size_per_vector*kryloc dimension
		device_values_xi_k: krylov dimension

	rhs = b - A*x0 (b becomes rhs, so b = b-A*x0)
	residual = rhs
	result = 0's
	norm_rhs = norm2(residual)
	rho_0 = norm_rhs
	rho = 1.0

for resCount = 0 to max_restarts
	1) if(resCount > 0)
		residual = rhs - A*result
		rho_0 = norm2(residual)
	2) test rho_0 <= abs_tolerance (trivial rhs)
		break
	
	3) residual /= rho_0
		
	4) rho = 1.0

	5) test if rho_0 / norm_rhs < tolerance or rho_0 < abs_tolerance
		break

	6) for k = 0 to krylov_dim
		6.1) if k == 0
			6.1.a) device_krylov_basis(:,0) = v_0 = cg_csr_prod(A, residual)
					device_inner_prod_buffer holds resulting inner product buffer
			6.1.b) gmres_normalize_vk(device_krylov_basis(:,0), residual, device_buffer_R,device_inner_prod_buffer)
					device_r_dot_vk_buffer holds resulting inner product (first stage) of <r, v_0>
			else
			6.1.c) compute v_k = A*v_k-1 and perform first reduction stage for ||v_k|| 
					i.e. v_0 = cg_csr_prod(A, residual) with device_inner_product_buffer storing partially reduced values

			6.1.d) pipelined_gmres_gram_schmidt_stage1(device_krylov_basis, rhs.size(), rhs.internal_size(), k, device_vi_in_vk_buffer, buffer_size_per_vector);
			6.1.e) pipelined_gmres_normalize_vk(vk, residual, device_buffer_R, k*tag.krylov_dim() + k,
                                                         device_inner_prod_buffer, device_r_dot_vk_buffer,
                                                         buffer_size_per_vector, k*buffer_size_per_vector);
	keep k iterated above (so should be = krylov_dim - 1 after for loop)
	7) complete reductions of host_r_dot_vk_buffer (10 total reductions)
		host_values_xi_k[i] = reduce(host_r_dot_vk_buffer[i,:]) (can be done with simple reduce kernel
	8) copy device_buffer_R to host		
		vcl_size_t full_krylov_dim = k;
		8.1) for i 0:k
				if(abs(device_buffer_R[i + i*krylov_dim]) < tolerance * device_buffer_R[0])
					k = i
					break
	9) compute error estimator:
		for i = 0:k
			increment iteration number
			if( host_values_xi_k[i] >= rho || host_values_xi_k[i] <= -rho)
				k = i
				break

			rho *= sin(acos(host_values_xi_k[i] / rho))

	10) solve minimization problem:
		set host_values_eta_k_buffer = host_values_xi_k
		for( i2 = k-1:-1:-1)
			i = i2
			for j = i+1:k
				host_values_eta_k_buffer[i] = host_buffer_R[i + j*full_krylov_dim] * host_values_eta_k_buffer[j]
						
			host_values_eta_k_buffer[i] /= host_buffer_R[i + i*full_krylov_dim]
		
	11) update x += rho * z, z = ...
		for i = 1:k
			host_update_coefficient[i] = rho_0*host_values_eta_k_buffer[i]

	12) copy host_update_coefficients to device_values_xi_k

	13) pipelined_gmres_update_result(result, residual,
                                                      device_krylov_basis, rhs.size(), rhs.internal_size(),
                                                      device_values_xi_k, k);

	14) can have additional convergence test and break if true here

	15)
*/		





class GMRESGenerator : public sourceGenerator
{
private:
	static GMRESGenerator* b_instance;
	GMRESGenerator(const std::string type_, const clEnv::programType ptype_)
	{
		type = type_;
		pType = ptype_;
	}

	~GMRESGenerator()
	{
	}

public:

	static GMRESGenerator* GMRESInstance()
	{
		if (!b_instance)
			b_instance = new GMRESGenerator("GMRESGenerator", clEnv::GMRESType);
		return b_instance;
	}

	static const std::string NormalizeResidualKernel;
	static const std::string NormalizeResidualKernelWithNorm2Base;
	static const std::string Norm2IntermediateKernel;
	static const std::string Norm2FinalKernel;
	static const std::string ReduceRVkKernel;
	static const std::string Norm2BaseKernel;
	static const std::string addVectors;
	static const std::string addIniGuess;

	static const cl_uint reduce_group_size;
	cl_uint reduceNextToLastSize;
	
	// Generates the string containing all kernels necessary for solving with BiCGS
	void ini(const int xsize_, const int xsizefull_, const int ysize_);

	std::string genFinalNorm2Kernel(int wgsize);
	std::string genIntermediateNorm2Kernel(int wgsize);
	std::string genIntermediateNorm2KernelNext2Last(int wgsize);

};


#endif