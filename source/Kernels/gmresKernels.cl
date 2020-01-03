#pragma OPENCL EXTENSION cl_amd_fp64 : enable



__kernel void normVCL( 
          __global const double * vec, 
          unsigned int size1, 
          __local double * tmp_buffer, 
          __global double * group_buffer) 
{ 
    unsigned int start1_ = (get_group_id(0)  * size1) / get_num_groups(0);
    unsigned int size1_ = ((1 + get_group_id(0)) * size1) / get_num_groups(0) -
                  (get_group_id(0)*size1) / get_num_groups(0); 
    double tmp = 0.0; 
    double vec_entry = 0; 
    for (unsigned int i = get_local_id(0); i < size1_; i += get_local_size(0)) 
    { 
      vec_entry = vec[i + start1_]; 
      tmp += vec_entry * vec_entry; 
    } 
    tmp_buffer[get_local_id(0)] = tmp; 
    for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2) 
    { 
      barrier(CLK_LOCAL_MEM_FENCE); 
      if (get_local_id(0) < stride) 
        tmp_buffer[get_local_id(0)] += tmp_buffer[get_local_id(0)+stride]; 
    } 

  if (get_local_id(0) == 0) 
    group_buffer[get_group_id(0)] = tmp_buffer[0]; 
} 

__kernel void cg_csr_prod( 
  __global const unsigned int * row_indices, 
  __global const unsigned int * column_indices, 
  __global const unsigned int * row_blocks, 
  __global const double * elements, 
  unsigned int num_blocks, 
  __global const double * p, 
  __global double * Ap, 
  unsigned int size, 
  __global double * inner_prod_buffer, 
  unsigned int buffer_size, 
  __local double * shared_array_ApAp, 
  __local double * shared_array_pAp, 
  __local double * shared_elements) 
{ 
  double inner_prod_ApAp = 0; 
  double inner_prod_pAp = 0; 
  for (unsigned int block_id = get_group_id(0); block_id < num_blocks; block_id += get_num_groups(0)) { 
    unsigned int row_start = row_blocks[block_id]; 
    unsigned int row_stop  = row_blocks[block_id + 1]; 
    unsigned int rows_to_process = row_stop - row_start; 
    unsigned int element_start = row_indices[row_start]; 
    unsigned int element_stop = row_indices[row_stop]; 
    if (rows_to_process > 1) { 
      for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0)) 
        shared_elements[i - element_start] = elements[i] * p[column_indices[i]]; 
      barrier(CLK_LOCAL_MEM_FENCE); 
      for (unsigned int row = row_start + get_local_id(0); row < row_stop; row += get_local_size(0)) { 
        double dot_prod = 0; 
        unsigned int thread_row_start = row_indices[row]     - element_start; 
        unsigned int thread_row_stop  = row_indices[row + 1] - element_start; 
        for (unsigned int i = thread_row_start; i < thread_row_stop; ++i) 
          dot_prod += shared_elements[i]; 
        Ap[row] = dot_prod; 
        inner_prod_ApAp += dot_prod * dot_prod; 
        inner_prod_pAp  +=   p[row] * dot_prod; 
      } 
    } 
    else  
    { 
      shared_elements[get_local_id(0)] = 0; 
      for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0)) 
        shared_elements[get_local_id(0)] += elements[i] * p[column_indices[i]]; 
      for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2) { 
        barrier(CLK_LOCAL_MEM_FENCE); 
        if (get_local_id(0) < stride) 
          shared_elements[get_local_id(0)] += shared_elements[get_local_id(0) + stride]; 
      } 
      if (get_local_id(0) == 0) { 
        Ap[row_start] = shared_elements[0]; 
        inner_prod_ApAp += shared_elements[0] * shared_elements[0]; 
        inner_prod_pAp  +=       p[row_start] * shared_elements[0]; 
      } 
    } 
    barrier(CLK_LOCAL_MEM_FENCE); 
  } 
  shared_array_ApAp[get_local_id(0)] = inner_prod_ApAp; 
  shared_array_pAp[get_local_id(0)]  = inner_prod_pAp; 
  for (uint stride=get_local_size(0)/2; stride > 0; stride /= 2) 
  { 
    barrier(CLK_LOCAL_MEM_FENCE); 
    if (get_local_id(0) < stride) { 
      shared_array_ApAp[get_local_id(0)] += shared_array_ApAp[get_local_id(0) + stride];  
      shared_array_pAp[get_local_id(0)]  += shared_array_pAp[get_local_id(0) + stride];  
    }   }   if (get_local_id(0) == 0) { 
     inner_prod_buffer[  buffer_size + get_group_id(0)] = shared_array_ApAp[0]; 
    inner_prod_buffer[2*buffer_size + get_group_id(0)] = shared_array_pAp[0]; 
  } 
} 

 __kernel void gmres_csr_prod( 
  __global const unsigned int * row_indices, //0
  __global const unsigned int * column_indices, //1 
  __global const unsigned int * row_blocks, //2
  __global const double * elements, //3
  unsigned int num_blocks, //4
  __global const double * p, //5
  unsigned int offset_p, //6
  __global double * Ap, //7
  unsigned int offset_Ap, //8
  unsigned int size, //9
  __global double * inner_prod_buffer, //10
  unsigned int buffer_size, //11
  __local double * shared_array_ApAp, //12 
  __local double * shared_array_pAp, //13
  __local double * shared_elements) //14 
{ 
  cg_csr_prod(row_indices, column_indices, row_blocks, elements, num_blocks, p + offset_p, Ap + offset_Ap, size, inner_prod_buffer, buffer_size, shared_array_ApAp, shared_array_pAp, shared_elements); 
} 


__kernel void gmres_gram_schmidt_1( 
          __global double const * krylov_basis, 
          unsigned int size, 
          unsigned int internal_size, 
          unsigned int k, 
          __global double * vi_in_vk_buffer, 
          unsigned int chunk_size) 
{ 
  __local double shared_array[7*128]; 
  double vi_in_vk[7]; 
  double value_vk = 0; 
  unsigned int k_base = 0;   
  while (k_base < k) {   
    unsigned int vecs_in_iteration = (k - k_base > 7) ? 7 : (k - k_base);   
    vi_in_vk[0] = 0;
    vi_in_vk[1] = 0;
    vi_in_vk[2] = 0;
    vi_in_vk[3] = 0;
    vi_in_vk[4] = 0;
    vi_in_vk[5] = 0;
    vi_in_vk[6] = 0;
    for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) { 
      value_vk = krylov_basis[i + k * internal_size]; 
       
      for (unsigned int j=0; j<vecs_in_iteration; ++j) 
        vi_in_vk[j] += value_vk * krylov_basis[i + (k_base + j) * internal_size]; 
    }  
    for (uint j=0; j<vecs_in_iteration; ++j) 
      shared_array[get_local_id(0) + j*chunk_size] = vi_in_vk[j]; 
    for (uint stride=get_local_size(0)/2; stride > 0; stride /= 2) 
    { 
      barrier(CLK_LOCAL_MEM_FENCE); 
      if (get_local_id(0) < stride) { 
        for (uint j=0; j<vecs_in_iteration; ++j) 
          shared_array[get_local_id(0) + j*chunk_size] += shared_array[get_local_id(0) + j*chunk_size + stride];  
      }     }     if (get_local_id(0) == 0) 
       for (unsigned int j=0; j<vecs_in_iteration; ++j) 
        vi_in_vk_buffer[get_group_id(0) + (k_base + j) * chunk_size] = shared_array[j*chunk_size];     k_base += vecs_in_iteration;   
  }  
} 



__kernel void gmres_gram_schmidt_2( 
          __global double * krylov_basis, 
          unsigned int size, 
          unsigned int internal_size, 
          unsigned int k, 
          __global double const * vi_in_vk_buffer, 
          unsigned int chunk_size, 
          __global double * R_buffer, 
          unsigned int krylov_dim, 
          __global double * inner_prod_buffer, 
         __local double * shared_array) 
{ 
  double vk_dot_vk = 0; 
  double value_vk = 0; 
  unsigned int k_base = 0;   
  while (k_base < k) {   
    unsigned int vecs_in_iteration = (k - k_base > 7) ? 7 : (k - k_base);   
    for (uint j=0; j<vecs_in_iteration; ++j) 
      shared_array[get_local_id(0) + j*chunk_size] = vi_in_vk_buffer[get_local_id(0) + (k_base + j) * chunk_size]; 
    for (uint stride=get_local_size(0)/2; stride > 0; stride /= 2) 
    { 
      barrier(CLK_LOCAL_MEM_FENCE); 
      if (get_local_id(0) < stride) { 
        for (uint j=0; j<vecs_in_iteration; ++j) 
          shared_array[get_local_id(0) + j*chunk_size] += shared_array[get_local_id(0) + j*chunk_size + stride];  
      }     }     barrier(CLK_LOCAL_MEM_FENCE); 
    for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) { 
      value_vk = krylov_basis[i + k * internal_size]; 
       
      for (unsigned int j=0; j<vecs_in_iteration; ++j) 
        value_vk -= shared_array[j*chunk_size] * krylov_basis[i + (k_base + j) * internal_size]; 
      vk_dot_vk += (k_base + vecs_in_iteration == k) ? (value_vk * value_vk) : 0;  
      krylov_basis[i + k * internal_size] = value_vk;  
    }  
    if (get_group_id(0) == 0) 
      for (unsigned int j=0; j<vecs_in_iteration; ++j) 
        R_buffer[(k_base + j) + k*krylov_dim] = shared_array[j*chunk_size];     barrier(CLK_LOCAL_MEM_FENCE); 
    k_base += vecs_in_iteration;   
  }  
  shared_array[get_local_id(0)] = vk_dot_vk; 
  for (uint stride=get_local_size(0)/2; stride > 0; stride /= 2) 
  { 
    barrier(CLK_LOCAL_MEM_FENCE); 
    if (get_local_id(0) < stride) 
      shared_array[get_local_id(0)] += shared_array[get_local_id(0) + stride];  
  }   if (get_local_id(0) == 0) 
     inner_prod_buffer[chunk_size+get_group_id(0)] = shared_array[0]; } 



__kernel void gmres_normalize_vk( 
          __global double * vk, //0
          unsigned int vk_offset, //1
          __global double const * residual, //2
          __global double * R_buffer, //3
          unsigned int R_offset, //4
          __global double const * inner_prod_buffer, //5 
          unsigned int chunk_size, //6
          __global double * r_dot_vk_buffer, //7 
          unsigned int chunk_offset, //8
          unsigned int size, //9
         __local double * shared_array) //10
{ 
  double norm_vk = 0; 
  shared_array[get_local_id(0)] = inner_prod_buffer[get_local_id(0) + chunk_size]; 
  for (uint stride=get_local_size(0)/2; stride > 0; stride /= 2) 
  { 
    barrier(CLK_LOCAL_MEM_FENCE); 
    if (get_local_id(0) < stride) 
      shared_array[get_local_id(0)]  += shared_array[get_local_id(0) + stride];  
  }   barrier(CLK_LOCAL_MEM_FENCE); 
  norm_vk = sqrt(shared_array[0]); 
  double inner_prod_contrib = 0; 
  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) { 
    double value_vk = vk[i + vk_offset] / norm_vk; 
     
    inner_prod_contrib += residual[i] * value_vk; 
     
    vk[i + vk_offset] = value_vk; 
  }  
  barrier(CLK_LOCAL_MEM_FENCE); 
  shared_array[get_local_id(0)] = inner_prod_contrib; 
  for (uint stride=get_local_size(0)/2; stride > 0; stride /= 2) 
  { 
    barrier(CLK_LOCAL_MEM_FENCE); 
    if (get_local_id(0) < stride)  
      shared_array[get_local_id(0)] += shared_array[get_local_id(0) + stride];  
  }   if (get_local_id(0) == 0) 
     r_dot_vk_buffer[get_group_id(0) + chunk_offset] = shared_array[0];   if (get_global_id(0) == 0) 
     R_buffer[R_offset] = norm_vk; 
} 



__kernel void gmres_update_result( 
          __global double * result, 
          __global double const * residual, 
          __global double const * krylov_basis, 
          unsigned int size, 
          unsigned int internal_size, 
          __global double const * coefficients, 
          unsigned int k) 
{ 
  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) 
  { 
    double value_result = result[i] + coefficients[0] * residual[i]; 
    for (unsigned int j = 1; j < k; ++j) 
    {
      value_result += coefficients[j] * krylov_basis[i + (j-1)*internal_size]; 
    }
    result[i] = value_result; 
  }  
} 


__kernel void vec_mul( 
  __global const unsigned int * row_indices, 
  __global const unsigned int * column_indices, 
  __global const unsigned int * row_blocks, 
  __global const double * elements, 
  unsigned int num_blocks, 
  __global const double * x, 
  uint4 layout_x, 
  __global double * result, 
  uint4 layout_result 
) { 
  __local double shared_elements[1024]; 
  unsigned int row_start = row_blocks[get_group_id(0)]; 
  unsigned int row_stop  = row_blocks[get_group_id(0) + 1]; 
  unsigned int rows_to_process = row_stop - row_start; 
  unsigned int element_start = row_indices[row_start]; 
  unsigned int element_stop = row_indices[row_stop]; 
  if (rows_to_process > 4) { 
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0)) 
      shared_elements[i - element_start] = elements[i] * x[column_indices[i] * layout_x.y + layout_x.x]; 
    barrier(CLK_LOCAL_MEM_FENCE); 
    for (unsigned int row = row_start + get_local_id(0); row < row_stop; row += get_local_size(0)) { 
      double dot_prod = 0; 
      unsigned int thread_row_start = row_indices[row]     - element_start; 
      unsigned int thread_row_stop  = row_indices[row + 1] - element_start; 
      for (unsigned int i = thread_row_start; i < thread_row_stop; ++i) 
        dot_prod += shared_elements[i]; 
      result[row * layout_result.y + layout_result.x] = dot_prod; 
    } 
  } 
  else if (rows_to_process > 1) 
  {
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0))
      shared_elements[i - element_start] = elements[i] * x[column_indices[i] * layout_x.y + layout_x.x];
    barrier(CLK_LOCAL_MEM_FENCE); 
    for (unsigned int row = row_start; row < row_stop; ++row)
    {
      unsigned int current_row_start = row_indices[row]     - element_start;
      unsigned int current_row_stop  = row_indices[row + 1] - element_start;
      unsigned int thread_base_id  = current_row_start + get_local_id(0);
      for (unsigned int j = thread_base_id + get_local_size(0); j < current_row_stop; j += get_local_size(0))
        shared_elements[thread_base_id] += shared_elements[j];
      for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)
      {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (get_local_id(0) < stride && thread_base_id < current_row_stop)
          shared_elements[thread_base_id] += (thread_base_id + stride < current_row_stop) ? shared_elements[thread_base_id+stride] : 0;
      }
      double row_result = 0; 
      if (current_row_stop > current_row_start)
        row_result = shared_elements[current_row_start]; 
      if (get_local_id(0) == 0)
        result[row * layout_result.y + layout_result.x] = row_result;
    }
  }
  else  
  { 
    shared_elements[get_local_id(0)] = 0; 
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0)) 
      shared_elements[get_local_id(0)] += elements[i] * x[column_indices[i] * layout_x.y + layout_x.x]; 
    for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2) { 
      barrier(CLK_LOCAL_MEM_FENCE); 
      if (get_local_id(0) < stride) 
        shared_elements[get_local_id(0)] += shared_elements[get_local_id(0) + stride]; 
    } 
    if (get_local_id(0) == 0) 
      result[row_start * layout_result.y + layout_result.x] = shared_elements[0]; 
  } 
} 


__kernel void updateRHS( 
  __global const unsigned int * row_indices, 
  __global const unsigned int * column_indices, 
  __global const unsigned int * row_blocks, 
  __global const double * elements, 
  unsigned int num_blocks, 
  __global const double * x, 
  uint4 layout_x, 
  __global double * result, 
  uint4 layout_result) { 
  __local double shared_elements[1024]; 
  unsigned int row_start = row_blocks[get_group_id(0)]; 
  unsigned int row_stop  = row_blocks[get_group_id(0) + 1]; 
  unsigned int rows_to_process = row_stop - row_start; 
  unsigned int element_start = row_indices[row_start]; 
  unsigned int element_stop = row_indices[row_stop]; 
  if (rows_to_process > 4) { 
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0)) 
      shared_elements[i - element_start] = elements[i] * x[column_indices[i] * layout_x.y + layout_x.x]; 
    barrier(CLK_LOCAL_MEM_FENCE); 
    for (unsigned int row = row_start + get_local_id(0); row < row_stop; row += get_local_size(0)) { 
      double dot_prod = 0; 
      unsigned int thread_row_start = row_indices[row]     - element_start; 
      unsigned int thread_row_stop  = row_indices[row + 1] - element_start; 
      for (unsigned int i = thread_row_start; i < thread_row_stop; ++i) 
        dot_prod += shared_elements[i]; 
      result[row * layout_result.y + layout_result.x] -= dot_prod; 
    } 
  } 
  else if (rows_to_process > 1) 
  {
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0))
      shared_elements[i - element_start] = elements[i] * x[column_indices[i] * layout_x.y + layout_x.x];
    barrier(CLK_LOCAL_MEM_FENCE); 
    for (unsigned int row = row_start; row < row_stop; ++row)
    {
      unsigned int current_row_start = row_indices[row]     - element_start;
      unsigned int current_row_stop  = row_indices[row + 1] - element_start;
      unsigned int thread_base_id  = current_row_start + get_local_id(0);
      for (unsigned int j = thread_base_id + get_local_size(0); j < current_row_stop; j += get_local_size(0))
        shared_elements[thread_base_id] += shared_elements[j];
      for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)
      {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (get_local_id(0) < stride && thread_base_id < current_row_stop)
          shared_elements[thread_base_id] += (thread_base_id + stride < current_row_stop) ? shared_elements[thread_base_id+stride] : 0;
      }
      double row_result = 0; 
      if (current_row_stop > current_row_start)
        row_result = shared_elements[current_row_start]; 
      if (get_local_id(0) == 0)
        result[row * layout_result.y + layout_result.x] -= row_result;
    }
  }
  else  
  { 
    shared_elements[get_local_id(0)] = 0; 
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0)) 
      shared_elements[get_local_id(0)] += elements[i] * x[column_indices[i] * layout_x.y + layout_x.x]; 
    for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2) { 
      barrier(CLK_LOCAL_MEM_FENCE); 
      if (get_local_id(0) < stride) 
        shared_elements[get_local_id(0)] += shared_elements[get_local_id(0) + stride]; 
    } 
    if (get_local_id(0) == 0) 
      result[row_start * layout_result.y + layout_result.x] -= shared_elements[0]; 
  } 
} 


__kernel void updateResidual( 
  __global const unsigned int * row_indices, 
  __global const unsigned int * column_indices, 
  __global const unsigned int * row_blocks, 
  __global const double * elements, 
  unsigned int num_blocks, 
  __global const double * x, 
  uint4 layout_x, 
  __global double * result, 
  uint4 layout_result,
  __global const double* rhs) { 
  __local double shared_elements[1024]; 
  unsigned int row_start = row_blocks[get_group_id(0)]; 
  unsigned int row_stop  = row_blocks[get_group_id(0) + 1]; 
  unsigned int rows_to_process = row_stop - row_start; 
  unsigned int element_start = row_indices[row_start]; 
  unsigned int element_stop = row_indices[row_stop]; 
  if (rows_to_process > 4) { 
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0)) 
      shared_elements[i - element_start] = elements[i] * x[column_indices[i] * layout_x.y + layout_x.x]; 
    barrier(CLK_LOCAL_MEM_FENCE); 
    for (unsigned int row = row_start + get_local_id(0); row < row_stop; row += get_local_size(0)) { 
      double dot_prod = 0; 
      unsigned int thread_row_start = row_indices[row]     - element_start; 
      unsigned int thread_row_stop  = row_indices[row + 1] - element_start; 
      for (unsigned int i = thread_row_start; i < thread_row_stop; ++i) 
        dot_prod += shared_elements[i]; 
      result[row * layout_result.y + layout_result.x] = rhs[row * layout_result.y + layout_result.x] - dot_prod; 
    } 
  } 
  else if (rows_to_process > 1) 
  {
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0))
      shared_elements[i - element_start] = elements[i] * x[column_indices[i] * layout_x.y + layout_x.x];
    barrier(CLK_LOCAL_MEM_FENCE); 
    for (unsigned int row = row_start; row < row_stop; ++row)
    {
      unsigned int current_row_start = row_indices[row]     - element_start;
      unsigned int current_row_stop  = row_indices[row + 1] - element_start;
      unsigned int thread_base_id  = current_row_start + get_local_id(0);
      for (unsigned int j = thread_base_id + get_local_size(0); j < current_row_stop; j += get_local_size(0))
        shared_elements[thread_base_id] += shared_elements[j];
      for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)
      {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (get_local_id(0) < stride && thread_base_id < current_row_stop)
          shared_elements[thread_base_id] += (thread_base_id + stride < current_row_stop) ? shared_elements[thread_base_id+stride] : 0;
      }
      double row_result = 0; 
      if (current_row_stop > current_row_start)
        row_result = shared_elements[current_row_start]; 
      if (get_local_id(0) == 0)
        result[row * layout_result.y + layout_result.x] = rhs[row * layout_result.y + layout_result.x] - row_result;
    }
  }
  else  
  { 
    shared_elements[get_local_id(0)] = 0; 
    for (unsigned int i = element_start + get_local_id(0); i < element_stop; i += get_local_size(0)) 
      shared_elements[get_local_id(0)] += elements[i] * x[column_indices[i] * layout_x.y + layout_x.x]; 
    for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2) { 
      barrier(CLK_LOCAL_MEM_FENCE); 
      if (get_local_id(0) < stride) 
        shared_elements[get_local_id(0)] += shared_elements[get_local_id(0) + stride]; 
    } 
    if (get_local_id(0) == 0) 
      result[row_start * layout_result.y + layout_result.x] = rhs[row_start * layout_result.y + layout_result.x] - shared_elements[0]; 
  } 
} 


