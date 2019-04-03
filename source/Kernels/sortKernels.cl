uint lowerBoundBinarylocal(local par* data, uint left, uint right, par searchVal)
{
	uint firstIndex = left;
	uint lastIndex = right;

	while (firstIndex < lastIndex)
	{
		uint midIndex = (firstIndex + lastIndex) / 2;
		par midValue = data[midIndex];

		if (midValue.loc < searchVal.loc)
		{
			firstIndex = midIndex + 1;
		}
		else
		{
			lastIndex = midIndex;
		}
	}
	return firstIndex;
}

uint upperBoundBinarylocal(local par* data, uint left, uint right, par searchVal)
{
	uint upperBound = lowerBoundBinarylocal(data, left, right, searchVal);

	if (upperBound != right)
	{
		uint mid = 0;
		par upperValue = data[upperBound];
		while ((upperValue.loc == searchVal.loc) && (upperBound < right))
		{
			mid = (upperBound + right) / 2;
			par midValue = data[mid];
			if (midValue.loc == searchVal.loc)
			{
				upperBound = mid + 1;
			}
			else
			{
				right = mid;
				upperBound++;
			}
			upperValue = data[upperBound];
		}
	}
	return upperBound;
}

uint lowerBoundLinear(global par* data, uint left, uint right, par searchVal)
{
	uint firstIndex = left;
	uint lastIndex = right;

	while (firstIndex < lastIndex)
	{
		par dataVal = data[firstIndex];

		if (dataVal.loc < searchVal.loc)
		{
			firstIndex = firstIndex + 1;
		}
		else
		{
			break;
		}
	}

	return firstIndex;
}

uint lowerBoundBinary(__global par *source_ptr, uint left, uint right, par searchVal)
{
	uint firstIndex = left;
	uint lastIndex = right;

	while (firstIndex < lastIndex)
	{
		uint midIndex = (firstIndex + lastIndex) / 2;
		par midValue = source_ptr[midIndex];

		if (midValue.loc < searchVal.loc)
		{
			firstIndex = midIndex + 1;
		}
		else
		{
			lastIndex = midIndex;
		}
	}
	return firstIndex;
}

uint upperBoundBinary(__global par *source_ptr, uint left, uint right, par searchVal)
{
	uint upperBound = lowerBoundBinary(source_ptr, left, right, searchVal);

	if (upperBound != right)
	{
		uint mid = 0;
		par upperValue = source_ptr[upperBound];
		while ((searchVal.loc == upperValue.loc) && (upperBound < right))
		{
			mid = (upperBound + right) / 2;
			par midValue = source_ptr[mid];
			if (midValue.loc == searchVal.loc)
			{
				upperBound = mid + 1;
			}
			else
			{
				right = mid;
				upperBound++;
			}
			upperValue = source_ptr[upperBound];
		}
	}
	return upperBound;
}


kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
	void Sort_merge_global(global par* source_ptr,
	global par* result_ptr,
	const uint srcVecSize,
	const uint srcLogicalBlockSize)
{
	size_t globalID = get_global_id(0);// * get_global_size(0) + get_global_id(0);
	size_t groupID = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
	size_t localID = get_local_id(0);// +get_local_size(0) + get_local_id(0);
	size_t wgSize = get_local_size(0);// * get_local_size(1);


	if (globalID >= srcVecSize)
		return;
	uint srcBlockNum = globalID / srcLogicalBlockSize;
	uint srcBlockIndex = globalID % srcLogicalBlockSize;


	uint dstLogicalBlockSize = srcLogicalBlockSize << 1;
	uint leftBlockIndex = globalID & ~(dstLogicalBlockSize - 1);

	leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
	leftBlockIndex = (((leftBlockIndex) < (srcVecSize)) ? (leftBlockIndex) : (srcVecSize));
	uint rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (srcVecSize)) ? (leftBlockIndex + srcLogicalBlockSize) : (srcVecSize));

	uint insertionIndex = 0;
	par search_val = source_ptr[globalID];
	if ((srcBlockNum & 0x1) == 0)
	{
		insertionIndex = lowerBoundBinary(source_ptr, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
	}
	else
	{
		insertionIndex = upperBoundBinary(source_ptr, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
	}

	uint dstBlockIndex = srcBlockIndex + insertionIndex;
	uint dstBlockNum = srcBlockNum / 2;

	result_ptr[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = source_ptr[globalID];


}

kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
	void Sort_merge_local(global par* data_ptr,
	const uint vecSize,
	local par* lds,
	local par* lds2)
{
	size_t gloId = get_global_id(0);// *get_global_size(0) + get_global_id(0);
	size_t groId = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
	size_t locId = get_local_id(0);// +get_local_size(0) + get_local_id(0);
	size_t wgSize = get_local_size(0);// *get_local_size(1);


	if (gloId < vecSize)
	{
		lds[locId] = data_ptr[gloId];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	uint end = wgSize;
	if ((groId + 1)*(wgSize) >= vecSize)
		end = vecSize - (groId*wgSize);

	uint numMerges = SORT_NUM_MERGES;
	uint pass;
	for (pass = 1; pass <= numMerges; ++pass)
	{
		uint srcLogicalBlockSize = 1 << (pass - 1);
		if (gloId < vecSize)
		{
			uint srcBlockNum = (locId) / srcLogicalBlockSize;
			uint srcBlockIndex = (locId) % srcLogicalBlockSize;

			uint dstLogicalBlockSize = srcLogicalBlockSize << 1;
			uint leftBlockIndex = (locId)& ~(dstLogicalBlockSize - 1);

			leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
			leftBlockIndex = (((leftBlockIndex) < (end)) ? (leftBlockIndex) : (end));
			uint rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (end)) ? (leftBlockIndex + srcLogicalBlockSize) : (end));

			uint insertionIndex = 0;
			if (pass % 2 != 0)
			{
				if ((srcBlockNum & 0x1) == 0)
				{
					insertionIndex = lowerBoundBinarylocal(lds, leftBlockIndex, rightBlockIndex, lds[locId]) - leftBlockIndex;
				}
				else
				{
					insertionIndex = upperBoundBinarylocal(lds, leftBlockIndex, rightBlockIndex, lds[locId]) - leftBlockIndex;
				}
			}
			else
			{
				if ((srcBlockNum & 0x1) == 0)
				{
					insertionIndex = lowerBoundBinarylocal(lds2, leftBlockIndex, rightBlockIndex, lds2[locId]) - leftBlockIndex;
				}
				else
				{
					insertionIndex = upperBoundBinarylocal(lds2, leftBlockIndex, rightBlockIndex, lds2[locId]) - leftBlockIndex;
				}
			}
			uint dstBlockIndex = srcBlockIndex + insertionIndex;
			uint dstBlockNum = srcBlockNum / 2;
			if (pass % 2 != 0)
				lds2[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = lds[locId];
			else
				lds[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = lds2[locId];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (gloId < vecSize)
	{
		data_ptr[gloId] = lds[locId];
	}
}

__kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_PLOC, 1, 1)))
	void Sort_update_loc(__global par *P,
	__global int *Ploc_array)
{
	int i = get_global_id(0);
	if (i >= TRC_NUM_TRACERS)
		return;
	int ploc = P[i].loc;
	
	if (i == TRC_NUM_TRACERS - 1)
	{
		Ploc_array[2*(ploc + 2)+1] = i+1;
		return;
	}
	
	int next_loc = P[i + 1].loc;

	if(i == 0)
	{
		Ploc_array[2*(ploc + 2)] = 0;
	}
	
	if (ploc != next_loc)
	{
		Ploc_array[2*(ploc + 2)+1] = i + 1;
		Ploc_array[2*(next_loc + 2)] = i + 1;
		for (int k = ploc + 1; k < next_loc; k++)
		{
			Ploc_array[2*(k + 2)] = -1;
		}
	}
}


//uint lowerBoundBinarylocal_BL(local par* data, uint left, uint right, par searchVal)
//{
//	uint firstIndex = left;
//	uint lastIndex = right;
//
//	while (firstIndex < lastIndex)
//	{
//		uint midIndex = (firstIndex + lastIndex) / 2;
//		par midValue = data[midIndex];
//
//		if (midValue.loc < searchVal.loc)
//		{
//			firstIndex = midIndex + 1;
//		}
//		else
//		{
//			lastIndex = midIndex;
//		}
//	}
//	return firstIndex;
//}
//
//uint upperBoundBinarylocal_BL(local par* data, uint left, uint right, par searchVal)
//{
//	uint upperBound = lowerBoundBinarylocal_BL(data, left, right, searchVal);
//
//	if (upperBound != right)
//	{
//		uint mid = 0;
//		par upperValue = data[upperBound];
//		while ((upperValue.loc == searchVal.loc) && (upperBound < right))
//		{
//			mid = (upperBound + right) / 2;
//			par midValue = data[mid];
//			if (midValue.loc == searchVal.loc)
//			{
//				upperBound = mid + 1;
//			}
//			else
//			{
//				right = mid;
//				upperBound++;
//			}
//			upperValue = data[upperBound];
//		}
//	}
//	return upperBound;
//}
//
//uint lowerBoundLinear_BL(global par* data, uint left, uint right, par searchVal)
//{
//	uint firstIndex = left;
//	uint lastIndex = right;
//
//	while (firstIndex < lastIndex)
//	{
//		par dataVal = data[firstIndex];
//
//		if (dataVal.loc < searchVal.loc)
//		{
//			firstIndex = firstIndex + 1;
//		}
//		else
//		{
//			break;
//		}
//	}
//
//	return firstIndex;
//}
//
//uint lowerBoundBinary_BL(__global par *source_ptr, uint left, uint right, par searchVal)
//{
//	uint firstIndex = left;
//	uint lastIndex = right;
//
//	while (firstIndex < lastIndex)
//	{
//		uint midIndex = (firstIndex + lastIndex) / 2;
//		par midValue = source_ptr[midIndex];
//
//		if (midValue.loc < searchVal.loc)
//		{
//			firstIndex = midIndex + 1;
//		}
//		else
//		{
//			lastIndex = midIndex;
//		}
//	}
//	return firstIndex;
//}
//
//uint upperBoundBinary_BL(__global par *source_ptr, uint left, uint right, par searchVal)
//{
//	uint upperBound = lowerBoundBinary_BL(source_ptr, left, right, searchVal);
//
//	if (upperBound != right)
//	{
//		uint mid = 0;
//		par upperValue = source_ptr[upperBound];
//		while ((searchVal.loc == upperValue.loc) && (upperBound < right))
//		{
//			mid = (upperBound + right) / 2;
//			par midValue = source_ptr[mid];
//			if (midValue.loc == searchVal.loc)
//			{
//				upperBound = mid + 1;
//			}
//			else
//			{
//				right = mid;
//				upperBound++;
//			}
//			upperValue = source_ptr[upperBound];
//		}
//	}
//	return upperBound;
//}
//
//
//kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
//void Sort_merge_global_BL(global par* source_ptr,
//global par* result_ptr,
//const uint srcVecSize,
//const uint srcLogicalBlockSize)
//{
//	size_t globalID = get_global_id(0);// * get_global_size(0) + get_global_id(0);
//	size_t groupID = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
//	size_t localID = get_local_id(0);// +get_local_size(0) + get_local_id(0);
//	size_t wgSize = get_local_size(0);// * get_local_size(1);
//
//
//	if (globalID >= srcVecSize)
//		return;
//	uint srcBlockNum = globalID / srcLogicalBlockSize;
//	uint srcBlockIndex = globalID % srcLogicalBlockSize;
//
//
//	uint dstLogicalBlockSize = srcLogicalBlockSize << 1;
//	uint leftBlockIndex = globalID & ~(dstLogicalBlockSize - 1);
//
//	leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
//	leftBlockIndex = (((leftBlockIndex) < (srcVecSize)) ? (leftBlockIndex) : (srcVecSize));
//	uint rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (srcVecSize)) ? (leftBlockIndex + srcLogicalBlockSize) : (srcVecSize));
//
//	uint insertionIndex = 0;
//	par search_val = source_ptr[globalID];
//	if ((srcBlockNum & 0x1) == 0)
//	{
//		insertionIndex = lowerBoundBinary_BL(source_ptr, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
//	}
//	else
//	{
//		insertionIndex = upperBoundBinary_BL(source_ptr, leftBlockIndex, rightBlockIndex, search_val) - leftBlockIndex;
//	}
//
//	uint dstBlockIndex = srcBlockIndex + insertionIndex;
//	uint dstBlockNum = srcBlockNum / 2;
//
//	result_ptr[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = source_ptr[globalID];
//
//
//}
//
//kernel __attribute__((reqd_work_group_size(WORKGROUPSIZE_SORT, 1, 1)))
//void Sort_merge_local_BL(global par* data_ptr,
//const uint vecSize,
//local par* lds,
//local par* lds2)
//{
//	size_t gloId = get_global_id(0);// *get_global_size(0) + get_global_id(0);
//	size_t groId = get_group_id(0);// *get_num_groups(0) + get_group_id(0);
//	size_t locId = get_local_id(0);// +get_local_size(0) + get_local_id(0);
//	size_t wgSize = get_local_size(0);// *get_local_size(1);
//
//
//	if (gloId < vecSize)
//	{
//		lds[locId] = data_ptr[gloId];
//	}
//	barrier(CLK_LOCAL_MEM_FENCE);
//	uint end = wgSize;
//	if ((groId + 1)*(wgSize) >= vecSize)
//		end = vecSize - (groId*wgSize);
//
//	uint numMerges = SORT_NUM_MERGES;
//	uint pass;
//	for (pass = 1; pass <= numMerges; ++pass)
//	{
//		uint srcLogicalBlockSize = 1 << (pass - 1);
//		if (gloId < vecSize)
//		{
//			uint srcBlockNum = (locId) / srcLogicalBlockSize;
//			uint srcBlockIndex = (locId) % srcLogicalBlockSize;
//
//			uint dstLogicalBlockSize = srcLogicalBlockSize << 1;
//			uint leftBlockIndex = (locId)& ~(dstLogicalBlockSize - 1);
//
//			leftBlockIndex += (srcBlockNum & 0x1) ? 0 : srcLogicalBlockSize;
//			leftBlockIndex = (((leftBlockIndex) < (end)) ? (leftBlockIndex) : (end));
//			uint rightBlockIndex = (((leftBlockIndex + srcLogicalBlockSize) < (end)) ? (leftBlockIndex + srcLogicalBlockSize) : (end));
//
//			uint insertionIndex = 0;
//			if (pass % 2 != 0)
//			{
//				if ((srcBlockNum & 0x1) == 0)
//				{
//					insertionIndex = lowerBoundBinarylocal_BL(lds, leftBlockIndex, rightBlockIndex, lds[locId]) - leftBlockIndex;
//				}
//				else
//				{
//					insertionIndex = upperBoundBinarylocal_BL(lds, leftBlockIndex, rightBlockIndex, lds[locId]) - leftBlockIndex;
//				}
//			}
//			else
//			{
//				if ((srcBlockNum & 0x1) == 0)
//				{
//					insertionIndex = lowerBoundBinarylocal_BL(lds2, leftBlockIndex, rightBlockIndex, lds2[locId]) - leftBlockIndex;
//				}
//				else
//				{
//					insertionIndex = upperBoundBinarylocal_BL(lds2, leftBlockIndex, rightBlockIndex, lds2[locId]) - leftBlockIndex;
//				}
//			}
//			uint dstBlockIndex = srcBlockIndex + insertionIndex;
//			uint dstBlockNum = srcBlockNum / 2;
//			if (pass % 2 != 0)
//				lds2[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = lds[locId];
//			else
//				lds[(dstBlockNum*dstLogicalBlockSize) + dstBlockIndex] = lds2[locId];
//		}
//		barrier(CLK_LOCAL_MEM_FENCE);
//	}
//	if (gloId < vecSize)
//	{
//		data_ptr[gloId] = lds[locId];
//	}
//}














