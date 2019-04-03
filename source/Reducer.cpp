#include "Reducer.h"

template<typename T>
ReduceGenerator::varType ReducerBase<T>::getVarType()
{
	return ReduceGenerator::NULL_T;
}


template<>
ReduceGenerator::varType ReducerBase<double>::getVarType()
{
	return ReduceGenerator::DOUBLE_T;
}

template<>
ReduceGenerator::varType ReducerBase<int>::getVarType()
{
	return ReduceGenerator::INT_T;
}

template<>
ReduceGenerator::varType ReducerBase<char>::getVarType()
{
	return ReduceGenerator::CHAR_T;
}

template<>
ReduceGenerator::varType ReducerBase<cl_uint>::getVarType()
{
	return ReduceGenerator::UINT_T;
}

template<>
ReduceGenerator::varType ReducerBase<bool>::getVarType()
{
	return ReduceGenerator::BOOL_T;
}














