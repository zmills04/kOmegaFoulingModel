#include "Reducer.h"

template<typename T, typename T2>
ReduceGenerator::varType ReducerBase<T, T2>::getVarType()
{
	return ReduceGenerator::NULL_T;
}


template<>
ReduceGenerator::varType ReducerBase<double,double>::getVarType()
{
	return ReduceGenerator::DOUBLE_T;
}

template<>
ReduceGenerator::varType ReducerBase<int,int>::getVarType()
{
	return ReduceGenerator::INT_T;
}

template<>
ReduceGenerator::varType ReducerBase<char,char>::getVarType()
{
	return ReduceGenerator::CHAR_T;
}

template<>
ReduceGenerator::varType ReducerBase<cl_uint,cl_uint>::getVarType()
{
	return ReduceGenerator::UINT_T;
}

template<>
ReduceGenerator::varType ReducerBase<bool,bool>::getVarType()
{
	return ReduceGenerator::BOOL_T;
}

template<>
ReduceGenerator::varType ReducerBase<int, cl_short>::getVarType()
{
	return ReduceGenerator::INT_T;
}















