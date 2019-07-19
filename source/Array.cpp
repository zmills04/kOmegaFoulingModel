// Array Classes which hold both Host and Device memory and provide
// convenience functions for transfering memory between host and device
// (c) Zachary Grant Mills, 2016 
//////////////////////////////////////////////////////////////////////

#include "Array.h"

cl_context* staticBaseVar::context = nullptr;

std::string staticBaseVar::curDumpSaveDir = "";





//template <typename T>
//cl_command_queue* ArrayBase<T>::ioQue = {};
//
//template <typename T>
//cl_context* ArrayBase<T>::context = {};
//

