#include "SourceGenerator.h"


cl_device_id* sourceGenerator::device = nullptr;
cl_context* sourceGenerator::context = nullptr;
cl_command_queue* sourceGenerator::ioQue = nullptr;
sourceGenerator* sourceGenerator::s_instance = nullptr;


void sourceGenerator::addString2Kernel(const std::string &kerStr)
{
	programSrc.append(kerStr);
	programSrc.append("\n");
}

void sourceGenerator::addString2Defines(const std::string &defstr_)
{
	defineStr.append(defstr_);
	defineStr.append("\n");
}

void sourceGenerator::addFile2Kernel(const std::string &fname_)
{
	std::string fullfname = "source\\Kernels\\" + fname_;
	std::ifstream in_file(fullfname);
	std::string kernel_string((std::istreambuf_iterator<char>(in_file)), std::istreambuf_iterator<char>());
	addString2Kernel(kernel_string);
}

void sourceGenerator::addDefine(std::string &kerstr, const std::string varname, const double value)
{
	char char_temp[100];
	sprintf(char_temp, "(%12.11g)\n", value);
	kerstr.append("#define " + varname + "\t\t" + char_temp);
}

void sourceGenerator::addDefine(std::string &kerstr, const std::string varname, int value)
{
	char char_temp[100];
	sprintf(char_temp, "(%d)\n", value);
	kerstr.append("#define " + varname + "\t\t" + char_temp);
}

void sourceGenerator::addDefine(std::string& kerstr, const std::string varname, unsigned int value)
{
	char char_temp[100];
	sprintf(char_temp, "(%d)\n", value);
	kerstr.append("#define " + varname + "\t\t" + char_temp);
}

void sourceGenerator::addDefine(std::string &kerstr, const std::string varname, std::string value)
{
	kerstr.append("#define " + varname + "\t\t" + value + "\n");
}

void sourceGenerator::addDefine(std::string &kerstr, const std::string varname)
{
	kerstr.append("#define " + varname + "\n");
}


sourceGenerator::hashResult sourceGenerator::checkHash(std::string &kerstr, std::string &kername)
{
	auto hash = rsHash(kerstr);

	auto kernel_iterator = kernel_map.find(hash);
	if (kernel_iterator != kernel_map.end())
	{
		kername = kernel_iterator->second;
		return oldKernel;
	}
	else //build program and compile the kernel;
	{
		kernel_map[hash] = kername;
		programSrc.append(kerstr + "\n\n");
		return newKernel;
	}
}


void sourceGenerator::callIniKernels()
{
	std::list<std::function<void(void)>>::iterator createKernelsIter = createKernelList.begin();

	while (createKernelsIter != createKernelList.end())
	{
		std::function<void(void)> createFunc_ = *createKernelsIter;

		if (createFunc_ != NULL)
		{
			createFunc_();
			++createKernelsIter;
		}
		else
		{
			// BTW, who is deleting pItem, a.k.a. (*iter)?
			createKernelsIter = createKernelList.erase(createKernelsIter);
		}
	}

	std::list<std::function<void(void)>>::iterator setKernelArgsIter = setArgumentList.begin();

	while (setKernelArgsIter != setArgumentList.end())
	{
		std::function<void(void)> setArgFunc_ = *setKernelArgsIter;

		if (setArgFunc_ != NULL)
		{
			setArgFunc_();
			++setKernelArgsIter;
		}
		else
		{
			// BTW, who is deleting pItem, a.k.a. (*iter)?
			setKernelArgsIter = setArgumentList.erase(setKernelArgsIter);
		}
	}
}


void sourceGenerator::findAndReplace(std::string & data, std::string toSearch, std::string replaceStr)
{
	// Get the first occurrence
	size_t pos = data.find(toSearch);
	data.replace(pos, toSearch.size(), replaceStr);
}

void sourceGenerator::Gen_Error_Msg(int Errcode, std::string mesg)
{
#ifdef LOG_ERROR_IN_FILE
	FILE *stream;
	stream = fopen("error.txt", "w+");
	printf("Error in %s class: Code: %d\nnMessage: %s\n", type.c_str(), Errcode, mesg.c_str());
	fclose(stream);
#endif
	printf("Error in %s class: Code: %d\nnMessage: %s\n", type.c_str(), Errcode, mesg.c_str());
	exit(Errcode);
}


unsigned int sourceGenerator::rsHash(const std::string &key)
{
	unsigned int b = 378551;
	unsigned int a = 63689;
	unsigned int hash = 0;
	unsigned int i = 0;

	auto len = key.size();

	for (i = 0; i < len; i++)
	{
		hash = hash * a + key[i];
		a = a * b;
	}

	return hash;
}