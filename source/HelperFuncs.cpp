#include "HelperFuncs.h"

template <class T>
void functionPointerWrapper(void* pt2Object, FuncPtrType fptype_)
{
	T* mySelf = (T*)pt2Object;
	if (fptype_ == ptr2CreateKernels)
		mySelf->createKernels();
	else if (fptype_ == ptr2SetKernelArgs)
		mySelf->setKernelArgs();
	else
		mySelf->loadParams();
}


void delay_func(double sec)
{
	time_t StartTime, CurrTime;
	time(&StartTime);
	do
	{
		time(&CurrTime);
	} while (difftime(CurrTime, StartTime) < sec);
}


bool fileopen(std::fstream &stream, std::string &FileName, FileMode mode_, const int att)
{
	for (int i = 0; i < att; i++)
	{
		switch (mode_)
		{
		case FileIn:
		{
			stream.open(FileName, std::ios_base::in);
			break;
		}
		case FileOut:
		{
			stream.open(FileName, std::ios_base::out);
			break;
		}
		case FileAppend:
		{
			stream.open(FileName, std::ios_base::out | std::ios_base::app);
			break;
		}
		case BinaryIn:
		{
			stream.open(FileName, std::ios_base::in | std::ios_base::binary);
			break;
		}
		case BinaryOut:
		{
			stream.open(FileName, std::ios_base::out | std::ios_base::binary);
			break;
		}
		}

		if (stream.is_open())
			return true;
		delay_func(0.1);
	}
	return false;
}

bool fileopen(std::fstream &stream, std::string &FileName, FileMode Mode)
{
	return fileopen(stream, FileName, Mode, 100);
}

void MakeDir(std::string& NewDir)
{
#if defined(_WIN32)
#pragma warning(suppress: 6031)
	_mkdir(NewDir.c_str());
#else 
	mkdir(NewDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif //_WIN32
}

void CopyFile(std::string NameSrc, std::string NameDest)
{
	
	std::fstream  hSrc, hDst;
	fileopen(hSrc, NameSrc, BinaryIn);
	fileopen(hDst, NameDest, BinaryOut);
	hDst << hSrc.rdbuf();
	hSrc.close();
	hDst.close();
}

void RenameFile(std::string SourceFileName, std::string DestFileName)
{
	int maxAttempt = 100;

	for (int i = 0; i < maxAttempt; i++)
	{
		if (std::rename(SourceFileName.c_str(), DestFileName.c_str()))
		{
			return;
		}
		delay_func(0.1);
	}
}

double boxMuller(double meanval, double stdval)
{
	static double y2;
	static int use_last = 0;
	double x1, x2, Wval, y1;
	if (use_last) {
		y1 = y2;
		use_last = 0;
	}
	else {
		do {
			x1 = 2.*(double)rand() / ((double)RAND_MAX) - 1.;
			x2 = 2.*(double)rand() / ((double)RAND_MAX) - 1.;
			Wval = x1*x1 + x2*x2;
		} while (Wval >= 1.);
		Wval = sqrt((-2. * log(Wval)) / Wval);
		y1 = x1*Wval; y2 = x2*Wval;
	}
	return (meanval + y1*stdval);
}


void operator >> (const YAML::Node& node, cl_uint4 vec_)
{
	node[0] >> vec_.x;
	node[1] >> vec_.y;
	node[2] >> vec_.z;
	node[3] >> vec_.w;
}

YAML::Emitter& operator << (YAML::Emitter& out, const cl_uint4 vec_)
{
	out << YAML::Flow;
	out << YAML::BeginSeq << vec_.x << vec_.y << vec_.z << vec_.w;
	out << YAML::EndSeq;
	return out;
}



