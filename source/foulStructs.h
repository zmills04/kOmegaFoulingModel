#pragma once


#include "Array.h"
#include "particleStructs.h"


typedef struct legacyRampInfo
{
	double Ybegin;
	double Coeff;
	cl_uint IOind; //Index of LS nodes at beginning and end of TR areatypename
	cl_uint Cind; //direction traversed from IO_end when creating ramp
} legacyRampInfo;

class rampI : public trStructBase
{
public:

	Array1Dd Ybegin, Coeff;
	Array1Du IOind, Cind;



	rampI(std::string name_ = "") : trStructBase(name_)
	{}

	~rampI() {}

void setStruct(legacyRampInfo& struct_, const int i)
{
	Ybegin(i) = struct_.Ybegin;
	Coeff(i) = struct_.Coeff;
	IOind(i) = struct_.IOind;
	Cind(i) = struct_.Cind;
}

legacyRampInfo getStruct(const int i)
{
	legacyRampInfo struct_;
	struct_.Ybegin = Ybegin(i);
	struct_.Coeff = Coeff(i);
	struct_.IOind = IOind(i);
	struct_.Cind = Cind(i);
	return struct_;
}


	void allocateArrays()
	{
		Ybegin.setName("RampI_Ybegin");
		Coeff.setName("RampI_Coeff");
		IOind.setName("RampI_IOind");
		Cind.setName("RampI_Cind");


		Ybegin.allocate(fullSize);
		Coeff.allocate(fullSize);
		IOind.allocate(fullSize);
		Cind.allocate(fullSize);
	}

	void allocateBuffers(int bufSize = -1)
	{
		if (bufSize == -1)
			bufSize = fullSize;

		bufFullSize = bufSize;

		Ybegin.allocate_buffer_size(bufFullSize);
		Coeff.allocate_buffer_size(bufFullSize);
		IOind.allocate_buffer_size(bufFullSize);
		Cind.allocate_buffer_size(bufFullSize);

	}

	void copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
	{
		if (writeSize == -1)
			writeSize = fullSize;
		Ybegin.copy_to_buffer_size(writeSize, que_, bFlag_);
		Coeff.copy_to_buffer_size(writeSize, que_, bFlag_);
		IOind.copy_to_buffer_size(writeSize, que_, bFlag_);
		Cind.copy_to_buffer_size(writeSize, que_, bFlag_);

	}

	void copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
	{
		if (readSize == -1)
			readSize = fullSize;

		Ybegin.read_from_buffer_size(readSize, que_, bFlag_);
		Coeff.read_from_buffer_size(readSize, que_, bFlag_);
		IOind.read_from_buffer_size(readSize, que_, bFlag_);
		Cind.read_from_buffer_size(readSize, que_, bFlag_);
	}

	//Ybegin
	//	Coeff
	//	IOind
	//	Cind

	void ini()
	{
		//initialization functions specific to class
	}

	bool load()
	{
		bool ret = true;
		ret &= Ybegin.load("load" SLASH + Ybegin.getName());
		ret &= Coeff.load("load" SLASH + Coeff.getName());
		ret &= IOind.load("load" SLASH + IOind.getName());
		ret &= Cind.load("load" SLASH + Cind.getName());
		return ret;
	}

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
	{
		bool ret = true;
		if(!saveAll)
		{
			std::cout << "Warning: to save RampI, set saveAll to true when calling save\n";
		}

		if (saveAll)
		{
			switch (saveFlag_)
			{
			case saveFl:
			{
				ret &= Ybegin.save();
				ret &= Coeff.save();
				ret &= IOind.save();
				ret &= Cind.save();
				break;
			}
			case saveTxtFl:
			{
				ret &= Ybegin.save2file();
				ret &= Coeff.save2file();
				ret &= IOind.save2file();
				ret &= Cind.save2file();
				break;
			}
			case saveBinFl:
			{
				ret &= Ybegin.savebin();
				ret &= Coeff.savebin();
				ret &= IOind.savebin();
				ret &= Cind.savebin();
				break;
			}
			}
		}
		return ret;
	}

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr)
	{
		bool ret = true;
		if(!saveAll)
		{
			std::cout << "Warning: to save foulI from device, set saveAll to true when calling saveFromDevice\n";
		}

		if (saveAll)
		{
			switch (saveFlag_)
			{
			case saveFl:
			{
				ret &= Ybegin.save_from_device();
				ret &= Coeff.save_from_device();
				ret &= IOind.save_from_device();
				ret &= Cind.save_from_device();
				break;
			}
			case saveTxtFl:
			{
				ret &= Ybegin.save_txt_from_device();
				ret &= Coeff.save_txt_from_device();
				ret &= IOind.save_txt_from_device();
				ret &= Cind.save_txt_from_device();
				break;
			}
			case saveBinFl:
			{
				ret &= Ybegin.save_bin_from_device();
				ret &= Coeff.save_bin_from_device();
				ret &= IOind.save_bin_from_device();
				ret &= Cind.save_bin_from_device();
				break;
			}
			}
		}
		return ret;
	}
};



typedef struct legacyFoulInfo
{
	cl_double4 WeightsL;	//Weights applied to BL deposits to left
	cl_double4 WeightsR;	//Weights applied to BL deposits to right
	cl_double2 vN;			//Normal vector of C node
	cl_double disp;			//Current displacement distance
	cl_uint BL_ind;			//index of BL to left of node
	cl_uint C_ind;			//index of wall node this struct corresponds to
} legacyFoulInfo;


class foulI : public trStructBase
{
public:

	//weightsL
	//weightsR
	//vN
	//disp
	//blInd
	//cInd

	Array1Dd weightsL, weightsR;
	Array1Dv2d vN;
	Array1Dd disp;
	Array1Du blInd, cInd;

	foulI(std::string name_ = "") : trStructBase(name_)
	{}

	~foulI() {}

	void setStruct(legacyFoulInfo& struct_, const int i)
	{
		weightsL(4 * i) = struct_.WeightsL.x;
		weightsL(4 * i+1) = struct_.WeightsL.y;
		weightsL(4 * i+2) = struct_.WeightsL.z;
		weightsL(4 * i+3) = struct_.WeightsL.w;
		weightsR(4 * i) = struct_.WeightsR.x;
		weightsR(4 * i + 1) = struct_.WeightsR.y;
		weightsR(4 * i + 2) = struct_.WeightsR.z;
		weightsR(4 * i + 3) = struct_.WeightsR.w;
		vN(i) = struct_.vN;
		disp(i) = struct_.disp;
		blInd(i) = struct_.BL_ind;
		cInd(i) = struct_.C_ind;
	}

	legacyFoulInfo getStruct(const int i)
	{
		legacyFoulInfo struct_;
		struct_.WeightsL = { { weightsL(4 * i), weightsL(4 * i + 1), weightsL(4 * i + 2), weightsL(4 * i + 3)} };
		struct_.WeightsR = { { weightsR(4 * i), weightsR(4 * i + 1), weightsR(4 * i + 2), weightsR(4 * i + 3)} };
		struct_.vN = vN(i);
		struct_.disp = disp(i);
		struct_.BL_ind = blInd(i);
		struct_.C_ind = cInd(i);

		return struct_;
	}


	void allocateArrays()
	{
		weightsL.setName("foulI_weightsL");
		weightsR.setName("foulI_weightsR");
		vN.setName("foulI_vN");
		disp.setName("foulI_disp");
		blInd.setName("foulI_blInd");
		cInd.setName("foulI_cInd");

		weightsL.allocate(fullSize*4);
		weightsR.allocate(fullSize*4);
		vN.allocate(fullSize);
		disp.allocate(fullSize);
		blInd.allocate(fullSize);
		cInd.allocate(fullSize);

	}

	void allocateBuffers(int bufSize = -1)
	{
		if (bufSize == -1)
			bufSize = fullSize;

		bufFullSize = bufSize;

		weightsL.allocate_buffer_size(bufFullSize * 4);
		weightsR.allocate_buffer_size(bufFullSize * 4);
		vN.allocate_buffer_size(bufFullSize);
		disp.allocate_buffer_size(bufFullSize);
		blInd.allocate_buffer_size(bufFullSize);
		cInd.allocate_buffer_size(bufFullSize);
	}

	void copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
	{
		if (writeSize == -1)
			writeSize = fullSize;
		weightsL.copy_to_buffer_size(writeSize*4, que_, bFlag_);
		weightsR.copy_to_buffer_size(writeSize*4, que_, bFlag_);
		vN.copy_to_buffer_size(writeSize, que_, bFlag_);
		disp.copy_to_buffer_size(writeSize, que_, bFlag_);
		blInd.copy_to_buffer_size(writeSize, que_, bFlag_);
		cInd.copy_to_buffer_size(writeSize, que_, bFlag_);
	}

	void copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
	{
		if (readSize == -1)
			readSize = fullSize;
		weightsL.read_from_buffer_size(readSize*4, que_, bFlag_);
		weightsR.read_from_buffer_size(readSize*4, que_, bFlag_);
		vN.read_from_buffer_size(readSize, que_, bFlag_);
		disp.read_from_buffer_size(readSize, que_, bFlag_);
		blInd.read_from_buffer_size(readSize, que_, bFlag_);
		cInd.read_from_buffer_size(readSize, que_, bFlag_);
	}

	void ini()
	{
		//initialization functions specific to class
	}

	cl_double4 getWeightR(const int ind)
	{
		cl_double4 ret = { { weightsR(ind * 4), 
			weightsR(ind * 4 + 1), weightsR(ind * 4 + 2),
			weightsR(ind * 4 + 3) } };

		return ret;
	}

	cl_double4 getWeightL(const int ind)
	{
		cl_double4 ret = { { weightsL(ind * 4),
			weightsL(ind * 4 + 1), weightsL(ind * 4 + 2),
			weightsL(ind * 4 + 3) } };

		return ret;
	}

	void setWeightR(const int ind, cl_double4 &weight_)
	{
		weightsR(ind * 4) = weight_.x;
		weightsR(ind * 4 + 1) = weight_.y;
		weightsR(ind * 4 + 2) = weight_.z;
		weightsR(ind * 4 + 3) = weight_.w;
	}

	void setWeightL(const int ind, cl_double4& weight_)
	{
		weightsL(ind * 4) = weight_.x;
		weightsL(ind * 4 + 1) = weight_.y;
		weightsL(ind * 4 + 2) = weight_.z;
		weightsL(ind * 4 + 3) = weight_.w;
	}

	bool load()
	{
		bool ret = true;
		ret &= weightsL.load("load" SLASH + weightsL.getName());
		ret &= weightsR.load("load" SLASH + weightsR.getName());
		ret &= vN.load("load" SLASH + vN.getName());
		ret &= disp.load("load" SLASH + disp.getName());
		ret &= blInd.load("load" SLASH + blInd.getName());
		ret &= cInd.load("load" SLASH + cInd.getName());
		return ret;
	}

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
	{
		bool ret = true;
		if(!saveAll)
		{
			std::cout << "Warning: to save foulI, set saveAll to true when calling save\n";
		}

		if (saveAll)
		{
			switch (saveFlag_)
			{
			case saveFl:
			{
				ret &= weightsL.save();
				ret &= weightsR.save();
				ret &= vN.save();
				ret &= disp.save();
				ret &= blInd.save();
				ret &= cInd.save();
				break;
			}
			case saveTxtFl:
			{
				ret &= weightsL.save2file();
				ret &= weightsR.save2file();
				ret &= vN.save2file();
				ret &= disp.save2file();
				ret &= blInd.save2file();
				ret &= cInd.save2file();
				break;
			}
			case saveBinFl:
			{
				ret &= weightsL.savebin();
				ret &= weightsR.savebin();
				ret &= vN.savebin();
				ret &= disp.savebin();
				ret &= blInd.savebin();
				ret &= cInd.savebin();
				break;
			}
			}
		}
		return ret;
	}

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr)
	{
		bool ret = true;
		if(!saveAll)
		{
			std::cout << "Warning: to save foulI from device, set saveAll to true when calling saveFromDevice\n";
		}

		if (saveAll)
		{
			switch (saveFlag_)
			{
			case saveFl:
			{
				ret &= weightsL.save_from_device();
				ret &= weightsR.save_from_device();
				ret &= vN.save_from_device();
				ret &= disp.save_from_device();
				ret &= blInd.save_from_device();
				ret &= cInd.save_from_device();
				break;
			}
			case saveTxtFl:
			{
				ret &= weightsL.save_txt_from_device();
				ret &= weightsR.save_txt_from_device();
				ret &= vN.save_txt_from_device();
				ret &= disp.save_txt_from_device();
				ret &= blInd.save_txt_from_device();
				ret &= cInd.save_txt_from_device();
				break;
			}
			case saveBinFl:
			{
				ret &= weightsL.save_bin_from_device();
				ret &= weightsR.save_bin_from_device();
				ret &= vN.save_bin_from_device();
				ret &= disp.save_bin_from_device();
				ret &= blInd.save_bin_from_device();
				ret &= cInd.save_bin_from_device();
				break;
			}
			}
		}
		return ret;
	}
};


//class  : public trStructBase
//{
//public:
//
//	Par(std::string name_ = "") {}
//
//	~Par() {}
//
//void setStruct(legacyPParam& struct_, const int i)
//{
//}
//
//legacyPParam getStruct(const int i)
//{
//	legacyPParam struct_;
//	
//	return struct_;
//}

//
//	void allocateArrays()
//	{
//		.setName("Par_Pos");
//
//		.allocate(fullSize);
//
//	}
//
//	void allocateBuffers(int bufSize = -1)
//	{
//		if (bufSize == -1)
//			bufSize = fullSize;
//
//		bufFullSize = bufSize;
//
//		.allocate_buffer_size(bufFullSize);
//
//	}
//
	//void copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
	//{
	//	if (writeSize == -1)
	//		writeSize = fullSize;
	//	.copy_to_buffer_size(writeSize, que_, bFlag_);

	//}
//
//	void copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
//	{
//		if (readSize == -1)
//			readSize = fullSize;
//		.read_from_buffer_size(readSize, que_, bFlag_);
//
//	}
//
//	void ini()
//	{
//		//initialization functions specific to class
//	}
//
//	bool load()
//	{
//		int ret = true;
//		ret &= .load("load" SLASH + pos.getName());
//
//		return ret;
//	}
//
//	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
//	{
//		bool ret = true;
//		switch (saveFlag_)
//		{
//		case saveFl:
//		{
//			ret &= .save();
//
//			break;
//		}
//		case saveTxtFl:
//		{
//			ret &= .save2file();
//
//			break;
//		}
//		case saveBinFl:
//		{
//			ret &= .savebin();
//
//			break;
//		}
//		}
//
//		if (saveAll)
//		{
//			switch (saveFlag_)
//			{
//			case saveFl:
//			{
//				ret &= .save();
//
//				break;
//			}
//			case saveTxtFl:
//			{
//				ret &= .save2file();
//
//				break;
//			}
//			case saveBinFl:
//			{
//				break;
//			}
//			}
//		}
//		return ret;
//	}
//
//	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
//		cl_command_queue* que_ = nullptr)
//	{
//		bool ret = true;
//		switch (saveFlag_)
//		{
//		case saveFl:
//		{
//			ret &= .save_from_device();
//
//			break;
//		}
//		case saveTxtFl:
//		{
//			ret &= .save_txt_from_device("", que_);
//
//			break;
//		}
//		case saveBinFl:
//		{
//			ret &= .save_bin_from_device("", que_);
//
//			break;
//		}
//		}
//
//		if (saveAll)
//		{
//			switch (saveFlag_)
//			{
//			case saveFl:
//			{
//				ret &= .save_from_device();
//
//				break;
//			}
//			case saveTxtFl:
//			{
//				ret &= .save_txt_from_device("", que_);
//
//				break;
//			}
//			case saveBinFl:
//			{
//				break;
//			}
//			}
//		}
//		return ret;
//	}
//};
