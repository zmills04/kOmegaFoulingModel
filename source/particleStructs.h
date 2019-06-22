// particleStructs.h: Structures holding various clVariablesTR data sets. Groups data for particle
// solver into structures each of which corresponds to a different part of solver.
// Legacy implementation used arrays of structs to hold data, but this was non-performant
// on GPUs due to it requiring large chunks of data to be read/written (especially bad for
// performance when reading/writing/passing full structure to use a single item from it.)
// Data can still be accessed through these structs on cpu side to ease transition
// from arrays of structs to structures of arrays. This should only be done during 
// initialization and maybe output steps since it will be very slow, and should be
// avoided in any frequently used functions (i.e. fouling layer update steps, etc)
// Note: classes are being used for structures to hold arrays in new implementation,
// just to maintain similar coding style as rest of project
//
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PARTICLESTRUCTS_H__INCLUDED_)
#define AFX_PARTICLESTRUCTS_H__INCLUDED_
#pragma once


#include "StdAfx.h"


////////////////////////////////
// Legacy Structures: data can be passed in/out to
// newly implemented structures of arrays using these
// old structures (all prepended with legacy in name)


class trStructBase
{
public:
	std::string structName;
	int xSize;
	int xSizeFull;
	int ySize;
	int fullSize;
	int bufFullSize;
	int restartFlag;
	enum saveFlags { saveFl, saveBinFl, saveTxtFl };
	
	trStructBase(std::string name_ = "")
	{
		ERROR_CHECKING((name_.length() == 0), "All trStructs must be initialized with a name", ERROR_INITIALIZING_VTR);
		structName = name_;
		restartFlag = false;
	}

	~trStructBase() {}

	virtual void setSizes(const int xsize, const int xsizefull, const int ysize)
	{
		xSize = xsize;
		ySize = ysize;
		xSizeFull = xsizefull;
		fullSize = ySize * xSizeFull;
	}

	virtual void allocateArrays() = 0;

	virtual void allocateBuffers(int bufSize = -1) = 0;

	virtual void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE) = 0;
	
	virtual void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE) = 0;

	virtual void ini() = 0;

	virtual bool load() = 0;

	virtual bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl) = 0;

	virtual bool saveFromDevice(bool saveAll = false, 
		saveFlags saveFlag_ = saveFl, cl_command_queue * que_ = nullptr) = 0;
};

typedef struct legacyPar
{
	cl_double2 pos;		///position vector
	cl_uint Num_rep;	///number of particles represented by this particle
	cl_short type;		///type of particle (cooresponds to Param element
	cl_short Dep_Flag; 	//-2 if waiting for re-release, -1 for not deposited, > 0 signifies the BL location it has deposited at.
	cl_ushort Dep_timer;	//timer set to specified value once particle deposits and decrements at 0, particle is deposited
	cl_ushort timer;		///time for use in re-releasing
	cl_int loc;			//Node number particle is located within

	legacyPar()
	{
	}

	//copy constructor
	legacyPar(const legacyPar& p_)
	{
		pos = p_.pos;		
		Num_rep = p_.Num_rep;
		type = p_.type;	
		Dep_Flag = p_.Dep_Flag; 
		Dep_timer = p_.Dep_timer;
		timer = p_.timer;
		loc = p_.loc;
	}

	legacyPar& addPar(legacyPar& p_)
	{
		this->Num_rep += p_.Num_rep;
		this->pos = Add2(this->pos, p_.pos);

//		this->timer = this->timer + p_.timer;
		return *this;
	}

	legacyPar operator+(legacyPar& p_)
	{
		legacyPar psum(*this);
		psum.Num_rep += p_.Num_rep;
		psum.pos.x += p_.pos.x;
		psum.pos.y += p_.pos.y;
		psum.timer += p_.timer;
		psum.timer = MIN(20000, psum.timer);
		return psum;
	}


} legacyPar;

class Par : public trStructBase
{
public:
	Array1Dv2d pos;		///position vector
	Array1Du Num_rep;	///number of particles represented by this particle
	Array1Ds type;		///type of particle (cooresponds to Param element
	Array1Di Dep_Flag; 	//-2 if waiting for re-release, -1 for not deposited, > 0 signifies the BL location it has deposited at.
	Array1D<cl_ushort> Dep_timer;	//timer set to specified value once particle deposits and decrements at 0, particle is deposited
	Array1D<cl_uint> timer;		///time for use in re-releasing
	Array1Di loc;			// probably easier to just calculate on demand than waste space with this

	enum arrName { posArr, numRepArr, typeArr, depFlagArr, depTimerArr, timerArr, locArr };

	Par(std::string name_ = "") {}

	~Par() {}

	void setStruct(legacyPar& struct_, const int i);

	legacyPar operator()(const int i) { return getStruct(i); }

	legacyPar getStruct(const int i);

	void setBuffers(Kernel& ker, int& curind, arrName arrList[], int numArrs);

	void setBuffersAll(Kernel& ker, int& curind)
	{
		arrName arrList_[] = { posArr, numRepArr, typeArr, depFlagArr, depTimerArr, timerArr, locArr };
		setBuffers(ker, curind, arrList_, 7);
	}

	void allocateArrays();

	void allocateBuffers(int bufSize = -1);

	void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	// Copies contents from Ptemp arrays on host to these arrays on device
	void writeParToBuffer(Par &Ptemp, int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	// Copies contents from Ptemp arrays on device to these arrays on device
	// pass copyLoc = false to skip copy of loc array
	void copyToParOnDevice(Par& Ptemp, bool copyLoc = true, int writeSize = -1,
		cl_command_queue* que_ = nullptr, int num_wait = 0, cl_event* wait = NULL,
		cl_event* evt = NULL);
	
	// Copies contents from Ptemp array specified with arrname_ on device to 
	// these arrays on device
	// pass copyLoc = false to skip copy of loc array
	void copyToArrayOnDevice(Par& Ptemp, arrName arrname_, int writeSize = -1,
		cl_command_queue* que_ = nullptr, int num_wait = 0, cl_event* wait = NULL,
		cl_event* evt = NULL);

	// Copies contents from Ptemp arrays on device to these arrays on device
	// using blocking copy call
	// pass copyLoc = false to skip copy of loc array
	void copyToParOnDeviceBlocking(Par& Ptemp, bool copyLoc = true, int writeSize = -1,
		cl_command_queue* que_ = nullptr, int num_wait = 0, cl_event* wait = NULL);

	// Copies contents from Ptemp array specified with arrname_ on device to 
	// these arrays on device using blocking copy call
	// pass copyLoc = false to skip copy of loc array
	void copyToArrayOnDeviceBlocking(Par& Ptemp, arrName arrname_, int writeSize = -1,
		cl_command_queue* que_ = nullptr, int num_wait = 0, cl_event* wait = NULL);



	void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void ini();

	bool load();

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl);

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr);
};


typedef struct legacyPParam
{
	cl_double2 Q_A_prime;	///Used in deposition/rebound calculation
	cl_double2 Q_A;			///used in dep/reb
	cl_double2 tau_crit;	///critical shear stress
	double Dp;			///particle diameter
	double Mp;			///particle mass
	double Kth;			///Thermophoretic Coeff
	double D_dist;		///% of particles distributed in bin
	double L_coeff;		///Lift coefficient (Lift force = L_coeff*Tau^1.5)
	double D_coeff;		///Drag coefficient (Drag force = D_coeff*Tau)
} legacyPParam;


class PParam : public trStructBase
{
public:
	Array1Dv2d Q_A_prime;	///Used in deposition/rebound calculation
	Array1Dv2d Q_A;			///used in dep/reb
	Array1Dv2d tau_crit;	///critical shear stress
	Array1Dd Dp;			///particle diameter
	Array1Dd Mp;			///particle mass
	Array1Dd Kth;			///Thermophoretic Coeff
	Array1Dd D_dist;		///% of particles distributed in bin
	Array1Dd L_coeff;		///Lift coefficient (Lift force = L_coeff*Tau^1.5)
	Array1Dd D_coeff;		///Drag coefficient (Drag force = D_coeff*Tau)

	enum arrName { qaPrimeArr, qaArr, tauCritArr, dpArr, mpArr, kthArr, distArr, dCoeffArr, lCoeffArr };

	PParam(std::string name_ = "") {}

	~PParam() {}

	void setBuffers(Kernel& ker, int& curind, arrName arrList[], int numArrs);

	void setBuffersAll(Kernel& ker, int& curind)
	{
		arrName arrList_[] = { qaPrimeArr, qaArr, tauCritArr, dpArr, mpArr, kthArr, distArr, dCoeffArr, lCoeffArr };
		setBuffers(ker, curind, arrList_, 9);
	}

	void setStruct(legacyPParam& struct_, const int i);

	legacyPParam getStruct(const int i);

	void allocateArrays();

	void allocateBuffers(int bufSize = -1);

	void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void ini();

	// all data created from parameters at initialization
	// nothing is stored in bin files
	bool load() { return true;  }

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl);

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr);
};


typedef struct legacyNodeI
{
	int BLind;	//Index of BL with center closest to center of trNode (i+0.5,j+0.5)
	cl_short wallFlag;  // -1 represents a full wall node,
						// 0 is for no walls nearby,
						// 1 if bottom wall nearby and
						// 2 if top wall nearby
} legacyNodeI;


class NodeI : public trStructBase
{
public:

	//Starting index of three closest boundary nodes (BLind, BLind+1, BLind+2 are the nodes). Note: must test for bounds to ensure no OOB accesses
	Array2Di BLind; 

	//stores flag for if wall is nearby 
	Array2D<cl_short> wallFlag;	// -1 represents a full wall node,
						// 0 is for no walls nearby,
						// 1 if bottom wall nearby and
						// 2 if top wall nearby


	NodeI(std::string name_ = "") {}

	~NodeI() {}

	void setStruct(legacyNodeI& struct_, const int i, const int j);

	legacyNodeI getStruct(const int i, const int j);

	void allocateArrays();

	void allocateBuffers(int bufSize = -1);

	void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void ini();

	// can easily be generated from other data, so no need to save/load checkpoint data
	bool load() { return true; }

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl);

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr);
};

typedef struct legacyNeighs
{
	cl_uint ii00;
	cl_uint ii10;
	cl_uint ii01;
	cl_uint ii11;

} legacyNeighs;


class Neighs : public trStructBase
{
public:

	// Old Implementation
	//typedef struct legacyNeighs
	//{
	//	cl_int2 ii00;
	//	cl_int2 ii10;
	//	cl_int2 ii01;
	//	cl_int2 ii11;

	//} legacyNeighs;

	Array2Di ii00, ii10, ii01, ii11;

	Neighs(std::string name_ = "") {}

	~Neighs() {}

void setStruct(legacyNeighs& struct_, const int i, const int j)
{
	ii00(i, j) = struct_.ii00;
	ii10(i, j) = struct_.ii10;
	ii01(i, j) = struct_.ii01;
	ii11(i, j) = struct_.ii11;
}

legacyNeighs getStruct(const int i, const int j)
{
	legacyNeighs struct_;
	struct_.ii00 = ii00(i, j);
	struct_.ii10 = ii10(i, j);
	struct_.ii01 = ii01(i, j);
	struct_.ii11 = ii11(i, j);
	return struct_;
}


	void allocateArrays()
	{
		ii00.setName("Neighs_ii00");
		ii10.setName("Neighs_ii10");
		ii01.setName("Neighs_ii01");
		ii11.setName("Neighs_ii11");

		ii00.zeros(xSize, xSizeFull, ySize, ySize);
		ii01.zeros(xSize, xSizeFull, ySize, ySize);
		ii10.zeros(xSize, xSizeFull, ySize, ySize);
		ii11.zeros(xSize, xSizeFull, ySize, ySize);
	}

	void allocateBuffers(int bufSize = -1)
	{
		if (bufSize == -1)
			bufSize = fullSize;

		bufFullSize = bufSize;

		ii00.allocate_buffer_size(bufFullSize);
		ii10.allocate_buffer_size(bufFullSize);
		ii01.allocate_buffer_size(bufFullSize);
		ii11.allocate_buffer_size(bufFullSize);

	}

	void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void ini();

	// can be reconstructed
	bool load() { return true; }

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl);

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr);
};


typedef struct legacyTrParam
{//Used for re-releasing particles
	cl_double Top_location;	//Y value of uppermost node
	cl_double Bottom_location;	//Y value of lowermost node
	cl_double umax_val;		//Max velocity at inlet
	cl_double bval;			//spacing between upper and lower wall at particle inlet
	cl_double offset_y;		//Location of wall at particle inlet
	
	// can be set using define cl_double X_release;	//X location of release
	// Can be set using define cl_uint BL_rel_bot;	//index of bottom BL at particle inlet (shifted to put bottom wall at zero) 
	// can be set using define cl_uint BL_rel_top;	//index of top BL at particle inlet (shifted to put bottom wall at zero)
	
	// can be set using define (shouldnt ever change since no dep immediatly)
	//cl_uint Uvals_start;	//location of starting point for Uvals used in re-distribution

} legacyTrParam;



class TrParam : public trStructBase
{
public:
	Array1Dd fVals;


	TrParam(std::string name_ = "") 
	{
		xSize = 5;
		ySize = 1;
		xSizeFull = 5;
		fullSize = ySize * xSizeFull;
	}

	~TrParam() {}

	legacyTrParam getStruct();

	void allocateArrays();

	void allocateBuffers(int bufSize = -1);

	void setSizes(const int xsize, const int xsizefull, const int ysize) override
	{
		// want to avoid letting sizes which were set in constructor be set to
		// different values
	}

	void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void ini();

	// can be reconstructed
	bool load() { return true; }

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl);

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr);

	void setBuffers(Kernel& ker, int& curind)
	{
		ker.set_argument(curind++, fVals.get_buf_add());
	}

	double& operator[](const int ind_)
	{
		return fVals(ind_);
	}

};

typedef struct legacyNodeC
{
	cl_double4 CoeffT00;
	cl_double4 CoeffT10;
	cl_double4 CoeffT01;
	cl_double4 CoeffT11;
	cl_double4 CoeffU00;
	cl_double4 CoeffU10;
	cl_double4 CoeffU01;
	cl_double4 CoeffU11;
} legacyNodeC;

// TODO: remove neigh if unnecessary

class NodeC : public trStructBase
{
public:
	Array2Dv4d CoeffT00, CoeffT10, CoeffT01, CoeffT11;
	Array2Dv4d CoeffU00, CoeffU10, CoeffU01, CoeffU11;

	void setTempBuffers(Kernel& ker, int& curind);
	void setVelBuffers(Kernel& ker, int& curind);

	NodeC(std::string name_ = "") {}

	~NodeC() {}

	void setStruct(legacyNodeC& struct_, const int i, const int j);

	legacyNodeC getStruct(const int i, const int j);

	void allocateArrays();

	void allocateBuffers(int bufSize = -1);

	void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void ini()
	{
		//initialization functions specific to class
	}

	// can re-construct without saving data
	bool load()	{ return true; }

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl);

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr);
};


typedef struct legacyNodeV
{
	cl_double4 Temps;
	cl_double2 U00;
	cl_double2 U10;
	cl_double2 U01;
	cl_double2 U11;
} legacyNodeV;



class NodeV : public trStructBase
{
public:
	Array2Dv4d Temps;
	Array2Dv2d U00;
	Array2Dv2d U10;
	Array2Dv2d U01;
	Array2Dv2d U11;
	NodeV(std::string name_ = "") {}

	~NodeV() {}

	void setStruct(legacyNodeV& struct_, const int i, const int j);

	legacyNodeV getStruct(const int i, const int j);

	void allocateArrays();

	void allocateBuffers(int bufSize = -1);

	void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void ini();

	bool load()	{ return true; }

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl);

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr);
};


typedef struct legacyBLinks
{
	//cl_double2 vP0;		//Location of left node
	//cl_double2 vP1;		//Location of right node
	//cl_double2 vTvec;	//tangential vector (points downstream) not needed as vT = (vN.y, -vN.x)
	cl_double2 vNvec;   //normal vector pointing into domain
	cl_double Tau;		//Shear stress at location (negative for directed upstream, pos downstream)
	cl_double blLen;	//length of BL
	cl_uint Node_loc;	//points to location BL is located in (or majority of BL when across two)
	cl_ushort2 P01ind; // basically a reduced version of vls.BL with only BLinks indicies
	cl_short int_type;  //0 if wall, 1 if soot (changes after sufficiently thick deposit forms)
	int colorInd; // tracks indicies of boundary layers for opengl display.
} legacyBLinks;


class BLinks : public trStructBase
{
public:
	//Array1Dv2d vTvec;	//tangential vector (points downstream)
	Array1Dv2d vNvec;   //normal vector pointing into domain
	Array1Dd Tau;		//Shear stress at location (negative for directed upstream, pos downstream)
	Array1Dd blLen;		//length of BL
	Array1Du Node_loc;	//points to location BL is located in (or majority of BL when across two)
						// note: this is the vtr node it is closest to, so it is based on vtr.trDomainSize
						//			not p.nX and p.nY
	Array1Dv2us P01ind;	//basically a reduced version of vls.BL with only BLinks indicies
	Array1Ds int_type;  //0 if wall, 1 if soot (changes after BLinks displaced certain FOUL_SIZE_SWITCH_SOOT2)
	Array1Di colorInds; //tracks indicies of boundary layers for opengl display.
	
	enum arrName { vNArr, tauArr, lenArr, nodLocArr, P01Arr, typeArr, colorArr };

	bool glFlag = false;
	
	BLinks(std::string name_ = "") {}

	~BLinks() {}

	void setBuffers(Kernel& ker, int& curind, arrName arrList[], int numArrs);

	void setBuffersAll(Kernel& ker, int& curind)
	{
		arrName arrList_[] = { vNArr, tauArr, lenArr, nodLocArr, P01Arr, typeArr, colorArr };
		setBuffers(ker, curind, arrList_, 7);
	}
	
	void setStruct(legacyBLinks& struct_, const int i);

	legacyBLinks getStruct(const int i);

	void allocateArrays();

	void allocateBuffers(int bufSize = -1);

	void copyToDevice(int writeSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);

	void copyToHost(int readSize = -1, cl_command_queue* que_ = nullptr,
		cl_bool bFlag_ = CL_TRUE);
	
	void createColorInds();

	void ini();
	
	bool load() { return true; }

	bool save(bool saveAll = false, saveFlags saveFlag_ = saveFl);

	bool saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
		cl_command_queue* que_ = nullptr);
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


















#endif