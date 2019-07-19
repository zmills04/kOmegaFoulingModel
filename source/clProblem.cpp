// clProblem.cpp: implementation of the clProblem class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "clProblem.h"
#include "BiCGStabGenerator.h"



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////             Functions for Outputting Data                  ////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// Save Files, Check for save bin, check for savedump
void clProblem::DumpStep()
{
	// save Arrays and Time Data
	// Functions being called will check flags to see if
	// they should save data (rarely called, so not calling
	// function that returns immediatly will not affect
	// performance
	vlb.save2file();
	vlb.saveTimeData();
	
	vfd.save2file();
	vfd.saveTimeData();
	
	vtr.save2file();
	vtr.saveTimeData();
	
	vls.save2file();
	vls.saveTimeData();

	vfl.save2file();
	vfl.saveTimeData();
		
	currentDumpStep++;
};

// Creates folder in results and copies data currently in
// main dir to new folder
void clProblem::RenameOutputFiles()
{
	if (flOutputDump == false)
		return;

	std::ostringstream streamObj;
	streamObj << std::scientific;
	streamObj << std::setprecision(6);
	streamObj << (double)Time;

	std::string curDumpSaveDir = outDir;
	curDumpSaveDir.append(SLASH);
	curDumpSaveDir.append(streamObj.str());

	staticBaseVar::setCurDumpDir(curDumpSaveDir);

	MakeDir(curDumpSaveDir);

	vls.renameSaveFiles();
	vlb.renameSaveFiles();
	vfd.renameSaveFiles();
	vtr.renameSaveFiles();
	vfl.renameSaveFiles();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////                Simulation Control Functions                ////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


void clProblem::ini()
{
	srand(1);

	vfd.thermalSolverFlag = false;
	vtr.trSolverFlag = false;
	vfl.flSolverFlag = false;

	loadRunParameters("load" SLASH YAML_PARAM_FILE);

	outDir = OUTPUT_DIR;
	MakeDir(outDir);
	setSourceDefines();
};


void clProblem::step()
{
	if (Next_Clump_Time < TimeN && vtr.trSolverFlag)
	{
		vtr.parSort.clumpParticles();
		Next_Clump_Time += Time_Btw_Clump;
	}

	if (nextSaveAvgTime < TimeN)
	{
		vlb.updateTimeData();
		nextSaveAvgTime += (int)dTlb * saveTimeStepAvg;
	}

	if (nextSaveNuTime < TimeN && vfd.calcNuFlag)
	{
		vfd.updateTimeData();
		nextSaveNuTime += (int)dTlb * saveTimeStepNu;
	}

	if (nextSaveSSTime < TimeN && vtr.wallShear.calcSSFlag)
	{
		vtr.wallShear.updateTimeData();
		nextSaveSSTime += (int)dTlb * saveTimeStepSS;
	}

	if (nextSaveIOTime < TimeN && vtr.calcIOFlag)
	{
		vtr.updateTimeData();
		nextSaveIOTime += (int)dTlb * saveTimeStepIO;
	}

	bool dumpBinFlag = false; 
	bool saveDumpFlag = false;
	bool dumpFlag = false;

	if (nextDumpStepTime < TimeN)
	{
		p.DumpStep();
		dumpFlag = true;
		clFinish(*clEnv::instance()->getIOqueue());
		nextDumpStepTime += dumpTimeStep;
	}

	// if time to save bin files, do so
	if (nextDumpBin < currentDumpStep)
	{
		dumpBinFlag = true;
		nextDumpBin += saveBinStep;
		vls.saveRestartFiles();
		vlb.saveRestartFiles();
		vfd.saveRestartFiles();
		vtr.saveRestartFiles();
		vfl.saveRestartFiles();
	}

	// If dump save step, copy files from main folder
	// to newly created results folder
	if (nextDumpSave < currentDumpStep)
	{
		saveDumpFlag = true;
		nextDumpSave += saveDumpStep;
		RenameOutputFiles();
	}

	if (dumpFlag)
	{
		saveParameters();
	}
	if (dumpBinFlag)
	{
		CopyFile(YAML_PARAM_FILE, "load" SLASH YAML_PARAM_FILE);
	}
	if (saveDumpFlag)
	{
		std::string saveYAMLName = staticBaseVar::curDumpSaveDir + SLASH + YAML_PARAM_FILE;
		RenameFile(YAML_PARAM_FILE, saveYAMLName);
	}
	
#ifdef MAINTAIN_FLOWRATE
	if (vlb.Next_Update_Time < TimeN)
	{
		vlb.Calculate_U_mean();
	}
#endif

	Time = TimeN;
	TimeN += (unsigned int)dTlb;
}


void clProblem::finish()
{
	DumpStep();
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////          Function to Add Defines to Kernel Src             ////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


// TODO: remove unncessary define and consolidate remaining defines
// TODO: make sure that all domain size definitions are used correctly
//			especially FULLSIZEY (seems like it should be nY, not channel height)
void clProblem::setSourceDefines()
{
	/////////////  Domain Sizes///////////////////////////////////////////
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "FULLSIZEX", static_cast<int>(nX));
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "FULLSIZEY", static_cast<int>(Channel_Height));
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "CHANNEL_LENGTH", static_cast<int>(nX));
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "CHANNEL_HEIGHT", static_cast<int>(nY));
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "CHANNEL_LENGTH_FULL", static_cast<int>(XsizeFull));
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "DIST_SIZE", static_cast<int>(XsizeFull*nY));
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "FULLSIZEXY", static_cast<int>(Channel_Height*nX));
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "FULLSIZEY_UT", static_cast<int>(Channel_Height)+1);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "DOMAIN_SIZE_Y", static_cast<int>(nY));

	/////////////  Workgroup Sizes ////////////////////////////////////////
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "BLOCK_SIZE", BLOCKSIZE);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_RED", WORKGROUPSIZE_RED);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZEX_LB", WORKGROUPSIZEX_LB);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZEY_LB", WORKGROUPSIZEY_LB);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_INLET", WORKGROUPSIZE_INLET);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_IBB", WORKGROUPSIZE_IBB);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZEX_FD", WORKGROUPSIZEX_FD);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZEY_FD", WORKGROUPSIZEY_FD);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_NU", WORKGROUPSIZE_NU);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_RERELEASE", (int)WORKGROUPSIZE_RERELEASE);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_TR", WORKGROUPSIZE_TR);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_TR_WALL", WORKGROUPSIZE_TR_WALL);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_TR_GL", WORKGROUPSIZE_TR_GL);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_TR_SHEAR", WORKGROUPSIZE_TR_SHEAR);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_SORT", WORKGROUPSIZE_SORT);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_PLOC", WORKGROUPSIZE_PLOC);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_FL_SHIFT", WORKGROUPSIZE_FL_SHIFT);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_FL_RAMP", NUM_INLET_OUTLET_NODES);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_UPDATEM", WORKGROUPSIZE_UPDATEM);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_UPDATEFD", WORKGROUPSIZE_UPDATEFD);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_UPDATEWALL", WORKGROUPSIZE_UPDATEWALL);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_UPDATE_GL", WORKGROUPSIZE_UPDATE_GL);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_SHIFT_PAR", WORKGROUPSIZE_SHIFT_PAR);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_TR_DISTS", WORKGROUPSIZE_TR_DISTS);
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "WORKGROUPSIZE_TR_WALL_REFLECT", WORKGROUPSIZE_TR_WALL_REFLECT);

	if(useOpenGL)
		SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "USE_OPENGL");

#ifdef OPENCL_VERSION_1_2
	SOURCEINSTANCE->addDefine(SOURCEINSTANCE->getDefineStr(), "OPENCL_VERSION_1_2");
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////               Parameter Saving  Functions                  ////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void clProblem::saveParameters()
{
	yamlOut = new YAML::Emitter();

	yamlOut->SetDoublePrecision(15);
	yamlOut->SetFloatPrecision(15);
	yamlOut->SetBoolFormat(YAML::TrueFalseBool);


	*p.yamlOut << YAML::BeginDoc;
	saveSystemParams();
	*p.yamlOut << YAML::EndDoc;
	
	*p.yamlOut << YAML::BeginDoc;
	vls.saveParams();
	*p.yamlOut << YAML::EndDoc;
	
	*p.yamlOut << YAML::BeginDoc;
	vlb.saveParams();
	*p.yamlOut << YAML::EndDoc;
	
	*p.yamlOut << YAML::BeginDoc;
	vfd.saveParams();
	*p.yamlOut << YAML::EndDoc;
	
	*p.yamlOut << YAML::BeginDoc;
	vtr.saveParams();
	*p.yamlOut << YAML::EndDoc;
	
	*p.yamlOut << YAML::BeginDoc;
	vfl.saveParams();
	*p.yamlOut << YAML::EndDoc;

	std::string fname_ = YAML_PARAM_FILE;

	std::fstream outstream;

	fileopen(outstream, fname_, FileOut);

	outstream << p.yamlOut->c_str();

	delete yamlOut;
}




template <typename T>
void clProblem::setParameter(std::string pname_, T val_)
{
	*yamlOut << YAML::BeginMap;
	*yamlOut << YAML::Key << pname_;
	*yamlOut << YAML::Value << val_;
	*yamlOut << YAML::EndMap;
}

template <>
void clProblem::setParameter(std::string pname_, std::string val_)
{
	*yamlOut << YAML::BeginMap;
	*yamlOut << YAML::Key << pname_;
	*yamlOut << YAML::Value << val_;
	*yamlOut << YAML::EndMap;
}


// THis should be called before incrementing Time,
// but after all other outputs have been done, and 
// their counters have been increased
void clProblem::saveSystemParams()
{

	p.setParameter("ParamsName", systemParamNum);
	if (p.Time > 0)
		p.setParameter("Restart Run", true);
	else
		p.setParameter("Restart Run", false);

	p.setParameter("Device ID", DeviceID);
	p.setParameter("nX", nX);
	p.setParameter("nY", nY);
	p.setParameter("Pipe radius", Pipe_radius);
	p.setParameter("Current Time", Time + 1);
	p.setParameter("Stop Time", StopTime);
	p.setParameter("TR Steps Per LB", trSteps);
	p.setParameter("TR Wall Steps Per LB", trSteps_wall);
	p.setParameter("T Steps Per LB", tfdSteps);
	p.setParameter("Clumping Time", Time_Btw_Clump);
	p.setParameter("Save Bin Step", saveBinStep);
	p.setParameter("Save Avgs Start Time", saveAvgStartTime);
	p.setParameter("Save Nu Start Time", saveNuStartTime);
	p.setParameter("Save Shear Start Time", saveSSStartTime);
	p.setParameter("Save IO Start Time", saveIOStartTime);
	p.setParameter("Dump Start Time", dumpStartTime);

	if (avgNumStepDef)
		p.setParameter("Save Avgs Num Steps", saveStepNumAvg);
	else
		p.setParameter("Save Avgs Timestep", saveTimeStepAvg);

	if (nuNumStepDef)
		p.setParameter("Save Nu Num Steps", saveStepNumNu);
	else
		p.setParameter("Save Nu Timestep", saveTimeStepNu);

	if (ssNumStepDef)
		p.setParameter("Save Shear Num Steps", saveStepNumSS);
	else
		p.setParameter("Save Shear Timestep", saveTimeStepSS);

	if (ioNumStepDef)
		p.setParameter("Save IO Num Steps", saveStepNumIO);
	else
		p.setParameter("Save IO Timestep", saveTimeStepIO);

	if (dumpNumStepDef)
		p.setParameter("Dump Num Steps", dumpStepNum);
	else
		p.setParameter("Dump Timestep", dumpTimeStep);

	p.setParameter("Dump Steps Per Save", saveDumpStep);
	p.setParameter("Start Dump Step Save", saveDumpStepStart);
	p.setParameter("Start Dump Step Save", saveDumpStepStart);
	p.setParameter("Save Dump Files", flOutputDump);
	p.setParameter("Display Signal Freq", displaySignalFreq);
	p.setParameter("currentDumpStep", currentDumpStep);
	p.setParameter("nextDumpSave", nextDumpSave);
	p.setParameter("nextDumpBin", nextDumpBin);
	p.setParameter("Next_Clump_Time", Next_Clump_Time);
	p.setParameter("nextSaveAvgTime", nextSaveAvgTime);
	p.setParameter("nextSaveNuTime", nextSaveNuTime);
	p.setParameter("nextSaveSSTime", nextSaveSSTime);
	p.setParameter("nextDumpStepTime", nextDumpStepTime);
	p.setParameter("nextSaveIOTime", nextSaveIOTime);
	p.setParameter("Use OpenGL", useOpenGL);

}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////                Parameter Loading Functions                 ////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
T clProblem::getParameter(std::string pname_)
{
	T retval_;
	if (const YAML::Node *pName = yamlIn.FindValue(pname_))
	{
		*pName >> retval_;
	}
	else
	{
		ERROR_CHECKING(true, "Parameter " + pname_ + " does not have "\
			"a default value and must be provided in the parameter file", ERROR_LOADING_PARAMS);
	}
	return retval_;
}



void clProblem::getSaveStepInfo(std::string varname_, unsigned int &starttime_,
	unsigned int &numsteps_, unsigned int &timestep_, bool &numdefined_,
	const unsigned int definedstart_, const unsigned int definednumsteps_)
{
	starttime_ = p.getParameter(varname_ + " Start Time", definedstart_);
	int numstepsFl = p.getParameter(varname_ + " Num Steps",
		varname_ + " Timestep", numsteps_, timestep_);

	if (numstepsFl == 1)
	{
		numdefined_ = true;
		timestep_ = (StopTime - starttime_) / numsteps_;
	}
	else if (numstepsFl == 2)
	{
		numdefined_ = false;
		numsteps_ = (StopTime - starttime_) / timestep_;
	}
	else
	{
		numdefined_ = true;
		numsteps_ = definednumsteps_;
		timestep_ = (StopTime - starttime_) / numsteps_;
	}
}

void clProblem::loadParams()
{
	DeviceID = p.getParameter<int>("Device ID", DEVICE_ID);
	nX = p.getParameter<int>("nX");
	nY = p.getParameter<int>("nY");
	XsizeFull = int(ceil((double)nX / (double)BLOCKSIZE)) * BLOCKSIZE;

	xWavyStart = p.getParameter("Wavy Start Loc", -1);
	wavyPeriodLen = p.getParameter("Period Length", -1);
	numWavyPeriods = p.getParameter("Num Wave Periods", -1);


	nn = { { nX, nY } };
	Pipe_radius = p.getParameter<double>("Pipe Radius");
	Channel_Height = 2 * (int)Pipe_radius;
	FullSize = 2;
	unsigned int domainsize_ = (unsigned int)(XsizeFull*nY);
	while (FullSize < domainsize_)
		FullSize *= 2;
	distSize = domainsize_;
	DELTA_L = ((2.*Pipe_radius) / PIPE_DIAMETER_REAL);

	useOpenGL = p.getParameter<bool>("Use OpenGL", OPENGL_DISPLAY);

	clEnv::instance()->iniOpenCL(DeviceID, useOpenGL);
	sourceGenerator::SourceInstance();
	BiCGStabGenerator::BiCGStabInstance()->ini(nX, XsizeFull, nY);
	ReduceGenerator::ReduceInstance()->ini(nX, nY);



#ifndef PROFILING
	StopTime = p.getParameter<int>("Stop Time", STOP_TIME);
#else
	StopTime = STOP_TIME;
#endif

	RestartTime = p.getParameter("Current Time", 0);
	if (RestartTime > 0)
		restartRunFlag = true;
	else
		restartRunFlag = false;

	trSteps = p.getParameter("TR Steps Per LB", TR_STEPS_PER_LB);
	trSteps_wall = p.getParameter("TR Wall Steps Per LB", TR_STEPS_PER_LB_WALL);
	tfdSteps = p.getParameter("T Steps Per LB", FD_STEPS_PER_LB);
	Time_Btw_Clump = p.getParameter("Clumping Time", CLUMP_TIME);

	dTlb = 1.;
	dTtfd = dTlb / (double)tfdSteps;
	dTtr = dTlb / (double)trSteps;
	dTtr_wall = dTlb / (double)trSteps_wall;

	// Number of dump steps between save bin steps (keep at 1)
	saveBinStep = p.getParameter("Save Bin Step", DUMP_STEPS_PER_DUMP_BIN);

	getSaveStepInfo("Save Avgs", saveAvgStartTime, saveStepNumAvg, saveTimeStepAvg,
		avgNumStepDef, CURRENT_SAVE_START_TIME, CURRENT_SAVE_STEP_NUM);

	getSaveStepInfo("Save Nu", saveNuStartTime, saveStepNumNu, saveTimeStepNu,
		nuNumStepDef, CURRENT_SAVE_START_TIME_NU, CURRENT_SAVE_STEP_NUM_NU);

	getSaveStepInfo("Save Shear", saveSSStartTime, saveStepNumSS, saveTimeStepSS,
		ssNumStepDef, CURRENT_SAVE_START_TIME_SS, CURRENT_SAVE_STEP_NUM_SS);

	getSaveStepInfo("Save IO", saveIOStartTime, saveStepNumIO, saveTimeStepIO,
		ioNumStepDef, CURRENT_SAVE_START_TIME_IO, CURRENT_SAVE_STEP_NUM_IO);

	getSaveStepInfo("Dump", dumpStartTime, dumpStepNum, dumpTimeStep,
		dumpNumStepDef, DUMP_STEP_START, DUMP_STEP_NUM);

	saveDumpStep = p.getParameter("Dump Steps Per Save", DUMP_STEPS_PER_SAVE);

	saveDumpStepStart = p.getParameter("Start Dump Step Save", DUMP_STEP_START_SAVE);

	saveDumpStepStart = p.getParameter("Start Dump Step Save", DUMP_STEP_START_SAVE);

	flOutputDump = p.getParameter("Save Dump Files", FLAG_DUMP);

	displaySignalFreq = p.getParameter("Display Signal Freq", StopTime / 10000);

	// These will always exist with a restart. If not, set nextDumpSave to saveDumpStepStart if > 0, saveDumpStep otherwise
	currentDumpStep = p.getParameter("currentDumpStep", 0);
	nextDumpSave = p.getParameter("nextDumpSave", saveDumpStepStart);
	nextDumpBin = p.getParameter("nextDumpBin", saveBinStep);
	Next_Clump_Time = p.getParameter("Next_Clump_Time", Time_Btw_Clump);
	nextSaveAvgTime = p.getParameter("nextSaveAvgTime", MAX(saveAvgStartTime, saveTimeStepAvg));
	nextSaveNuTime = p.getParameter("nextSaveNuTime", MAX(saveNuStartTime, saveTimeStepNu));
	nextSaveSSTime = p.getParameter("nextSaveSSTime", MAX(saveSSStartTime, saveTimeStepSS));
	nextDumpStepTime = p.getParameter("nextDumpStepTime", dumpStartTime);
	nextSaveIOTime = p.getParameter("nextSaveIOTime", MAX(saveIOStartTime, saveTimeStepIO));
	if (displaySignalFreq < 1)
		displaySignalFreq = 1;
}

void clProblem::parseAndLoadParams(std::string &runparams_, int paramnum_, 
	std::function<void(void)> &loadparamptr_)
{
	std::ifstream fin(runparams_);

	parser.Load(fin);
	//yamlIn.Clear();
	
	while (parser.GetNextDocument(yamlIn))
	{
		int ParamNum;
		yamlIn["ParamsName"] >> ParamNum;

		if (ParamNum == paramnum_)
		{
			loadparamptr_();
			
			return;
		}
	}

	parseEmptyFile(paramnum_);
	loadparamptr_();
	remove("empty.yaml");
}

void clProblem::parseEmptyFile(int addnum_)
{
	std::ofstream  out_;
	out_.open("empty.yaml");
	out_ << "---\n";
	out_ << "ParamsName: " << addnum_ << "\n";
	out_.close();
	std::ifstream fintemp;
	fintemp.open("empty.yaml");
	parser.Load(fintemp);
	parser.GetNextDocument(yamlIn);
}

void clProblem::loadRunParameters(std::string runparams_)
{
	using std::string;
	//std::ifstream fin2("test.yaml");

	//YAML::Parser yamnode(fin2);
	//YAML::Node doc2;
	//parser.GetNextDocument(doc2);

	// Note: from what I can determine, the YAML-cpp library is limited to
	// keys with less than 16 characters. Also, the release version is having
	// issues that I could not resolve and decided that since its negligible
	// additional time using debug version, just using that.

	parseAndLoadParams(runparams_, systemParamNum, p.loadParamsPtr);
	parseAndLoadParams(runparams_, lsParamNum, vls.loadParamsPtr);
	parseAndLoadParams(runparams_, fluidParamNum, vlb.loadParamsPtr);
	parseAndLoadParams(runparams_, thermalParamNum, vfd.loadParamsPtr);
	parseAndLoadParams(runparams_, trParamNum, vtr.loadParamsPtr);
	parseAndLoadParams(runparams_, flParamNum, vfl.loadParamsPtr);

	restartRunFlag &= vlb.restartRunFlag;
	restartRunFlag &= vls.restartRunFlag;

	if (vfd.thermalSolverFlag)
		restartRunFlag &= vfd.restartRunFlag;
	if (vtr.trSolverFlag)
		restartRunFlag &= vtr.restartRunFlag;
	if (vfl.flSolverFlag)
		restartRunFlag &= vfl.restartRunFlag;

	if (!restartRunFlag)
	{
		restartRunFlag = false;
		RestartTime = 0;
		currentDumpStep = 0;
		nextDumpSave = saveDumpStepStart;
		nextDumpBin = 1; // dont want to save bin at Time=0

		// Only the next dumpstep should start at 0 (if dumpStartTime = 0)
		Next_Clump_Time = Time_Btw_Clump;
		nextSaveAvgTime = MAX(saveAvgStartTime, saveTimeStepAvg);
		nextSaveNuTime = MAX(saveNuStartTime, saveTimeStepNu);
		nextSaveSSTime = MAX(saveSSStartTime, saveTimeStepSS);
		nextDumpStepTime = dumpStartTime;
		nextSaveIOTime = MAX(saveIOStartTime, saveTimeStepIO);
	}
	Time = RestartTime;
	TimeN = RestartTime + 1;
	if (Next_Clump_Time <= Time)
		Next_Clump_Time = Time + 1;

	// Pass restart flag to clEnv so the generator can access the
	// information
	clEnv::instance()->setRestartFlag(restartRunFlag);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////             Initialization of Static Variables             ////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

const double clProblem::R = 287.;
const double clProblem::Pi = PI_NUMBER;
const double clProblem::Kb = 1.38066e-23;
const double clProblem::eps = 1.e-5;

double clProblem::DELTA_L = 0.;
double clProblem::DELTA_T = 0.;
double clProblem::DELTA_M = 0.;
double clProblem::DELTA_F = 0.;
double clProblem::DELTA_P = 0.;

double clProblem::Pipe_radius;
int clProblem::Channel_Height = 0;
int clProblem::nX = 0;
int clProblem::nY = 0;
int clProblem::XsizeFull = 0;
unsigned int clProblem::FullSize = 0;
unsigned int clProblem::distSize = 0;
cl_int2 clProblem::nn = { { 0, 0 } };

bool clProblem::restartRunFlag = 0;
bool clProblem::flOutputDump = FLAG_DUMP;

unsigned int clProblem::StopTime = STOP_TIME;

double clProblem::dTlb = LB_STEPS;
int clProblem::tfdSteps = FD_STEPS_PER_LB;
int clProblem::trSteps = TR_STEPS_PER_LB;
int clProblem::trSteps_wall = TR_STEPS_PER_LB_WALL;
cl_short clProblem::IBB_Flag[8] = { 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80 };

int clProblem::DeviceID = DEVICE_ID;
double clProblem::dTtr = 0.;
double clProblem::dTtr_wall = 0.;
double clProblem::dTtfd = 0.;
unsigned int clProblem::Time = 0;
unsigned int clProblem::RestartTime = 0;
unsigned int clProblem::TimeN = 0;
unsigned int clProblem::StartTime = 0;
unsigned int clProblem::Time_Btw_Clump = 0;
unsigned int clProblem::Next_Clump_Time = 0;
unsigned int clProblem::saveBinStep = 0;

unsigned int clProblem::saveAvgStartTime = 0;
unsigned int clProblem::saveTimeStepAvg = 0;
unsigned int clProblem::saveStepNumAvg = 0;
unsigned int clProblem::nextSaveAvgTime = 0;
unsigned int clProblem::saveNuStartTime = 0;
unsigned int clProblem::saveTimeStepNu = 0;
unsigned int clProblem::saveStepNumNu = 0;
unsigned int clProblem::nextSaveNuTime = 0;
unsigned int clProblem::saveSSStartTime = 0;
unsigned int clProblem::saveTimeStepSS = 0;
unsigned int clProblem::saveStepNumSS = 0;
unsigned int clProblem::nextSaveSSTime = 0;
unsigned int clProblem::saveIOStartTime = 0;
unsigned int clProblem::saveTimeStepIO = 0;
unsigned int clProblem::saveStepNumIO = 0;
unsigned int clProblem::nextSaveIOTime = 0;
unsigned int clProblem::dumpStartTime = 0;
unsigned int clProblem::dumpTimeStep = 0;
unsigned int clProblem::dumpStepNum = 0;
unsigned int clProblem::nextDumpStepTime = 0;
unsigned int clProblem::saveDumpStep = 0;
unsigned int clProblem::saveDumpStepStart = 0;
unsigned int clProblem::displaySignalFreq = 0;
int clProblem::currentDumpStep = 0;
int clProblem::nextDumpSave = 0;
int clProblem::nextDumpBin = 0;

bool clProblem::avgNumStepDef = false;
bool clProblem::nuNumStepDef = false;
bool clProblem::ssNumStepDef = false;
bool clProblem::ioNumStepDef = false;
bool clProblem::dumpNumStepDef = false;
bool clProblem::useOpenGL = false;



///// Not updated and checked for errors
//void clProblem::SaveNu()
//{
//	FILE *hFile;
//	char FileName[] = NUSSELT_OUTPUT_FILE;
//	int maxAttempt = 100;
//	for (int i = 0; i < maxAttempt; i++)
//	{
//		hFile = fopen(FileName, "a");
//		if (hFile != NULL)
//			break;
//		printf("cannot open output file %s, waiting", FileName);
//		delay_func(0.1);
//	}
//
//	if (hFile == NULL)
//	{
//		printf("cannot open output file %s, data lost", FileName);
//		return;
//	}
//
//	clFinish(*clEnv::instance()->getFDqueue());
//	if (vfd.Save_loc_Nu == 0)
//		return;
//	vfd.Nu.read_from_buffer_size(vfd.Save_loc_Nu*vfd.Nu.getSizeY());
//
//	for (int i = 0; i < vfd.Save_loc_Nu; i++)
//	{
//		for (int j = 0; j < vfd.Nu.getSizeY(); j++)
//		{
//			fprintf(hFile, "%g\t", vfd.Nu(i, j));
//		}
//		fprintf(hFile, "\n");
//
//	}
//
//	fclose(hFile);
//
//	return;
//}
//
//
///// Not updated and checked for errors
//void clProblem::SaveSS()
//{
//	FILE *hFile;
//	char FileName[] = STRESS_OUTPUT_FILE;
//	int maxAttempt = 100;
//	for (int i = 0; i < maxAttempt; i++)
//	{
//		hFile = fopen(FileName, "a");
//		if (hFile != NULL)
//			break;
//		printf("cannot open output file %s, waiting", FileName);
//		delay_func(0.1);
//	}
//
//	if (hFile == NULL)
//	{
//		printf("cannot open output file %s, data lost", FileName);
//		return;
//	}
//	clFinish(*clEnv::instance()->getTRqueue());
//	if (vtr.Save_loc_SS == 0)
//		return;
//	vtr.SS_output.read_from_buffer_size(vtr.Save_loc_SS*vtr.numbl_bounds);
//
//	for (int i = 0; i < vtr.Save_loc_SS; i++)
//	{
//		for (int j = 0; j < vtr.numbl_bounds; j++)
//		{
//			fprintf(hFile, "%g\t", vtr.SS_output(i, j));
//		}
//		fprintf(hFile, "\n");
//	}
//
//	fclose(hFile);
//
//	return;
//}
//
///// Not updated and checked for errors
//void clProblem::SaveIO()
//{
//	FILE *hFile;
//	char FileName[] = IO_OUTPUT_FILE;
//	int maxAttempt = 100;
//	for (int i = 0; i < maxAttempt; i++)
//	{
//		hFile = fopen(FileName, "a");
//		if (hFile != NULL)
//			break;
//		printf("cannot open output file %s, waiting", FileName);
//		delay_func(0.1);
//	}
//
//	if (hFile == NULL)
//	{
//		printf("cannot open output file %s, data lost", FileName);
//		return;
//	}
//	clFinish(*clEnv::instance()->getTRqueue());
//	if (vtr.Save_IO_Loc == 0)
//		return;
//	vtr.IO_dists_save.read_from_buffer_size(vtr.Save_IO_Loc * 2 * vtr.Nd);
//
//	for (int i = 0; i < vtr.Save_IO_Loc; i++)
//	{
//		for (int j = 0; j < vtr.Nd * 2; j++)
//		{
//			fprintf(hFile, "%d\t", vtr.IO_dists_save(i * 2 * vtr.Nd + j));
//		}
//		fprintf(hFile, "\n");
//	}
//
//	fclose(hFile);
//
//	return;
//}

/// Not updated and checked for errors
//void clProblem::SaveAvgs()
//{
//	//FILE *hFile;
//	//char FileName[] = AVG_OUTPUT_FILE;
//	//int maxAttempt = 100;
//	//for (int i = 0; i < maxAttempt; i++)
//	//{
//	//	hFile = fopen(FileName, "a");
//	//	if (hFile != NULL)
//	//		break;
//	//	printf("cannot open output file %s, waiting", FileName);
//	//	delay_func(0.1);
//	//}
//
//	//if (hFile == NULL)
//	//{
//	//	printf("cannot open output file %s, data lost", FileName);
//	//	curLine = 0;
//	//	return;
//	//}
//	//clFinish(*clEnv::instance()->getLBqueue());
//	//if (vlb.Save_loc == 0)
//	//	return;
//	//vlb.RoJout.read_from_buffer_size(vlb.Save_loc * 5);
//	//vtr.Umean_Current = vlb.RoJout((vlb.Save_loc - 1) * 5 + 1) / vls.Masses(0);
//	//for (int i = 0; i < vlb.Save_loc; i++)
//	//{
//	//	double Umean = vlb.RoJout[i * 5 + 1] / vls.Masses(0);
//	//	double Tmin = vlb.RoJout[i * 5 + 2];
//	//	double Tmout1 = vlb.RoJout[i * 5 + 3];
//	//	double Tmout2 = vlb.RoJout[i * 5 + 4];
//
//	//	double Nuss1 = 4.*vlb.Pipe_radius*vlb.Pipe_radius*Umean / (double)(X_STOP_POS - X_RELEASE_POS) / vfd.Alpha_fluid*log(Tmin / Tmout1);
//	//	double Nuss2 = 4.*vlb.Pipe_radius*vlb.Pipe_radius*Umean / (double)(CHANNEL_LENGTH - X_RELEASE_POS) / vfd.Alpha_fluid*log(Tmin / Tmout2);
//
//	//	double Fluid_tot = vls.Masses(0);
//	//	double Foul_tot = vls.Masses(1);
//	//	double Ro = vlb.RoJout[i * 5 + 0];
//	//	double Roavg = Ro / vls.Masses(0);
//	//	fprintf(hFile, "%g\t", Fluid_tot);
//	//	fprintf(hFile, "%g\t", Foul_tot);
//	//	fprintf(hFile, "%g\t", Ro);
//	//	fprintf(hFile, "%g\t", Roavg);
//	//	fprintf(hFile, "%g\t", vlb.RoJout(i * 5 + 1));
//	//	fprintf(hFile, "%g\t", Umean);
//	//	fprintf(hFile, "%g\t", Tmin);
//	//	fprintf(hFile, "%g\t", Tmout1);
//	//	fprintf(hFile, "%g\t", Tmout2);
//	//	fprintf(hFile, "%g\t", Nuss1);
//	//	fprintf(hFile, "%g\n", Nuss2);
//	//}
//
//	//fclose(hFile);
//
//	return;
//}
//
////


//void clProblem::updateAvgs()
//{
//	if (Time < saveAvgStartTime)
//		return;
//	vlb.updateOutputs();
//	if (vlb.Save_loc == OUTPUT_MAX_LINES)
//	{
//		p.SaveAvgs();
//		vlb.resetOutputs();
//	}
//}
//
//void clProblem::updateNu()
//{
//	if (Time < saveNuStartTime)
//		return;
//	vfd.UpdateNu();
//	if (vfd.Save_loc_Nu == OUTPUT_MAX_LINES_NU)
//	{
//		p.SaveNu();
//		vfd.reset_Nu_out();
//	}
//}
//
//void clProblem::updateSS()
//{
//	if (Time < saveSSStartTime)
//		return;
//	vtr.UpdateSS();
//	if (vtr.Save_loc_SS == OUTPUT_MAX_LINES_SS)
//	{
//		p.SaveSS();
//		vtr.reset_SS_out();
//	}
//}
//
//void clProblem::updateIODist()
//{
//	if (Time < saveIOStartTime)
//		return;
//	vtr.Update_IO_Dists();
//	if (vtr.Save_IO_Loc == OUTPUT_MAX_LINES_IO)
//	{
//		p.SaveIO();
//		vtr.Reset_IO_Dists();
//	}
//}