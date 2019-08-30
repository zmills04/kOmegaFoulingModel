// Logger singleton class for logging warnings, since this code has
// been written to be executed on a windows machine, and I'm not sure
// how to pipe that to a file.
// (c) Zachary Grant Mills, 2019 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LOGGER_H__INCLUDED_)
#define AFX_LOGGER_H__INCLUDED_

#pragma once


#include "StdAfx.h"
#include <fstream>

class Logger
{
public:

	static Logger* Instance()
	{
		if (!m_pInstance)   // Only allow one instance of class to be generated.
			m_pInstance = new Logger();

		return m_pInstance;
	}
		
	//static Logger* SetupInstance(bool restartFlag);
	void iniLogger(bool restartFlag);
	void writeToLogFile(std::string logMessage);
	void closeLogFile(bool finishedSim = true);
	void writePartialToLogFile(std::string logMessage);
	//void iniLogFile(bool restartFlag)

private:
	std::string logName = "logFile.txt";
	std::fstream logStream;
	Logger() 
	{
		//openLogFile(restartFlag);
	};
	~Logger()
	{
		closeLogFile(false);
	}

	Logger(Logger const&) {};             // copy constructor is private
	Logger& operator=(Logger const&) {};  // assignment operator is private
	static Logger* m_pInstance;
};









#endif //AFX_LOGGER_H__INCLUDED_