#include "logger.h"
#include <chrono>
#include <ctime>  

// Global static pointer used to ensure a single instance of the class.
Logger* Logger::m_pInstance = nullptr;

Logger* Logger::SetupInstance(bool restartFlag)
{
	if (!m_pInstance)   // Only allow one instance of class to be generated.
		m_pInstance = new Logger(restartFlag);

	return m_pInstance;
}

void Logger::openLogFile(bool restartFlag)
{
	using std::chrono::system_clock;
		
	if(!restartFlag)
		logStream.open(logName, std::ios_base::out);
	else
		logStream.open(logName, std::ios_base::app);

	ERROR_CHECKING(!logStream.is_open(), "Unable to open logfile",
		ERROR_IN_LOGGER);

	auto timenow = system_clock::now();
	std::time_t timenowformatted = system_clock::to_time_t(timenow);

	if (!restartFlag)
	{
		logStream << "Started new simulation at " 
			<< std::ctime(&timenowformatted) << "\n";
	}
	else
	{
		logStream << "Restarted simulation at "
			<< std::ctime(&timenowformatted) << "\n";
	}

	logStream.flush();
}


// This should be called directly at end of simulation 
// with finishedSim = true to show simulation was
// successfully completed. If not, the destructor
// will call with finishedSim = false to write that
// simulation was not completed.
void Logger::closeLogFile(bool finishedSim)
{
	using std::chrono::system_clock;

	// if log file already closed closeLogFile was called
	// directly and no need to close file using destructor
	if (!logStream.is_open()) { return; }

	auto timenow = system_clock::now();
	std::time_t timenowformatted = system_clock::to_time_t(timenow);

	if(finishedSim)
	{
		logStream << "Simulation finished at "
			<< std::ctime(&timenowformatted) << "\n";
	}
	else
	{
		logStream << "Simulation ended without finishing at "
			<< std::ctime(&timenowformatted) << "\n";
	}

	logStream.close();
}

void Logger::writeToLogFile(std::string logMessage)
{
	ERROR_CHECKING(!logStream.is_open(), "Trying to write to "
		"unopened log file", ERROR_IN_LOGGER);

	logStream << logMessage;
	logStream << "\n";
	logStream.flush();
}



