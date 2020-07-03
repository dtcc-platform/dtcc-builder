// Copyright (C) 2020 ReSpace AB
// Licensed under the MIT License

#ifndef DTCC_LOGGING_H
#define DTCC_LOGGING_H

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace DTCC
{
  // Log levels
  enum LogLevels
  {
    PROGRESS = 0,
    INFO = 1,
    WARNING = 2,
    ERROR = 3
  };

  // Global log level
  LogLevels LogLevel = INFO;

  // Global log file
  std::ofstream logFile;

  // Interface for printable objects
  class Printable
  {
  public:
    virtual std::string __str__() const = 0;
  };

  // Print message to stdout
  void __print__(const std::string& message)
  {
    std::cout << message << std::endl;
  }

  // Print message
  void __print_logfile__(const std::string& message, bool closeLogFile=false)
  {
    if (logFile)
    {
      logFile << message << std::endl;
      if (closeLogFile) logFile.close();
    }
  }

  /// Return current time (as a string)
  ///
  /// @return Current time
  std::string CurrentTime()
  {
    time_t tt{};
    time(&tt);
    struct tm *ti = localtime(&tt);
    return asctime(ti);
  }

  // Print information message
  void Info(const std::string& message)
  {
    if (LogLevel <= INFO)
      __print__(message);
    __print_logfile__(message);
  }

  // Print information about printable object
  void Info(const Printable &printable)
  {
    Info(printable.__str__());
  }

  // Print progress message
  void Progress(const std::string& message)
  {
    if (LogLevel <= PROGRESS)
      __print__(message);
    __print_logfile__(message);
  }

  // Print warning message
  void Warning(const std::string& message)
  {
    if (LogLevel <= WARNING)
      __print__(message);
    __print_logfile__(message);
  }

  // Print error message and throw exception
  void Error(const std::string& message)
  {
    if (LogLevel <= ERROR)
      __print__(message);
    __print_logfile__(message, true);
    throw std::runtime_error(message);
  }

  // Set log level
  void SetLogLevel(LogLevels logLevel)
  {
    LogLevel = logLevel;
  }

  // Set log file
  void SetLogFile(std::string fileName, bool append=false)
  {
    if (append)
      logFile.open(fileName, std::ofstream::out | std::ofstream::app);
    else
      logFile.open(fileName, std::ofstream::out);
    if (!logFile)
      Error("Unable to write to logfile " + fileName);
    logFile << std::endl << "Time: " << CurrentTime();
  }

  // Convert printable object to string
  std::string str(const Printable &x) { return x.__str__(); }

  // FIXME: Functions below should be possible to handle with a common
  // template but does not seem to work (overloaded str for Printable
  // seems to be hidden).

  // Convert const char* to string
  std::string str(const char *x)
  {
    std::string s(x);
    return s;
  }

  // Convert unsigned integer to string
  std::string str(int x) { return std::to_string(x); }

  // Convert unsigned integer to string
  std::string str(size_t x) { return std::to_string(x); }

  // Convert double to string
  std::string str(double x)
  {
    std::ostringstream out;
    out.precision(6);
    out << std::scientific << x;
    return out.str();
  }

} // namespace DTCC

#endif
