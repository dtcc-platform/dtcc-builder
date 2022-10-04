// Copyright (C) 2020 ReSpace AB
// Licensed under the MIT License
//
// Simple logging system for use in DTCC C++ code modeled
// after the standard Python logging module.

#ifndef DTCC_LOGGING_H
#define DTCC_LOGGING_H

#include <ctime>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace DTCC
{
  // Log levels
  enum LogLevels
  {
    DEBUG = 10,
    INFO = 20,
    WARNING = 30,
    ERROR = 40,
    PROGRESS = 25
  };

  // Global log level
  LogLevels __log_level__ = INFO;

  // Interface for printable objects
  class Printable
  {
  public:

    // Convert to string (pretty-print)
    virtual std::string __str__() const = 0;

    // Conversion operator
    operator std::string() const { return __str__(); }

  };

  // Format message
  std::string __format__(const std::string& message)
  {
    return message;
  }

  // Print message to stdout
  void __print__(const std::string& message)
  {
    std::cout << __format__(message) << std::endl;
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

  //--- Public interface of logging system ---

  // Set log level
  void set_log_level(LogLevels log_level)
  {
    __log_level__ = log_level;
  }

  // Print message at given log level
  void log(int log_level, const std::string &message = "")
  {
    if (log_level >= __log_level__)
      __print__(message);
  }

  // Print debug message
  void Debug(const std::string& message)
  {
    log(DEBUG, message);
  }

  // Print information message (string)
  void Info(const std::string &message = "")
  {
    log(INFO, message);
  }

  // Print warning message
  void Warning(const std::string& message)
  {
    log(WARNING, message);
  }

  // Print error message and throw exception
  void Error(const std::string& message)
  {
    log(ERROR, message);
    throw std::runtime_error(message);
  }

  //--- Utility functions for string conversion ---

  // Convert printable object to string
  std::string str(const Printable &x) { return x.__str__(); }

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

  // Convert unsigned integer to string
  std::string str(uint x) {return std::to_string(x); }

  // Convert double to string
  std::string str(double x, std::streamsize precision = 6)
  {
    std::ostringstream out;
    out.precision(precision);
    out << std::defaultfloat << x;
    return out.str();
  }

} // namespace DTCC

#endif
