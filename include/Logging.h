// Copyright (C) 2020 ReSpace AB
// Licensed under the MIT License
//
// Simple logging system for use in DTCC C++ code modeled after the standard
// standard Python logging module with the following configuration:
//
// import logging
// format = "%(asctime)s [%(name)s] [%(levelname)s] %(message)s"
// logging.basicConfig(format=format)
// logging.addLevelName(25, 'PROGRESS')

#ifndef DTCC_LOGGING_H
#define DTCC_LOGGING_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace DTCC_BUILDER
{
  // Log levels
enum LogLevel
{
  DEBUG = 10,
  INFO = 20,
  WARNING = 30,
  ERROR = 40,
  PROGRESS = 25
};

// Global log level
LogLevel __log_level__ = INFO;

// Interface for printable objects
class Printable
{
public:
  // Convert to string (pretty-print)
  virtual std::string __str__() const = 0;

  // Conversion operator
  operator std::string() const { return __str__(); }
  };

  // Return current time
  std::string current_time()
  {
    // Stackoverflow: get-current-time-in-milliseconds-or-hhmmssmmm-format
    using namespace std::chrono;
    auto now = system_clock::now();
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
    auto timer = system_clock::to_time_t(now);
    std::tm bt = *std::localtime(&timer);
    std::ostringstream oss;
    oss << std::put_time(&bt, "%Y-%m-%d %H:%M:%S");
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    return oss.str();
  }

  // Format message
  std::string __format__(LogLevel log_level, const std::string &message)
  {
    // Set component (hard-coded for now)
    std::string component{"[dtcc-builder]"};

    // Format level
    std::string level{};
    switch (log_level)
    {
    case DEBUG:
      level = "[DEBUG]";
      break;
    case INFO:
      level = "[INFO]";
      break;
    case WARNING:
      level = "[WARNING]";
      break;
    case ERROR:
      level = "[ERROR]";
      break;
    case PROGRESS:
      level = "[PROGRESS]";
      break;
    default:
      level = "[UNKNOWN]";
    };

    return current_time() + " " + component + " " + level + " " + message;
  }

  // Print message to stdout
  void __print__(const std::string &message)
  {
    std::cout << message << std::endl;
  }

  //--- Public interface of logging system ---

  // Set log level
  void set_log_level(LogLevel log_level) { __log_level__ = log_level; }

  // Print message at given log level
  void log(LogLevel log_level, const std::string &message = "")
  {
    // Skip if below log level threshold
    if (log_level < __log_level__)
      return;

    // Format and print
    std::string formatted_message = __format__(log_level, message);
    __print__(formatted_message);
  }

  // Print debug message
  void debug(const std::string &message) { log(DEBUG, message); }

  // Print information message (string)
  void info(const std::string &message = "") { log(INFO, message); }

  // Print warning message
  void warning(const std::string &message) { log(WARNING, message); }

  // Print error message and throw exception
  void error(const std::string &message)
  {
    log(ERROR, message);
    throw std::runtime_error(message);
  }

  // Report progress (a number between 0 and 1)
  void progress(double x)
  {
    x = std::max(0.0, std::min(1.0, x));
    std::ostringstream ss{};
    ss << std::setprecision(2) << std::fixed << 100.0 * x << "%";
    log(PROGRESS, ss.str());
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
    std::ostringstream ss{};
    ss << std::setprecision(precision) << std::defaultfloat << x;
    return ss.str();
  }

  } // namespace DTCC_BUILDER

#endif
