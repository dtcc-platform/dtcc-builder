// Anders Logg 2019
// Licensed under the MIT License

#ifndef DTCC_COMMAND_LINE_H
#define DTCC_COMMAND_LINE_H

#include <dirent.h>
#include <experimental/filesystem>
#include <string>
#include <sys/types.h>
#include <vector>

#include "Filesystem.h"
#include "Utils.h"

namespace DTCC_BUILDER
{

/// Simple utilities for command-line parsing
class CommandLine
{
public:
  static bool HasOption(const std::string& option, int argc, char *argv[])
  {
    for (int i = 1; i < argc; i++)
      if (option == argv[i])
        return true;
    return false;
  }

  static std::string GetOption(const std::string& option, int argc, char *argv[])
  {
    for (int i = 1; i < argc; i++)
      if (option == argv[i])
        return argv[i + 1];
    return "";
  }

  static int GetIntOption(const std::string& option, int argc, char *argv[])
  {
    return atoi(GetOption(option, argc, argv).c_str());
  }
};

} // namespace DTCC_BUILDER

#endif
