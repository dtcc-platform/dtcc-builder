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

namespace DTCC
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

  static std::pair<std::string, std::string> GetDataParameters(int argc,
                                                               char *argv[])
  {
    std::string dataDirectory = "";
    std::string parameterFile = "";
    if (argc == 1)
    {
      // no arguments
      dataDirectory = std::experimental::filesystem::current_path().string();
    }
    else if (argc == 2)
    {
      std::string arg = argv[1];

      if (Filesystem::IsDirectory(arg))
      {
        dataDirectory = arg;
      }
      else if (Filesystem::IsFile(arg))
      {
        parameterFile = arg;
      }
    }
    if (dataDirectory.size() > 0 && dataDirectory.back() != '/')
      dataDirectory += "/";
    if (parameterFile.size() == 0)
    {
      if (Filesystem::IsFile(dataDirectory + "Parameters.json"))
      {
        parameterFile = dataDirectory + "Parameters.json";
      }
    }
    return std::make_pair(dataDirectory, parameterFile);
  }
};

} // namespace DTCC

#endif
