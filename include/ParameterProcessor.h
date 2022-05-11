// Dag Wästberg 2022
// Licensed under the MIT License

#ifndef DTCC_PARAMETE_PROCESSOR_H
#define DTCC_PARAMETE_PROCESSOR_H

#include <experimental/filesystem>
#include <stdlib.h>

#include "Logging.h"
#include "Parameters.h"
namespace DTCC
{
class ParameterProcessor
{
public:
  static Parameters ProcessArgs(int argc, char *argv[])
  {
    Parameters parameters;
    std::string dataDirectory = "";
    std::string parameterFile = "";
    std::string outputDirectory = "";
    if (argc == 1)
    {
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
      else
      {
        Error("Unknown argument: " + arg);
        exit(1);
      }
    }
    if (dataDirectory.size() > 0 && dataDirectory.back() != '/')
    {
      dataDirectory += "/";
    }
    if (parameterFile.empty())
    {
      if (Filesystem::IsFile(dataDirectory + "Parameters.json"))
      {
        parameterFile = dataDirectory + "Parameters.json";
      }
    }

    if (!parameterFile.empty())
    {
      JSON::Read(parameters, parameterFile);
    }
    if (!dataDirectory.empty())
    {
      parameters["DataDirectory"] = dataDirectory;
    }
    outputDirectory = (std::string)parameters["OutputDirectory"];
    if (outputDirectory.empty())
      outputDirectory = dataDirectory;
    if (outputDirectory.back() != '/')
      outputDirectory += "/";
    std::experimental::filesystem::create_directories(outputDirectory);
    parameters["OutputDirectory"] = outputDirectory;

    return parameters;
  }
};
} // namespace DTCC

#endif