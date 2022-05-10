// Dag Wästberg 2022
// Licensed under the MIT License

#ifndef DTCC_PARAMETE_PROCESSOR_H
#define DTCC_PARAMETE_PROCESSOR_H

#include <experimental/filesystem>

#include "Parameters.h"
namespace DTCC
{
class ParameterProcessor
{
public:
  /// get (and create if necessary) the data and output path
  static std::pair<std::string, std::string>
  getDataAndOutputPath(const Parameters &p)
  {
    // Get data directory
    std::string dataDirectory = p["DataDirectory"];
    if (dataDirectory.empty())
      dataDirectory = "./";

    if (dataDirectory.back() != '/')
      dataDirectory += '/';

    std::string outputDirectory = p["OutputDirectory"];
    if (outputDirectory.empty())
      outputDirectory = dataDirectory;
    if (outputDirectory.back() != '/')
      outputDirectory += "/";
    std::experimental::filesystem::create_directories(outputDirectory);
    return std::pair<std::string, std::string>(dataDirectory, outputDirectory);
  }
};
} // namespace DTCC

#endif