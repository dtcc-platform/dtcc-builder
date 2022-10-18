// Dag Wästberg 2022
// Licensed under the MIT License

#ifndef DTCC_FILESYSTEM_H
#define DTCC_FILESYSTEM_H

#include <experimental/filesystem>
#include <string>
#include <vector>

#include "Utils.h"

namespace DTCC_BUILDER
{
class Filesystem
{
public:
  static bool IsDirectory(const std::string &path)
  {
    return (std::experimental::filesystem::exists(path) &&
            std::experimental::filesystem::is_directory(path));
  }

  static bool IsFile(const std::string &path)
  {
    return (std::experimental::filesystem::exists(path) &&
            std::experimental::filesystem::is_regular_file(path));
  }

  static std::vector<std::string>
  ListDirectory(const std::string &directory, const std::string &extension = "")
  {
    std::vector<std::string> fileNames;
    bool extension_filter = extension.size() > 0;

    for (auto &p : std::experimental::filesystem::directory_iterator(directory))
    {
      auto path_string = p.path().string();
      if (extension_filter)
      {
        if (Utils::EndsWith(path_string, extension))
          fileNames.push_back(path_string);
      }
      else
      {
        fileNames.push_back(path_string);
      }
    }
    return fileNames;
  }
};
} // namespace DTCC_BUILDER

#endif