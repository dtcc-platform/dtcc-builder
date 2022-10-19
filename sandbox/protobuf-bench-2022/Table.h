// Copyright (C) 2021 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_TABLE_H
#define DTCC_TABLE_H

#include <string>
#include <vector>

#include "Logging.h"

namespace DTCC_BUILDER
{

/// Table is used for easy printing of tabular data.

class Table : public Printable
{
public:
  /// Title
  std::string Title{};

  /// Rows of data (number of columns may vary)
  std::vector<std::vector<std::string>> Rows{};

  /// Create table with given title
  Table(const std::string &title) : Title(title) {}

  /// Pretty-print
  std::string __str__() const override
  {
    // Padding between columns
    const size_t h = 2;

    // Compute number of columns
    size_t numCols{0};
    for (const auto &row : Rows)
      numCols = std::max(numCols, row.size());

    // Return if no data
    if (numCols == 0)
      return "";

    // Compute column sizes
    std::vector<size_t> colSizes(numCols);
    std::fill(colSizes.begin(), colSizes.end(), 0);
    for (size_t i = 0; i < Rows.size(); i++)
      for (size_t j = 0; j < Rows[i].size(); j++)
        colSizes[j] = std::max(colSizes[j], Rows[i][j].size() + h);

    // Compute table width
    size_t width{0};
    for (size_t j = 0; j < numCols; j++)
      width += colSizes[j];
    width -= h;

    // Print title
    std::string s{};
    s += "\n";
    for (size_t j = 0; j < width; j++)
      s += "=";
    s += "\n";
    s += Title + "\n";
    for (size_t j = 0; j < width; j++)
      s += "=";
    s += "\n";

    // Print rows
    for (size_t i = 0; i < Rows.size(); i++)
    {
      for (size_t j = 0; j < Rows[i].size(); j++)
      {
        s += Rows[i][j];
        if (j + 1 < Rows[i].size())
        {
          const size_t padding = colSizes[j] - Rows[i][j].size();
          for (size_t k = 0; k < padding; k++)
            s += " ";
        }
      }
      s += "\n";
      if (i == 0)
      {
        for (size_t j = 0; j < width; j++)
          s += "-";
        s += "\n";
      }
    }
    for (size_t j = 0; j < width; j++)
      s += "=";

    return s;
  }
};

} // namespace DTCC_BUILDER

#endif
