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
  /// title
  std::string title{};

  /// rows of data (number of columns may vary)
  std::vector<std::vector<std::string>> rows{};

  /// Create table with given title
  Table(const std::string &title) : title(title) {}

  /// Pretty-print
  std::string __str__() const override
  {
    // Padding between columns
    const size_t h = 2;

    // Compute number of columns
    size_t num_cols{0};
    for (const auto &row : rows)
      num_cols = std::max(num_cols, row.size());

    // Return if no data
    if (num_cols == 0)
      return "";

    // Compute column sizes
    std::vector<size_t> col_sizes(num_cols);
    std::fill(col_sizes.begin(), col_sizes.end(), 0);
    for (size_t i = 0; i < rows.size(); i++)
      for (size_t j = 0; j < rows[i].size(); j++)
        col_sizes[j] = std::max(col_sizes[j], rows[i][j].size() + h);

    // Compute table width
    size_t width{0};
    for (size_t j = 0; j < num_cols; j++)
      width += col_sizes[j];
    width -= h;

    // print title
    std::string s{};
    s += "\n";
    for (size_t j = 0; j < width; j++)
      s += "=";
    s += "\n";
    s += title + "\n";
    for (size_t j = 0; j < width; j++)
      s += "=";
    s += "\n";

    // print rows
    for (size_t i = 0; i < rows.size(); i++)
    {
      for (size_t j = 0; j < rows[i].size(); j++)
      {
        s += rows[i][j];
        if (j + 1 < rows[i].size())
        {
          const size_t padding = col_sizes[j] - rows[i][j].size();
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
