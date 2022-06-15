// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_PARAMETERS_H
#define DTCC_PARAMETERS_H

#include <map>
#include <string>

#include "Logging.h"
#include "Parameter.h"
#include "Table.h"

namespace DTCC
{

/// Parameters holds a database of parameters (key-value pairs).

class Parameters : public Printable
{
public:
  /// Set default parameters (and clear all previous parameters)
  void SetDefaultParameters()
  {
    Map.clear();

    Add("ModelName", "");
    Add("AutoDomain", true);
    Add("GenerateSurfaceMeshes", true);
    Add("GenerateVolumeMeshes", true);
    Add("WriteJSON", true);
    Add("WriteVTK", true);
    Add("Debug", false);

    Add("GroundSmoothing", 5);
    Add("NumRandomBuildings", 25);

    Add("DomainMargin", 10.0);
    Add("X0", 0.0);
    Add("Y0", 0.0);
    Add("XMin", 0.0);
    Add("YMin", 0.0);
    Add("XMax", 0.0);
    Add("YMax", 0.0);
    Add("ElevationModelResolution", 1.0);
    Add("MinBuildingDistance", 1.0);
    Add("MinVertexDistance", 1.0);
    Add("GroundMargin", 1.0);
    Add("MeshResolution", 10.0);
    Add("DomainHeight", 100.0);
    Add("GroundPercentile", 0.1);
    Add("RoofPercentile", 0.9);
    Add("OutlierMargin", 2.0);
    Add("MinBuildingSize", 15.0);
    Add("BuildingLOD", 1);

    Add("StatisticalOutlierRemover", true);
    Add("OutlierNeighbors", 5);
    Add("OutlierSTD", 1.5);

    Add("RANSACOutlierRemover", true);
    Add("RANSACOutlierMargin", 3.0);
    Add("RANSACIterations", 250);

    Add("NaiveVegitationFilter", true);
    Add("DataDirectory", "");
    Add("OutputDirectory", "");
  }

  // Map (dictionary) of parameters
  std::map<std::string, Parameter> Map;

  /// Create default parameter set
  explicit Parameters() { SetDefaultParameters(); }

  /// Access parameter
  Parameter &operator[](const std::string &key)
  {
    auto it = Map.find(key);
    if (it == Map.end())
      Error("Unknown parameter: \"" + key + "\"");
    return it->second;
  }

  /// Access parameter (const version)
  const Parameter &operator[](const std::string &key) const
  {
    auto it = Map.find(key);
    if (it == Map.end())
      Error("Unknown parameter: \"" + key + "\"");
    return it->second;
  }

  /// Check whether parameter with given key exists
  bool HasKey(const std::string &key) const
  {
    return Map.find(key) != Map.end();
  }

  /// Add bool parameter
  void Add(const std::string &key, bool value)
  {
    if (HasKey(key))
      Error("Unable to add parameter; key \"" + key + "\" already exists");
    Parameter p(ParameterType::Bool, key);
    p = value;
    Map[key] = p;
  }

  /// Add int parameter
  void Add(const std::string &key, int value)
  {
    if (HasKey(key))
      Error("Unable to add parameter; key \"" + key + "\" already exists");
    Parameter p(ParameterType::Int, key);
    p = value;
    Map[key] = p;
  }

  /// Add float parameter
  void Add(const std::string &key, double value)
  {
    if (HasKey(key))
      Error("Unable to add parameter; key \"" + key + "\" already exists");
    Parameter p(ParameterType::Float, key);
    p = value;
    Map[key] = p;
  }

  /// Add string parameter
  void Add(const std::string &key, const std::string &value)
  {
    if (HasKey(key))
      Error("Unable to add parameter; key \"" + key + "\" already exists");
    Parameter p(ParameterType::String, key);
    p = value;
    Map[key] = p;
  }

  /// Add string parameter (handle string literals, will otherwise go to bool)
  void Add(const std::string &key, const char *value)
  {
    if (HasKey(key))
      Error("Unable to add parameter; key \"" + key + "\" already exists");
    Parameter p(ParameterType::String, key);
    p = std::string(value);
    Map[key] = p;
  }

  // FIXME: Consider making the following proper parameters

  // Tolerance for geometric tests
  static constexpr double Epsilon = 1e-6;

  // Precision for output and printing
  static constexpr double Precision = 16;

  // Threshold for filtering points with small angles in building footprints
  static constexpr double FootprintAngleThreshold = 0.01;

  // Threshold for filtering outliers (clouds?) from point cloud
  static constexpr double PointCloudOutlierThreshold = 150.0;

  // Number of digits of precision used when writing files
  static constexpr double OutputPrecision = 3;

  /// Pretty-print
  std::string __str__() const override
  {
    // Build table
    Table table("Parameters");
    table.Rows.push_back({"Key", "Value", "Type", "Access"});
    for (const auto &it : Map)
    {
      std::vector<std::string> row;
      row.push_back(it.first);
      row.push_back(it.second.ValueString());
      row.push_back(it.second.TypeString());
      row.push_back(str(it.second.AccessCount));
      table.Rows.push_back(row);
    }

    // Return table string
    return str(table);
  }
};

} // namespace DTCC

#endif
