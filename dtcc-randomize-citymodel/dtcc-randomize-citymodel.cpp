// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <string>
#include <vector>

#include "CityModelGenerator.h"
#include "CommandLine.h"
#include "ElevationModelGenerator.h"
#include "JSON.h"
#include "Logging.h"
#include "Parameters.h"
#include "VTK.h"

using namespace DTCC;

void Help() { Error("Usage: dtcc-randomize-citymodel Parameters.json"); }

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Read parameters
  Parameters p;
  JSON::Read(p, argv[1]);
  Info(p);
  const std::string modelName = Utils::GetFilename(argv[1], true);

  // Set data directory
  const std::string dataDirectory{p.DataDirectory + "/"};

  // Set bounding box
  Point2D O{0.0, 0.0};
  const Point2D P{p.XMin, p.YMin};
  const Point2D Q{p.XMax, p.YMax};
  BoundingBox2D bbox{P, Q};

  // Check size of bounding box
  Info("Bounding box of city model: " + str(bbox));
  if (bbox.Area() < 100.0)
  {
    Error("Domain too small to generate a city model");
    return 1;
  }

  // Randomize elevation model
  GridField2D dtm;
  ElevationModelGenerator::RandomizeElevationModel(dtm, bbox,
                                                   p.ElevationModelResolution);

  // Randomize city model
  CityModel cityModel;
  cityModel.Name = modelName;
  CityModelGenerator::RandomizeCityModel(cityModel, dtm, p.NumRandomBuildings);
  Info(cityModel);

  // Clean city model
  CityModelGenerator::CleanCityModel(cityModel, p.MinVertexDistance);

  // Write JSON
  if (p.WriteJSON)
  {
    JSON::Write(dtm, dataDirectory + "DTM.json", O);
    JSON::Write(cityModel, dataDirectory + "CityModel.json", O);
  }

  // Write VTK
  if (p.WriteVTK)
  {
    VTK::Write(dtm, dataDirectory + "DTM.vts");
  }

  return 0;
}
