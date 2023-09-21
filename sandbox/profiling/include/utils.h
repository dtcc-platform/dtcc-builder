
#ifndef SNDBX_UTILS_H
#define SNDBX_UTILS_H

#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>
#include <vector>

#include "model/City.h"
#include "model/GridField.h"
#include "model/VolumeMesh.h"

namespace DTCC_BUILDER
{

class JSON
{
public:
  // Utility functions for JSON input

  /// Read JSON data from file
  static void Read(nlohmann::json &json, std::string fileName)
  {
    info("JSON: Reading from file " + fileName + "...");
    std::ifstream f(fileName);
    if (!f)
      error("Unable to read from file " + fileName);
    f >> json;
  }

  /// Read object from file
  template <class T> static void Read(T &t, std::string fileName)
  {
    nlohmann::json json{};
    Read(json, fileName);
    Deserialize(t, json);
  }

private:
  /// Deserialize City
  // Bounds:
  // {"xmax":102999.99,"xmin":102000.0,"ymax":6214000.0,"ymin":6213004.0}
  // georef: {"crs":"epsg:3008","epsg":0,"x0":0.0,"y0":0.0}
  // landuse: []
  //  - Name (string)
  //  - Buildings (std::vector<Building>)
  //    |- AttachedUUIDs (std::vector<std::string>)
  //    |- UUID (std::string)
  //    |- PropertyFNR (size_t)
  //    |- PropertyUUID (std::string)
  //    |- Footprint (Polygon)
  //    |- BaseAreaID (Polygon)
  //    |- Height (double)
  //    |- error (size_t)
  //    |- GroundHeight (double)
  //    |- GroundPoints (std::vector<Vector3D>)
  //    |- RoofPoints (std::vector<Vector3D>)
  //    |- RoofSegments (std::vector<std::vector<size_t>>)
  //  - Origin (Vector2D)
  static void Deserialize(City &city, const nlohmann::json &json)
  {
    // CheckType("City", json);
    //  name
    city.Name = json["name"];

    // Buildings
    auto jsonBuildings = json["buildings"];
    city.Buildings.resize(jsonBuildings.size());
    for (size_t i = 0; i < jsonBuildings.size(); i++)
    {
      auto jsonBuilding = jsonBuildings[i];

      std::string error_string = jsonBuilding["error"];
      city.Buildings[i].error = std::stoi(error_string);
      city.Buildings[i].Height = jsonBuilding["height"];
      city.Buildings[i].GroundHeight = jsonBuilding["groundHeight"];
      city.Buildings[i].UUID = jsonBuilding["uuid"];
      // city.Buildings[i].SHPFileID = jsonBuilding["SHPFileID"];

      auto jsonFootprint = jsonBuilding["footprint"]["shell"]["vertices"];
      city.Buildings[i].Footprint.Vertices.resize(jsonFootprint.size());
      for (size_t j = 0; j < jsonFootprint.size(); j++)
      {
        city.Buildings[i].Footprint.Vertices[j].x = jsonFootprint[j]["x"];
        city.Buildings[i].Footprint.Vertices[j].y = jsonFootprint[j]["y"];
      }

      auto jsonRoofPoints = jsonBuilding["roofpoints"]["points"];
      size_t num_of_roofpoints = jsonRoofPoints.size() / 3;
      city.Buildings[i].RoofPoints.resize(num_of_roofpoints);
      for (size_t j = 0; j < num_of_roofpoints; j++)
      {
        city.Buildings[i].RoofPoints[j].x = jsonRoofPoints[3 * j];
        city.Buildings[i].RoofPoints[j].y = jsonRoofPoints[3 * j + 1];
        city.Buildings[i].RoofPoints[j].z = jsonRoofPoints[3 * j + 2];
      }
    }
  }

  static void Deserialize(GridField &gridfield, const nlohmann::json &json)
  {
    auto bounds = json["bounds"];

    double px = bounds["xmin"];
    double py = bounds["ymin"];
    double qx = bounds["xmax"];
    double qy = bounds["ymax"];
    auto bbox = BoundingBox2D(Vector2D(px, py), Vector2D(qx, qy));
    std::cout << "BBOX P x: " << bbox.P.x << " y: " << bbox.P.y << std::endl;
    std::cout << "BBOX Q x: " << bbox.Q.x << " y: " << bbox.Q.y << std::endl;

    auto terrain = json["terrain"];
    size_t XSize = terrain["width"];
    size_t YSize = terrain["height"];
    std::cout << "Terrain Width: " << XSize << " Terrain Height: " << YSize
              << std::endl;
    std::cout << "Terrain transform:" << terrain["transform"] << std::endl;

    Grid grid(bbox, XSize, YSize);
    gridfield.grid = grid;

    auto terrain_values = terrain["values"];
    size_t num_terrain_values = static_cast<size_t>(terrain_values.size());
    gridfield.Values.resize(num_terrain_values);
    for (size_t i = 0; i < num_terrain_values; i++)
    {
      gridfield.Values[i] = static_cast<double>(terrain_values[i]);
    }
  }

  static void Deserialize(VolumeMesh &volume_mesh, const nlohmann::json &json)
  {

    auto vertices = json["vertices"];
    volume_mesh.Vertices.reserve(vertices.size());
    for (size_t i = 0; i < vertices.size(); i = i + 3)
    {
      Vector3D pt(vertices[i], vertices[i + 1], vertices[i + 2]);
      volume_mesh.Vertices.push_back(pt);
    }

    auto cells = json["cells"];
    volume_mesh.Cells.reserve(cells.size());
    for (size_t i = 0; i < cells.size(); i = i + 4)
    {
      size_t v0 = cells[i];
      size_t v1 = cells[i + 1];
      size_t v2 = cells[i + 2];
      size_t v3 = cells[i + 3];
      Simplex3D cell(v0, v1, v2, v3);

      volume_mesh.Cells.push_back(cell);
    }

    auto markers = json["markers"];

    volume_mesh.Markers.reserve(markers.size());
    for (size_t i = 0; i < markers.size(); i++)
    {
      const int marker = markers[i];
      volume_mesh.Markers.push_back(marker);
    }
    info(volume_mesh.__str__());
  }
};

class Parameters
{
public:
  std::string data_directory = "";

  std::string output_directory = "";

  uint num_of_smoothings = 0;

  std::vector<std::string> city_filenames;

  std::vector<std::string> volume_mesh_filenames;

  std::vector<size_t> max_iter;

  std::vector<double> domain_height;

  std::vector<double> rel_tol;

  std::vector<bool> fix_buildings;

  Parameters(const std::string parameters_filename)
  {
    nlohmann::json _json;
    DTCC_BUILDER::JSON::Read(_json, parameters_filename);

    this->data_directory = _json["data_directory"];
    this->output_directory = _json["output_directory"];

    auto smoothing_cases = _json["cases"];
    this->num_of_smoothings = smoothing_cases.size();

    this->city_filenames.resize(num_of_smoothings);
    this->volume_mesh_filenames.resize(num_of_smoothings);
    this->max_iter.resize(num_of_smoothings);
    this->domain_height.resize(num_of_smoothings);
    this->rel_tol.resize(num_of_smoothings);
    this->fix_buildings.resize(num_of_smoothings);

    for (size_t i = 0; i < num_of_smoothings; i++)
    {
      city_filenames[i] = smoothing_cases[i]["city"];
      city_filenames[i] = data_directory + "/" + city_filenames[i];
      
      volume_mesh_filenames[i] = smoothing_cases[i]["volume_mesh"];
      volume_mesh_filenames[i] =
          data_directory + "/" + volume_mesh_filenames[i];
      
      max_iter[i] = smoothing_cases[i]["max_iter"];
      rel_tol[i] = smoothing_cases[i]["rel_tol"];
      domain_height[i] = smoothing_cases[i]["domain_height"];
      fix_buildings[i] = smoothing_cases[i]["fix_buildings"];
    }
  }

  inline void print()
  {

    std::cout << this->data_directory << std::endl;
    std::cout << this->output_directory << std::endl;
    for (size_t i = 0; i < this->num_of_smoothings; i++)
    {
      std::cout << "Case: [" << i << "]" << std::endl;
      std::cout << "City: " << this->city_filenames[i] << std::endl;
      std::cout << "mesh: " << this->volume_mesh_filenames[i] << std::endl;
      std::cout << "max_iter: " << this->max_iter[i] << std::endl;
      std::cout << "Rel_tol: " << this->rel_tol[i] << std::endl;
      std::cout << "Domain_height: " << this->domain_height[i] << std::endl;
      std::cout << "fix Buildings: " << this->fix_buildings[i] << std::endl;
    }
  }
};

} // namespace DTCC_BUILDER
#endif