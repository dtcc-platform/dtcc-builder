// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_ROAD_H
#define DTCC_ROAD_H

#include <vector>
#include <unordered_map>

#include "Point.h"

namespace DTCC
{
class RoadNetwork : public Printable
{
public:

  /** Holds the network's vertex positions. */
  std::vector<Point2D> Vertices;

  /** Holds the vertex indices indicating the road edges. */
  std::vector<std::pair<size_t, size_t>> Edges;

  /** Holds vectors of additional vertex values.  */
  std::unordered_map<std::string, std::vector<double>> VertexValues;

  /** Holds vectors of additional edge values. */
  std::unordered_map<std::string, std::vector<double>> EdgeValues;

  std::unordered_map<std::string, std::vector<std::string>> VertexProperties;

  std::unordered_map<std::string, std::vector<std::string>> EdgeProperties;

  std::string __str__() const override
  {
    return "RoadNetwork with " + str(Edges.size()) + " edges, " +
           str(VertexValues.size()) + " vertex value arrays, and " +
           str(EdgeValues.size()) + " edge value arrays";
  }
  };
} // namespace DTCC

#endif // DTCC_ROAD_H
