// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_ROAD_H
#define DTCC_ROAD_H

#include <vector>
#include <unordered_map>

#include "Point.h"

namespace DTCC
{
  class Road : public Printable
  {
  public:
    /** Holds the road vertex positions. */
    std::vector<Point2D> Vertices;

    /** Holds the vertex indices indicating the road edges. */
    std::vector<unsigned int> Edges;

    /** Holds vectors of additional vertex values.  */
    std::unordered_map<std::string, std::vector<double>> VertexValues;

    /** Holds vectors of additional edge values. */
    std::unordered_map<std::string, std::vector<double>> EdgeValues;

    std::string __str__() const override
    {
      return "Road with " +
             str(Vertices.size()) + " vertices, " +
             str(Edges.size() / 2) + " edges, " +
             str(VertexValues.size()) + " vertex value arrays, and " +
             str(EdgeValues.size()) + " edge value arrays";
    }
  };
} // namespace DTCC

#endif // DTCC_ROAD_H
