// Representation of a CityJSON
// Copyright (C) 2020 Anders Logg, Vasilis Naserentin.

#ifndef DTCC_CITYJSON_H
#define DTCC_CITYJSON_H

#include <vector>

#include "JSON.h"

namespace DTCC
{

struct Extensions {
};

struct Appearance {
};

struct Geometry-templates
};

struct Transform {
    std::vector<nlohmann::json> Scale;
    std::vector<nlohmann::json> Translate;
};

struct Metadata {
};

struct Attributes {
};

struct Semantics {
    std::vector<nlohmann::json> Surfaces;
    std::vector<nlohmann::json> Values;
};

struct CJGeometry { // temp till we merge with Geometry.h
    std::vector<nlohmann::json> Boundaries;
    int64_t Lod;
    Semantics semantics;
    std::string Type;
};

struct CityObjects {
    Attributes attributes;
    std::vector<CJGeometry> geometry;
    std::string type;
};

class CityJSON
{
public:

    CityObjects cityObjects;
    std::string Type;
    std::string Version;
    std::vector<nlohmann::json> Vertices;

  // Create empty CityJSON
  CityJSON() : Version("1.0") {}
};


} // namespace DTCC
#endif

// Template CityJson.json
/*
{
 "CityObjects": {
         "attributes": {},
         "geometry": [
           { "boundaries": [], "lod": 1, "semantics":{"surfaces":[], "values":[]},
             "type": "Solid"
             }
           
           ],
           "type": "Building"
 },
  "type": "CityJSON",
  "version": "1.0",
  "vertices": [
  ]
}
------------------------------------------
{
  "type": "CityJSON",
  "version": "1.0",
  "extensions": {},
  "metadata": {},
  "transform": {
    "scale": [],
    "translate": []
  },
  "CityObjects": {},
  "vertices": [],
  "appearance": {},
  "geometry-templates": {}
}
*/
