// Representation of a CityJSON
// Copyright (C) 2020 Anders Logg, Vasilis Naserentin.

#ifndef DTCC_CITYJSON_H
#define DTCC_CITYJSON_H

#include <vector>
#include "JSON.h"
#include "Point.h"

/*
 {
  "type": "CityJSON",
  "version": "1.0",
  "CityObjects": {},
  "vertices": []
}

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

*/
namespace DTCC
{
struct Attributes {
    double MeasuredHeight;
};
struct Semantics {
    std::vector<nlohmann::json> Surfaces;
    std::vector<nlohmann::json> Values;
};
struct _Geometry {
    std::vector<nlohmann::json> Boundaries;
    int64_t Lod;
    Semantics _Semantics;
    std::string Type;
};
struct CityObjects {
    Attributes _Attributes;
    std::vector<_Geometry> Geometry;
    std::string Type;
};
class CityJSON{
public:

    std::vector<std::string> UUID;
    CityObjects _CityObjects;
    std::string Type;
    std::string Version;
    nlohmann::json Vertices;
    // Create empty CityJSON
    CityJSON() : Type("CityJSON"), Version("1.0"), Vertices{} {}

    void printInfo()
    {
        Info("Information of CityJSON object");
        Info(Type);
        Info(Version);
        Info("Included UUIDs: ");
        for (size_t i = 0; i < UUID.size(); i++)
        {
            Info(UUID[i]);
        }
    }
};
}
#endif
