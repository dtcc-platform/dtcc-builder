// NetCDF4 Class.
// Copyright (C) Vasilis Naserentin
// Licensed under the MIT License

#ifndef DTCC_NETCDF4_H
#define DTCC_NETCDF4_H

#include <netcdf>
#include <tuple>

#include "Logging.h"

namespace DTCC_BUILDER
{
using namespace netCDF;

// template <class T>
class NetCDF4
{

public:
  std::string Name;
  std::string FileName;
  std::vector<double> Origin;
  std::vector<NcVar> Coordinates;
  // std::tuple<std::string, NcVar> CoordinatesTuple;
  size_t Size;
  NcVar mNcVar;
  std::vector<size_t> CoordinateDimensions;
  // std::unique_ptr<T[]> Pointer;
  // T *Vector;
  explicit NetCDF4(size_t size)
  {
    info("Creating new NetCDF4 with known size " + str(size));
    // Vector = new T[size];
    // Pointer.reset(new T[Size]);
    // Size = size;
  }
  NetCDF4()
  {
    info("Creating new NetCDF4 without known size");

    Size = 0;
  }
  ~NetCDF4()
  {
    info("Destroying NetCDF4");
    info("Size is " + str(Size));
    /*if (Size!=0)
        {
        delete Vector;
    Vector = NULL;
    }*/
  }
  template <typename T> void allocateArray(NcVar &iData)
  {
    std::unique_ptr<T[]> Pointer(new T[Size]);
    T *Vector = new T[Size];
    iData.getVar(Pointer.get());
    iData.getVar(Vector);
    toJson(Pointer);
    /*if (Size != 0)
    {
      this->Vector = new T[Size];
          this->Pointer.reset(new T[Size]);
          iData.getVar(this->Vector);
          iData.getVar(this->Pointer.get());
    }
    else
    {
      std::cout << "Size is 0. Define Size first!" << std::endl;
    }*/
  }
  // void allocateArray(T *&p, size_t size, NcVar &iData) { iData.getVar(p); }
  template <typename T> void initialize(NcVar &iVar, NcFile &iFile)
  {
    NcDim dim;
    std::vector<size_t> Dimensions;
    std::vector<std::string> CoordinateNames;
    int numDims = iVar.getDimCount();
    for (int i = 0; i < numDims; i++)
    {
      dim = iVar.getDim(i);
      size_t dims = dim.getSize();
      Dimensions.push_back(dims);

      debug("Dim is " + str(dims) + " of " + str(numDims) + " with name " +
            dim.getName());
      CoordinateNames.push_back(dim.getName());

      // NcVar tempvar=iFile.getVar(dim.getName());
    }
    Size = std::accumulate(begin(Dimensions), end(Dimensions), 1,
                           std::multiplies<double>());
    for (size_t i = 0; i < CoordinateNames.size(); i++)
    {
      NcVar tempvar = iFile.getVar(CoordinateNames[i]);
      Coordinates.push_back(tempvar);
      // double *tempvector = NULL;
      // tempvector = new double[Dimensions[i]];
      // tempvar.getVar(tempvector);
      // CoordinatesTuple = std::make_tuple(CoordinateNames[i], tempvar);
    }
    CoordinateDimensions = Dimensions;
    getOrigin(iVar);
    allocateArray<T>(iVar);
  }
  void extractCoordinates(int i)
  {

    auto tempvector = new double[CoordinateDimensions[i]];

    debug("into extractCoordinates");
    debug(Coordinates[i].getName());

    Coordinates[i].getVar(tempvector);

    debug("Got em" + str(CoordinateDimensions[i]));
    debug(str(tempvector[0]));
    // iVector.push_back(tempvector);

    debug("Got em all");
    delete[] tempvector;
    tempvector = nullptr;
  }
  void getOrigin(NcVar &iData)
  {
    NcGroup group = iData.getParentGroup();
    debug(str(group.getAttCount()));
    std::multimap<std::string, NcGroupAtt> myMap;
    std::multimap<std::string, NcGroupAtt>::iterator itr;
    std::map<std::string, NcGroup> myCoordMap;
    std::map<std::string, NcGroup>::iterator citr;
    myCoordMap = group.getCoordVars();
    for (citr = myCoordMap.begin(); citr != myCoordMap.end(); ++citr)
    {
      debug("Coordinates");
      debug(citr->first);
    }
    myMap = group.getAtts();
    int counter = 0;
    for (itr = myMap.begin(); itr != myMap.end(); ++itr)
    {
      debug(str(counter) + " " + itr->first + " " + itr->second.getName());
      counter++;
      NcType type = itr->second.getType();
      debug(type.getName());
      if (type.getName() == "char")
        debug(" char found");
      if (itr->first == "origin_x")
      {
        debug(" X origin found");
        double *originx = new double[itr->second.getAttLength()];
        itr->second.getValues(originx);
        Origin.push_back(originx[0]);
        debug(str(originx[0]));
      }
      if (itr->first == "origin_y")
      {
        debug(" y origin found");
        double *originy = new double[itr->second.getAttLength()];
        itr->second.getValues(originy);
        Origin.push_back(originy[0]);
        debug(str(originy[0]));
      }
      if (itr->first == "origin_z")
      {
        debug(" z origin found");
        double *originz = new double[itr->second.getAttLength()];
        itr->second.getValues(originz);
        Origin.push_back(originz[0]);
        debug(str(originz[0]));
      }
    }
  }

  void getCoordinates(NcVar &iVar)
  {
    // Coordinates[0].getVar(CoordinatesName[0]);
    for (size_t i = 0; i < Coordinates.size(); i++)
    {
      info("Coordinate Name " + Coordinates[i].getName());
    }
  }

  template <typename T> void toJson(std::unique_ptr<T[]> &iPtr)
  {
    std::unique_ptr<double[]> tempPtr;
    nlohmann::json json;
    json["Type"] = "NetCDF4Parsed";
    json["Variable"] = Name;
    // for (size_t i=0;i<Size;i++)
    /*{
    json["VariableValues"].push_back(iPtr[i]);
}*/
    for (size_t i = 0; i < CoordinateDimensions.size(); i++)
    {
      json["VariableDimension"].push_back(CoordinateDimensions[i]);
    }
    for (size_t i = 0; i < Origin.size(); i++)
    {
      json["Origin"].push_back(Origin[i]);
    }
    for (size_t i = 0; i < CoordinateDimensions.size(); i++)
    {
      if (Coordinates[i].getName() != Name)
      {
        debug(Coordinates[i].getName());
        ;
        tempPtr.reset(new double[CoordinateDimensions[i]]);
        Coordinates[i].getVar(tempPtr.get());
        for (size_t j = 0; j < CoordinateDimensions[i]; j++)
        {
          json[Coordinates[i].getName()].push_back(tempPtr[j]);
        }
      }
    }

    // VariableDimension 1 130 512 512
    // Origin x0,y0,z0
    // AssociatedCoordinates xu zu_3d Values
    std::ofstream o(FileName);
    o << json << std::endl;
    o.close();
  }
};
#if 0
static void readGlobalAtts(NcVar &datau, std::vector<double> origins)
{
  NcGroup group = datau.getParentGroup();
  std::cout << group.getAttCount() << std::endl;
  std::multimap<std::string, NcGroupAtt> myMap;
  std::multimap<std::string, NcGroupAtt>::iterator itr;
  std::map<std::string, NcGroup> myCoordMap;
  std::map<std::string, NcGroup>::iterator citr;
  myCoordMap = group.getCoordVars();
  for (citr = myCoordMap.begin(); citr != myCoordMap.end(); ++citr)
  {
    std::cout << "Coordinates" << std::endl;
    std::cout << citr->first << std::endl;
  }

  myCoordMap = group.getCoordVars();
  myMap = group.getAtts();
  int counter = 0;
  for (itr = myMap.begin(); itr != myMap.end(); ++itr)
  {

    std::cout << counter << " " << itr->first << " " << itr->second.getName()
              << std::endl;
    counter++;
    NcType type = itr->second.getType();
    std::cout << type.getName() << std::endl;
    if (type.getName() == "char")
      std::cout << " char found" << std::endl;
    if (itr->first == "origin_x")
    {
      std::cout << " X origin found" << std::endl;
      double *originx = new double[itr->second.getAttLength()];
      itr->second.getValues(originx);
      origins.push_back(originx[0]);
      std::cout << originx[0] << std::endl;
    }
    if (itr->first == "origin_y")
    {
      std::cout << " y origin found" << std::endl;
      double *originy = new double[itr->second.getAttLength()];
      itr->second.getValues(originy);
      origins.push_back(originy[0]);
      std::cout << originy[0] << std::endl;
    }
    if (itr->first == "origin_z")
    {
      std::cout << " z origin found" << std::endl;
      double *originz = new double[itr->second.getAttLength()];
      itr->second.getValues(originz);
      origins.push_back(originz[0]);
      std::cout << originz[0] << std::endl;
    }

    if (counter == 16)
    {
      std::cout << "getting values" << std::endl;
      // create
      char *values = new char[itr->second.getAttLength()];
      std::cout << itr->second.getAttLength() << std::endl;
      itr->second.getValues(values);
      std::string test(values);
      std::cout << values[0] << std::endl;
      std::cout << test << std::endl;
    }
  }
}
static void readVarAtts(NcVar &datau)
{
  std::cout << "Var atts" << std::endl;
  std::map<std::string, NcVarAtt> myMap;
  std::map<std::string, NcVarAtt>::iterator itr;
  myMap = datau.getAtts();
  NcGroup mygroup = datau.getParentGroup();
  std::map<std::string, NcGroup> myCoordMap;
  std::map<std::string, NcGroup>::iterator citr;
  myCoordMap = mygroup.getCoordVars(NcGroup::Location::Current);
  for (citr = myCoordMap.begin(); citr != myCoordMap.end(); ++citr)
  {
    std::cout << "Coordinates" << std::endl;
    std::cout << citr->first << std::endl;
  }
  for (itr = myMap.begin(); itr != myMap.end(); ++itr)
  {
    std::cout << itr->first << "\n";
    std::cout << itr->second.getName() << "\n";
    std::cout << itr->second.getAttLength() << "\n";
    char *values = new char[itr->second.getAttLength()];
    std::string test(values);
    std::cout << values << std::endl;
  }
}
#endif
} // namespace DTCC_BUILDER

#endif
