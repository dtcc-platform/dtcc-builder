#include "JSON.h"
#include "NetCDF4.h"
#include "Parameters.h"
#include "PostProcessParameters.h"
#include <iostream>
#include <memory>
#include <netcdf>
#include <numeric>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace DTCC;
void Help()
{
  std::cerr << "Usage: vc-read-netcdf filein.nc fileout.json variable"
            << std::endl;
}

void unflatted(int height, int width, int depth, int index)
{

  int width_index =
      index / (height * depth); // Note the integer division . This is x
  int height_index = (index - width_index * height * depth) / depth; // This is
                                                                     // y
  int depth_index =
      index - width_index * height * depth - height_index * depth; // This is z
  std::cout << height_index << std::endl;
  std::cout << width_index << std::endl;
  std::cout << depth_index << std::endl;
}

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    Help();
    return 1;
  }
  std::string fileout = argv[2];
  std::string filein = argv[1];
  nlohmann::json json;
  std::cout << "I am here" << std::endl;
  // PostProcessParameters ppparameters;
  json["Type"] = "NetCDF4Parsed";
  // JSON::Read(ppparameters, argv[3]);
  // json["Variable"] = ppparameters.Variable;
  json["Variable"] = argv[3];
  // std::cout << ppparameters << std::endl;
  NcDim dim;
  std::cout << "I am here 2" << std::endl;
  NcFile dataFile(filein, NcFile::read);
  std::cout << "I am here 3" << std::endl;
  // NcVar datau = dataFile.getVar(ppparameters.Variable.c_str());
  NcVar datau = dataFile.getVar(argv[3]);
  NcType type = datau.getType();

  std::vector<size_t> dimensions;
  if (type == NC_FLOAT)
  {
    // NetCDF4<float> myNet;
    NetCDF4 myNet;
    // myNet.Name = ppparameters.Variable.c_str();
    myNet.Name = argv[3];
    myNet.FileName = fileout;
    myNet.initialize<float>(datau, dataFile);
    // myNet.toJson("test.in");
    // dimensions=NetCDF4<float>::createDimensions(datau, dataFile);
    // NetCDF4<float> myNet(std::accumulate(begin(dimensions), end(dimensions),
    // 1, 	  std::multiplies<double>()));
    // myNet.Name = ppparameters.Variable.c_str();
    std::unique_ptr<double[]> p;
    p.reset(new double[myNet.CoordinateDimensions[2]]);
    NcVar datav = dataFile.getVar(myNet.Coordinates[2].getName());
    datav.getVar(p.get());
    // myNet.CoordinateDimensions = dimensions;
    // myNet.allocateArray(myNet.Vector, myNet.Size, datau);
    // myNet.allocateArray(datau);

    // myNet.getOrigin(datau);
    // myNet.getCoordinates(datau);
    // std::cout<<myNet.Vector[0]<<std::endl;

    // myNet.Coordinates[1].getVar(dataIn);
    myNet.extractCoordinates(1);
    myNet.extractCoordinates(2);
    // syNet.extractCoordinates(5);

    // readGlobalAtts(datau, myNet.Origin);

    // myNet.Size = std::accumulatea(begin(dimensions), end(dimensions), 1,
    // std::multiplies<double>());
    // myNet(5);
  }
  return 0;

  // if (type == NC_FLOAT)
  {

    // std::unique_ptr< NetCDF4 < float> > uptr(new NetCDF4<float>(multi));
  }
  //       {
  //		   //NetCDF<float> myvariable;
  //		   size_t multi = std::accumulate(begin(dimensions),
  // end(dimensions), 1, 			   std::multiplies<double>());
  //
  //		   //NetCDF4<float> myvariable(multi);
  //		   //NetCDF4<float>* floatPtr= new NetCDF4<float>(multi);
  //		   //std::unique_ptr<NetCDF4 < float> > uptr(floatPtr);
  //		   std::unique_ptr< NetCDF4 < float> > uptr(new
  // NetCDF4<float>(multi)); 		   uptr->Name = ppparameters.Variable;
  //		   //myvariable.Name = ppparameters.Variable;
  //
  //		   std::cout << "float !";
  //         //float *arr2 = NULL;
}

//   try
//   {
//    NcFile dataFile(filein, NcFile::read);
//    NcVar datau=dataFile.getVar(ppparameters.Variable.c_str());
//    std::vector<size_t> dimensions;
//	std::vector<NcVar> coordinates;
//	ppparameters.Dimension = datau.getDimCount();
//	json["Dimension"] = ppparameters.Dimension;
//	for (int i=0;i<ppparameters.Dimension;i++)
//    {
//      dim = datau.getDim(i);
//      size_t dims = dim.getSize();
//      dimensions.push_back(dims);
//      std::cout << "Dim is "<<dims <<" of "<< ppparameters.Dimension<< " with
//      name " <<dim.getName()<<std::endl;
//	  NcVar tempvar=dataFile.getVar(dim.getName());
//	  coordinates.push_back(tempvar);
//      json["Boundaries"].push_back(dims);
//
//    }
//   NcType type = datau.getType();
//   std::cout << type.getName() << std::endl;
//
//
//
//       if (type == NC_FLOAT)
//       {
//		   //NetCDF<float> myvariable;
//		   size_t multi = std::accumulate(begin(dimensions),
// end(dimensions), 1, 			   std::multiplies<double>());
//
//		   //NetCDF4<float> myvariable(multi);
//		   //NetCDF4<float>* floatPtr= new NetCDF4<float>(multi);
//		   //std::unique_ptr<NetCDF4 < float> > uptr(floatPtr);
//		   std::unique_ptr< NetCDF4 < float> > uptr(new
// NetCDF4<float>(multi)); 		   uptr->Name = ppparameters.Variable;
//		   //myvariable.Name = ppparameters.Variable;
//
//		   std::cout << "float !";
//         //float *arr2 = NULL;
//         uptr->allocateArray(uptr->Vector, dimensions,
//         ppparameters.Variable,datau, json);
//		 unflatted( dimensions[1], dimensions[2], dimensions[3],529);
//
//       }
//       else if (type == NC_DOUBLE)
//       {
//		   size_t multi = std::accumulate(begin(dimensions),
// end(dimensions), 1, 			   std::multiplies<double>());
//         std::cout << "DOUBLE !";
//         //double *arr2 = NULL;
//		 NetCDF4<double> myvariable(multi);
//		 myvariable.Name = ppparameters.Variable;
//		 myvariable.allocateArray(myvariable.Vector, dimensions,
// ppparameters.Variable,datau,json);
//	   }
//	   else
//	   {
//		   std::cout << "Only variables of type int, double or float
// currently supported." << std::endl;; 		   return 0;
//	   }
//	   std::cout << "Writing json" << std::endl;
//	   std::vector<double> origins;
//	   readGlobalAtts(datau,origins);
//	   for (int i=0;i<3;i++)
//	   {
//	   json["Origin"].push_back(origins[i]);
//       }
//	   for (size_t i = 0; i < coordinates.size(); i++)
//	   {
//		   NcType temp=coordinates[i].getType();
//		   std::cout << "Coordinate type is " << temp.getName() <<" with
// dimension "<< dimensions[i]<< std::endl;
//		   //double *arr2 = NULL;
//		   NetCDF4<double> myvariable(dimensions[i]);
//		   myvariable.Name = ppparameters.Variable;
//		   std::cout << myvariable.Name << " vs " <<
// coordinates[i].getName(); 		   if
// (coordinates[i].getName()!=myvariable.Name)
//		   {myvariable.allocateArray(myvariable.Vector, dimensions[i],
// coordinates[i]);
// json["Coordinates"].push_back(coordinates[i].getName()); 		   for
//(size_t j = 0; j < dimensions[i]; j++)
//		   {
//			   json[coordinates[i].getName()].push_back(myvariable.Vector[j]);
//		   }
//		   }
//	   }
//	     std::cout  << std::endl;
//
//
//
//
//
//
//   //static float *arr2;
//       //= new float[NX * NY * NZ];
//
//       //allocateArray(&arr2, NX, NY, NZ);
//
//
//   // if(datau.isNull()) return NC_ERR;
//
//   }catch(NcException& e)
//     {
//       e.what();
//       cout<<"FAILURE*************************************"<<endl;
//       std::cout<< e.what() << std::endl;
//	   return 1;
//     }
// std::ofstream o(fileout);
//  o << json << std::endl;
//
//}
