#include <iostream>
#include <netcdf>
#include "JSON.h"
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

void Help()
{
  std::cerr << "Usage: vc-read-netcdf filein.nc fileout.json var"
            << std::endl;
}


// We are reading 2D data, a 6 x 12 grid. 
static const int NX = 40;
static const int NY = 40;
static const int NZ = 42;
static const int TIME = 4;

// Return this in event of a problem.
static const int NC_ERR = 2;

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    Help();
    return 1;
  }
  std::string varname = argv[3];
  std::string fileout = argv[2];
  std::string filein = argv[1];
  nlohmann::json json;
  NcFile dataFile(filein, NcFile::read);

  json["Type"]="VectorFieldStructured";
  json["Dimension"]="3";
//json ["Grid"]={0,2,40,0,2,40,0,2,80};

   try
   {
   // This is the array we will read.
   float dataInu[NX][NY][NZ][TIME]; 
   //float dataInv[NX][NY][NZ][TIME];
   //float dataInw[NX][NY][NZ][TIME];
   // Open the file for read access

   // Retrieve the variable named "data"
   NcVar datau=dataFile.getVar(varname);
   //NcVar datav=dataFile.getVar("v");
   //NcVar dataw=dataFile.getVar("w");
   if(datau.isNull()) return NC_ERR;
   //if(datav.isNull()) return NC_ERR;
   //if(dataw.isNull()) return NC_ERR; 
   datau.getVar(dataInu);
   //datav.getVar(dataInv);
   //dataw.getVar(dataInw);

   // Pivoting x then y then z
   for (int k = 0; k < NZ; k++)
      for (int j = 0; j < NY; j++)
	for (int i = 0; i < NX; i++)
	{
        
	std::cout<<"U "<<i<<" "<<j<<" "<<k<<" "<<dataInu[i][j][k][0]<<std::endl;
	//std::cout<<"V "<<i<<" "<<j<<" "<<k<<" "<<dataInv[i][j][k][0]<<std::endl;
    //std::cout<<"W "<<i<<" "<<j<<" "<<k<<" "<<dataInw[i][j][k][0]<<std::endl;
json["Values"].push_back(dataInu[i][j][k][0]);
//json["Values"].push_back(dataInv[i][j][k][0]);
//json["Values"].push_back(dataInw[i][j][k][0]);

	}

   }catch(NcException& e)
     {
       e.what();
       cout<<"FAILURE*************************************"<<endl;
       return NC_ERR;
     }
std::ofstream o(fileout);
  o << json << std::endl;

}
