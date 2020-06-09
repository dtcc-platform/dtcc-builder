// vc-cityjson-hello-world
// Vasilis Naserentin 2020
// Orfeas Eleftheriou 2020

#include <iostream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "HeightMap.h"
#include "HeightMapGenerator.h"
#include "JSON.h"
#include "LAS.h"
#include "Parameters.h"
#include "CityJSON.h"
#include "Logging.h"

using namespace DTCC;

void Help()
{
  std::cerr << "Usage: vc-generate-heightmap Parameters.json" << std::endl;
}

void assign(const nlohmann::json iJson, const std::string iUUID, const std::string iKey)
{
if (!iJson[iUUID][iKey].is_null())
{
std::cout<<"Found it!"<<std::endl;
}
}


int main(int argc, char *argv[])
{
  // Check command-line arguments
  //if (argc != 2)
 // {
 //   Help();
 //   return 1;
//  }

  // Read parameters
//  Parameters parameters;
  nlohmann::json inputjson;
//  JSON::Read(parameters, argv[1]);
//  std::cout << parameters << std::endl;
  inputjson=JSON::Read("../data/sampleCityJSON/testBox.json");

//for (nlohmann::json::iterator it = j.begin(); it != j.end(); ++it) {
//  std::cout << *it << '\n';
//}

//CityJSON cityObj;
//////for (size_t i=0;i<inputjson["Vertices"].size();i++)
//////{
//////cityObj.Vertices.push_back(inputjson["vertices"]);
//////}
//cityObj.Vertices = inputjson["vertices"];
//std::cout << cityObj.Vertices << std::endl;
//cityObj.Type=inputjson["type"];
//nlohmann::json cityjson=inputjson["CityObjects"];
//////cityObj.Version=inputjson["version"];
//////std::cout<<j<<std::endl;
//////nlohmann::json::iterator it = j.begin(); 
//////std::cout<<it.key()<<it.value()<<std::endl;
//
//auto obj = cityjson.get<nlohmann::json::object_t>();
//for (auto kvp = obj.begin(); kvp!=obj.end();kvp++)
//////for (nlohmann::json::iterator it = obj.begin(); it != obj.end(); ++it) {
//        {
//           
//            std::cout << kvp->first << std::endl;
//            cityObj.UUID.push_back(kvp->first);
//        }
//std::cout<<cityjson[cityObj.UUID[0]]["attributes"]<<std::endl;
//if (!cityjson["attributes"]["measuredHeight"].is_null())    
//{
//cityObj._CityObjects._Attributes.MeasuredHeight=cityjson["attributes"]["measuredHeight"];
//std::cout<<"Found it"<<std::endl;
//}
//cityObj.printInfo();

//assign(cityjson, cityObj.UUID[i], "attributes");

//std::map<std::string, std::string> m = j.at("CityObjects").get<std::map<std::string, std::string>>();



  CityJSON* cityJson = new CityJSON(inputjson);
  std::cout << "Printing city objs.."<<std::endl;
  std::vector<CityObject> cityObjects = cityJson->getCityObjects();
  std::cout << "city objects size:"<<cityObjects.size()<<std::endl;

  for (auto it = cityObjects.begin(); it != cityObjects.end(); ++it)
  {
    std::cout<<*it<<std::endl;
  }

delete cityJson;



  return 0;
}

