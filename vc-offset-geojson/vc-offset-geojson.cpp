// VirtualCity@Chalmers: vc-offset-geojson
// Anders Logg 2019
// Vasilis Naserentin 2019

#include <string>
#include <iostream>
#include "JSON.h"

void Help()
{
    std::cerr << "Usage: vc-offset-geojson filein.geojson fileout.geojson x y " << std::endl;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 5)
    {
      Help();
      return 1;
    }
	int x, y;
	try
	{
	x = std::stoi(argv[3]);
	y = std::stoi(argv[4]);
	}
	catch (...)
	{
	std::cout<<"Error; x and y must be integers."<<std::endl;
	return 1;
	}
//	y = std::stoi(argv[4]);
	std::string fileout = argv[2];
	std::string filein = argv[1];
	std::ifstream in(filein); 
	if (in.fail()) {
		std::cout << "Input file issue; check input file exists." << std::endl;
		return 1;
	}
	nlohmann::json j;
	in >> j;
	double tempx, tempy;
	int offsetxx = x;
	int offsetyy = y;
	std::cout << "JSON stuff" << std::endl;
	std::cout << j["type"]<<std::endl;
	std::cout << j["features"][0] << std::endl;
	std::cout << j["features"].size()<< std::endl;
	//std::cout << j["features"][0]["geometry"] << std::endl;
	//show whole polygon, this returns [x,y] and change the coordinates
	for (int k = 0; k < j["features"].size(); k++)
	{
		for (int i = 0; i < j["features"][k]["geometry"]["coordinates"][0][0].size(); i++)
		{
			//std::cout << j["features"][k]["geometry"]["coordinates"][0][0][i] << std::endl;
			tempx = j["features"][k]["geometry"]["coordinates"][0][0][i][0];
			tempy = j["features"][k]["geometry"]["coordinates"][0][0][i][1];

			j["features"][k]["geometry"]["coordinates"][0][0][i][0] = tempx + offsetxx;
			j["features"][k]["geometry"]["coordinates"][0][0][i][1] = tempy + offsetyy;

		}
		
			std::cout << j["features"][k]["geometry"]["coordinates"][0][0] << std::endl;


		
	}
	
	
	std::cout << j["features"][0]["geometry"]["coordinates"][0][0][0].size() << std::endl;//first x
	std::cout << j["features"][0]["geometry"]["coordinates"][0][0].size() << std::endl;//first x
	//std::cout << j["features"][0]["geometry"]["coordinates"][0].size() << std::endl;//first x
	std::ofstream o(fileout);
	//o << std::setw(4) << j << std::endl;
	o << j << std::endl;
	
	return 0;
}
