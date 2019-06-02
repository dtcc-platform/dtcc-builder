// VirtualCity@Chalmers: vc-offset-geojson
// Anders Logg 2019
// Vasilis Naserentin 2019
// Offset in EPSG:3006 used in UE -148000 -6398600 (meters)
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
	double curx, cury;
	int offsetx = x;
	int offsety = y;
	std::cout << "JSON stuff" << std::endl;
	std::cout << "JSON has "<<j["features"].size()<<" of type "<<j["features"][0]["geometry"]["type"] << std::endl;

	//std::cout << j["features"][0]["geometry"] << std::endl;
	//show whole polygon, this returns [x,y] and change the coordinates
	std::cout<<"Offsetting"<<std::endl;
	for (int k = 0; k < j["features"].size(); k++)
	{
		for (int i = 0; i < j["features"][k]["geometry"]["coordinates"][0].size(); i++)
		{
			//tempx = j["features"][k]["geometry"]["coordinates"][0][0][i][0];
			//tempy = j["features"][k]["geometry"]["coordinates"][0][0][i][1];
			//std::cout<< j["features"][k]["geometry"]["coordinates"][0][i][0]<<std::endl;
			//std::cout<< j["features"][k]["geometry"]["coordinates"][0][i][1]<<std::endl;
			std::cout<<"Feature "<<k<<" point " <<i<<" out of "<<j["features"][k]["geometry"]["coordinates"][0].size()<<" has x "
			<<j["features"][k]["geometry"]["coordinates"][0][i][0]<<"and y "<<j["features"][k]["geometry"]["coordinates"][0][i][1]<<std::endl;
	                double originalx=j["features"][k]["geometry"]["coordinates"][0][i][0];
			double originaly=j["features"][k]["geometry"]["coordinates"][0][i][1];
			double newx=originalx+offsetx;
			double newy=originaly+offsety;
			j["features"][k]["geometry"]["coordinates"][0][i][0]=newx;
			j["features"][k]["geometry"]["coordinates"][0][i][1]=newy;
 			std::cout<<"Feature "<<k<<" point " <<i<<" out of "<<j["features"][k]["geometry"]["coordinates"][0].size()<<" has x "
                        <<j["features"][k]["geometry"]["coordinates"][0][i][0]<<"and y "<<j["features"][k]["geometry"]["coordinates"][0][i][1]<<std::endl;			
			double validateoffsetx=newx-originalx;
			double validateoffsety=newy-originaly;
			std::cout<<"Offsetx is "<<validateoffsetx<<std::endl;
   			std::cout<<"Offsety is "<<validateoffsety<<std::endl;
			//j["features"][k]["geometry"]["coordinates"][0][0][i][0] = tempx + offsetxx;
			//j["features"][k]["geometry"]["coordinates"][0][0][i][1] = tempy + offsetyy;

		}
		
			//std::cout << j["features"][k]["geometry"]["coordinates"][0][0] << std::endl;


		
	}
	std::cout<<"Offset complete"<<std::endl;	
	
	std::cout << j["features"][0]["geometry"]["coordinates"][0][0][0].size() << std::endl;//first x
	std::cout << j["features"][0]["geometry"]["coordinates"][0][0].size() << std::endl;//first x
	//std::cout << j["features"][0]["geometry"]["coordinates"][0].size() << std::endl;//first x
	std::ofstream o(fileout);
	//o << std::setw(4) << j << std::endl;
	o << j << std::endl;
	
	return 0;
}
