// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <random>

#include "CityModel.h"
#include "Geometry.h"
#include "JSON.h"

using namespace VirtualCity;

void Help()
{
    std::cerr << "Usage: vc-generate-citymodel"
              << std::endl;
}

// Generate building at given position with side length a and height h
void GenerateBuilding(Building& building, const Point2D& p, double a, double h)
{
    building.Footprint.push_back(Point2D(p.x - 0.5 * a, p.y - 0.5 * a));
    building.Footprint.push_back(Point2D(p.x + 0.5 * a, p.y - 0.5 * a));
    building.Footprint.push_back(Point2D(p.x + 0.5 * a, p.y + 0.5 * a));
    building.Footprint.push_back(Point2D(p.x - 0.5 * a, p.y + 0.5 * a));
    building.Height = h;
}

// Generate maximal city model (fill out the height map)
void GenerateMaximalCityModel(CityModel& cityModel)
{
    // Geometry of height map (hard coded)
    const Point2D C(-1081.75, -1418.75);
    const double W = 29504.0;
    const double H = 29584.0;

    // Parameters for building size
    const double a = 100.0;
    const double h = 1.0;
    const double margin = 0.9;

    // Displacement of buildings relative to center
    const double d = margin*(0.5 * std::min(W, H) - a);
    //const double d = 1000;

    // Add a building in each corner
    Building b0, b1, b2, b3;
    GenerateBuilding(b0, C + Point2D(-d, 0), a, h);
    GenerateBuilding(b1, C + Point2D(+d, 0), a, h);
    GenerateBuilding(b2, C + Point2D(0, -d), a, h);
    GenerateBuilding(b3, C + Point2D(0, +d), a, h);
    cityModel.Buildings.push_back(b0);
    cityModel.Buildings.push_back(b1);
    cityModel.Buildings.push_back(b2);
    cityModel.Buildings.push_back(b3);
}

// Generate random city model
void GenerateRandomCityModel(CityModel& cityModel, size_t numBuildings)
{
    // City center
    const Point2D C(-1081.75, -1418.75);

    // Parameters for building sizes and locations
    const double R = 100.0; // city radius
    const double a = 10.0;  // building side length
    const double h = 100.0; // maximum building height

    // Generate the buildings
    std::vector<Point2D> centers;
    for (size_t i = 0; i < numBuildings; i++)
    {
        while (true)
        {
            // Randomize point [-1, 1]
            double x = 2.0 * (rand() / double(RAND_MAX) - 0.5);
            double y = 2.0 * (rand() / double(RAND_MAX) - 0.5);
            Point2D p(C.x + R * x, C.y + R * y);
            std::cout << "Trying to add building at p = " << p << std::endl;

            // Check that we are not too close to other buildings
            bool ok = true;
            for (auto const & q : centers)
            {
                const double d = Geometry::Distance2D(p, q);
                if (d < 2.0 * a)
                {
                    ok = false;
                    break;
                }
            }

            // Add building
            if (ok)
            {
                // Randomize height
                double height = a + (rand() / double(RAND_MAX)) * h;

                // Generate building
                Building building;
                GenerateBuilding(building, p, a, h);

                // Add building
                cityModel.Buildings.push_back(building);
                centers.push_back(p);

                std::cout << "Building added" << std::endl;

                break;
            }
        }
    }
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 1)
    {
        Help();
        return 1;
    }

    // Generate maximal city cityModel
    CityModel cityModel;
    GenerateMaximalCityModel(cityModel);

    // Generate random city model
    // CityModel cityModel;
    //GenerateRandomCityModel(cityModel, 16);

    // Write to file
    JSON::Write(cityModel, "CityModel.json");

    return 0;
}
