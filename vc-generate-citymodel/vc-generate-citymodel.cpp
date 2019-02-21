// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <iomanip>
#include <random>

#include "CityModel.h"
#include "Geometry.h"
#include "JSON.h"

using namespace VirtualCity;

// FIXME: Consider moving to Parameters.json

// Get parameters for building sizes and locations
const double CITY_RADIUS = 100.0;     // city radius
const double BUILDING_SIZE = 10.0;    // building side length
const double BUILDING_HEIGHT = 100.0; // maximum building height
const double VELOCITY = 1.0;          // velocity including time step

void Help()
{
    std::cerr << "Usage: vc-generate-citymodel"
              << std::endl;
}

// Return random number between 0 and 1
double Random()
{
    return std::rand() / double(RAND_MAX);
}

// Generate building at given position with side length a and height h
Building GenerateBuilding(const Point2D& p, double a, double h)
{
    Building building;
    building.Footprint.push_back(Point2D(p.x - 0.5 * a, p.y - 0.5 * a));
    building.Footprint.push_back(Point2D(p.x + 0.5 * a, p.y - 0.5 * a));
    building.Footprint.push_back(Point2D(p.x + 0.5 * a, p.y + 0.5 * a));
    building.Footprint.push_back(Point2D(p.x - 0.5 * a, p.y + 0.5 * a));
    building.Height = h;
    return building;
}

// Generate maximal city model (fill out the height map)
std::vector<CityModel> GenerateMaximalCityModel()
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
    const double d = margin * (0.5 * std::min(W, H) - a);
    //const double d = 1000;

    // Create empty city model
    CityModel cityModel;

    // Add a building in each corner
    Building b0 = GenerateBuilding(C + Point2D(-d, 0), a, h);
    Building b1 = GenerateBuilding(C + Point2D(+d, 0), a, h);
    Building b2 = GenerateBuilding(C + Point2D(0, -d), a, h);
    Building b3 = GenerateBuilding(C + Point2D(0, +d), a, h);
    cityModel.Buildings.push_back(b0);
    cityModel.Buildings.push_back(b1);
    cityModel.Buildings.push_back(b2);
    cityModel.Buildings.push_back(b3);

    // Pack single model in list
    std::vector<CityModel> cityModels;
    cityModels.push_back(cityModel);

    return cityModels;
}

// Generate random city model
std::vector<CityModel> GenerateRandomCityModel(size_t numBuildings)
{
    // City center (approximate location of Johanneberg)
    const Point2D C(-3000.0, -4000.0);

    // Get parameters for building sizes and locations
    const double R = CITY_RADIUS;
    const double a = BUILDING_SIZE;
    const double h = BUILDING_HEIGHT;

    // Create empty city model
    CityModel cityModel;

    // Generate the buildings
    std::vector<Point2D> centers;
    for (size_t i = 0; i < numBuildings; i++)
    {
        while (true)
        {
            // Randomize point [-1, 1]
            double x = 2.0 * (Random() - 0.5);
            double y = 2.0 * (Random() - 0.5);
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
                // Randomize dimension
                double width = (0.5 + 0.5 * Random()) * a;
                double height = (0.2 + 0.8 * Random()) * h;

                // Generate building
                Building building = GenerateBuilding(p, width, height);

                // Add building
                cityModel.Buildings.push_back(building);
                centers.push_back(p);

                std::cout << "Building added" << std::endl;

                break;
            }
        }
    }

    // Pack single model in list
    std::vector<CityModel> cityModels;
    cityModels.push_back(cityModel);

    return cityModels;
}

// Generate random city model
std::vector<CityModel> GenerateEvolvingCityModel(size_t numBuildings,
        size_t numFrames)
{
    // Generate a random city model
    std::vector<CityModel> cityModels;
    cityModels = GenerateRandomCityModel(numBuildings);

    // Initialize building positions
    std::vector<Point2D> x(numBuildings);
    for (size_t i = 0; i < numBuildings; i++)
    {
        std::vector<Point2D>& footPrint = cityModels[0].Buildings[i].Footprint;
        for (auto const & p : footPrint)
            x[i] += p;
        x[i] /= footPrint.size();
    }

    // Compute center of gravity
    Point2D c;
    for (size_t i = 0; i < numBuildings; i++)
        c += x[i];
    c /= numBuildings;

    // Compute radius for perimeter and minimal building distance
    const double R = 2.0 * CITY_RADIUS;
    const double D = BUILDING_SIZE;

    // Initialize building velocities
    std::vector<Point2D> v(numBuildings);
    for (size_t i = 0; i < numBuildings; i++)
    {
        v[i].x = Random() - 0.5;
        v[i].y = Random() - 0.5;
    }

    // Bounce buildings around
    for (size_t n = 1; n < numFrames; n++)
    {
        // Check collisions (flip velocities)
        for (size_t i = 0; i < numBuildings; i++)
        {
            bool flip = false;

            // Check collision with other buildings
            for (size_t j = 0; j < numBuildings; j++)
            {
                if (i == j)
                    continue;

                const double d = Geometry::Distance2D(x[i], x[j]);
                if (d <= D)
                    flip = true;
            }

            // Check collision with perimeter
            const double r = Geometry::Distance2D(x[i], c);
            if (r >= R)
                flip = true;

            // Flip velocity
            if (flip)
            {
                v[i].x = -v[i].x;
                v[i].y = -v[i].y;
            }
        }

        // Update positions
        for (size_t i = 0; i < numBuildings; i++)
            x[i] += v[i];

        // Create new frame
        CityModel cityModel;
        cityModel = cityModels[n - 1];
        for (size_t i = 0; i < numBuildings; i++)
        {
            for (auto & p : cityModel.Buildings[i].Footprint)
                p += v[i];
        }
        cityModels.push_back(cityModel);
    }

    return cityModels;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 1)
    {
        Help();
        return 1;
    }

    // FIXME: Handle command-line and selection of city model

    // Generate empty list of city models
    std::vector<CityModel> cityModels;

    // Generate city model
    //cityModels = GenerateMaximalCityModel();
    //cityModels = GenerateRandomCityModel(64);
    cityModels = GenerateEvolvingCityModel(4, 100);

    // Write to file
    if (cityModels.size() == 1)
    {
        JSON::Write(cityModels[0], "CityModel.json");
    }
    else
    {
        for (size_t i = 0; i < cityModels.size(); i++)
        {
            std::stringstream s;
            s << "CityModel";
            s << std::setw(5) << std::setfill('0') << i;
            s << ".json";
            JSON::Write(cityModels[i], s.str());
        }
    }

    return 0;
}
