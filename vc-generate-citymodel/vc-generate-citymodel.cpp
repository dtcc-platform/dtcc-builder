// VirtualCity@Chalmers: vc-generate-citymodel
// Anders Logg 2019

#include <iostream>
#include <random>

#include "CityModel.h"
#include "Geometry.h"
#include "JSON.h"

using namespace VirtualCity;

void help()
{
    std::cerr << "Usage: vc-generate-citymodel"
              << std::endl;
}

int main(int argc, char* argv[])
{
    // Check command-line arguments
    if (argc != 1)
    {
        help();
        return 1;
    }

    // FIXME: This is just a test script for now.
    // We generate some random buildings that fit on the height map.
    // Note that the height map is wrong due to rescaling of the image.

    // Geometry of height map
    const Point2D C(-14756.8, 12294.2);
    const double W = 2154.0;
    const double H = 2158.0;
    const double w = 0.01*W;
    const double h = 0.01*H;

    // Generate the buildings
    const size_t n = 16;
    const double a = 0.05 * h;
    CityModel cityModel;
    std::vector<Point2D> centers;
    for (size_t i = 0; i < n; i++)
    {
        while (true)
        {
            // Randomize point
            double x = rand() / double(RAND_MAX);
            double y = rand() / double(RAND_MAX);
            Point2D p(C.x + (x - 0.5) * 0.3 * w, C.y + (y - 0.5) * 0.3 * h);

            std::cout << x << " " << y << std::endl;
            std::cout << "p = " << p << std::endl;

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
                double height = rand() / double(RAND_MAX);
                height *= 10.0 * a;
                height += a;

                // Create building
                Building b;
                b.Footprint.push_back(Point2D(p.x - 0.5 * a, p.y - 0.5 * a));
                b.Footprint.push_back(Point2D(p.x + 0.5 * a, p.y - 0.5 * a));
                b.Footprint.push_back(Point2D(p.x + 0.5 * a, p.y + 0.5 * a));
                b.Footprint.push_back(Point2D(p.x - 0.5 * a, p.y + 0.5 * a));
                b.Height = height;

                // Add building
                cityModel.Buildings.push_back(b);

                break;
            }
        }

    }

    // Write to file
    JSON::Write(cityModel, "CityModel.json");

    return 0;
}
