// VirtualCity@Chalmers: vc-randomize-citymodel
// Anders Logg 2019

#include <iomanip>
#include <iostream>
#include <random>

#include "CityModel.h"
#include "Geometry.h"
#include "JSON.h"

using namespace VirtualCity;

// FIXME: Consider moving to Parameters.json
const double BUILDING_SIZE = 10.0;    // building side length
const double BUILDING_HEIGHT = 100.0; // maximum building height
const double VELOCITY = 2.0;          // velocity including time step
const double BOUNDARY_MARGIN = 0.9;   // margin to domain boundary

void Help()
{
  std::cerr << "Usage: vc-randomize-citymodel Parameters.json" << std::endl;
}

// Return random number between 0 and 1
double Random() { return std::rand() / double(RAND_MAX); }

// Generate building at given position with side length a and height h
Building GenerateBuilding(const Point2D &p, double a, double h)
{
  Building building;
  building.Footprint.push_back(Point2D(p.x - 0.5 * a, p.y - 0.5 * a));
  building.Footprint.push_back(Point2D(p.x + 0.5 * a, p.y - 0.5 * a));
  building.Footprint.push_back(Point2D(p.x + 0.5 * a, p.y + 0.5 * a));
  building.Footprint.push_back(Point2D(p.x - 0.5 * a, p.y + 0.5 * a));
  building.Height = h;
  return building;
}

// Generate random city model
std::vector<CityModel>
GenerateRandomCityModel(size_t numBuildings, const Point2D &C, double R)
{
  // Get parameters for building dimensions
  const double a = BUILDING_SIZE;
  const double H = BUILDING_HEIGHT;
  const double m = BOUNDARY_MARGIN;

  // Create empty city model
  CityModel cityModel;

  // Generate the buildings
  std::vector<Point2D> centers;
  for (size_t i = 0; i < numBuildings; i++)
  {
    while (true)
    {
      // Randomize point [-1, 1]
      const double v = Random() * 2.0 * M_PI;
      const double r = Random() * m * R;
      const double x = C.x + r * std::cos(v);
      const double y = C.y + r * std::sin(v);
      Point2D p(x, y);
      // std::cout << "Trying to add building at p = " << p << std::endl;

      // Check that we are not too close to other buildings
      bool ok = true;
      for (auto const &q : centers)
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
        double width = (0.4 + 0.6 * Random()) * a;
        double height = (0.2 + 0.7 * Random()) * H;

        // Generate building
        Building building = GenerateBuilding(p, width, height);

        // Add building
        cityModel.Buildings.push_back(building);
        centers.push_back(p);

        std::cout << "Creating building at p = " << p << std::endl;

        break;
      }
    }
  }

  // Pack single model in list
  std::vector<CityModel> cityModels;
  cityModels.push_back(cityModel);

  return cityModels;
}

// Generate evolviing city model
std::vector<CityModel> GenerateEvolvingCityModel(size_t numBuildings,
                                                 size_t numFrames,
                                                 const Point2D &C,
                                                 double R)
{
  // Generate a random city model
  std::vector<CityModel> cityModels;
  cityModels = GenerateRandomCityModel(numBuildings, C, R);

  // Initialize building positions
  std::vector<Point2D> x(numBuildings);
  for (size_t i = 0; i < numBuildings; i++)
  {
    std::vector<Point2D> &footPrint = cityModels[0].Buildings[i].Footprint;
    for (auto const &p : footPrint)
      x[i] += p;
    x[i] /= footPrint.size();
  }

  // Compute center of gravity
  Point2D c;
  for (size_t i = 0; i < numBuildings; i++)
    c += x[i];
  c /= numBuildings;

  // Set radius for perimeter and minimal building distance
  const double D = BUILDING_SIZE;
  R -= BUILDING_SIZE;

  // Initialize building velocities
  std::vector<Point2D> v(numBuildings);
  for (size_t i = 0; i < numBuildings; i++)
  {
    v[i].x = VELOCITY * (Random() - 0.5);
    v[i].y = VELOCITY * (Random() - 0.5);
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
      const double dot = Geometry::Dot2D(x[i] - c, v[i]);
      if (r >= R && dot > 0.0)
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
      for (auto &p : cityModel.Buildings[i].Footprint)
        p += v[i];
    }
    cityModels.push_back(cityModel);
  }

  return cityModels;
}

int main(int argc, char *argv[])
{
  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }

  // Get filenames
  const std::string fileNameParameters(argv[1]);

  // Read parameters from file
  Parameters parameters;
  JSON::Read(parameters, fileNameParameters);
  std::cout << parameters << std::endl;

  // Report used parameters
  const Point2D C(parameters.DomainCenterX, parameters.DomainCenterY);
  const double R(parameters.DomainRadius);
  std::cout << "vc-randomize-citymodel: DomainCenter = " << C << std::endl;
  std::cout << "vc-randomize-citymodel: DomainRadius = " << R << std::endl;

  // FIXME: Handle command-line and selection of city model

  // Generate empty list of city models
  std::vector<CityModel> cityModels;

  // Generate city model
  cityModels = GenerateRandomCityModel(32, C, R);
  // cityModels = GenerateEvolvingCityModel(32, 1000, C, R);

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
