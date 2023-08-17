/*  Smoother Profiling
*   This is a simple program to isolate and test Smoother code,
*   so we can profile and further optimize it.
*   - Notes: 
*     - The input used in this program are the volume meshes that are 
*     produced by dtcc-builder build_volume_mesh() in steps 3.2 and 3.4. They are saved 
*     in JSON format using a temp function I wrote in dtcc-builder 
*     builders.py
*
*     - I checked if the results are correct by comparing the residuals for the same input used in builder.
*       The residuals in both smoothings are equal to the ones we get when running whole builder. 
*
*     - When compiling without the -O3 flag the solver runs for about 5 minutes for 1000 iterations.
*       When we use the -O3 flag ita takes 30-40 seconds. So I think it would help optimization
*       to see what changes the compiler makes in Assembly.
*     
*/


#include <iostream>

#include "VolumeMesh.h"
#include "model/City.h"

#include "Smoother.h"

#include "include/utils.h"

using namespace DTCC_BUILDER;

void Help() { error("Usage: dtcc_builder_profiling Parameters.json"); }

int main(int argc, char *argv[])
{
  std::cout << "DTCC-Builder Smoother Test!" << std::endl;

  // Check command-line arguments
  if (argc != 2)
  {
    Help();
    return 1;
  }
  const std::string parameters_path = argv[1];
  Parameters p(parameters_path);
  p.print();

  // Reading city from JSON
  City city;
  JSON::Read(city, p.city_filenames[0]);

  // Reading DEM from JSON
  GridField dem;
  JSON::Read(dem, p.city_filenames[0]);

  // 1) Step 3.3 volume mesh smoothing
  {
    // Reading Step 3.4 Volume Mesh from JSON
    VolumeMesh volume_mesh{};

    JSON::Read(volume_mesh, p.volume_mesh_filenames[0]);
    const double top_height = p.domain_height[0] + dem.Mean();
    Smoother::smooth_volume_mesh(volume_mesh, city, dem, top_height,
                                 p.fix_buildings[0], p.max_iter[0],
                                 p.rel_tol[0]);
  }

  // 2) Step 3.5 volume mesh smoothing
  {
    // Reading Step 3.4 Volume Mesh from JSON
    VolumeMesh volume_mesh{};

    JSON::Read(volume_mesh, p.volume_mesh_filenames[1]);
    const double top_height = p.domain_height[1] + dem.Mean();
    Smoother::smooth_volume_mesh(volume_mesh, city, dem, top_height,
                                 p.fix_buildings[1], p.max_iter[1],
                                 p.rel_tol[1]);
  }
  return 0;
}
