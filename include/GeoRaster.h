#ifndef DTCC_GEO_RASTER_H
#define DTCC_GEO_RASTER_H

#include <gdal/cpl_conv.h>
#include <gdal/gdal_priv.h>

#include "BoundingBox.h"
#include "GeoReference.h"
#include "GridField.h"
#include "Point.h"
#include "Logging.h"
namespace DTCC
{
// template<class T>
class GeoRaster
{
public:
  size_t XSize, YSize, Bands = 0;
  BoundingBox2D Bounds{};
  std::vector<GridField2D> Values{};
  // GeoReference geoReference = GeoReference();

  GeoRaster() {}

  GeoRaster(std::string fileName, std::vector<size_t> bands = {})
  {
    GDALDataset *rasterDataset;
    GDALRasterBand *rBand;
    size_t band_count;
    GDALAllRegister();
    rasterDataset = (GDALDataset *)GDALOpen(fileName.c_str(), GA_ReadOnly);
    if (rasterDataset == NULL)
    {
      throw std::runtime_error("Could not open file " + fileName);
    }
    XSize = rasterDataset->GetRasterXSize();
    YSize = rasterDataset->GetRasterYSize();
    band_count = rasterDataset->GetRasterCount();

    double gt[6];
    if (rasterDataset->GetGeoTransform(gt) != CE_None) {
      // file not georeferences
      gt[0] = 0;
      gt[1] = 1;
      gt[2] = 0;
      gt[3] = YSize;
      gt[4] = 0;
      gt[5] = -1;
    }


    Bounds.P.x = gt[0];
    Bounds.P.y = gt[3]+(YSize*gt[5]);
    Bounds.Q.x = gt[0]+(XSize*gt[1]);
    Bounds.Q.y = gt[3];
    
    

    if (bands.size() == 0) // load all bands in order
    {
      // bands 1-indexed in GDAL
      for (size_t b = 1; b <= band_count; b++)
      {
        bands.push_back(b);
      }
    }

    std::cout << "GeoRaster: Reading (" << XSize << "," << YSize << "," << bands.size() << ") raster" << std::endl;
    std::cout << "GeoRaster: Bounded by " << str(Bounds) << std::endl;

    for (auto band : bands)
    {
      std::cout << "loading band " << band << std::endl;
      GridField2D data(Grid2D(Bounds,XSize,YSize));
      rBand = rasterDataset->GetRasterBand(band);
      CPLErr err = rBand->RasterIO(GF_Read, 0, 0, XSize, YSize, &data.Values[0], XSize,
                                   YSize, GDT_Float64, 0, 0);
      if (err == CE_Failure)
      {
        throw std::runtime_error("Could not load data from " + fileName);
      }
      Values.push_back(data);
    }

    Bands = bands.size()+1; //Bands 1-indexed
    
    // accepted method of closing a GDAL dataset according to docs
    rasterDataset->~GDALDataset();
  }

  double operator()(const Point2D& p, size_t band = 1) const {
    if (band > Bands) 
    {
      throw std::runtime_error("Raster only has " + str(Bands) + " bands");
    }
    return Values[band-1].Nearest(p);
  }

  double Interpolate(const Point2D& p, size_t band = 1) {
    if (band > Bands) 
    {
      throw std::runtime_error("Raster only has " + str(Bands) + " bands");
    }
    return Values[band-1](p);
  }



};
} // namespace DTCC

#endif