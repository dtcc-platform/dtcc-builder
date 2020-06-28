#ifndef DTCC_GEO_RASTER_H
#define DTCC_GEO_RASTER_H

#include <gdal/cpl_conv.h>
#include <gdal/gdal_priv.h>

#include "GeoReference.h"
#include "Point.h"
namespace DTCC
{
// template<class T>
class GeoRaster
{
public:
  size_t XSize, YSize, Bands;
  std::vector<std::vector<double>> Values{};
  GeoReference geoReference = GeoReference();


  GeoRaster(std::string fileName, std::vector<size_t> bands = {}) 
  {
    GDALDataset *rasterDataset;
    GDALRasterBand  *rBand;
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
    if (rasterDataset->GetGeoTransform(gt) == CE_None)
    {
      geoReference.A = gt[1];
      geoReference.B = gt[2];
      geoReference.C = gt[0];
      geoReference.D = gt[4];
      geoReference.E = gt[5];
      geoReference.F = gt[3];
    }

    if (bands.size() == 0) // load all bands in order
    {
      //bands 1-indexed in GDAL
      for (size_t b = 1; b<=band_count; b++) 
      {
        bands.push_back(b);
      }
    }

    for (auto band: bands) 
    {
      std::cout << "loading band " << band << std::endl;
      std::vector<double> data(XSize * YSize);
      rBand = rasterDataset->GetRasterBand(band);
      CPLErr err = rBand->RasterIO( GF_Read, 0,0,XSize,YSize, &data[0], XSize, YSize, GDT_Float64,0,0   );
      if (err == CE_Failure)
      {
        throw std::runtime_error("Could not load data from " + fileName);
      }
      Values.push_back(data);
    }

    Bands = bands.size();


  }


};
} // namespace DTCC

#endif