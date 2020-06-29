// Copyright (C) 2020 Dag Wästberg
// Licensed under the MIT License

#include <gdal/cpl_conv.h>
#include <gdal/gdal_priv.h>

#include "GeoRaster.h"

#ifndef DTCC_RASTER_H
#define DTCC_RASTER_H
namespace DTCC
{
class RASTER
{
public:
  static void Read(GeoRaster &geoRaster, std::string fileName, std::vector<size_t> bands = {})
  {
    // we don't support adding data to a GeoRaster
    geoRaster.Values.clear();
    
    GDALDataset *rasterDataset;
    GDALRasterBand *rBand;
    size_t XSize, YSize, band_count;
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


    geoRaster.Bounds.P.x = gt[0];
    geoRaster.Bounds.P.y = gt[3]+(YSize*gt[5]);
    geoRaster.Bounds.Q.x = gt[0]+(XSize*gt[1]);
    geoRaster.Bounds.Q.y = gt[3];
    
    

    if (bands.size() == 0) // load all bands in order
    {
      // bands 1-indexed in GDAL
      for (size_t b = 1; b <= band_count; b++)
      {
        bands.push_back(b);
      }
    }

    std::cout << "GeoRaster: Reading (" << XSize << "," << YSize << "," << bands.size() << ") raster" << std::endl;
    std::cout << "GeoRaster: Bounded by " << str(geoRaster.Bounds) << std::endl;

    for (auto band : bands)
    {
      std::cout << "loading band " << band << std::endl;
      GridField2D data(Grid2D(geoRaster.Bounds,XSize,YSize));
      rBand = rasterDataset->GetRasterBand(band);
      CPLErr err = rBand->RasterIO(GF_Read, 0, 0, XSize, YSize, &data.Values[0], XSize,
                                   YSize, GDT_Float64, 0, 0);
      if (err == CE_Failure)
      {
        throw std::runtime_error("Could not load data from " + fileName);
      }
      geoRaster.Values.push_back(data);
    }

    geoRaster.Bands = bands.size(); 
    geoRaster.XSize = XSize;
    geoRaster.YSize = YSize;

    
    // accepted method of closing a GDAL dataset according to docs
    rasterDataset->~GDALDataset();

  }
};

} // namespace DTCC

#endif