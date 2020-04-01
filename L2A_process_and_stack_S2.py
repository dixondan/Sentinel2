#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 10:50:08 2019

@author: danieldixon
"""


import os
import zipfile
import subprocess
import fiona
import rasterio
import rasterio.mask
import geopandas as gpd
import sys

import numpy as np





#path = r"D:\#DATA\imagery\sentinel2\dandaragan\S2A_MSIL1C_20180102T022601_N0206_R060_T50JLL_20180102T085428\S2A_MSIL1C_20180102T022601_N0206_R060_T50JLL_20180102T085428.SAFE"
#path_L1C = r"D:\#DATA\imagery\sentinel2\dandaragan\S2A_MSIL1C_20180909T021341_N0206_R060_T50JLL_20180909T064151"
#cmd = "L2A_Process --resolution {} {}".format(resolution,path)
#os.system(cmd)

#zips = r"D:\#DATA\imagery\sentinel2\dandaragan"
#
#for i in os.listdir(zips):
#    if i.endswith(".zip"):
#        fullpath = os.path.join(zips, i)
#        outpath = os.path.join(zips, i.split(".")[0])
#        with zipfile.ZipFile(fullpath, 'r') as zip_ref:
#            zip_ref.extractall(outpath)
#    
#  
  




    
    
#safe = r"/Users/danieldixon/Desktop/dandaragan/S2_raw/S2A_MSIL1C_20181009T021341_N0206_R060_T50JLL_20181009T063846/S2A_MSIL1C_20181009T021341_N0206_R060_T50JLL_20181009T063846.SAFE"
#L2A = r"/Users/danieldixon/Sen2Cor-02.08.00-Darwin64/bin/L2A_Process"
#subprocess.call([L2A,  safe])
#
#






import sys
s2_raw = r"D:\#DATA\imagery\sentinel2\essentials"
all_jp2s = os.listdir(s2_raw)
for jp2 in all_jp2s:
    if ".DS" in jp2:
        pass
    else:
        #print(jp2)
        s2s = os.listdir(os.path.join(s2_raw, jp2))
        for s2 in s2s:
            pass
            if "MSIL2A" in s2:
                #sys.exit()
                get_base = os.listdir(os.path.join(s2_raw, jp2, s2, "GRANULE"))
                
                for base in get_base:
                    if ".DS" in base:
                        pass
                    else:
                        get_granules_path = os.path.join(s2_raw, jp2, s2, "GRANULE", base, "IMG_DATA", "R10m")
                        get_granules = os.listdir(get_granules_path)        
                        for granule in get_granules:
                            if "B0" in granule:
                                jpg2 = os.path.join(get_granules_path, granule)
                                print(os.path.exists(jpg2))                     
                                translate = r"C:\Users\22631228\AppData\Local\Continuum\anaconda3\Library\bin\gdal_translate.exe" 
                                out_tif = os.path.join(get_granules_path, granule.split(".")[0] + ".tif")
                                #sys.exit()
                                subprocess.call([translate, jpg2, out_tif])
            


band_dict = {"B02":0, 
             "B03":1,
             "B04":2,
             "B08":3}                   


shp1 = gpd.read_file(r"D:\#DATA\dandaragan\dandaragan\layers\POLYGON.shp")
shp1 = shp1.to_crs({'init': 'epsg:32750'})

shp1.to_file(r"D:\#DATA\dandaragan\dandaragan\layers\POLYGON_50s.shp")

shp = r"D:\#DATA\dandaragan\dandaragan\layers\POLYGON_50s.shp"
with fiona.open(shp, "r") as shapefile:
    features = [feature["geometry"] for feature in shapefile]

s2_raw = r"D:\#DATA\imagery\sentinel2\essentials"
all_jp2s = os.listdir(s2_raw)
for jp2 in all_jp2s:
    if ".DS" in jp2:
        pass
    else:
        #print(jp2)
        s2s = os.listdir(os.path.join(s2_raw, jp2))
        for s2 in s2s:
            pass
            if "MSIL2A" in s2:
                l=[]
                for root, dirs, files in os.walk(os.path.join(s2_raw, jp2)):
                    print(files)
                    for f in files:
                        if "10m.tif" in f:
                            l.append(os.path.join(root, f))
                            raster_2_fill = np.zeros((4, 182, 229))

                            for i in l:
                                pass
                                #create empty raster to fill
                                with rasterio.open(i) as src:
                                    out_image, out_transform = rasterio.mask.mask(src, features, crop=True)
                                    out_meta = src.meta.copy()
                                    out_image = out_image[0,:,:]
                                    print("outimage max = {}".format(out_image.max()))
                                    get_band = i.split("/")[-1].split("_")[-2]
                                    get_band = band_dict[get_band]
                                    
                                    print("raster band = {}".format(get_band))
                                    raster_2_fill[get_band,:,:] = out_image
                                
                                    out_meta.update({"driver": "GTiff",
                                                     "height": raster_2_fill.shape[1],
                                                     "width": raster_2_fill.shape[2],
                                                     "transform": out_transform, 
                                                     "count" : 4})
    
                            raster_name = i.split("/")[-1].split("_")[1][0:8]       
                            outpath = r"D:\#DATA\imagery\sentinel2\processed\tifs_4bands_dandaragan"                            
                            output_stacked = os.path.join(outpath, raster_name + ".tif")
                            
                            raster_2_fill = raster_2_fill.astype("uint16")
                            
                            with rasterio.open(output_stacked, "w", **out_meta) as dest:
                                dest.write(raster_2_fill)
                            
                               


ras = r"/Users/danieldixon/Desktop/dandaragan/S2_processed/20181128.tif"
    
with rasterio.open(ras) as src:
    #print(dir(src))
    array = src.read()
    print(src.crs)
    
    
    


r = r"/Users/danieldixon/Desktop/dandaragan/S2_processed/20181218.tif"



import geopandas as gpd
import os
import pandas as pd
import rasterio
import fiona
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from osgeo import gdal, ogr
from collections import OrderedDict
#import time
import itertools
import rasterio.mask
def makemydir(newfile):
  try:
    os.makedirs(newfile)
  except OSError:
    pass
  os.chdir(newfile)
import seaborn as sns
from osgeo import gdal
import numpy as np
from shapely.geometry import Point
from shapely.geometry import Polygon
import scipy
import warnings
import functools

warnings.filterwarnings("ignore")

#converting polygon centroid to square geometry for output
def centroid_2_polygon(row, x_size, y_size):
    x1 = row['x']
    y1 = row['y']
    halfx = np.abs(x_size  / 2)
    halfy = np.abs(y_size / 2) 
    topleft =  (x1 - halfx, y1 + halfy)
    topright = (x1 + halfx, y1 + halfy)
    botright = (x1 + halfx, y1 - halfy)
    botleft = (x1 - halfx, y1 - halfy)    
    poly = Polygon(np.array([topleft, topright, botright, botleft, topleft]))
    return (poly)

#converting raster to polygon fishnet
def ras_2_fish(ras, outname):
    r = gdal.Open(ras)
    band = r.GetRasterBand(1) #bands start at one
    a = band.ReadAsArray().astype(np.float)
    (y_index, x_index) = np.nonzero(a > 0)
    (upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = r.GetGeoTransform()
    x_coords = x_index * x_size + upper_left_x + (x_size / 2) #add half the cell size
    y_coords = y_index * y_size + upper_left_y + (y_size / 2)  
    df = pd.DataFrame()
    df['x'] = x_coords
    df['y'] = y_coords
    geometry = [Point(xy) for xy in zip(df.x, df.y)]
    crs = {'init': 'epsg:32750'}
    gdf = gpd.GeoDataFrame(df, geometry = geometry, crs = crs)
    gdf['geometry'] = gdf.apply(lambda row: centroid_2_polygon(row, x_size, y_size),axis=1)
    gdf['ID1'] = [str(i) + "_ID" for i in range(0, len(gdf))]
    outname = outname
    gdf = gdf.drop(['x', 'y'], axis=1)    
    #add the final intersection only include those squares that intersect  study_range_inside.shp
    #study_range_inside = gpd.read_file(r"/Users/danieldixon/Desktop/dandaragan/DRONE_IMAGE_COREGISTRATION/study_range_inside.shp")
    #gdf_subset = gpd.overlay(gdf,study_range_inside, how='intersection')
    gdf.to_file(outname)




r = r"/Users/danieldixon/Desktop/dandaragan/S2_processed/20181218.tif"
outname = r"/Users/danieldixon/Documents/research/PhD/paper_1/raw_data/S2/shapefiles/s2_fishnet.shp"
ras_2_fish(r, outname)
                
                




#
#    
#gdal_translate = r"C:\Users\22631228\AppData\Local\Continuum\miniconda3\pkgs\libgdal-2.3.3-h10f50ba_0\Library\bin\gdal_translate.exe"
#    
#cmd = "gdal_translate.exe -a_srs EPSG:32750 -of Gtiff" 
#
#
#jp2 = r"D:\#DATA\imagery\sentinel2\dandaragan\S2A_MSIL1C_20181108T021341_N0207_R060_T50JLL_20181108T052735\S2A_MSIL1C_20181108T021341_N0207_R060_T50JLL_20181108T052735.SAFE\GRANULE\L1C_T50JLL_A017646_20181108T022545\IMG_DATA\T50JLL_20181108T021341_B04.jp2"
#tif = r"D:\#DATA\imagery\sentinel2\dandaragan\S2A_MSIL1C_20181108T021341_N0207_R060_T50JLL_20181108T052735\S2A_MSIL1C_20181108T021341_N0207_R060_T50JLL_20181108T052735.SAFE\GRANULE\L1C_T50JLL_A017646_20181108T022545\IMG_DATA\T50JLL_20181108T021341_B04.tif"
#
#
#command = cmd + " " + jp2 + " " + tif
#


#
#
#r = r"D:\#DATA\imagery\sentinel2\dandaragan\S2B_MSIL1C_20180825T021339_N0206_R060_T50JLL_20180825T061521\S2B_MSIL1C_20180825T021339_N0206_R060_T50JLL_20180825T061521.SAFE"
#    
#import os
#
#os.environ.copy()
#
#
#import qgis
#os.environ = environment_variables
##have to run the gdal_translate from jp2 to tif and get the correct epsg output
##may have to export to wgs84, or try 50s, either way that needs a projection
##before stacking
#                    
                    
    
#
#jp2 = r"D:\#DATA\imagery\sentinel2\essentials\S2A_MSIL1C_20181009T021341_N0206_R060_T50JLL_20181009T063846\S2A_MSIL2A_20181009T021341_N0206_R060_T50JLL_20181009T063846.SAFE\GRANULE\L2A_T50JLL_A017217_20181009T022112\IMG_DATA\R10m\T50JLL_20181009T021341_B03_10m.jp2"
#
#gtiff = "-a_srs EPSG:32750 -of Gtiff"                                
#                    out_tif = os.path.join(get_granules_path, granule.split(".")[0] + ".tif")             
#                    command = gdal_translate + " " + gtiff + " " + jpg2 + " " + out_tif           
#
#    
#    
#tif = r"D:\#DATA\imagery\sentinel2\essentials\S2B_MSIL1C_20181213T021339_N0207_R060_T50JLL_20181213T052802\S2B_MSIL2A_20181213T021339_N0207_R060_T50JLL_20181213T052802.SAFE\GRANULE\L2A_T50JLL_A009238_20181213T022159\IMG_DATA\R10m\T50JLL_20181213T021339_B02_10m.tif"
#
#
#
#    

#crs = rasterio.crs.CRS({"init": "epsg:32750"})    # or whatever CRS you know the image is in    
    
    
    
#    
#
#
#import rasterio 
# 
#ras = r"D:\#DATA\imagery\sentinel2\dandaragan\S2A_MSIL1C_20181108T021341_N0207_R060_T50JLL_20181108T052735\S2A_MSIL1C_20181108T021341_N0207_R060_T50JLL_20181108T052735.SAFE\GRANULE\L1C_T50JLL_A017646_20181108T022545\IMG_DATA\T50JLL_20181108T021341_B03.jp2"
#    
#with rasterio.open(ras) as src:
#    #print(dir(src))
#    array = src.read()
#    print(src.crs)
#    
#    
    
    
    
    
    
    
    

shp2  = gpd.read_file(r"D:\#DATA\shapefiles\dandaragan_broad\layers\POLYGON.shp")                   
shp2.head()




shp = r"D:\#DATA\shapefiles\dandaragan_broad\layers\POLYGON.shp"
with fiona.open(shp, "r") as shapefile:
    features = [feature["geometry"] for feature in shapefile]

#want to get 4 band stacked tifs clip them later
base = r"D:\#DATA\imagery\sentinel2\essentials"
all_jp2s = os.listdir(base)
for jp2 in all_jp2s:
    #print(jp2)
    s2s = os.listdir(os.path.join(base, jp2))
    for s2 in s2s:
        pass
        if "MSIL2A" in s2:
            #sys.exit()
            get_base = os.listdir(os.path.join(base, jp2, s2, "GRANULE"))[0]
            get_granules_path = os.path.join(base, jp2, s2, "GRANULE", get_base, "IMG_DATA", "R10m")
            get_granules = os.listdir(get_granules_path)        
            temp_arrays = []
            for granule in get_granules:
                if "B0" in granule:
                    if granule.endswith(".tif"):
                        print(granule)
                        tif = os.path.join(get_granules_path, granule)
                        band = granule.split(".")[0].split("_")[-2]
                        print(band)

                       
                        with rasterio.open(tif) as src:
                            out_image, out_transform = rasterio.mask.mask(src, features, crop=True)
                            out_meta = src.meta.copy()
    
    
                        out_meta.update({"driver": "GTiff",
                                         "height": out_image.shape[1],
                                         "width": out_image.shape[2],
                                         "transform": out_transform})
    
    
    
                        outpath = r"D:\#DATA\imagery\sentinel2\processed\tifs_4bands_dandaragan"

                        output_stacked = os.path.join(outpath, )
                        with rasterio.open("RGB.byte.masked.tif", "w", **out_meta) as dest:
                            dest.write(out_image)
                        
                        
      






#mport gdal, gdalconst, osr
#from osgeo.gdalconst import *
#
#
#driver = gdal.GetDriverByName("GTiff")
#
#
#
#srs = osr.SpatialReference()
#srs.ImportFromEPSG(32750)
#
#
#
#f = r"D:\#DATA\imagery\sentinel2\essentials\S2B_MSIL1C_20181213T021339_N0207_R060_T50JLL_20181213T052802\S2B_MSIL2A_20181213T021339_N0207_R060_T50JLL_20181213T052802.SAFE\GRANULE\L2A_T50JLL_A009238_20181213T022159\IMG_DATA\R10m\T50JLL_20181213T021339_B03_10m.tif"
#ds = gdal.Open(f, gdal.GA_Update)
#
#
#
#
#ds.SetProjection(srs.ExportToWkt())
#
#
#
#22



