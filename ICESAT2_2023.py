#!/usr/bin/env python
"""
pySL4Land - this file provides a class to process ICESAT2 data.

See other source files for details
"""
# This file is part of 'pySL4Land'
# A set of tools to process spaceborne lidar (GEDI and ICESAT2) for land (pySL4Land) applications
#
# Copyright 2020 Pete Bunting
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# Purpose: class to process ICESAT2 data.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 25/06/2020
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import h5py
import numpy
import pandas as pd
import geopandas
import math
import os
import arcpy
from arcpy.sa import *
from shapely.geometry import Polygon
import logging
import pysl4land.pysl4land_utils
import subprocess
from whitebox_tools import WhiteboxTools
import laspy
wbt = WhiteboxTools()

arcpy.env.overwriteOutput = True

logger = logging.getLogger(__name__)


def get_beam_lst(input_file, strong_only, weak_only):
    """
    A function which returns a list of beam names.

    :param input_file: input file path.
    :return: list of strings

    """
    icesat2_h5_file = h5py.File(input_file, 'r')

    orientation = icesat2_h5_file['/orbit_info/sc_orient'][0]

    if strong_only == True:
        icesat2_keys = list(icesat2_h5_file.keys())
        strongOrientDict = {0: 'l', 1: 'r', 21: 'error'}
        icesat2_beams = ['gt1' + strongOrientDict[orientation], 'gt2' + strongOrientDict[orientation],
                         'gt3' + strongOrientDict[orientation]]
        icesat2_beams_lst = []
        for icesat2_beam_name in icesat2_keys:
            if icesat2_beam_name in icesat2_beams:
                icesat2_beams_lst.append(icesat2_beam_name)
        icesat2_h5_file.close()

    elif weak_only == True:
        icesat2_keys = list(icesat2_h5_file.keys())
        weakOrientDict = {0: 'r', 1: 'l', 21: 'error'}
        icesat2_beams = ['gt1' + weakOrientDict[orientation], 'gt2' + weakOrientDict[orientation],
                         'gt3' + weakOrientDict[orientation]]
        icesat2_beams_lst = []
        for icesat2_beam_name in icesat2_keys:
            if icesat2_beam_name in icesat2_beams:
                icesat2_beams_lst.append(icesat2_beam_name)
        icesat2_h5_file.close()

    else:
        icesat2_keys = list(icesat2_h5_file.keys())
        icesat2_beams = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
        icesat2_beams_lst = []
        for icesat2_beam_name in icesat2_keys:
            if icesat2_beam_name in icesat2_beams:
                icesat2_beams_lst.append(icesat2_beam_name)
        icesat2_h5_file.close()

    return icesat2_beams_lst


def get_segment_polygons(latitude, longitude, along_size=20.0, across_size=13.0):
    """
    Get polygon segments.

    """
    n_rows = latitude.shape[0]

    utm_zone = pysl4land.pysl4land_utils.latlon_to_mode_utm_zone_number(latitude, longitude)
    logger.debug("UTM Zone: {}".format(utm_zone))

    epsg_code = 32600 + utm_zone
    logger.debug("UTM Zone EPSG: {}".format(epsg_code))

    logger.debug("Creating pandas dataframe")
    icesat2_pts_df = pd.DataFrame({'latitude': latitude, 'longitude': longitude})
    logger.debug("Creating geopandas dataframe")
    icesat2_pts_wgs84_gdf = geopandas.GeoDataFrame(icesat2_pts_df, crs='EPSG:4326',
                                                   geometry=geopandas.points_from_xy(icesat2_pts_df.longitude,
                                                                                     icesat2_pts_df.latitude))
    logger.debug("Reprojecting geopandas dataframe to UTM")
    icesat2_pts_utm_gdf = icesat2_pts_wgs84_gdf.to_crs("EPSG:{}".format(epsg_code))

    x = numpy.array(icesat2_pts_utm_gdf.geometry.x)
    y = numpy.array(icesat2_pts_utm_gdf.geometry.y)

    # Step Forward
    x_fwd = numpy.zeros_like(x)
    x_fwd[1:] = x[0:n_rows - 1]
    x_fwd[0] = x[0]

    y_fwd = numpy.zeros_like(y)
    y_fwd[1:] = y[0:n_rows - 1]
    y_fwd[0] = y[0]

    # Step Backwards
    x_bck = numpy.zeros_like(x)
    x_bck[0:n_rows - 1] = x[1:n_rows]
    x_bck[-1] = x[-1]

    y_bck = numpy.zeros_like(y)
    y_bck[0:n_rows - 1] = y[1:n_rows]
    y_bck[-1] = y[-1]

    logger.debug("Calculate the angles.")
    sgl_flight_angle = math.atan((x[-1] - x[0]) / (y[-1] - y[0]))

    flight_angle = numpy.arctan((x_fwd - x_bck) / (y_fwd - y_bck))

    flight_angle_dif = numpy.absolute(flight_angle - sgl_flight_angle)
    flight_angle[flight_angle_dif > 0.1] = sgl_flight_angle

    theta_cos = numpy.cos(flight_angle)
    theta_sin = numpy.sin(flight_angle)

    logger.debug("Calculated the angles.")

    along_size_h = along_size / 2.0
    across_size_h = across_size / 2.0
    logger.debug("along_size_h: {} \t\t across_size_h: {}".format(along_size_h, across_size_h))

    icesat2_pts_utm_gdf['ul_x'] = icesat2_pts_utm_gdf.geometry.x - across_size_h
    icesat2_pts_utm_gdf['ul_y'] = icesat2_pts_utm_gdf.geometry.y + along_size_h
    icesat2_pts_utm_gdf['ul_x_rot'] = icesat2_pts_utm_gdf.geometry.x + theta_cos * (
            icesat2_pts_utm_gdf['ul_x'] - icesat2_pts_utm_gdf.geometry.x) + theta_sin * (
                                              icesat2_pts_utm_gdf['ul_y'] - icesat2_pts_utm_gdf.geometry.y)
    icesat2_pts_utm_gdf['ul_y_rot'] = icesat2_pts_utm_gdf.geometry.y + -theta_sin * (
            icesat2_pts_utm_gdf['ul_x'] - icesat2_pts_utm_gdf.geometry.x) + theta_cos * (
                                              icesat2_pts_utm_gdf['ul_y'] - icesat2_pts_utm_gdf.geometry.y)

    icesat2_pts_utm_gdf['ur_x'] = icesat2_pts_utm_gdf.geometry.x + across_size_h
    icesat2_pts_utm_gdf['ur_y'] = icesat2_pts_utm_gdf.geometry.y + along_size_h
    icesat2_pts_utm_gdf['ur_x_rot'] = icesat2_pts_utm_gdf.geometry.x + theta_cos * (
            icesat2_pts_utm_gdf['ur_x'] - icesat2_pts_utm_gdf.geometry.x) + theta_sin * (
                                              icesat2_pts_utm_gdf['ur_y'] - icesat2_pts_utm_gdf.geometry.y)
    icesat2_pts_utm_gdf['ur_y_rot'] = icesat2_pts_utm_gdf.geometry.y + -theta_sin * (
            icesat2_pts_utm_gdf['ur_x'] - icesat2_pts_utm_gdf.geometry.x) + theta_cos * (
                                              icesat2_pts_utm_gdf['ur_y'] - icesat2_pts_utm_gdf.geometry.y)

    icesat2_pts_utm_gdf['lr_x'] = icesat2_pts_utm_gdf.geometry.x + across_size_h
    icesat2_pts_utm_gdf['lr_y'] = icesat2_pts_utm_gdf.geometry.y - along_size_h
    icesat2_pts_utm_gdf['lr_x_rot'] = icesat2_pts_utm_gdf.geometry.x + theta_cos * (
            icesat2_pts_utm_gdf['lr_x'] - icesat2_pts_utm_gdf.geometry.x) + theta_sin * (
                                              icesat2_pts_utm_gdf['lr_y'] - icesat2_pts_utm_gdf.geometry.y)
    icesat2_pts_utm_gdf['lr_y_rot'] = icesat2_pts_utm_gdf.geometry.y + -theta_sin * (
            icesat2_pts_utm_gdf['lr_x'] - icesat2_pts_utm_gdf.geometry.x) + theta_cos * (
                                              icesat2_pts_utm_gdf['lr_y'] - icesat2_pts_utm_gdf.geometry.y)

    icesat2_pts_utm_gdf['ll_x'] = icesat2_pts_utm_gdf.geometry.x - across_size_h
    icesat2_pts_utm_gdf['ll_y'] = icesat2_pts_utm_gdf.geometry.y - along_size_h
    icesat2_pts_utm_gdf['ll_x_rot'] = icesat2_pts_utm_gdf.geometry.x + theta_cos * (
            icesat2_pts_utm_gdf['ll_x'] - icesat2_pts_utm_gdf.geometry.x) + theta_sin * (
                                              icesat2_pts_utm_gdf['ll_y'] - icesat2_pts_utm_gdf.geometry.y)
    icesat2_pts_utm_gdf['ll_y_rot'] = icesat2_pts_utm_gdf.geometry.y + -theta_sin * (
            icesat2_pts_utm_gdf['ll_x'] - icesat2_pts_utm_gdf.geometry.x) + theta_cos * (
                                              icesat2_pts_utm_gdf['ll_y'] - icesat2_pts_utm_gdf.geometry.y)

    def _polygonise_2Dcells(df_row):
        return Polygon([(df_row.ul_x_rot, df_row.ul_y_rot), (df_row.ur_x_rot, df_row.ur_y_rot),
                        (df_row.lr_x_rot, df_row.lr_y_rot), (df_row.ll_x_rot, df_row.ll_y_rot)])

    polys = icesat2_pts_utm_gdf.apply(_polygonise_2Dcells, axis=1)

    logger.debug("Creating geopandas polygons UTM dataframe")
    icesat2_polys_utm_gdf = geopandas.GeoDataFrame(icesat2_pts_df, crs="EPSG:{}".format(epsg_code), geometry=polys)
    logger.debug("Creating geopandas polygons WGS84 dataframe")
    icesat2_polys_wgs84_gdf = icesat2_polys_utm_gdf.to_crs("EPSG:4326".format(epsg_code))

    return numpy.asarray(icesat2_polys_wgs84_gdf.geometry.values)


def get_icesat2_alt08_beam_as_gdf(input_file, icesat2_beam_name, use_seg_polys=True, out_epsg_code=4326,
                                  strong_only=False, weak_only=False):
    """

    :param input_file:
    :param icesat2_beam_name:
    :param use_seg_polys:
    :param out_epsg_code:
    :return:
    """
    icesat2_beams = get_beam_lst(input_file, strong_only, weak_only)
    if icesat2_beam_name not in icesat2_beams:
        raise Exception("Beam '{}' is not available within the file: {}".format(icesat2_beam_name, input_file))

    icesat2_h5_file = h5py.File(input_file, 'r')
    if icesat2_h5_file is None:
        raise Exception("Could not open the input ICESAT2 file.")

    icesat2_beam = icesat2_h5_file[icesat2_beam_name]
    icesat2_beam_keys = list(icesat2_beam.keys())
    if 'land_segments' not in icesat2_beam_keys:
        raise Exception("Could not find land segments information.")

    icesat2_land_beam = icesat2_beam['land_segments']
    icesat2_land_beam_keys = list(icesat2_land_beam.keys())

    # Get canopy data.
    if 'canopy' not in icesat2_land_beam_keys:
        raise Exception("Could not find canopy information.")
    icesat2_beam_canopy = icesat2_land_beam['canopy']

    # Get terrain data.
    if 'terrain' not in icesat2_land_beam_keys:
        raise Exception("Could not find terrain information.")
    icesat2_beam_terrain = icesat2_land_beam['terrain']

    icesat2_beam_df = pd.DataFrame({
        'h_te_uncertainty': icesat2_beam_terrain['h_te_uncertainty'],
        'h_te_best_fit': icesat2_beam_terrain['h_te_best_fit'],
        'h_te_median': icesat2_beam_terrain['h_te_median'],
        'h_canopy': icesat2_beam_canopy['h_canopy'],
        'h_canopy_uncertainty': icesat2_beam_canopy['h_canopy_uncertainty'],
        'n_ca_photons': icesat2_beam_canopy['n_ca_photons'],
        'night_flag': icesat2_land_beam['night_flag']
    })

    df_list = []
    df_concat = pd.DataFrame(df_list)
    for i in range(5):
        icesat2_beam_df_20 = pd.DataFrame({
            'latitude_20m': icesat2_land_beam['latitude_20m'][:, i],
            'longitude_20m': icesat2_land_beam['longitude_20m'][:, i],
            'canopy_20m': icesat2_beam_canopy['h_canopy_20m'][:, i],
            'terrian_20m': icesat2_beam_terrain['h_te_best_fit_20m'][:, i]})
        df_vertical = pd.concat([icesat2_beam_df, icesat2_beam_df_20], axis=1)
        df_concat = pd.concat([df_concat, df_vertical])

    df_concat[df_concat > 100000] = -999
    df_concat = df_concat.loc[(df_concat["latitude_20m"] != -999) & (df_concat["longitude_20m"] != -999)]
    if use_seg_polys:
        latitude_arr = numpy.array(df_concat['latitude_20m'])
        longitude_arr = numpy.array(df_concat['longitude_20m'])

        polys = get_segment_polygons(latitude_arr, longitude_arr, along_size=20.0, across_size=13.0)
        icesat2_beam_gdf = geopandas.GeoDataFrame(df_concat, crs='EPSG:4326', geometry=polys)
    else:
        icesat2_beam_gdf = geopandas.GeoDataFrame(df_concat, crs='EPSG:4326',
                                                  geometry=geopandas.points_from_xy(df_concat.longitude_20m,
                                                                                    df_concat.latitude_20m))

    if out_epsg_code != 4326:
        icesat2_beam_gdf = icesat2_beam_gdf.to_crs("EPSG:{}".format(out_epsg_code))
    icesat2_h5_file.close()

    return icesat2_beam_gdf


def icesat2_alt08_beams_gpkg(input_file, out_vec_dir, use_seg_polys=True, out_epsg_code=4326, strong_only=False,
                             weak_only=False):
    icesat2_beams = get_beam_lst(input_file, strong_only, weak_only)
    print(icesat2_beams)
    for icesat2_beam_name in icesat2_beams:
        out_vec_file = out_vec_dir + icesat2_beam_name + ".shp"
        logger.info("Processing beam '{}'".format(icesat2_beam_name))
        icesat2_beam_gdf = get_icesat2_alt08_beam_as_gdf(input_file, icesat2_beam_name, use_seg_polys, out_epsg_code)
        icesat2_beam_gdf.to_file(out_vec_file, layer=icesat2_beam_name)
        logger.info("Finished processing beam '{}'".format(icesat2_beam_name))




def clip_shp():
    sites = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia', 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    dir = r"E:\ICESat2\Results\SHP_Misha\Individual/"
    out = r"E:\ICESat2\Results\SHP_Misha\Clip/"
    for site in sites:
        site_dir = "E:\GEDI\Boundingbox/" + site + ".shp"
        dir_site = dir + site + "/"
        out_site = out + site + "/"
        arcpy.env.workspace = dir_site
        fcs = arcpy.ListFeatureClasses()
        for fc in fcs:
            try:
                print("Clipping shapefile..." + fc)
                arcpy.MakeFeatureLayer_management(dir_site + fc, 'fc_lyr')
                arcpy.MakeFeatureLayer_management(site_dir, 'site_lyr')
                arcpy.SelectLayerByLocation_management('fc_lyr', 'intersect', 'site_lyr')
                arcpy.FeatureClassToShapefile_conversion('fc_lyr', out_site)
                arcpy.Rename_management(out_site + "fc_lyr.shp", out_site + fc)
            except:
                print("Fail to clip, skipped")


def to_shp_merge():
    # find track IDs
    dir = "E:\ICESat2\Results\SHP_Misha\Clip/"
    out = r"E:\ICESat2\Results\SHP_Misha\Merge/"
    sites = os.listdir(dir)
    for site in sites:
        print("processing..." + site)
        arcpy.env.workspace = dir + site
        output = out + "ICESat_" + site + ".shp"
        fcs = arcpy.ListFeatureClasses()
        for fc in fcs:
            print("processing..." + fc)
            arcpy.AddField_management(fc, 'filename', "TEXT")
            with arcpy.da.UpdateCursor(fc, ["filename"]) as cursor:
                for row in cursor:
                    row[0] = fc.split(".")[0]
                    cursor.updateRow(row)
            arcpy.AddField_management(fc, 'site', "TEXT")
            with arcpy.da.UpdateCursor(fc, ["site"]) as cursor:
                for row in cursor:
                    row[0] = site
                    cursor.updateRow(row)
        print("Output + " + output)
        arcpy.Merge_management(fcs, output)


def merge_SJ():
    dir = r"E:\ICESat2\Results\SHP_Misha\Merge/"
    tile_shp = r"E:\GEDI\Boundingbox/"
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    for site in studysite:
        print(site)
        targetFeatures = dir + "ICESat_" + site + ".shp"
        joinFeatures = tile_shp + site + ".shp"
        outfc = dir + "ICESat_" + site + "_SJ.shp"
        arcpy.analysis.SpatialJoin(targetFeatures, joinFeatures, outfc)
        outfc2 = dir + "ICESat_" + site + "_SJ_P.shp"
        if "Addo" in site or "Venetia" in site or "DNyala" in site:
            Proj = tile_shp + "UTM35S.prj"
        elif "Agulhas" in site:
            Proj = tile_shp + "UTM34S.prj"
        else:
            Proj = tile_shp + "UTM36S.prj"
        arcpy.Project_management(outfc, outfc2, Proj)


def basic_conversion():
    dir = "E:\ICESat2\Results\SHP_Misha\Merge/"
    out = "E:\ICESat2\Results\Basic_ICESAT2_Misha/"
    files = os.listdir(dir)
    for file in files:
        if "_SJ_P.shp" in file and "xml" not in file and "lock" not in file:
            print(file)
            xls = out + file.split("_")[1] + "_ICESat2.xlsx"
            fc = dir + file
            arcpy.TableToExcel_conversion(fc, xls)



# If geopandas not working (pyproj), please uninstall and reinstall pyproj, then Fiona/pyogrio
def basic_main():

    dir = r"E:\ICESat2\ICESAT2_archive\ATL08/"
    out = r"E:\ICESat2\Results\SHP_Misha\Individual/"
    sites = os.listdir(dir)
    for site in sites:
        path = dir + site + "/"
        files = os.listdir(path)
        for file in files:
            print(site)
            print(file)
            if not os.path.isdir(path + file.split(".")[0] + "/"):
                print("Processing + " + file)
                icesat2_alt08_beams_gpkg(path + file, out,
                                         use_seg_polys=True, out_epsg_code=4326,
                                         strong_only=True, weak_only=False)
                matches = []
                gt_files = os.listdir(out)
                for gt_file in gt_files:
                    if "gt" in gt_file and ".shp" in gt_file:
                        matches.append(out + gt_file)
                merge_file = file.split(".")[0] + ".shp"
                arcpy.Merge_management(matches, out + site + "/" + merge_file)
                arcpy.env.workspace = out
                for f in matches:
                    arcpy.management.Delete(f)
                # arcpy.management.Delete("Output.shp")
            else:
                print("File exists, skipped: " + file)

    #Clip shapefile to study site boundary
    clip_shp()
    #merge tracks to sites
    to_shp_merge()
    #Add lidar tiles to merge files
    merge_SJ()
    #Convert shapefile to csv
    basic_conversion()



def P98_lid_clip():
    dir = r"E:\ICESat2\Results\SHP\Merge/"
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    for site in studysite:
        dir_shp = dir + "ICESat_" + site + "_SJ_P.shp"
        arcpy.MakeFeatureLayer_management(dir_shp, 'SHP')
        with arcpy.da.SearchCursor('SHP', ['LAS_NM', 'site', 'TARGET_FID']) as cursor:
            for row in cursor:
                if row[1] == site:
                    if row[1] == 'Welverdiendt' and row[2] < 11000:
                        print("Skipped")
                    else:
                        dir_lidar = r"E:\ALS_archive\LiDAR_UTM/" + row[1] + "/" + row[0].replace(" ", "") + ".las"
                        output_las = r"E:\ICESat2\Results\ALS/" + row[1] + "/" + str(row[2]) + ".las"
                        query = """ "TARGET_FID" = {0}""".format(row[2])
                        arcpy.SelectLayerByAttribute_management("SHP", 'NEW_SELECTION', query)
                        dir_shp_temp = r'E:\ICESat2\Results\ALS/'
                        arcpy.FeatureClassToShapefile_conversion("SHP", dir_shp_temp)
                        wbt.clip_lidar_to_polygon(i=dir_lidar, polygons=dir_shp_temp + "SHP.shp", output=output_las)


def las_processing():
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    for site in studysite:
        list_shot = []
        list_P98 = []
        dir = r"E:\ICESat2\Results\ALS/" + site + "/"
        lastool_ground = r"D:\lastool\LAStools\bin\lasground.exe "
        lastool_height = r"D:\lastool\LAStools\bin\lasheight.exe "
        lastool_las = r"D:\lastool\LAStools\bin\las2las.exe "
        lastool_thin = r"D:\lastool\LAStools\bin\lasthin.exe "
        output_ground = r"E:\ICESat2\Results\P98\ground.las"
        output_height = r"E:\ICESat2\Results\P98\normalized.las"
        output_las = r"E:\ICESat2\Results\P98\filtered.las"
        files = os.listdir(dir)
        for file in files:
            if int(file.split(".")[0]) < 4000:
                print("skipped")
            else:
                print(site + " " + file)
                file_dir = dir + file
                output_thin = r"E:\ICESat2\Results\P98/" + site + "/" + file
                subprocess.call(lastool_ground + " -i " + file_dir + " -step 1 -o " + output_ground)
                subprocess.call(lastool_height + " -i " + output_ground + " -replace_z -o " + output_height)
                subprocess.call(lastool_las + " -i " + output_height + " -drop_class 2 -drop_z_below 0.01 -o " + output_las)
                subprocess.call(
                    lastool_thin + " -i " + output_las + " -step 1000 -percentile 98 -ignore_class 2 -o " + output_thin)
                try:
                    las = laspy.read(output_thin)
                    list_P98.append(max(las.Z / 1000))
                    list_shot.append(file.split(".")[0])
                except:
                    print("Error, skipped")

        df_shot = pd.DataFrame(list_shot)
        df_p98 = pd.DataFrame(list_P98)

        df = pd.concat([df_shot, df_p98], axis=1)
        df.columns = ["shot", "p98"]
        output_csv = r"E:\ICESat2\Results\Basic_ICESAT2_Misha/" + site + "_P98.csv"
        df.to_csv(output_csv, index=None)


def combine_basic_p98():
    list_testcase = []
    df = pd.DataFrame(list_testcase)
    out = r"E:\ICESat2\Results\Result\All_ICESat2.csv"
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    for site in studysite:
        print(site)
        dir1 = r"E:\ICESat2\Results\Basic_ICESAT2/" + site + "_P98.csv"
        dir2 = r"E:\ICESat2\Results\Basic_ICESAT2/" + site + "_ICESat2.xlsx"
        df1 = pd.read_csv(dir1)
        df2 = pd.read_excel(dir2)

        df_merge = pd.merge(df2, df1, left_on='TARGET_FID', right_on='shot', how='left')

        df = pd.concat([df, df_merge])
    df.to_csv(out, index=None)


def add_phenology():
    csv = r"E:\ICESat2\Results\Result/All_ICESat2_P_Misha.csv"
    out = r"E:\ICESat2\Results\Result/All_ICESat2_P_Misha_P.csv"
    df = pd.read_csv(csv,header=0)
    df['status'] = 0

    df.loc[(df['filename'].str[11:14].astype(int) >= 1) &
           (df['filename'].str[11:14].astype(int) <= 450), 'status'] = 'Leaf-on'
    df.loc[(df['filename'].str[11:14].astype(int) >= 1100) &
           (df['filename'].str[11:14].astype(int) <= 1250), 'status'] = 'Leaf-on'

    df.loc[(df['filename'].str[11:14].astype(int) >= 600) &
           (df['filename'].str[11:14].astype(int) <= 950), 'status'] = 'Leaf-off'

    df.loc[(df['filename'].str[11:14].astype(int) >= 500) &
           (df['filename'].str[11:14].astype(int) <= 550), 'status'] = 'Transition'
    df.loc[(df['filename'].str[11:14].astype(int) >= 1000) &
           (df['filename'].str[11:14].astype(int) <= 1050), 'status'] = 'Transition'
    print("All database output...")
    df.to_csv(out, index=None)



def merge_csv_final():
    Input_dir = r"E:\ICESat2\Results\Basic_ICESAT2_Misha/"
    Result_dir = r"E:\ICESat2\Results\Result/"
    out_csv = Result_dir + "All_ICESat2_Misha.csv"
    list_testcase = []
    df = pd.DataFrame(list_testcase)
    csv_folder = os.listdir(Input_dir)

    for i in csv_folder:
        print("Merging..." + i)
        in_csv = Input_dir + i
        df_item = pd.read_excel(in_csv)
        df = pd.concat([df, df_item])
    df.to_csv(out_csv, index=None)


def join_csv_final():
    list_testcase = []
    df = pd.DataFrame(list_testcase)
    out = r"E:\ICESat2\Results\Result\All_ICESat2_P_Misha.csv"
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    for site in studysite:
        print(site)
        dir1 = r"E:\ICESat2\Results\Basic_ICESAT2_2/" + site + "_P98.csv"
        dir2 = r"E:\ICESat2\Results\Basic_ICESAT2_2/" + site + "_ICESat2.xlsx"
        dir3 = r"E:\ICESat2\Results\Basic_ICESAT2_Misha/" + site + "_ICESat2.xlsx"

        df1 = pd.read_csv(dir1)
        df2 = pd.read_excel(dir2)
        df3 = pd.read_excel(dir3)

        df_merge1 = pd.merge(df2, df1, left_on='TARGET_FID', right_on='shot', how='left')
        df_merge2 = pd.merge(df_merge1, df3, left_on='TARGET_FID', right_on='TARGET_FID', how='left')

        df = pd.concat([df, df_merge2])
    df.to_csv(out, index=None)


def las_DTM_I2():
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    for site in studysite:
        dir = r"E:\ICESat2\Results\ALS/" + site + "/"
        output = r"E:\ICESat2\Results\DTM/"
        lastool_lasnoise = r"D:\lastool\LAStools\bin\lasnoise.exe "
        lastool_ground = r"D:\lastool\LAStools\bin\lasground.exe "
        lastool_lasclassify = r"D:\lastool\LAStools\bin\lasclassify.exe "
        lastool_height = r"D:\lastool\LAStools\bin\lasheight.exe "
        lastool_lascanopy = r"D:\lastool\LAStools\bin\lascanopy.exe "
        lastool_blast2dem = r"D:\lastool\LAStools\bin\blast2dem.exe "

        subprocess.call(lastool_lasnoise + " -i " + dir + "*.las" +
                        " -set_classification 0 -step_xy 1 -step_z 1  -odix _noise -olaz -odir "+
                        output + "noise/" + site)

        subprocess.call(lastool_ground + " -i " + output + "noise/" + site + "/*.laz" +
                        " -odix _ground -olaz -step 10 -all_returns -odir " +
                        output + "ground/" + site)
        subprocess.call(lastool_lasclassify + " -i " + output + "ground/" + site + "/*.laz" +
                        " -olaz -odix _class -odir " +
                        output + "class/" + site)
        subprocess.call(lastool_height + " -i " + output + "class/" + site + "/*.laz" +
                        " -olaz  -replace_z -odix _height -odir " +
                        output + "height/" + site)
        subprocess.call(lastool_lascanopy + " -i " + output + "height/" + site + "/*.laz" +
                        " -step 1 -cover_cutoff 0.5 -dns -p 98 -height_cutoff 0.0 -odix _1m -otif -odir " +
                        output + "canopy/" + site)
        subprocess.call(lastool_blast2dem + " -i " + output + "ground/" + site + "/*.laz" +
                        " -odix _class_2 -otif -keep_class 2 -step 1 -odir " +
                        output + "DTM/" + site)

        #subprocess.call(lastool_lasnoise + " -i *.las -set_classification 0 -step_xy 1 -step_z 1 -odir -odix _noise -olaz")
        #subprocess.call(lastool_ground + " -i *_noise.laz -odix _ground_all_returns_step_10 -olaz -step 10 -all_returns")
        #subprocess.call(lastool_lasclassify+" -i *_ground_all_returns_step_10.laz -odix _class –olaz")
        #subprocess.call(lastool_height+" -i *_ground_all_returns_step_10_class.laz -replace_z -odix _height –olaz")
        #subprocess.call(lastool_lascanopy+" -i * _ground_all_returns_step_10_class_height.laz -step 1 -cover_cutoff 0.5 -dns -p 98 -height_cutoff 0.0 - odix _1m -otif")
        #subprocess.call(lastool_blast2dem+" -i * _ground_all_returns_step_10.laz -odix _class_2 -otif -keep_class 2 -step 1")



def las_DTM_GEDI():
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    for site in studysite:
        dir = r"E:\GEDI\Result\ALS\P98/" + site + "/"
        output = r"E:\GEDI\Result\ALS\DTM/"
        lastool_lasnoise = r"D:\lastool\LAStools\bin\lasnoise.exe "
        lastool_ground = r"D:\lastool\LAStools\bin\lasground.exe "
        lastool_lasclassify = r"D:\lastool\LAStools\bin\lasclassify.exe "
        lastool_height = r"D:\lastool\LAStools\bin\lasheight.exe "
        lastool_lascanopy = r"D:\lastool\LAStools\bin\lascanopy.exe "
        lastool_blast2dem = r"D:\lastool\LAStools\bin\blast2dem.exe "

        subprocess.call(lastool_lasnoise + " -i " + dir + "*.las" +
                        " -set_classification 0 -step_xy 1 -step_z 1  -odix _noise -olaz -odir "+
                        output + "noise/" + site)
        subprocess.call(lastool_ground + " -i " + output + "noise/" + site + "/*.laz" +
                        " -odix _ground -olaz -step 10 -all_returns -odir " +
                        output + "ground/" + site)
        subprocess.call(lastool_lasclassify + " -i " + output + "ground/" + site + "/*.laz" +
                        " -olaz -odix _class -odir " +
                        output + "class/" + site)
        subprocess.call(lastool_height + " -i " + output + "class/" + site + "/*.laz" +
                        " -olaz  -replace_z -odix _height -odir " +
                        output + "height/" + site)
        subprocess.call(lastool_lascanopy + " -i " + output + "height/" + site + "/*.laz" +
                        " -step 1 -cover_cutoff 0.5 -dns -p 98 -height_cutoff 0.0 -odix _1m -otif -odir " +
                        output + "canopy/" + site)
        subprocess.call(lastool_blast2dem + " -i " + output + "ground/" + site + "/*.laz" +
                        " -odix _class_2 -otif -keep_class 2 -step 1 -odir " +
                        output + "DTM/" + site)

        #subprocess.call(lastool_lasnoise + " -i *.las -set_classification 0 -step_xy 1 -step_z 1 -odir -odix _noise -olaz")
        #subprocess.call(lastool_ground + " -i *_noise.laz -odix _ground_all_returns_step_10 -olaz -step 10 -all_returns")
        #subprocess.call(lastool_lasclassify+" -i *_ground_all_returns_step_10.laz -odix _class –olaz")
        #subprocess.call(lastool_height+" -i *_ground_all_returns_step_10_class.laz -replace_z -odix _height –olaz")
        #subprocess.call(lastool_lascanopy+" -i * _ground_all_returns_step_10_class_height.laz -step 1 -cover_cutoff 0.5 -dns -p 98 -height_cutoff 0.0 - odix _1m -otif")
        #subprocess.call(lastool_blast2dem+" -i * _ground_all_returns_step_10.laz -odix _class_2 -otif -keep_class 2 -step 1")


def las_DTM_I2_merge():
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']

    for site in studysite:
        dir = r"E:\ICESat2\Results\DTM/"
        output_csv = dir + "Results/"
        output_shp = r"E:\ICESat2\Results\SHP_Misha\Merge/"
        lastool_lascanopy = r"D:\lastool\LAStools\bin\lascanopy.exe "

        print("Extract P98..." + site)
        subprocess.call(lastool_lascanopy + " -i " + dir + "height/" + site + "/*.laz" +
                        " -merged -drop_class 7 -lop " + output_shp + "ICESat_" + site + "_SJ_P.shp " +
                        " TARGET_FID -cover_cutoff 0.5 -dns -height_cutoff 0.0 -c 0.0 0.5 2 60 -max -p 50 90 98 -o " +
                        output_csv + site + "_P98.csv")

        print("Extract DTM..." + site)
        subprocess.call(lastool_lascanopy + " -i " + dir + "ground/" + site + "/*.laz" +
                        " -merged -keep_class 2 -lop " + output_shp + "ICESat_" + site + "_SJ_P.shp " +
                        " TARGET_FID -avg -min -max -p 50 -o " +
                        output_csv + site + "_DTM.csv")



def las_DTM_GEDI_merge():
    studysite = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']

    for site in studysite:
        dir = r"E:\GEDI\Result\ALS\DTM/"
        output_csv = dir + "Results/"
        output_shp = r"E:\GEDI\Result_Misha\SHP/"
        lastool_lascanopy = r"D:\lastool\LAStools\bin\lascanopy.exe "

        print("Extract P98..." + site)
        subprocess.call(lastool_lascanopy + " -i " + dir + "height/" + site + "/*.laz" +
                        " -merged -drop_class 7 -lop " + output_shp + site + "_p_buffer.shp " +
                        " shot_numbe -cover_cutoff 0.5 -dns -height_cutoff 0.0 -c 0.0 0.5 2 60 -max -p 50 90 98 -o " +
                        output_csv + site + "_P98.csv")

        print("Extract DTM..." + site)
        subprocess.call(lastool_lascanopy + " -i " + dir + "ground/" + site + "/*.laz" +
                        " -merged -keep_class 2 -lop " + output_shp + site + "_p_buffer.shp " +
                        " shot_numbe -avg -min -max -p 50 -o " +
                        output_csv + site + "_DTM.csv")
