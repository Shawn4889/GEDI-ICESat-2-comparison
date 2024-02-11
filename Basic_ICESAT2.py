#Candidacy exam script (April 2022)
#Author: Xiaoxuan Li
import pandas as pd
import numpy as np
import os
import h5py
import shutil
import arcpy
from arcpy.sa import *
arcpy.env.overwriteOutput = True


def filter_h5():
    dir = r"E:\ICESat2\ATL08\Welverdiendt/"
    folders = os.listdir(dir)
    for folder in folders:
        files = os.listdir(dir + folder)
        for file in files:
            if "h5" in file:
                print("moving..." + file)
                shutil.move(dir + folder + "/" + file, dir + "/" + file)


def seg_100_20_atl08():
    ICESAT2_loc = r"E:\ICESat2\ATL08\Duku/"
    Result_dir = r"E:\ICESat2\Results\Basic_ICESAT2/"
    output = r"E:\ICESat2\Results\Basic_ICESAT2\Duku.csv"
    beam = ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]
    folders = os.listdir(ICESAT2_loc)
    # Welverdiendt
    # bounding_box = '31.1573,-24.6478,31.3902,-24.544'
    # Duku
    # bounding_box = '32.1716,-28.4874,32.356,-28.322'

    lat_max = -28.322
    lon_min = 32.1716
    lat_min = -28.4874
    lon_max = 32.356

    ini = []
    df_ini = pd.DataFrame(ini)
    for folder in folders:
        HDFs = os.listdir(ICESAT2_loc+folder)
        for file in HDFs:
            if ".h5" in file:
                print(file)
                # read HDF file
                ICESAT2_ATL08 = ICESAT2_loc + folder + "/" + file
                f_ATL08 = h5py.File(ICESAT2_ATL08, 'r')
                # beam list
                for b in beam:
                    print(b)
                    try:
                        #lat and long: 100m segemnts, 20m segments
                        df_lat_100 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/latitude'][:]), columns=['ICESAT2_lat'])
                        df_long_100 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/longitude'][:]), columns=['ICESAT2_long'])
                        df_lat_20 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/latitude_20m'][:]),
                                                 columns=['lat20_1','lat20_2','lat20_3','lat20_4','lat20_5'])
                        df_long_20 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/longitude_20m'][:]),
                                                  columns=['long20_1','long20_2','long20_3','long20_4','long20_5'])
                        df_Coor = pd.concat([df_lat_100, df_long_100], axis=1)
                        df_interaction = df_Coor[(df_Coor["ICESAT2_lat"] > lat_min) &
                                                 (df_Coor["ICESAT2_lat"] < lat_max) &
                                                 (df_Coor["ICESAT2_long"] > lon_min) &
                                                 (df_Coor["ICESAT2_long"] < lon_max)]
                        index_result = np.asarray(df_interaction.index.values)
                        #ID: segment_id_beg
                        df_seg = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/segment_id_beg'][:]), columns=['ICESAT2_ID'])
                        #quality flag 1: subset_can_flag, subset_te_flag
                        df_pnc = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/subset_can_flag'][:]),
                                                 columns=['pnc_1','pnc_2','pnc_3','pnc_4','pnc_5'])
                        df_pnt = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/subset_te_flag'][:]),
                                                  columns=['pnt_1','pnt_2','pnt_3','pnt_4','pnt_5'])
                        #quality flag 2: photon_rate_can, photon_rate_te
                        df_prc = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/photon_rate_can'][:]), columns=['prc'])
                        df_prt = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/photon_rate_te'][:]), columns=['prt'])
                        #quality flag 3: segment_cover
                        df_cover = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/segment_cover'][:]),columns=['ICESAT2_cover'])
                        #quality flag 4: terrain_slope
                        df_slope = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/terrain_slope'][:]),columns=['ICESAT2_slope'])
                        #quality flag 5: night_flag
                        df_night = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/night_flag'][:]), columns=['night_flag'])
                        #quality_flag 6: canopy_rh_conf
                        df_conf = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/canopy_rh_conf'][:]),columns=['ICESAT2_conf'])
                        #canopy 100m, 20m: h_canopy, h_canopy_20m
                        df_rh98 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/h_canopy'][:]),columns=['ICESAT2_RH98'])
                        df_rh98_20 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/h_canopy_20m'][:]),
                                                 columns=['RH98_1','RH98_2','RH98_3','RH98_4','RH98_5'])
                        #terrain 100m, 20m: h_te_best_fit, h_te_mean, h_te_best_fit_20m
                        df_ground = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/h_te_best_fit'][:]),columns=['ICESAT2_ground'])
                        df_mground = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/h_te_mean'][:]),columns=['ICESAT2_mground'])
                        df_ground_20 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/h_te_best_fit_20m'][:]),
                                                 columns=['ground_1','ground_2','ground_3','ground_4','ground_5'])
                        #merge all dfs
                        df = pd.concat([df_seg, df_Coor, df_lat_20, df_long_20,
                                        df_pnc, df_pnt,df_prc,df_prt, df_cover, df_slope, df_night,
                                        df_conf, df_rh98, df_rh98_20, df_ground, df_mground, df_ground_20], axis=1)
                        df = df.loc[index_result, :]
                        df['beam'] = b
                        df['orient'] = np.array(f_ATL08['orbit_info/sc_orient'][:])[0]
                        df['orbit'] = np.array(f_ATL08['orbit_info/orbit_number'][:])[0]
                        df['filename'] = file
                        #append to loop dataframe
                        df_ini = pd.concat([df_ini, df], sort=False)
                        print(df_ini)
                    except:
                        print("No beam ATL08 info, skip...")

    #define power and weak

    df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt1l'), 'Type'] = 'Strong'

    df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt1l'), 'Type'] = 'Weak'

    df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt2l'), 'Type'] = 'Strong'

    df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt2l'), 'Type'] = 'Weak'

    df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt3l'), 'Type'] = 'Strong'

    df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt3l'), 'Type'] = 'Weak'

    df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt1r'), 'Type'] = 'Weak'

    df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt1r'), 'Type'] = 'Strong'

    df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt2r'), 'Type'] = 'Weak'

    df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt2r'), 'Type'] = 'Strong'

    df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt3r'), 'Type'] = 'Weak'

    df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt3r'), 'Type'] = 'Strong'

    df_ini.to_csv(output, index=None)




#use 20m data instead
def quickstart_atl08():
    sites = ['Addo', 'Agincourt', 'Agulhas', 'DNyala', 'Duku', 'Welverdiendt',
            'Venetia', 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    ICESAT2_loc = r"E:\ICESat2\Archive\ATL08/"
    beam = ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]
    CSV_dir = r"E:\ICESat2/SA.csv"
    df_geo = pd.read_csv(CSV_dir, sep=",")
    ini = []
    df_ini = pd.DataFrame(ini)
    for site in sites:
        print(site)
        lat_max = df_geo.loc[df_geo['Location'] == site, 'ul_lat'].iloc[0]
        lon_min = df_geo.loc[df_geo['Location'] == site, 'ul_lon'].iloc[0]
        lat_min = df_geo.loc[df_geo['Location'] == site, 'lr_lat'].iloc[0]
        lon_max = df_geo.loc[df_geo['Location'] == site, 'lr_lon'].iloc[0]
        track_sa = df_geo.loc[df_geo['Location'] == site, 'Track_combine'].iloc[0]
        track_list = track_sa.split(' ')
        for track in track_list:
            if int(track) < 1000:
                track = "0" + str(track)
            if int(track) < 100:
                track = "0" + str(track)
            ICESAT2_loc_t = ICESAT2_loc + track
            print(ICESAT2_loc_t)
            folder2 = os.listdir(ICESAT2_loc_t)
            for folder1 in folder2:
                if "zip" not in folder1:
                    folders = os.listdir(ICESAT2_loc_t + "/" + folder1)
                    for folder in folders:
                        files = os.listdir(ICESAT2_loc_t + "/" + folder1 + "/" + folder)
                        for file in files:
                            if "h5" in file:
                                ICESAT2_ATL08 = ICESAT2_loc_t + "/" + folder1 + "/" + folder + "/" + file
                                f_ATL08 = h5py.File(ICESAT2_ATL08, 'r')
                                # beam list
                                for b in beam:
                                    try:
                                        # lat and long: 100m segemnts, 20m segments
                                        df_lat_100 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/latitude'][:]),
                                                                  columns=['ICESAT2_lat'])
                                        df_long_100 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/longitude'][:]),
                                                                   columns=['ICESAT2_long'])
                                        df_Coor = pd.concat([df_lat_100, df_long_100], axis=1)
                                        df_interaction = df_Coor[(df_Coor["ICESAT2_lat"] > lat_min) &
                                                                 (df_Coor["ICESAT2_lat"] < lat_max) &
                                                                 (df_Coor["ICESAT2_long"] > lon_min) &
                                                                 (df_Coor["ICESAT2_long"] < lon_max)]
                                        index_result = np.asarray(df_interaction.index.values)
                                        # ID: segment_id_beg
                                        df_seg = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/segment_id_beg'][:]),
                                                              columns=['ICESAT2_ID'])
                                        # quality flag 2: photon_rate_can, photon_rate_te
                                        df_prc = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/photon_rate_can'][:]),
                                                              columns=['prc'])
                                        df_prt = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/photon_rate_te'][:]),
                                                              columns=['prt'])
                                        # quality flag 3: segment_cover
                                        df_cover = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/segment_cover'][:]),
                                                                columns=['ICESAT2_cover'])
                                        # quality flag 4: terrain_slope
                                        df_slope = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/terrain_slope'][:]),
                                                                columns=['ICESAT2_slope'])
                                        # quality flag 5: night_flag
                                        df_night = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/night_flag'][:]),
                                                                columns=['night_flag'])
                                        # quality_flag 6: canopy_rh_conf
                                        df_conf = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/canopy_rh_conf'][:]),
                                                               columns=['ICESAT2_conf'])
                                        # canopy 100m, 20m: h_canopy, h_canopy_20m
                                        df_rh98 = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/canopy/h_canopy'][:]),
                                                               columns=['ICESAT2_RH98'])
                                        # terrain 100m, 20m: h_te_best_fit, h_te_mean, h_te_best_fit_20m
                                        df_ground = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/h_te_best_fit'][:]),
                                                                 columns=['ICESAT2_ground'])
                                        df_mground = pd.DataFrame(np.array(f_ATL08[b + '/land_segments/terrain/h_te_mean'][:]),
                                                                  columns=['ICESAT2_mground'])
                                        # merge all dfs
                                        df = pd.concat([df_seg, df_Coor,df_prc, df_prt, df_cover, df_slope,
                                                        df_night, df_conf, df_rh98, df_ground, df_mground], axis=1)
                                        df = df.loc[index_result, :]
                                        df['beam'] = b
                                        df['orient'] = np.array(f_ATL08['orbit_info/sc_orient'][:])[0]
                                        df['filename'] = str(file.split(".")[0])
                                        track_num = str(file.split(".")[0].split("_")[2][0:4])
                                        df['track'] = track_num
                                        df['date'] = file.split(".")[0].split("_")[1][0:8]
                                        '''
                                        #Wel 0517 right dir
                                        if track_num == "0517":
                                            df['azimuth'] = 354.2
                                        else:
                                            df['azimuth'] = 5.8
                                        '''
                                        # append to loop dataframe
                                        df_ini = pd.concat([df_ini, df], sort=False)
                                    except:
                                        print("No beam ATL08 info, skip..."+ICESAT2_ATL08)
                                if len(df) != 0:
                                    shutil.copy2(ICESAT2_ATL08, "E:\ICESat2\ATL08/" + site + "/" + file)

        # define power and weak
        df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt1l'), 'Type'] = 'Strong'
        df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt1l'), 'Type'] = 'Weak'
        df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt2l'), 'Type'] = 'Strong'
        df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt2l'), 'Type'] = 'Weak'
        df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt3l'), 'Type'] = 'Strong'
        df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt3l'), 'Type'] = 'Weak'
        df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt1r'), 'Type'] = 'Weak'
        df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt1r'), 'Type'] = 'Strong'
        df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt2r'), 'Type'] = 'Weak'
        df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt2r'), 'Type'] = 'Strong'
        df_ini.loc[(df_ini['orient'] == 0) & (df_ini['beam'] == 'gt3r'), 'Type'] = 'Weak'
        df_ini.loc[(df_ini['orient'] == 1) & (df_ini['beam'] == 'gt3r'), 'Type'] = 'Strong'
        output = r"E:\ICESat2\Results\Basic_ICESAT2/" + site + ".csv"
        df_ini.to_csv(output, index=None)





# calculate P98 using lastool first, then run the following codes

#this method is similar to pysl4land, use pysl4land instead
def xy_to_ICESAT(site):
    if site == "Addo" in site or site == "Venetia" or site == "DNyala":
        Proj = r'E:\GEDI\Boundingbox/UTM35S.prj'
    elif site == "Agulhas":
        Proj = r'E:\GEDI\Boundingbox/UTM34S.prj'
    else:
        Proj = r'E:\GEDI\Boundingbox/UTM36S.prj'
    table = r"E:\ICESat2\Results\Basic_ICESAT2/" + site + ".csv"
    out = r"E:\ICESat2\Results\SHP/" + site + "_ICESAT2_point.shp"
    print("XY table to shp...")
    arcpy.management.XYTableToPoint(table, out, "ICESAT2_long", "ICESAT2_lat", '', "")
    print("Shp to ellipse...")
    out_featureclass = r"E:\ICESat2\Results\SHP/" + site + "_ICESAT2_ellipse.shp"
    arcpy.AddField_management(out, "Elong", "Short")
    arcpy.AddField_management(out, "Eshort", "Short")
    arcpy.management.CalculateField(out, "Elong", 100, "PYTHON")
    arcpy.management.CalculateField(out, "Eshort", 13, "PYTHON")
    arcpy.management.TableToEllipse(out, out_featureclass, "ICESAT2_lo","ICESAT2_la", "Elong", "Eshort","METERS",
                                    "azimuth", "DEGREES", "ICESAT2_ID", "")
    print("Ellipse to ICESAT2 rectangle...")
    out_featureclass = r"E:\ICESat2\Results\SHP/" + site + "_ICESAT2_ellipse.shp"
    project_featureclass = r"E:\ICESat2\Results\SHP/" + site + "_ICESAT2_ellipse_p.shp"
    arcpy.management.Project(out_featureclass, project_featureclass, Proj)
    out_feature_class_2 = r"E:\ICESat2\Results\SHP/" + site + "_ICESAT2_footprint.shp"
    arcpy.management.GenerateRectanglesAlongLines(project_featureclass, out_feature_class_2, 100, 13, "")


def join_output(site):
    Point = r"E:\ICESat2\Results\SHP/" + site + "_ICESAT2_point.shp"
    Footprint = r"E:\ICESat2\Results\SHP/" + site + "_ICESAT2_footprint.shp"
    arcpy.MakeFeatureLayer_management(Footprint, 'footprint')
    arcpy.SelectLayerByAttribute_management('footprint', '','"Angle" < 100')
    out_feature_class = r"E:\ICESat2\Results\SHP/" + site + "_ICESAT2_sj.shp"
    print("Spatial join...")
    arcpy.SpatialJoin_analysis('footprint', Point, out_feature_class,"","KEEP_COMMON","","CLOSEST",5)
    print("Zonal stats: sj footprint + CHM...")
    arcpy.MakeFeatureLayer_management(out_feature_class, "temp_layer")
    CHM = "E:\ICESat2\Results\ALS\CHM/" + site + "/Mosaic.tif"
    table_CHM = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/table_chm_" + site
    ZonalStatisticsAsTable("temp_layer", "FID", CHM, table_CHM,"","Maximum")
    arcpy.AddJoin_management("temp_layer", "FID", table_CHM, "FID", "KEEP_COMMON")
    print("Zonal stats: sj footprint + DTM...")
    DTM = "E:\ICESat2\Results\ALS\DTM/" + site + "/Mosaic.tif"
    table_DTM = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/table_DTM_" + site
    ZonalStatisticsAsTable("temp_layer", "FID", DTM, table_DTM,"","MEAN")
    arcpy.AddJoin_management("temp_layer", "FID", table_DTM, "VALUE", "KEEP_COMMON")
    xlsx_DTM = "E:\ICESat2\Results\CSV/" + site + "_ICESAT2.xlsx"
    arcpy.TableToExcel_conversion("temp_layer", xlsx_DTM)


def join_output_GEDI(site):
    GEDI_footprint = r"E:\ICESat2\Results\SHP/" + site + "_GEDI_footprint.shp"
    print("Zonal stats: sj footprint + CHM...")
    arcpy.MakeFeatureLayer_management(GEDI_footprint, "temp_layer")
    CHM = "E:\ICESat2\Results\ALS\CHM/" + site + "/Mosaic.tif"
    table_chm = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/table_chm_G_" + site
    ZonalStatisticsAsTable("temp_layer", "shot_numbe", CHM, table_chm,"","Maximum")
    arcpy.AddJoin_management("temp_layer", "shot_numbe", table_chm, "shot_numbe", "KEEP_COMMON")
    print("Zonal stats: sj footprint + DTM...")
    DTM = "E:\ICESat2\Results\ALS\DTM/" + site + "/Mosaic.tif"
    table_DTM = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/table_DTM_G_" + site
    ZonalStatisticsAsTable("temp_layer", "shot_numbe", DTM, table_DTM,"","MEAN")
    arcpy.AddJoin_management("temp_layer", "shot_numbe", table_DTM, "VALUE", "KEEP_COMMON")
    xlsx_DTM = "E:\ICESat2\Results\CSV/" + site + "_GEDI.xlsx"
    arcpy.TableToExcel_conversion("temp_layer", xlsx_DTM)


