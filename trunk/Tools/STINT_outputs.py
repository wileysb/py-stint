#  Output completed intersections to CSV or other product formats
#   part of Spatial/Temporal Intersection Toolset
#
#  (c) Copyright Wiley Bogren 2014
#  Authors:      Wiley Bogren
#  Department: Industrial Ecology
#              Norges Teknisk og Naturvitenskapelige Universitet
#
#  Email: wiley dot bogren at gmail dot com
#
##################################################################
#
#  This MODIS Python file is licensed under the terms of GNU GPL 2.
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.
#
##################################################################
__author__ = 'wiley'

import os,sys,math,cPickle,csv
import rtree,ogr,osr,gdal,gdalconst,h5py
import numpy as np


def Unique_values(src_dsn,field):
    '''
    value_list = uniqueValues(src_dsn, fieldname)
    ### This program returns a list of all unique values within a specified
    ### field. It takes two arguments, 'f' for the filename (including extension)
    ### and 'field' for the desired field. Both arguments passed to the function
    ### should be strings.
    from: Chandler Sterling
    https://github.com/csterling/ogrTools/blob/master/uniqueValues.py
    '''
    #driver = ogr.GetDriverByName("GeoJSON")
    #inFile = driver.Open(f,0)
    inFile = ogr.Open(src_dsn+'.shp',gdalconst.GA_ReadOnly)
    layer  = inFile.GetLayer()
    src_fn = os.path.split(src_dsn)[-1]

    uniqueValues = "select distinct " + field + " from " + src_fn

    result = inFile.ExecuteSQL(uniqueValues)
    resultFeat = result.GetNextFeature()

    uniqueFieldList = []

    while resultFeat:
        field = resultFeat.GetField(0)

        uniqueFieldList.append(field)

        resultFeat = result.GetNextFeature()

    return uniqueFieldList


def Veclc2csv( project ):
    '''Workflow for projects with shapefile aoi (project['lc_type']=='shp')
    Open prj_lcm.shp

    *Collect unique values for field era_id
    Write output file with era timeseries, 1 line/cell, *scale+off

    Count unique values for attribute field mod_id
    Determine number of output tables numfiles,
     such that each file reports max n mod_id (n=50?)

    Write output files for MODIS and landcover:
    for i in range(numfiles):
        for mod_id in that outnum:
            write fileA: forestry/lc attrs + mod_id,era_id: 1 line/feat
            write fileB1, B2, ...: modis timeseries: 1 line/cell, *scale+off
    '''
    lcm_dsn = project['prj_name']+'_lcm'
    lcm_path = os.path.join(project['shp_dir'],lcm_dsn)
    lcm_ds  = ogr.Open(lcm_path+'.shp',gdalconst.GA_ReadOnly)

    era_id_list = Unique_values( lcm_path,"era_id")
    era_ind_list= []
    ind_sql_ = 'SELECT era_x_ind,era_y_ind FROM %s WHERE era_id=%i'
    for era_id in era_id_list:
        ind_sql = ind_sql_ % (lcm_dsn, era_id)
        ind_lyr = lcm_ds.ExecuteSQL(ind_sql)
        ind_feat= ind_lyr.GetNextFeature()

        era_x_ind = ind_feat.GetField(0)
        era_y_ind = ind_feat.GetField(1)

        era_ind_list.append((era_id,era_x_ind,era_y_ind))

    for era_sds in project['era'].keys():
        Era2csv(project, era_sds, era_ind_list)


    mod_id_list = Unique_values( lcm_path,"mod_id")
    # Split mod_ids into groups of 50
    # for each group of 50, open new outfile_number
       # for each modis id within the group of 50:
            # layer.SetAttributeFilter("mod_id = '%s'" % mod_id_val)
            # for feature in layer:
            #     find modis sum
            #     write climate data for that feat


def Raslc2csv( project, ras_lcm):
    # ras_lcm will consist of ras_fn,poly_dsn,dst_p
    isect  = cPickle.load(open(ras_lcm['p'],'r'))
    hmm = 'not really sure how this will work'
    hmm+= 'modis and climate get intersected in vectorland?'
    hmm+= 'Then poly_dsn = mod+era, ras_fn = landcover, dst_p = (mod+era)X(landcover)'
    print hmm
    return hmm


def Era2csv(project, era_sds, era_ind_list):
    '''
    Get matching era_x_ind and era_y_ind
    for each era sds:
        Open hdf stack, get scale and offset
        Open outfile
        for each era id, write climate timeseries
    '''
    era_csv_fn = os.path.join(project['csv_dir'],
                              project['prj_name']+'_'+era_sds+'.csv')
    era_hdf_fn = os.path.join(project['hdf_dir'],
                              project['prj_name']+'_'+era_sds+'.hdf5')

    era_hdf    = h5py.File(era_hdf_fn,'r')

    if len(project['modis_days']) != len(era_hdf['time']):
        print '[ERROR] LENGTH MISMATCH: hdf time dimension and modis_days'
        print era_hdf_fn

    hdr = project['modis_days'][:]
    for col_name in ['era_y_ind','era_x_ind','era_id']:
        hdr.insert(0,col_name)

    era_csv_f = open(era_csv_fn,'wt')
    era_csv    = csv.writer(era_csv_f)
    era_csv.writerow(hdr)

    scale = era_hdf[era_sds].attrs['scale_factor']
    offset= era_hdf[era_sds].attrs['add_offset']
    era_nan = era_hdf[era_sds].attrs['fill_value']

    for era_line in era_ind_list:
        era_id    = era_line[0]
        era_x_ind = era_line[1]
        era_y_ind = era_line[2]
        out = era_hdf[era_sds][:,era_y_ind,era_x_ind]
        out = list(np.where(out==era_nan,np.nan,out)*scale+offset)
        for val in [era_y_ind,era_x_ind,era_id]:
            out.insert(0,val)
        # Should we reduce the precision of the csv output?
        era_csv.writerow(out)

    era_csv_f.close()
    era_hdf.close()