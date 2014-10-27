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

import os
import sys
import csv
import math
import ogr
import gdalconst
import h5py
import sqlite3
import numpy as np
from ORG_tools import Countdown
from SPATIAL_tools import Ogr_open


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
    src_fn = os.path.split(src_dsn)[-1]

    uniqueValues = "select distinct " + field + " from " + src_fn

    result = inFile.ExecuteSQL(uniqueValues)
    resultFeat = result.GetNextFeature()

    uniqueFieldList = []

    while resultFeat:
        field = resultFeat.GetField(0)

        uniqueFieldList.append(field)

        resultFeat = result.GetNextFeature()


    uniqueFieldList.sort()
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
    ## Define some paths
    lcm_dsn  = project['prj_name']+'_lcm'
    lcm_path = os.path.join(project['shp_dir'],lcm_dsn)

    mod_dsn = project['prj_name']+'_modis_reprj'
    mod_path= os.path.join(project['shp_dir'],mod_dsn)
    ## Get ERA id and hdf5 indices from lcm shapefile
    lcm_ds  = ogr.Open(lcm_path+'.shp',gdalconst.GA_ReadOnly)

    era_id_list = Unique_values( lcm_path,"era_id")
    era_ind_list= []
    ind_sql_ = 'SELECT era_x_ind,era_y_ind FROM %s WHERE era_id=%i'
    print 'Summarizing %s ERA cells from lcm' % len(era_id_list)
    for era_id in era_id_list:
        ind_sql = ind_sql_ % (lcm_dsn, era_id)
        ind_lyr = lcm_ds.ExecuteSQL(ind_sql)
        ind_feat= ind_lyr.GetNextFeature()

        era_x_ind = ind_feat.GetField(0)
        era_y_ind = ind_feat.GetField(1)

        era_ind_list.append((era_id,era_x_ind,era_y_ind))

    del lcm_ds,ind_lyr

    ## Write each ERA dataset to CSV, one row per cell within aoi
    for era_sds in project['era'].keys():
        Era2csv(project, era_sds, era_ind_list)


    ## Get MODIS id and hdf5 indices from lcm shapefile
    lcm_ds  = ogr.Open(lcm_path+'.shp',gdalconst.GA_ReadOnly)
    mod_ds  = ogr.Open(mod_path+'.shp',gdalconst.GA_ReadOnly)
    mod_id_list = Unique_values( lcm_path,"mod_id")
    mod_ind_list= []
    ind_sql_ = 'SELECT mod_x_ind,mod_y_ind FROM %s WHERE mod_id=%i'
    mod_sql_ = 'SELECT ctr_x,ctr_y FROM %s WHERE id=%i'
    print 'Summarizing %s MODIS cells from lcm' % len(mod_id_list)
    progress_bar = Countdown(len(mod_id_list))
    for i in range(len(mod_id_list)):
        mod_id = mod_id_list[i]
        ind_sql = ind_sql_ % (lcm_dsn, mod_id)
        mod_sql = mod_sql_ % (mod_dsn, mod_id)
        ind_lyr = lcm_ds.ExecuteSQL(ind_sql)
        ind_feat= ind_lyr.GetNextFeature()

        mod_x_ind = ind_feat.GetField(0)
        mod_y_ind = ind_feat.GetField(1)

        mod_lyr   = mod_ds.ExecuteSQL(mod_sql)
        mod_feat  = mod_lyr.GetNextFeature()
        geom      = mod_feat.GetGeometryRef()
        mod_area  = geom.Area()

        mod_ind_list.append((mod_id,mod_x_ind,mod_y_ind,mod_area))
        progress_bar.check(i)

    progress_bar.flush()

    del lcm_ds,ind_lyr

    ## Check for number of landcover features
    lcm_ds  = ogr.Open(lcm_path+'.shp',gdalconst.GA_ReadOnly)
    lcm_lyr = lcm_ds.GetLayer(0)
    lcm_sz  = lcm_lyr.GetFeatureCount()
    del lcm_ds,lcm_lyr

    # Export Landcover/Modis datasets to csv:
    # use lcm_sz to control regionification
    Lcmod_manager( project, lcm_sz, mod_ind_list)


def Lcmod_manager( project, lcm_sz, mod_ind_list ):
    if 'lcm_thresh' in project.keys():
        try:
            lcm_thresh = int(project['lcm_thresh'])
        except:
            'lcm_thresh not an integer (INPUT.txt)'
            lcm_thresh = 20000
    else:
        lcm_thresh = 20000

    if 'mod_per_region' in project.keys():
        mod_per_candidate = project['mod_per_region']
        try:
            mod_per_reg = float(mod_per_candidate)
        except:
            mod_default = 50.
            print 'mod_per_region not an integer (INPUT.txt)'
            print 'assigning default %s modis IDs per region' % mod_default
            mod_per_reg = float(mod_default)
    else:
        mod_default = 50.
        mod_per_reg = float(mod_default)

    if lcm_sz > lcm_thresh:
        # How many times does mod_per_reg go into len(mod_id) ?
        num_reg = int(math.ceil(len(mod_ind_list)/mod_per_reg))
        modr = (lcm_sz,num_reg,int(mod_per_reg))
        print 'Splitting %s lcm features into %s regions, %s MODIS cells each' % modr

        # For each region, export one set of CSV files
        # for the modis ids within that region
        progress_bar = Countdown(num_reg)
        for region in range(num_reg):
            s = int(region * mod_per_reg)
            e = int((region+1) * mod_per_reg)
            if e>=len(mod_ind_list):
                e=-1
            region_mod_ind = mod_ind_list[s:e]
            numrows = Land2csv(project, region_mod_ind, region=region)
            if numrows:
                print 'rows in %s, rows out %s' % (lcm_sz,numrows)
            for mod_type in project['modis'].keys():
                for modis_sds in project['modis'][mod_type].values():
                    Mod2csv(project, modis_sds, region_mod_ind, region=region)
            progress_bar.check(region)
        progress_bar.flush()


    else:
        print 'Writing %s lcm features to 1 CSV per dataset:' % lcm_sz
        print 'lc.csv: '
        numrows = Land2csv(project, mod_ind_list, region=None)
        if numrows:
            toprint = 'rows in %s, rows out %s' % (lcm_sz,numrows)
        for mod_type in project['modis'].keys():
            for modis_sds in project['modis'][mod_type].values():
                sys.stdout.write("%s.csv . . " % modis_sds)
                sys.stdout.flush()
                Mod2csv(project, modis_sds, mod_ind_list)

        sys.stdout.write("Done \n")
        sys.stdout.flush()


def Get_unique_ids(mc_dsn, fid_list ):
    '''WORDS!

    # MODIS: id, area, x/y_ind, x/y_sin, dates/attrs
    # ERA  : id, x/y_ind, x/y_geog, dates/attrs
    # LC   : id, area, lc_id, mod_id, era_id
    #   report x/y_ind in lc.csv! also x/y proj?
    # mod_id, era_id come from mc[fid]
    # fid, px, py (or px/py/area) constitute unique identifying characteristics for lc.db:inside
    # px,py alone constitute unique identifying characteristics for lc.db:border

    :param mc_dsn:
    :param fid_list:
    :return:
    '''

    mc_ds, mc_lyr = Ogr_open(mc_dsn)

    mod_id_list = []
    era_id_list = []

    mod_ind_list = []
    era_ind_list = []

    for fid in fid_list:
        src_feat = mc_lyr.GetFeature(fid)
        mod_id = src_feat.GetField('mod_id')
        era_id = src_feat.GetField('era_id')
        if mod_id not in mod_id_list:
            mod_id_list.append(mod_id)
            mod_area  = src_feat.GetField('mod_area')
            mod_x_ind = src_feat.GetField('mod_x_ind')
            mod_y_ind = src_feat.GetField('mod_y_ind')
            mod_ind_list.append((mod_id, mod_x_ind, mod_y_ind, mod_area))
        if era_id not in era_id_list:
            era_id_list.append(era_id)
            era_x_ind = src_feat.GetField('era_x_ind')
            era_y_ind = src_feat.GetField('era_y_ind')
            era_ind_list.append((era_id, era_x_ind, era_y_ind))

    # sort both lists by mod_id and era_id:
    era_ind_list = sorted(era_ind_list, key=lambda feat: feat[0])
    mod_ind_list = sorted(mod_ind_list, key=lambda feat: feat[0])

    del era_id_list, mod_id_list, mc_ds, mc_lyr, src_feat
    return era_ind_list, mod_ind_list


def Get_mc_sz(c):
    sz_qry = '''SELECT COUNT(*) FROM inside
                UNION ALL
                SELECT COUNT(*) FROM border'''

    c.execute(sz_qry)
    sz_inside, sz_border = c.fetchall()
    mc_sz = sz_inside[0] + sz_border[0]

    return mc_sz


def Get_mcfid_list(c):
    fid_qry = '''SELECT DISTINCT fid FROM
                (SELECT fid FROM inside
                UNION ALL
                SELECT fid FROM border)'''
    c.execute(fid_qry)

    fid_list = c.fetchall()

    return [fid[0] for fid in fid_list]


def Raslc2csv( project ):
    mc_name = project['prj_name']+'_mc'
    mc_dsn  = os.path.join(project['shp_dir'], mc_name)

    isect_fn = os.path.join(project['prj_directory'],'lcmc.db')

    conn = sqlite3.connect(isect_fn)
    c    = conn.cursor()

    mc_sz = Get_mc_sz(c)

    fid_list = Get_mcfid_list(c)

    conn.close()
    del c

    # IF csv output is too slow, implement indexed sqlite conversion
    # for Get_unique_ids and Rlc2csv
    era_ind_list, mod_ind_list =Get_unique_ids(mc_dsn, fid_list)

    ## Write each ERA dataset to CSV, one row per cell within aoi
    for era_sds in project['era'].keys():
        Era2csv(project, era_sds, era_ind_list)

    ## Export MOD and LC to csv
    # mc_sz used to determine regionisation
    Lcmod_manager( project, mc_sz, mod_ind_list )


def Land2csv(project, mod_ind_list, region=None):
    if project['lc_type']=='shp':
        Vlc2csv(project, mod_ind_list, region)
    elif project['lc_type']=='tif_dir':
        num_rows = Rlc2csv(project, mod_ind_list, region)
        return num_rows


def Rlc2csv(project, mod_ind_list, region=None):

    ### FIRST BIG SECTION: IDENTIFY FILE PATHS AND OPEN FILES
    # Open lc.hdf5 as source for landcover parameters and cell size
    hdf_name = project['prj_name']+'_lc.hdf5'
    lc_path = os.path.join( project['hdf_dir'], hdf_name )
    lc_stack= h5py.File(lc_path,'r')

    # Open mc.shape to link modis_id to intersect rows lcmc.db rows
    # IF csv output is too slow, implement indexed sqlite conversion
    # for mod_id lookups
    mc_dsn  = project['prj_name']+'_mc'
    mc_path = os.path.join(project['shp_dir'],mc_dsn)
    mc_ds  = ogr.Open(mc_path+'.shp',gdalconst.GA_ReadOnly)

    # Couldn't era_id and mod_id be incorporated in lcmc.db?
    # then this pain in the ass would go MUCH quicker, without shapefile queries...

    # Open lcmc.db as climate-modis X landcover intersect
    isect_fn = os.path.join(project['prj_directory'],'lcmc.db')
    conn = sqlite3.connect(isect_fn)
    c    = conn.cursor()

    # Define output file, depending on regionalization
    if region==None:
        lc_csv_fn = os.path.join(project['csv_dir'],
                                 project['prj_name']+'_lc.csv')
    else:
        out_dir =  os.path.join(project['csv_dir'],'lc')
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        lc_csv_fn = os.path.join(out_dir,
                                 project['prj_name']+'_lc_'+str(region)+'.csv')

    lc_csv_f = open(lc_csv_fn,'wt')
    lc_csv   = csv.writer(lc_csv_f)


    ### DEFINE SOME VARIABLES
    lc_nan =lc_stack['lc'].attrs['fill_value']
    scale_factor = 'None?'
    add_offset = 'None?'
    cell_size = np.abs(lc_stack['lc'].attrs['dx'] * lc_stack['lc'].attrs['dx'])
    mc_sql_ = 'SELECT id, era_id FROM '+mc_dsn+' WHERE mod_id=%i'
    hdr = ['id','area','modis_id','era_id']
    for attrib in lc_stack['fields']:
        hdr.append(attrib)

    lc_csv.writerow(hdr)
    ### START THE LOOP!
    ### IF this is taking too long per cycle, implement sqlite in place of
    ### ogr modis_id -> attributes linking
    row_id=0
    progress_bar = Countdown(len(mod_ind_list), update_interval=.01)
    for i in range(len(mod_ind_list)):
        mod_id = mod_ind_list[i][0]
        mc_sql = mc_sql_ % mod_id
        mc_lyr = mc_ds.ExecuteSQL(mc_sql)
        mc_feat= mc_lyr.GetNextFeature()
        while mc_feat:
            fid = mc_feat.GetField('id')
            era_id   = mc_feat.GetField('era_id')
            # Get lcmc features within modis/era cells
            c.execute('SELECT px,py FROM inside WHERE fid=%s' % fid)
            area = cell_size
            lcmc_feat = c.fetchone()
            mc_feat= mc_lyr.GetNextFeature()
            while lcmc_feat:
                row_id += 1
                px = lcmc_feat[0]
                py = lcmc_feat[1]
                # out = [row_id,area,mod_id,era_id]
                out = lc_stack['lc'][:,py,px]
                out = list(np.where(out==lc_nan,np.nan,out)) #*scale+offset)
                for val in [era_id, mod_id, area, row_id]:
                    out.insert(0,val)

                lc_csv.writerow(out)
                lcmc_feat = c.fetchone()
            c.execute('SELECT px,py,area FROM border WHERE fid=%s' % fid)
            lcmc_feat = c.fetchone()
            while lcmc_feat:
                row_id += 1
                px = lcmc_feat[0]
                py = lcmc_feat[1]
                area = lcmc_feat[2]
                # out = [row_id,area,mod_id,era_id]
                out = lc_stack['lc'][:,py,px]
                out = list(np.where(out==lc_nan,np.nan,out)) #*scale+offset)
                for val in [era_id, mod_id, area, row_id]:
                    out.insert(0,val)

                lc_csv.writerow(out)
                lc_csv.writerow(out)
                lcmc_feat = c.fetchone()
        progress_bar.check(i)

    conn.close()
    del c
    progress_bar.flush()
    lc_csv_f.close()

    return row_id


def Vlc2csv(project, mod_ind_list, region=None):


    lcm_dsn  = project['prj_name']+'_lcm'
    lcm_path = os.path.join(project['shp_dir'],lcm_dsn)
    lcm_ds  = ogr.Open(lcm_path+'.shp',gdalconst.GA_ReadOnly)
    if region==None:
        lc_csv_fn = os.path.join(project['csv_dir'],
                                 project['prj_name']+'_lc.csv')
    else:
        out_dir =  os.path.join(project['csv_dir'],'lc')
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        lc_csv_fn = os.path.join(out_dir,
                                 project['prj_name']+'_lc_'+str(region)+'.csv')

# Do we want cell-center coordinates for modis? era? lc?
# Which projections?
    hdr = ['id','area','lc_id','modis_id','era_id']
    lc_attrs = ''
    for attrib in project['lc'].keys():
        hdr.append(attrib)
        lc_attrs += ',lc_'+project['lc'][attrib]

    lc_csv_f = open(lc_csv_fn,'wt')
    lc_csv   = csv.writer(lc_csv_f)
    lc_csv.writerow(hdr)

    lc_sql_ = 'SELECT id,area,lc_id,mod_id,era_id'+lc_attrs+' FROM '+lcm_dsn+\
              ' WHERE mod_id=%i'

    if region==None:
        progress_bar = Countdown(len(mod_ind_list))
    for i in range(len(mod_ind_list)):
        mod_id = mod_ind_list[i][0]
        lc_sql = lc_sql_ % mod_id
        lc_lyr = lcm_ds.ExecuteSQL(lc_sql)
        lc_feat= lc_lyr.GetNextFeature()
        while lc_feat:
            lcm_id   = lc_feat.GetField('id')
            lcm_area = lc_feat.GetField('area')
            lc_id    = lc_feat.GetField('lc_id')
            mod_id   = lc_feat.GetField('mod_id')
            era_id   = lc_feat.GetField('era_id')
            out = [lcm_id,lcm_area,lc_id,mod_id,era_id]
            for attrib in lc_attrs.split(','):
                if len(attrib)>0:
                    out.append(lc_feat.GetField(attrib))
            lc_csv.writerow(out)
            lc_feat = lc_lyr.GetNextFeature()
        if region==None:
            progress_bar.check(i)

    del lcm_ds
    if region==None:
        progress_bar.flush()
    lc_csv_f.close()


def Mod2csv(project,mod_sds,mod_ind_list, region = None):
    '''
    if not region:
        lc_csv_fn = os.path.join(project['csv_dir'],
                                 project['prj_name']+'_lc.csv')
    else:
        out_dir =  os.path.join(project['csv_dir'],'lc')
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        lc_csv_fn = os.path.join(out_dir,
                                 project['prj_name']+'_lc_'+str(region)+'.csv')
    '''
    if region==None:
        mod_csv_fn = os.path.join(project['csv_dir'],
                     project['prj_name']+'_'+mod_sds+'.csv')
    else:
        out_dir =  os.path.join(project['csv_dir'],mod_sds)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        mod_csv_fn = os.path.join(out_dir,
                     project['prj_name']+'_'+mod_sds+'_'+str(region)+'.csv')

    mod_hdf_fn = os.path.join(project['hdf_dir'],
                              project['prj_name']+'_'+mod_sds+'.hdf5')

    mod_hdf    = h5py.File(mod_hdf_fn,'r')
    mod_x  = mod_hdf['x'][:]
    mod_y  = mod_hdf['y'][:]

    if len(project['modis_days']) != len(mod_hdf['time']):
        print '[ERROR] LENGTH MISMATCH: hdf time dimension and modis_days'
        print mod_hdf_fn

    hdr = project['modis_days'][:]
    for col_name in ['mod_y_sin','mod_x_sin','mod_y_ind','mod_x_ind','mod_area','mod_id']:
        hdr.insert(0,col_name)

    mod_csv_f = open(mod_csv_fn,'wt')
    mod_csv    = csv.writer(mod_csv_f)
    mod_csv.writerow(hdr)

    scale = mod_hdf[mod_sds].attrs['scale_factor']
    offset= mod_hdf[mod_sds].attrs['add_offset']
    mod_nan = mod_hdf[mod_sds].attrs['fill_value']

    for mod_line in mod_ind_list:
        mod_id    = mod_line[0]
        mod_x_ind = mod_line[1]
        mod_y_ind = mod_line[2]
        mod_area  = mod_line[3]
        mod_x_sin = mod_x[mod_x_ind]
        mod_y_sin = mod_y[mod_y_ind]
        out = mod_hdf[mod_sds][:,mod_y_ind,mod_x_ind]
        out = list(np.where(out==mod_nan,np.nan,out)*scale+offset)
        for val in [mod_y_sin,mod_x_sin,mod_y_ind,mod_x_ind,mod_area,mod_id]:
            out.insert(0,val)
        # Should we reduce the precision of the csv output?
        mod_csv.writerow(out)

    mod_csv_f.close()
    mod_hdf.close()


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
    # Get axes coordinates, hdf grid
    era_x  = era_hdf['x'][:]
    era_y  = era_hdf['y'][:]

    if len(project['modis_days']) != len(era_hdf['time']):
        print '[ERROR] LENGTH MISMATCH: hdf time dimension and modis_days'
        print era_hdf_fn

    hdr = project['modis_days'][:]
    for col_name in ['era_y_wgs84','era_x_wgs84','era_y_ind','era_x_ind','era_id']:
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
        era_x_wgs84 = era_x[era_x_ind]
        era_y_wgs84 = era_y[era_y_ind]
        out = era_hdf[era_sds][:,era_y_ind,era_x_ind]
        out = list(np.where(out==era_nan,np.nan,out)*scale+offset)
        for val in [era_y_wgs84,era_x_wgs84,era_y_ind,era_x_ind,era_id]:
            out.insert(0,val)
        era_csv.writerow(out)

    era_csv_f.close()
    era_hdf.close()