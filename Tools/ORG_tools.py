#  Assorted support tools
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

import os,ogr,osr,gdal,gdalconst
import datetime as dt
import numpy as np

datestr_fmt = '%Y-%m-%d'

def Parse_input(fn):
    inp_params = {}
    inp_params['modis'] = {}
    inp_params['era']   = {}
    inp_params['lc']    = {} # tif paths or shapefile fields
    inp = open(fn,'r')
    lines = inp.readlines()

    int_list = ['start_year','end_year','NODATA',
                'lc_fine','lc_exclude_nans']

    for line in lines:
        if line[0]!='#' and line[:2]!='\n':
            ln_pair = line.split('=')[:2] # take only two splits
            k = ln_pair[0].strip()
            val = ln_pair[1].split('#')[0].strip()
            if k == 'aoi':
                if os.path.isdir(val):
                    inp_params['lc_type'] = 'tif_dir'
                    inp_params['lc_dir']  = val
                elif os.path.splitext(val)[-1] == '.shp':
                    inp_params['lc_type'] = 'shp'
                    inp_params['lc_src'] = val
            if k in int_list:
                inp_params[k] = int(val)
            elif ',' in val:
                if ':' in val:
                    inp_params['modis'][k] = {}
                    mod_set = val.split(',')
                    for mod_pair in mod_set:
                        if ':' in mod_pair:
                            dnum,dset = mod_pair.split(':')
                            inp_params['modis'][k][int(dnum)] = dset
                else:
                    vals  = val.split(',')
                    nc_list = filter(bool,vals) # must preserve order!
                    nc_files = [nc.strip() for nc in nc_list]
                    prj_dir  = inp_params['prj_directory']
                    era_path = inp_params['ERA_dir']
                    era_dir  = os.path.join(prj_dir,era_path)
                    nc_paths = [os.path.join(era_dir,ncf) 
                                for ncf in nc_files]
                    inp_params['era'][k] = nc_paths
            elif k[:3] == 'lc_':
                if os.path.splitext(val)[-1]=='.tif':
                    val_path = os.path.join(inp_params['lc_dir'],val)
                    inp_params['lc'][k[3:]] = val_path
                else:
                    inp_params['lc'][k[3:]] = val
            
            else:
                inp_params[k] = val
        
    if inp_params['lc_type']=='tif_dir':
        # define lc_src:
        src_k  = inp_params['rlc_init'][3:]
        lc_fn = inp_params['lc'][src_k]
        inp_params['lc_src'] = lc_fn
        
        # define spatial reference system:
        lc = gdal.Open(lc_fn, gdalconst.GA_ReadOnly)
        src_projection = lc.GetProjectionRef()
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_projection)
        inp_params['srs'] = srs # ogr/osr object
        del lc
    elif inp_params['lc_type']=='shp':
        lc_fn = inp_params['lc_src']
        lc  = ogr.Open(lc_fn)
        lyr = lc.GetLayer()
        srs = lyr.GetSpatialRef()
        inp_params['srs'] = srs
        del lc
    
    jnk = inp_params.pop('aoi')
    
    
    return inp_params


def Check_input(inp_params):
    err = []
    
    # Collect project directories
    dir_keys = ['ERA_dir','MODIS_dir','prj_directory']
    if inp_params['lc_type']=='tif_dir':
        dir_keys.append('lc_dir')
    
    # Check that directories exist
    for path_key in dir_keys:
        path = inp_params[path_key]
        if not os.path.isdir(path):
            err.append('[ERROR] Bad path to %s: %s' % (path_key,path))
    
    # Check paths to climate and landcover files
    fn_key = 'lc_src'
    fn     = inp_params[fn_key]
    if not os.path.isfile(fn):
        err.append('[ERROR] Bad path to %s: %s' % (fn_key,fn))
    for era_dset in inp_params['era'].keys():
        for nc_file in inp_params['era'][era_dset]:
            if not os.path.isfile(nc_file):
                 err.append('[ERROR] Bad path to ERA.nc: %s' % nc_file)
    
    # Check paths to landcover tifs
    if inp_params['lc_type']=='tif_dir':
        for lc_key in inp_params['lc'].keys():
            lc_fn = inp_params['lc'][lc_key]
            if not os.path.isfile(lc_fn):
                err_pair = (lc_key, lc_fn)
                err.append('[ERROR] Bad path to lc_%s: %s' % err_pair)
    
    # Check landcover projection
    prj_error = 'Landcover file should have projected reference system'
    shp_error = 'Landcover shapefile lacks lc_%s attribute: %s'
    tif_error = 'Landcover raster does not match initial extents: %s'
    tif_prjerr= 'Landcover raster does not match initial projection: %s'
    lc_fn = inp_params['lc_src']
    short_lc_fn = os.path.split(lc_fn)[-1]
    
    # Check shapefile projection and fields
    if inp_params['lc_type']=='shp':
        lc  = ogr.Open(lc_fn)
        lyr = lc.GetLayer()
        srs = lyr.GetSpatialRef()
        if not srs.IsProjected():
            err.append('[ERROR] %s: %s' % (prj_error,short_lc_fn))
        lyr_defn = lyr.GetLayerDefn()
        fcnt = range(lyr_defn.GetFieldCount())
        field_names = [lyr_defn.GetFieldDefn(i).GetName() for i in fcnt]
        for attr_key in inp_params['lc'].keys():
            lc_att = inp_params['lc'][attr_key]
            if not lc_att in field_names:
                err.append(shp_error % (attr_key,lc_att))
    
    # Check that raster has projected reference system
    if inp_params['lc_type']=='tif_dir':
        lc = gdal.Open(lc_fn, gdalconst.GA_ReadOnly)
        src_projection = lc.GetProjectionRef()
        srs = osr.SpatialReference()
        srs.ImportFromWkt(src_projection)
        if not srs.IsProjected():
            err.append('[ERROR] %s: %s' % (prj_error,short_lc_fn))
        
        src_extent = lc.GetGeoTransform()
        
        # Check that all rasters have same extent and projection
        for attr_key in inp_params['lc'].keys():
            fn = inp_params['lc'][attr_key]
            fn_short = os.path.split(fn)[-1]
            lc_att = gdal.Open(fn, gdalconst.GA_ReadOnly)
            lc_projection = lc_att.GetProjectionRef()
            lc_extent = lc_att.GetGeoTransform()
            if lc_projection != src_projection:
                err.append(tif_prjerr % fn_short)
            if lc_extent != src_extent:
                err.append(tif_error % fn_short)
    
    return err


def Get_modis_days(start_year, end_year):
    modis_dates = []
    for yr in range(start_year, end_year+1):
        for day in range(1, 369, 8):
            modis_dates.append('20%02d%03d' % (yr, day))

    if start_year == 0:
        modis_dates = modis_dates[7:]

    return modis_dates


def Yearday2datetime(yearday):
    ydate = dt.datetime.strptime(str(yearday), '%Y%j')
    return ydate


def Yearday2hrnum(basedate, yearday):
    ydate = dt.datetime.strptime(str(yearday), '%Y%j')
    tdelta = ydate - basedate
    hrnum = tdelta.total_seconds()/3600
    return np.int32(hrnum)


def Daynum2date(daynum, basedatestr):
    '''modis_start_date = Daynum2date(daynum, basedate)
    dt.datetime.strptime(str(mday), '%Y%j').date()-
      basedate).total_seconds()/86400.
    '''
    basedate = dt.datetime.strptime(basedatestr, datestr_fmt)
    out = basedate + dt.timedelta(days=daynum)
    out_string = dt.datetime.strftime(out, datestr_fmt)
    return out_string
