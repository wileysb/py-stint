#  NetCDF to HDF5 conversion for ERA climate archive 
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

import os,sys,osr,h5py,multiprocessing
import numpy as np
import datetime as dt
#from Scientific.IO.NetCDF import NetCDFFile
from scipy.io.netcdf import netcdf_file as NetCDFFile
#Can this simply drop in like this?
from MODIS_aoi import Mk_bbox
from ORG_tools import Yearday2hrnum

### DEFINE GEOGRAPHIC PROJECTION, WGS84
# sin_prj_string directly from MODIS hdf files:
wgs84 = osr.SpatialReference()
wgs84.ImportFromEPSG(4326)
wgs84_string = wgs84.ExportToWkt()

# ERA time_var standard
basedate= dt.datetime(1900,01,01,00)

### Builder functions
def Start_era_hdf( project, hdfp):
    print 'Creating',hdfp['sds'],hdfp['h5f']
    src_nc = Build_src_nc( project, hdfp)
    missing_val  = src_nc[0].variables[hdfp['sds']].missing_value
    nclon   = src_nc[0].variables['longitude'][:]
    nclat   = src_nc[0].variables['latitude'][:]

    dlon = nclon[1]-nclon[0]
    dlat = nclat[1]-nclat[0]

    nclon = nclon - 0.5*dlon
    nclat = nclat + 0.5*dlat

    x_s,x_e,y_s,y_e = Get_spatial_indexes(nclon, nclat, project['aoi'])
    hdlon = nclon[x_s:x_e]
    hdlat = nclat[y_s:y_e]

    with h5py.File(hdfp['h5f'],"w") as hdf:
        dshape=(len(hdfp['hdtime']),len(hdlat),len(hdlon))
        arr_out = hdf.create_dataset(hdfp['sds'],dshape,
                      dtype='int16',chunks=True,compression='lzf')
        arr_out.attrs['scale_factor'] = src_nc[0].variables[hdfp['sds']].scale_factor
        arr_out.attrs['add_offset']   = src_nc[0].variables[hdfp['sds']].add_offset
        arr_out.attrs['units']        = src_nc[0].variables[hdfp['sds']].units
        arr_out.attrs['projection']   = wgs84_string
        arr_out.attrs['fill_value']   = missing_val
        arr_out.attrs['dx']           = dlon
        arr_out.attrs['dy']           = dlat
        x_out = hdf.create_dataset("x",data=hdlon)
        y_out = hdf.create_dataset("y",data=hdlat)
        t_out = hdf.create_dataset("time", (len(hdfp['hdtime']),),
                                    dtype='int32')
        t_out.attrs['units'] = 'hours since 1900-01-01 00:00:0.0'
        yr_out= hdf.create_dataset("year",data=hdfp['years'])
        day_out=hdf.create_dataset("yday",data=hdfp['ydays'])
    del src_nc


def Build_src_nc( project, hdfp ):
    '''Can I make this a class, so I can include a close/destroy?'''
    src_nc = []
    for nc_fn in project['era'][hdfp['sds']]:
        src_nc.append(NetCDFFile(nc_fn,'r'))
    return src_nc


def Gen_appendexes( project, hdfp ):
    '''Find times already done
    split times-to-do into groups of $appendnum
    return groups'''
    with h5py.File(hdfp['h5f'],"r") as hdf:
        progress = hdf['time'][:]
        diff_p = np.diff(progress)
        if (diff_p<=0).any(): # start or mid
            if (diff_p<0).any(): # mid
                i_st = np.where(np.diff(progress)<0)[0][0]
            elif (diff_p==0).all():
                i_st = 0
            appind = []
            while (i_st+hdfp['appendnum'])<len(project['modis_days']):
                appind.append((i_st,i_st+hdfp['appendnum']))
                i_st+=hdfp['appendnum']
            appind.append((i_st,len(project['modis_days'])))
        else:
            print hdfp['sds'],'already done??'
            appind=False
    return appind


def ID_intervals(start_time,end_time,src_nc):
    t_ind = []
    for nc in src_nc:
        t_interval = np.where((nc.variables['time'][:]<=end_time)&(nc.variables['time'][:]>=start_time))[0]
        if len(t_interval)>0:
            t_ind.append([t_interval[0],t_interval[-1]])
        else:
            t_ind.append(False)
    return t_ind


def Append_to_hdf( project, hdfp, st_i,end_i):
    """thread worker function"""
    src_nc = Build_src_nc( project, hdfp )
    with h5py.File(hdfp['h5f'],"a") as hdf:
        for itime in range(st_i,end_i):
            try:
                a = ERA_to_mdays( project,
                                  hdfp['sds'], 
                                  hdfp['hdtime'][itime],
                                  src_nc)
                hdf[hdfp['sds']][itime,:,:]=a
                del a
                toprint = project['modis_days'][itime]
                hdf['time'][itime] = hdfp['hdtime'][itime]
            except:
                print project['modis_days'][itime],'NOVALUE'
                hdf['time'][itime] = hdfp['hdtime'][itime]
    del src_nc


def ERA_to_mdays(project,dset,start_time,src_nc):
    '''era_series = Avg_to_mdays(x,y,dset=['sd' OR 't2m']'''
    
    nclon = src_nc[0].variables['longitude'][:]
    nclat = src_nc[0].variables['latitude'][:]
    dlon  = nclon[1]-nclon[0]
    dlat  = nclat[1]-nclat[0]

    nclon = nclon - 0.5*dlon
    nclat = nclat + 0.5*dlat

    x_s,x_e,y_s,y_e = Get_spatial_indexes(nclon,nclat,project['aoi'])
    # start_time assigned
    end_time   = start_time+16*24
    intervals = ID_intervals(start_time,end_time,src_nc)
    era_vals = False
    for k in range(len(intervals)):
        if intervals[k]:
            st = int(intervals[k][0])
            fi = int(intervals[k][1])
            tot= len(src_nc[k].variables['time'][:])
            if fi==(tot-1):
                fi = None
            
            data = src_nc[k].variables[dset][st:fi,y_s:y_e,x_s:x_e]
            if era_vals:
                era_series = np.concatenate([era_series,data],axis=0) 
            else:
                era_series = data
                era_vals = True
    if (era_series.shape[0]==16) or (era_series.shape[0]==64):
        mean16day = np.mean(era_series,axis=0)
    else:
        print 'PROBLEMS WITH ERA SAMPLE!'
        mean16day = 'PROBLEMS!!'
    return mean16day


def Get_wgs84_aoi( aoi ):
    # Load aoi extents
    aoi_bbox = Mk_bbox(aoi['xmin'],aoi['ymin'],aoi['xmax'],aoi['ymax'])

    # transform aoi from native prj to wgs84
    prj_transform = osr.CoordinateTransformation(aoi['srs'],wgs84) #src,dst
    # if CoordinateTransformation fails, it will return null:
    if prj_transform == None:
        print '[ERROR] Could not reproject AOI box to MODIS sinusoidal'
        sys.exit( 1 )

    aoi_bbox.Transform(prj_transform)  
    
    xmin,xmax,ymin,ymax = aoi_bbox.GetEnvelope() 
    return xmin,ymin,xmax,ymax


### Spatial functions
def Get_spatial_indexes(x_var,y_var,aoi):
    '''x_s,x_e,y_s,y_e = Get_spatial_indexes(x_var,y_var)
    x_sub = x_var[x_s:x_e]
    y_sub = y_var[y_s:y_e]
    
    era coords are cell ctr.
    '''
    aoi_xmin,aoi_ymin,aoi_xmax,aoi_ymax = Get_wgs84_aoi( aoi )
    
    try:
        x_s = int(np.where(x_var<aoi_xmin)[0][-1])
    except:
        x_s = None
    try:
        x_e = int(np.where(x_var>aoi_xmax)[0][0])
    except:
        x_e = None
    try:
        y_s = int(np.where(y_var>aoi_ymax)[0][-1])
    except:
        y_s = None
    try:
        y_e = int(np.where(y_var<aoi_ymin)[0][0])
    except:
        y_e = None

    return x_s,x_e,y_s,y_e


def Continue_era_hdf( project, hdfp ):
    print 'Continuing', hdfp['sds'], hdfp['h5f']
    start_end = Gen_appendexes( project, hdfp )
    if start_end:
        for i,s_e in enumerate(start_end):
            toprint = project['modis_days'][s_e[0]]
            p = multiprocessing.Process(target=Append_to_hdf, \
                                  args=(project, hdfp, s_e[0], s_e[1]))
            p.start()
            p.join() # main script waits for this child to grow up
            
            prog = int( float(i) / len(start_end) * 100 )
            report = '%s%% . ' % prog
            sys.stdout.write( report )
            sys.stdout.flush()
    
    sys.stdout.write("\n") 
    print hdfp['sds'],'FINISHED!'


def Era2hdf( project, sds ):
    hdfp = {}
    hdfp['sds'] = sds
    hdfp['hdtime'] = [Yearday2hrnum(basedate,mday) for \
                      mday in project['modis_days']]
    hdfp['years']  = [int(yearday[:4]) for \
                      yearday in project['modis_days']]
    hdfp['ydays']  = [int(yearday[4:]) for \
                      yearday in project['modis_days']]
    hdfp['appendnum'] = 5 # records to append per process
    hdfp['h5f'] = os.path.join(project['hdf_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')
    if not os.path.isfile(hdfp['h5f']):
        Start_era_hdf(project, hdfp )
    
    Continue_era_hdf( project, hdfp )

    
    print 'Bing!'
    
