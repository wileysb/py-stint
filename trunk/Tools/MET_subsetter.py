## This is a simpler subsetter script, specifically for metno tam and rr grids, initially in daily netcdf files
## two steps: aggregate/time-average

import os
import sys
import osr
import h5py
import multiprocessing
import numpy as np
import datetime as dt
from scipy.io.netcdf import netcdf_file as NetCDFFile
from MODIS_aoi import Mk_bbox
from ORG_tools import Yearday2hrnum
from ORG_tools import Countdown

### DEFINE PROJECTION, UTM33N
# sin_prj_string directly from MODIS hdf files:
utm33n = osr.SpatialReference()
utm33n.ImportFromEPSG(32633)
utm33n_string = utm33n.ExportToWkt()

dx = dy = 1000 # (m)
# The southwestern corner of the grid has the coordinates
llX = -75000  # and
llY = 6450000 # in the UTM33 system if you prefer this coordinate system.
nx = 1195
ny = 1550

# hdf creation and
#from Scientific.IO.NetCDF import NetCDFFile

def Aggregate_metno_grids(project):



def Start_metno_hdf( project, hdfp):
    print 'Creating',hdfp['sds'],hdfp['h5f']
    src_nc = Build_src_nc( project, hdfp)
    missing_val  = src_nc[0].variables[hdfp['sds']].missing_value
    ncdtype = src_nc[0].variables[hdfp['sds']].data.dtype
    nclon   = src_nc[0].variables['longitude'][:]
    nclat   = src_nc[0].variables['latitude'][:]

    dx = 1000
    dy = 1000

    nclon = nclon - 0.5*dlon
    nclat = nclat + 0.5*dlat

    hdlon = nclon[x_s:x_e]
    hdlat = nclat[y_s:y_e]

    with h5py.File(hdfp['h5f'],"w") as hdf:
        dshape=(len(hdfp['hdtime']),len(hdlat),len(hdlon))
        #arr_out = hdf.create_dataset(hdfp['sds'],dshape,dtype='int16',chunks=True,compression='lzf')
        arr_out = hdf.create_dataset(hdfp['sds'],dshape,dtype=ncdtype,chunks=True,compression='lzf')
        arr_out[:] = missing_val
        arr_out.attrs['scale_factor'] = src_nc[0].variables[hdfp['sds']].scale_factor
        arr_out.attrs['add_offset']   = src_nc[0].variables[hdfp['sds']].add_offset
        arr_out.attrs['units']        = src_nc[0].variables[hdfp['sds']].units
        arr_out.attrs['projection']   = utm33n_string
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



def Continue_era_hdf( project, hdfp ):
    print 'Continuing', hdfp['sds'], hdfp['h5f']
    start_end = Gen_appendexes( project, hdfp )
    if start_end:
        progress_bar = Countdown(len(start_end))
        for i,s_e in enumerate(start_end):
            toprint = project['modis_days'][s_e[0]]
            p = multiprocessing.Process(target=Append_to_hdf, \
                                  args=(project, hdfp, s_e[0], s_e[1]))
            p.start()
            p.join() # main script waits for this child to grow up

            progress_bar.check(i)
        progress_bar.flush()
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

