## This is a simpler subsetter script, specifically for metno tam and rr grids, initially in daily netcdf files
## two steps: aggregate/time-average (nc2hdf)
## spatial subset

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

x_var = None # todo
y_var = None # todo

# hdf creation and
#from Scientific.IO.NetCDF import NetCDFFile

def Aggregate_metno_grids(project):
    '''metnc2hdf'''
    # tam

    hdfp['sds'] = 'tam'

    modis_days = project['modis_days']

    hdfp = {} # hdf parameters
    hdfp['appendnum']   = 5
    hdfp['metno_dir']   = project['metno_dir']
    hdfp['modis_days']  = project['modis_days']
    hdfp['NODATA']      = np.nan
    hdfp['years'] = [int(str(modis_days[i])[:4]) for i in \
                       range(len(modis_days))]
    hdfp['ydays'] = [int(str(modis_days[i])[4:]) for i in \
                       range(len(modis_days))]

    hdfp['basedate'] = dt.datetime.strptime(str(modis_days[0]),\
                                                '%Y%j').date()

    hdfp['time_var'] = np.array([(dt.datetime.strptime(str(mday), \
      '%Y%j').date()-hdfp['basedate']).total_seconds()/86400. for \
      mday in modis_days],dtype='int16')

    hdfp['h5f'] = os.path.join(project['hdf_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    metno_md = Get_metno_md('/space/wib_data/CLIMATE/METNO/tam/tam24hNOgrd1957on_2006_07_02.nc')
    hdfp['fill_value'] = -999



def Mk_hdf( hdfp, x_var, y_var, metno_md ):
    '''Creates the hdf5 file, and defines the dimensions and datasets.

    :param project: (dict) py-stint project parameters
    :param hdfp: (dict) conversion parameters
    :param x_var: (np.array) [x] cell center coordinates bounding the AOI
    :param y_var: (np.array) [y] cell center coordinates bounding the AOI
    :param modis_md: (dict) keys to unpack the MODIS dataset stored values.
    :return:
    '''
    print 'Creating',hdfp['sds'],hdfp['h5f']

    with h5py.File(hdfp['h5f'],"w") as hdf:
        dshape=(len(hdfp['time_var']),len(y_var),len(x_var))
        #arr_out = hdf.create_dataset(hdfp['sds'],dshape,dtype='int16', \
        arr_out = hdf.create_dataset(hdfp['sds'],dshape,dtype=modis_md['dtype'], \
          chunks=True,compression='lzf') #compression='gzip' or 'szip'
        arr_out[:] = modis_md['fill_value']
        x_out = hdf.create_dataset("x",data=x_var)
        y_out = hdf.create_dataset("y",data=y_var)
        t_out = hdf.create_dataset("time", (len(hdfp['time_var']),), \
                                   dtype='int16')
        yr_out= hdf.create_dataset("year",data=hdfp['years'])
        day_out=hdf.create_dataset("yday",data=hdfp['ydays'])
        # Add metadata
        # access later: hdf[sds].attrs['scale_factor']
        t_out.attrs.create('time_format','Days since %s' \
                                        % hdfp['basedate'])
        t_out.attrs.create('basedate',str(hdfp['basedate']))
        arr_out.attrs.create('projection',sin_prj)
        arr_out.attrs.create('scale_factor',modis_md['scale_factor'])
        arr_out.attrs.create('add_offset',modis_md['add_offset'])
        arr_out.attrs.create('fill_value',modis_md['fill_value'])
        arr_out.attrs.create('dx',modis_md['dx'])
        arr_out.attrs.create('dy',modis_md['dy'])


def Get_metno_md(metno_fn):
    '''Return a dict giving the keys to unpack the optimized
    MODIS dataset values.

    metadata should be applied as follows:
    (1) REPLACE fill_value with NAN
    (2) v = c*z + b, where:
      v = measured value
      c = stored cell value
      z = scale_factor
      b = add_offset

    returns a dict with the following metadata parameters (all integer):

    'scale_factor' : z
    'add_offset'   : b
    'fill_value'   : dataset NAN value

    :param modis_dir: (string) path to MODIS_archive
    :param dset: (string) MODIS product, ie 'MCD43A3'
    :param dnum: (integer) subdataset number, ie 7 -> BSA_vis
    :param mtile: (string) MODIS tile number, ie 'h19v03'
    :return: (dict) MODIS metadata
    '''
    out      = {}

    metno = NetCDFFile(metno_fn,'r')

    out['dtype'] = np.dtype('int16') # or np.dtype('>i2') ?
    out['scale_factor'] = 0.1
    out['add_offset']   = 0
    out['fill_value']   = int(-999)
    out['dx'] = 1000
    out['dy'] = 1000
    return out


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

