## This is a simpler aggregator script, no temporal averaging, meant for metno (tam & rr) & nve (sd, swe, fsw) datasets
## Output is an hdf5 file for each dataset, basically just representing the entirety of each dataset's directory
## in a single file with a single format for metadata. We can see whether this gets subset later...

import os
import osr
import h5py
import gdal
import gdalconst
import multiprocessing
import numpy as np
import datetime as dt
from scipy.io.netcdf import netcdf_file as NetCDFFile
from ORG_tools import Countdown


### DEFINE PROJECTION, UTM33N
# sin_prj_string directly from MODIS hdf files:
utm33n = osr.SpatialReference()
utm33n.ImportFromEPSG(32633)
utm33n_string = utm33n.ExportToWkt()

dx = dy = 1000 # (m)
# The southwestern corner of the grid has the coordinates
llx = -75000  # and
lly = 6450000 # in the UTM33 system if you prefer this coordinate system.

nx = 1195
ny = 1550

ulx = llx
uly = lly + (ny * dy) + dy

x_var = ulx + dx*np.arange(nx)
y_var = uly - dy*np.arange(ny)

# hdf creation and
#from Scientific.IO.NetCDF import NetCDFFile

def Aggregate_climate_grids(project):
    '''metnc2hdf'''
    ## COMMON VARIABLES ######

    hdfp = {} # hdf parameters

    #  Get Date range from modis_days
    start_date = dt.datetime.strptime(str(project['modis_days'][0]), '%Y%j')
    end_date   = dt.datetime.strptime(str(project['modis_days'][-1]), '%Y%j')

    numdays = (end_date-start_date).days + 1
    daterange = [start_date + dt.timedelta(days=x) for x in range(0, numdays)]

    hdfp['start_date']  = start_date
    hdfp['daterange']   = daterange
    hdfp['appendnum']   = 5
    hdfp['climate_dir'] = '/space/wib_data/CLIMATE'
    hdfp['metno_fmt']   = os.path.join(hdfp['climate_dir'],'METNO','{0]}/{0}24hNOgrd1957on_{1}_{2}_{3}.nc') # .format (hdfp['sds'], date.year, '{:02d}'.format(date.month), '{:02d}'.format(date.day))
    hdfp['nve_fmt']   = os.path.join(hdfp['climate_dir'],'METNO','{0}/{0}_{1}_{2}_{3}.asc') # .format (hdfp['sds'], date.year, '{:02d}'.format(date.month), '{:02d}'.format(date.day))

    ### NVE ###
    # sd
    hdfp['dtype'] = np.dtype('int32')
    hdfp['scale_factor'] = 1
    hdfp['add_offset']   = 0
    hdfp['fill_value']   = int(65535)

    hdfp['sds'] = 'sd'
    hdfp['units'] = 'mm snowdepth'
    hdfp['h5f'] = os.path.join(hdfp['climate_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf(hdfp)

    Continue_hdf( hdfp )

    # swe
    hdfp['dtype'] = np.dtype('int32')
    hdfp['scale_factor'] = 0.1
    hdfp['add_offset']   = 0
    hdfp['fill_value']   = int(65535)

    hdfp['sds'] = 'swe'
    hdfp['units'] = 'mm water equivalent'
    hdfp['h5f'] = os.path.join(hdfp['climate_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf(hdfp)

    Continue_hdf( hdfp )

    # fsw
    hdfp['dtype'] = np.dtype('int16')
    hdfp['scale_factor'] = 1
    hdfp['add_offset']   = 0
    hdfp['fill_value']   = int(255) # 254 is also code, for 'bare ground'

    hdfp['sds'] = 'fsw'
    hdfp['units'] = 'mm water equivalent'
    hdfp['h5f'] = os.path.join(hdfp['climate_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf(hdfp)

    Continue_hdf( hdfp )

    ## MET.NO ###
    hdfp['dtype'] = np.dtype('int16') # or np.dtype('>i2') ?
    hdfp['scale_factor'] = 0.1
    hdfp['add_offset']   = 0
    hdfp['fill_value']   = int(-999)

    # tam
    hdfp['sds'] = 'tam'
    hdfp['units'] = 'degrees Kelvin'
    hdfp['h5f'] = os.path.join(hdfp['climate_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf(hdfp)

    Continue_hdf( hdfp )

    # rr
    hdfp['sds'] = 'rr'
    hdfp['units'] = 'mm/day'
    hdfp['h5f'] = os.path.join(hdfp['climate_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf(hdfp)

    Continue_hdf( hdfp )


    print 'BING!!'


def Mk_hdf( hdfp ):
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
        arr_out = hdf.create_dataset(hdfp['sds'],dshape,dtype=hdfp['dtype'], \
          chunks=True,compression='lzf') #compression='gzip' or 'szip'
        # arr_out[:] = hdfp['fill_value'] # Too slow!!
        x_out = hdf.create_dataset("x",data=x_var)
        y_out = hdf.create_dataset("y",data=y_var)
        t_out = hdf.create_dataset("time", (len(hdfp['time_var']),), \
                                   dtype='int16')
        date_out = hdf.create_dataset("date",data = [d.isoformat() for d in hdfp['daterange']])
        yr_out= hdf.create_dataset("year",data=hdfp['years'])
        # Add metadata
        # access later: hdf[sds].attrs['scale_factor']
        t_out.attrs.create('time_format','Days since %s' \
                                        % hdfp['start_date'])
        t_out.attrs.create('basedate',str(hdfp['start_date']))
        arr_out.attrs.create('projection',utm33n_string)
        arr_out.attrs.create('scale_factor',hdfp['scale_factor'])
        arr_out.attrs.create('add_offset',hdfp['add_offset'])
        arr_out.attrs.create('fill_value',hdfp['fill_value'])
        arr_out.attrs.create('dx',dx)
        arr_out.attrs.create('dy',dy)
        arr_out.attrs.create('units',hdfp['units'])


def Gen_appendexes( hdfp ):
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
            while (i_st+hdfp['appendnum'])<len(hdfp['modis_days']):
                appind.append((i_st,i_st+hdfp['appendnum']))
                i_st+=hdfp['appendnum']
            appind.append((i_st,len(hdfp['modis_days'])))
        else:
            print hdfp['sds'],'already done??'
            appind=False
    return appind


def Load_metno_arr(hdfp, date):
    metno_fn_fmt = hdfp['metno_fmt']
    metno_fn = metno_fn_fmt.format (hdfp['sds'], date.year, '{:02d}'.format(date.month), '{:02d}'.format(date.day))
    nc = NetCDFFile(metno_fn,'r')
    arr = nc.variables[hdfp['sds']][:]
    return np.flipud(arr)[:,:,0] # flipud??


def Load_nve_arr(hdfp, date):
    nve_fn_fmt = hdfp['nve_fmt']
    nve_fn = nve_fn_fmt.format (hdfp['sds'], date.year, '{:02d}'.format(date.month), '{:02d}'.format(date.day))
    ascii  = gdal.Open(nve_fn,gdalconst.GA_ReadOnly)
    band   = ascii.GetRasterBand(1)
    arr = band.ReadAsArray()
    del ascii, band
    return arr # np.flipud??


def Load_arr(hdfp, date):
    if hdfp['sds'] in ['sd', 'fsw', 'swe']:
        a = Load_nve_arr(hdfp, date)
    elif hdfp['sds'] in ['tam','rr']:
        a = Load_metno_arr(hdfp, date)
    else:
        print 'dataset not recognized'
        a = None

    return a


def Append_to_hdf(  hdfp, st_i,end_i):
    """thread worker function"""

    with h5py.File(hdfp['h5f'],"a") as hdf:
        for itime in range(st_i,end_i):
            try:
                a = Load_arr( hdfp, itime)
                hdf[hdfp['sds']][itime,:,:]=a
                del a
                toprint = hdfp['modis_days'][itime]
                hdf['time'][itime] = hdfp['time_var'][itime]
            except:
                hdf[hdfp['sds']][itime,:,:]=hdfp['fill_value']
                print hdfp['modis_days'][itime],'NODATA'

                hdf['time'][itime] = hdfp['time_var'][itime]


def Continue_hdf( hdfp ):
    print 'Continuing', hdfp['sds'], hdfp['h5f']
    start_end = Gen_appendexes( hdfp )
    if start_end:
        progress_bar = Countdown(len(start_end))
        for i,s_e in enumerate(start_end):
            p = multiprocessing.Process(target=Append_to_hdf, \
                                  args=( hdfp, s_e[0], s_e[1]))
            p.start()
            p.join() # main script waits for this child to grow up

            progress_bar.check(i)
        progress_bar.flush()
    print hdfp['sds'],'FINISHED!'
