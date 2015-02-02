## This is a simpler subsetter script, specifically for metno tam and rr grids, initially in daily netcdf files
## two steps: aggregate/time-average (nc2hdf)
## spatial subset

import os
import osr
import h5py
import multiprocessing
import gdalconst
import glob
import gdal
import numpy as np
import datetime as dt
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

def Aggregate_nve_grids(project):
    '''metnc2hdf'''
    ## COMMON VARIABLES ######

    modis_days = project['modis_days']
    hdfp = {} # hdf parameters
    hdfp['appendnum']   = 5
    hdfp['nve_dir']   = '/space/wib_data/CLIMATE/NVE' # project['metno_dir']
    hdfp['nve_fn_fmt']  = {}
    fsw_fmt =  'fsw/fsw_{0}_{1}_{2}.nc' # .format(YYYY, MM, DD)
    swe_fmt =  'swe/swe_{0}_{1}_{2}.nc' # .format(YYYY, MM, DD)
    sd_fmt =  'sd/sd_{0}_{1}_{2}.nc' # .format(YYYY, MM, DD)
    hdfp['nve_fn_fmt']['fsw'] = os.path.join(hdfp['nve_dir'], fsw_fmt)
    hdfp['nve_fn_fmt']['swe'] = os.path.join(hdfp['nve_dir'], swe_fmt)
    hdfp['nve_fn_fmt']['sd'] = os.path.join(hdfp['nve_dir'], sd_fmt)
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



    ##### fsw #######  mm?  #############
    hdfp['sds'] = 'fsw'
    fn = glob.glob(os.path.join(hdfp['nve_dir'], hdfp['sds'], '*.asc'))[0]
    nve_md = Get_nve_md(fn)
    hdfp['missing_value'] = nve_md['missing_value']
    hdfp['h5f'] = os.path.join(project['nve_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf(hdfp, nve_md)

    Continue_nve_hdf( project, hdfp )


    ##### swe #######  mm?  #############
    hdfp['sds'] = 'swe'
    fn = glob.glob(os.path.join(hdfp['nve_dir'], hdfp['sds'], '*.asc'))[0]
    nve_md = Get_nve_md(fn)
    hdfp['missing_value'] = nve_md['missing_value']
    hdfp['h5f'] = os.path.join(project['nve_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf(hdfp, nve_md)

    Continue_nve_hdf( project, hdfp )



    ##### sd  #######  cm?  #############
    hdfp['sds'] = 'sd'
    fn = glob.glob(os.path.join(hdfp['nve_dir'], hdfp['sds'], '*.asc'))[0]
    nve_md = Get_nve_md(fn)
    hdfp['missing_value'] = nve_md['missing_value']
    hdfp['h5f'] = os.path.join(project['nve_dir'],project['prj_name']+ \
                               '_'+hdfp['sds']+'.hdf5')

    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf(hdfp, nve_md)

    Continue_nve_hdf( project, hdfp )



def Mk_hdf( hdfp, nve_md ):
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
        arr_out = hdf.create_dataset(hdfp['sds'],dshape,dtype=nve_md['dtype'], \
          chunks=True,compression='lzf') #compression='gzip' or 'szip'
        # arr_out[:] = metno_md['fill_value'] # Too slow!!
        x_out = hdf.create_dataset("x",data=x_var)
        y_out = hdf.create_dataset("y",data=y_var)
        t_out = hdf.create_dataset("time", (len(hdfp['time_var']),), \
                                   dtype='int32') # unsigned int16 for sd & swe, unsigned int8 for fsw
        yr_out= hdf.create_dataset("year",data=hdfp['years'])
        day_out=hdf.create_dataset("yday",data=hdfp['ydays'])
        # Add metadata
        # access later: hdf[sds].attrs['scale_factor']
        t_out.attrs.create('time_format','Days since %s' \
                                        % hdfp['basedate'])
        t_out.attrs.create('basedate',str(hdfp['basedate']))
        arr_out.attrs.create('projection',utm33n_string)
        arr_out.attrs.create('scale_factor',nve_md['scale_factor'])
        arr_out.attrs.create('add_offset',nve_md['add_offset'])
        arr_out.attrs.create('fill_value',nve_md['fill_value'])
        arr_out.attrs.create('dx',nve_md['dx'])
        arr_out.attrs.create('dy',nve_md['dy'])


def Get_nve_md(fn):
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

    file = gdal.Open(fn, gdalconst.GA_ReadOnly)
    band = file.GetRasterBand(1)
    missing_val = band.GetNoDataValue()

    del file, band
    out      = {}

    out['dtype'] = np.dtype('int16') # int32 for 65532?
    out['scale_factor'] = 1 # todo confirm
    out['add_offset']   = 0 # todo confirm
    out['fill_value']   = int(missing_val)
    out['dx'] = 1000
    out['dy'] = 1000
    return out


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


def Load_nve_arr(hdfp, date):
    nve_fn_fmt = hdfp['nve_fn_fmt'][hdfp['sds']]
    nve_fn = nve_fn_fmt.format(date.year, '{:02d}'.format(date.month), '{:02d}'.format(date.day))
    ascii  = gdal.Open(nve_fn,gdalconst.GA_ReadONly)
    band   = ascii.GetRasterBand(1)
    arr = band.ReadAsArray()
    del ascii, band
    return arr[:,:,0] # np.flipud??


def NVE_to_mdays(hdfp, itime):
    # TODO: replace string placeholders with actual code

    modis_start = hdfp['modis_days'][itime]
    modis_start = dt.datetime.strptime(str(modis_start), '%Y%j')
    numdays = 16
    modis_range = [modis_start + dt.timedelta(days=x) for x in range(0, numdays)]
    src_range = np.ones((len(modis_range), ny, nx)) * hdfp['missing_value'] # missing value here would screw up interpolation if some dates failed to load.  However, the functions here don't have any error tolerance, so it should fail at 'append to hdf' rather than load_arr or _to_mdays
    for i in range(len(modis_range)):
        date = modis_range[i]
        src_range[i] = Load_nve_arr(hdfp, date)

    out = np.mean(src_range, axis=0) # todo nans? Since these are modelled / interpolated datasets, nans should be spatial constants?

    return out


def Append_to_hdf(  hdfp, st_i,end_i):
    """thread worker function"""

    with h5py.File(hdfp['h5f'],"a") as hdf:
        for itime in range(st_i,end_i):
            try:
                a = NVE_to_mdays( hdfp, itime)
                hdf[hdfp['sds']][itime,:,:]=a
                del a
                toprint = hdfp['modis_days'][itime]
                hdf['time'][itime] = hdfp['time_var'][itime]
            except:
                print hdfp['modis_days'][itime],'NOVALUE'
                hdf['time'][itime] = hdfp['time_var'][itime]


def Continue_nve_hdf( project, hdfp ):
    print 'Continuing', hdfp['sds'], hdfp['h5f']
    start_end = Gen_appendexes( project, hdfp )
    if start_end:
        progress_bar = Countdown(len(start_end))
        for i,s_e in enumerate(start_end):
            toprint = project['modis_days'][s_e[0]]
            p = multiprocessing.Process(target=Append_to_hdf, \
                                  args=( hdfp, s_e[0], s_e[1]))
            p.start()
            p.join() # main script waits for this child to grow up

            progress_bar.check(i)
        progress_bar.flush()
    print hdfp['sds'],'FINISHED!'
