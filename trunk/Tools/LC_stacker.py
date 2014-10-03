#  TIF directory to HDF5 conversion for MODIS archive
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
'''
This module contains functions for:
* Reading a specified set of geoTIFF files, each representing a data field
* Creating an hdf5 array file with dimensions [field,y,x]
from Tools.ORG_tools import Parse_input
project = Parse_input('../STINT/k3tif/k3tif_INPUT.txt')
project['hdf_dir'] = '/home/wiley/wrk/STINT/k3tif/HDF'
project['prj_name'] = 'k3tif'
fn = project['lc'].values()[0]
'''

__author__ = 'wiley'

import os,gdal,h5py
import numpy as np
try:
    from ORG_tools import Countdown
except:
    from Tools.ORG_tools import Countdown

### BUILDER FUNCTIONS
def Mk_hdf( hdfp ):
    '''Creates the hdf5 file, and defines the dimensions and datasets.

    :param project: (dict) py-stint project parameters
    :param hdfp: (dict) conversion parameters
    :param x_var: (np.array) [x] cell center coordinates bounding the AOI
    :param y_var: (np.array) [y] cell center coordinates bounding the AOI
    :param modis_md: (dict) keys to unpack the MODIS dataset stored values.
    :return:
    '''
    print 'Creating',hdfp['h5f']

    # field names

    with h5py.File(hdfp['h5f'],"w") as hdf:
        dshape=(len(hdfp['fields']),len(hdfp['y_var']),len(hdfp['x_var']))
        arr_out = hdf.create_dataset('lc',dshape,dtype=hdfp['datatype'], \
          chunks=True,compression='lzf') #compression='gzip' or 'szip'
        #arr_out[:] = hdfp['fill_vals'][0]
        x_out = hdf.create_dataset("x",data=hdfp['x_var'])
        y_out = hdf.create_dataset("y",data=hdfp['y_var'])
        fillv = hdf.create_dataset("fill_vals",data=hdfp['fill_vals'])
        field_names =  hdf.create_dataset("fields",data=hdfp['fields'].keys())
        arr_out.attrs.create('projection',hdfp['projection'])
        arr_out.attrs.create('fill_value',hdfp['fill_vals'][0])
        arr_out.attrs.create('dx',hdfp['dx'])
        arr_out.attrs.create('dy',hdfp['dy'])


def Fill_lc_hdf( hdfp ):
    '''Populate a Landcover hdf5 array file.

    There is no mosaicing or subsetting process here, since the landcover IS the AOI.
    Additionally, there are probably not as many landcover fields as there are
    MODIS or ERA days to aggregate. Because of this simplicity, the hdf5 population
    for landcover doesn't have the same multiprocessing complexity.  It just loops
    through the landcover tifs, reads in the array, and copies it to the appropriate
    part of the hdf5 array file.

    :param hdfp: (dict) array parameters
    :return: None
    '''
    print 'Continuing lc',hdfp['h5f']

    numfields = len(hdfp['fields'])
    progress_bar = Countdown(numfields)
    for i in range(numfields):
        with h5py.File(hdfp['h5f'], "a") as hdf:
            ras = gdal.Open(hdfp['fields'].values()[i])
            r = ras.GetRasterBand(1)
            a = r.ReadAsArray()
            hdf['lc'][i,:,:] = a

        progress_bar.check(i)
    progress_bar.flush()
    print 'Landcover hdf5 FINISHED!'


# Landcover Tif Functions
def Get_lc_params( fn ):
    '''Return Extents from a specified geoTIFF.


    'dx'   : raster cell size, x direction
    'dy'   : raster cell size, y direction
    'xmin' : minimum x extent (left edge)
    'ymin' : minimum y extent (bottom edge)
    'xmax' : maximum x extent (right edge)
    'ymax' : maximum y extent (top edge)

    :param modis_dir: (string) path to MODIS archive
    :param dset: (string) MODIS product name, ie 'MCD43A3'
    :param mtile: (string) MODIS tile number, ie 'h19v03'
    :return: (dict)
    '''
    out = {}
    ras      = gdal.Open( fn )
    gt       = ras.GetGeoTransform()
    r        = ras.GetRasterBand(1)
    out['dx']   = gt[1]
    out['dy']   = gt[-1] # abs if this should be negatives
    out['xmin'] = gt[0]
    out['ymax'] = gt[3]
    out['xmax'] = out['xmin'] + ras.RasterXSize*out['dx']
    out['ymin'] = out['ymax'] + ras.RasterYSize*out['dy']

    x_var = np.arange(out['xmin'],out['xmax'],out['dx']) # xmax+.5 to include last val
    y_var = np.arange(out['ymax'],out['ymin'],out['dy'])

    out['x_var'] = x_var
    out['y_var'] = y_var
    out['projection'] = ras.GetProjection()
    out['datatype']   = gdal.GetDataTypeName(r.DataType)
    out['fill_value'] = r.GetNoDataValue()
    if 'int' in out['datatype'] or 'Int' in out['datatype']:
        out['fill_value'] = int(out['fill_value'])
    # if datatype is int, fill_value should probably be converted to int?

    return out


def Load_lc_params( hdfp ):
    fill_vals = []
    # 'fill_vals'

    # Open first tif, fill common params:
    lcp = Get_lc_params( hdfp['fields'].values()[0])

    for k in ['dx','dy','x_var','y_var','projection','datatype']:
        hdfp[k] = lcp[k]

    # fill up fill vals, make sure datatypes are the same
    fill_vals.append(lcp['fill_value'])
    for i in range(1,len(hdfp['fields'].values())):
        lcp = Get_lc_params(hdfp['fields'].values()[i])
        if lcp['fill_value']!=fill_vals[-1]:
            print 'fill values vary between landcover datasets; this is ok!'
        if lcp['datatype']!=hdfp['datatype']:
            print '[ERROR] datatypes vary between landcover datasets; NOT ok!'
            print hdfp['fields'].keys()[i], '--> check if hdf5 transcriptions are faithful!'
        fill_vals.append(lcp['fill_value'])

    hdfp['fill_vals'] = fill_vals

    return hdfp


def Gen_lc_hdf( project ):
    hdfp = {}
    hdfp['fields'] = project['lc']
    hdfp['h5f'] = os.path.join(project['hdf_dir'], project['prj_name']+'_lc.hdf5')

    hdfp = Load_lc_params( hdfp )
    if not os.path.isfile(hdfp['h5f']):
        Mk_hdf( hdfp )

    Fill_lc_hdf( hdfp )
    print 'Bing!!'