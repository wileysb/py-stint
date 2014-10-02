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
'''

__author__ = 'wiley'

import os,sys,osr,gdal,h5py,multiprocessing
import numpy as np
import datetime as dt
from ORG_tools import Countdown


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


def Continue_modis_hdf( project, hdfp ):
    '''Continue population of a MODIS hdf5 array file.

    If the file was just created with Mk_hdf, it will be
    populated from the beginning.

    If the process was previously interrupted at timestep <i>, it will be
    populated from <i-1>, ensuring that all resulting output timesteps
    are complete.

    This function uses multiprocessing, it is not recommended to exploit
    this by calling more than one child at a time.  Every child process
    would try to modify the same file, which might adversely affect the
    integrity of the resulting file.

    Multiprocessing is used to make sure the RAM buffers are cleared
    periodically while looping through the potentially memory-intensive
    task of aggregating all MODIS tiles to a single hdf5 array file.

    :param project: (dict) py-stint project parameters
    :param hdfp: (dict) conversion parameters
    :return: None
    '''
    print 'Continuing',hdfp['sds'],hdfp['h5f']
    start_end = Gen_appendexes( hdfp )
    if start_end:
        progress_bar = Countdown(len(start_end))
        for i,s_e in enumerate(start_end):
            arrrrrgs = ( project, hdfp, s_e[0], s_e[1] )
            p = multiprocessing.Process(target=Append_to_hdf,
                                          args=arrrrrgs)
            p.start()
            p.join() # main script waits for this child to grow up

            progress_bar.check(i)

        progress_bar.flush()
        print hdfp['sds'],'FINISHED!'


def Gen_appendexes( hdfp ):
    '''Generates the timestep intervals for hdf5 population.

     The intervals are used to coordinate buffer-clearing child processes,
     as well as identifying where to pick up if the aggregation process
     was previously interrupted.

    :param hdfp: (dict) conversion parameters
    :return: None
    '''

    with h5py.File(hdfp['h5f'], "r") as hdf:
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


def Append_to_hdf( project, hdfp, st_i, end_i):
    '''Populate the hdf5 file by aggregating timesteps identified.

    st_i and end_i are start and end indices, aggregating
    MODIS intervals beginning with st_i and ending with the record
    before end_i.

    if, for example:
     modis_days = ['2001049','2001057','2001065','2001073']
     st_i  = 0
     end_i = 2
    then:
     MODIS intervals '2001049' and '2001057','2001065' will be aggregated
     and written out to the hdf5 array file


    :param project: (dict) py-stint project parameters
    :param hdfp: (dict) conversion parameters
    :param st_i: (int) timestep interval, start
    :param end_i: (int) timestep interval, end
    :return: None
    '''
    toprint = 'Continuing'#,hdfp['sds'],hdfp['h5f']
    with h5py.File(hdfp['h5f'], "a") as hdf:
        for itime in range(st_i,end_i):
            try:
                a=Build_modis_mosaic( project, hdfp,
                                      hdfp['modis_days'][itime] )
                hdf[hdfp['sds']][itime,:,:] = a
                hdf['time'][itime]=hdfp['time_var'][itime]
                toprint = hdfp['modis_days'][itime]
            except:
                print hdfp['modis_days'][itime],'NOVALUE'
                hdf['time'][itime]=hdfp['time_var'][itime]








def Gen_lc_hdf( project ):
    hdfp = {}
    hdfp{'fields'} = project['lc']
    hdfp['h5f'] = os.path.join(project['hdf_dir'],'lc.hdf5')
    # get x and y arrays
    # 