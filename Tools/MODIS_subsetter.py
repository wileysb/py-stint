#  HDF to HDF5 conversion for MODIS archive 
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
* Reading a specified set of MODIS product, dataset, tiles, and dates
* Reading a shapefile or raster defining an Area of Interest
* Creating an hdf5 array file with dimensions [time,y,x]
  where the y and x dimensions are subset to the minimum dimensions
  completely containing the Area of Interest
'''
import os,sys,osr,gdal,h5py,multiprocessing,glob
import numpy as np
import datetime as dt
from MODIS_aoi import Mk_bbox
from MODIS_librarian import Get_modis_fn
from ORG_tools import Countdown

### DEFINE MODIS SINUSOIDAL PROJECTION
# sin_prj_string directly from MODIS hdf files:
sin_prj = 'PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
sin_srs = osr.SpatialReference()
sin_srs.ImportFromWkt(sin_prj)

### READER FUNCTIONS
def Get_modis_extent(modis_dir, dset, mtile):
    '''Return Extents from a specified MODIS product and tile.

    Open the first MODIS file in the appropriate directory
    (MODIS_ARCHIVE/dset/mtile/*.hdf) and return a dictionary
    with the following parameters, specifiying the MODIS raster
    grid extents for that tile:

    'dx'   : raster cell size, x direction, sinusoidal projection
    'dy'   : raster cell size, y direction, sinusoidal projection
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
    src_pattern = os.path.join(modis_dir,dset,mtile,'*.hdf')
    src_fn = glob.glob(src_pattern)[0] # take first file in that folder!
    ras      = gdal.Open(src_fn)
    dset0    = ras.GetSubDatasets()[0][0] # first listed subdataset
    del ras
    ras = gdal.Open(dset0)
    gt          = ras.GetGeoTransform()
    out['dx']   = gt[1]
    out['dy']   = abs(gt[-1]) # abs if this should be negatives
    out['xmin'] = gt[0]
    out['ymax'] = gt[3]
    out['xmax'] = out['xmin'] + ras.RasterXSize*out['dx']
    out['ymin'] = out['ymax'] - ras.RasterYSize*out['dy']
    
    return out


def Get_mosaic_extent( project, dset ):
    '''Return x and y cell center coordinates of a
    MODIS mosaic bounding the project AOI.

    :param project: (dict) py-stint project parameters
    :param dset: (string) MODIS product, ie 'MCD43A3'
    :return: (np.array,np.array) x_var,y_var: 1-D arrays, cell center coords
    '''
    modis_tiles = project['modis_tiles']
    extents = Get_modis_extent(project['modis_dir'], 
                               dset,modis_tiles[0])
    for i in range(len(modis_tiles)): # fine that this repeats first tile
        nextent = Get_modis_extent(project['modis_dir'], 
                                   dset,modis_tiles[i])
        extents['xmax'] = max(extents['xmax'],nextent['xmax'])
        extents['ymax'] = max(extents['ymax'],nextent['ymax'])
        extents['xmin'] = min(extents['xmin'],nextent['xmin'])
        extents['ymin'] = min(extents['ymin'],nextent['ymin'])
    
    xmin,xmax,dx = extents['xmin'],extents['xmax'],extents['dx']
    ymin,ymax    = extents['ymin'],extents['ymax']
    dy           = -1*np.abs(extents['dy'])
    x_var = np.arange(xmin,xmax,dx) # xmax+.5 to include last val
    y_var = np.arange(ymax,ymin,dy) # ymin-.5 to include last val
    return x_var,y_var


def Get_modis_md(modis_dir, dset, dnum, mtile):
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
    src_pattern = os.path.join(modis_dir,dset,mtile,'*.hdf')
    src_fn = glob.glob(src_pattern)[0] # take first file in that folder!
    ras      = gdal.Open(src_fn)
    dname    = ras.GetSubDatasets()[dnum][0] # desired subdataset
    del ras
    ras  = gdal.Open(dname)
    out['dtype'] = ras.ReadAsArray().dtype
    gt    = ras.GetGeoTransform()
    rasmd = ras.GetMetadata_Dict()
    del ras
    if 'scale_factor' in rasmd.keys():
        out['scale_factor'] = float(rasmd['scale_factor'])
    else:
        out['scale_factor'] = 1
    if 'add_offset' in rasmd.keys():
        out['add_offset']   = float(rasmd['add_offset'])
    else:
        out['add_offset']   = 0
    
    # could fill value ever be float?
    out['fill_value']   = int(rasmd['_FillValue']) 
    
    out['dx'] = gt[1]
    out['dy'] = gt[5]
    return out


def Get_modis_tile( project, hdfp, modis_day, mtile ):
    '''

    a,x_var,y_var = Get_modis_tile(dset,dnum,modis_day,mtile)

    :param project: (dict) py-stint project parameters
    :param hdfp: (dict) conversion parameters
    :param modis_day: (string) YYYYDDD date, ie '2000049'
    :param mtile: (string) MODIS tile, ie 'h19v03'
    :return: (np.array,np.array,np.array)
              a,x_var,y_var
              a: [y,x] np.array giving stored MODIS values
              x_var: [x] x coordinates
              y_var: [y] y coordinates
    '''
    out      = {}
    src_fn   = Get_modis_fn(hdfp['modis_dir'], hdfp['dset'], modis_day, mtile)
    ras      = gdal.Open(src_fn)
    dname      = ras.GetSubDatasets()[hdfp['dnum']][0] # desired subdataset
    del ras
    ras = gdal.Open(dname)

    rasmd = ras.GetMetadata_Dict()
    # fill_value   = float(rasmd['_FillValue'])

    a = ras.ReadAsArray()
    #a = np.where(a==fill_value,hdfp['NODATA'],a)

    # define x and y coords for the array:
    gt          = ras.GetGeoTransform()
    out['dx']   = gt[1]
    out['dy']   = abs(gt[-1]) # abs if this should be negatives
    out['xmin'] = gt[0]
    out['ymax'] = gt[3]
    out['xmax'] = out['xmin'] + ras.RasterXSize*out['dx']
    out['ymin'] = out['ymax'] - ras.RasterYSize*out['dy']
    del ras

    x_var = np.arange(out['xmin'],out['xmax'],out['dx'])
    y_var = np.arange(out['ymax'],out['ymin'],-1*np.abs(out['dy']))
    
    # subset to aoi
    x_s,x_e,y_s,y_e = Get_spatial_indexes(x_var,y_var, project['aoi'] )
    
    x_sub = x_var[x_s:x_e]
    y_sub = y_var[y_s:y_e]
    a_sub = a[y_s:y_e,x_s:x_e]
    
    del a,x_var,y_var,out,gt,rasmd
    # array and the x/y coordinates that define it
    return a_sub,x_sub,y_sub


def Build_modis_mosaic( project, hdfp, modis_day):
    '''Return numpy array, a MODIS mosaic bounding the project AOI, populated
    with values for specified MODIS dset, tiles, and day.

    Gaps in coverage (ie missing or corrupt modis file) are filled
    with the nan/fill_value native to this modis dataset.

    :param project: (dict) py-stint project parameters
    :param hdfp: (dict) conversion parameters
    :param modis_day: (string) YYYYDDD date, ie '2000049'
    :return: (np.array) [y,x] MODIS mosaic bounding project AOI
    '''
    # Get extents of each tile
    x_var,y_var = Get_mosaic_extent( project, hdfp['dset'] )
    x_s,x_e,y_s,y_e = Get_spatial_indexes(x_var,y_var, project['aoi'])
    x_sub = x_var[x_s:x_e]
    y_sub = y_var[y_s:y_e]

    a = np.empty((len(y_sub),len(x_sub))) # or reverse x,y order?
    a[:] = hdfp['fill_value']
    for mtile in hdfp['modis_tiles']:
        try:
            a_,x_,y_ = Get_modis_tile( project, hdfp, modis_day, mtile )
            x_s,x_e,y_s,y_e = Get_subset_indexes(x_sub, y_sub ,x_, y_)
            a[y_s:y_e,x_s:x_e] = a_[:] 
        except:
            arrgs = (hdfp['dset'],modis_day,mtile)
            print 'Problems reading %s %s %s' % arrgs
            # Cells in this case stay NODATA
    
    return a #,x_sub,y_sub
    

def Get_subset_indexes(x_sub,y_sub,x_tile,y_tile):
    '''Return start and end indices for x and y coordinates.

    indices identify which portion of a MODIS tile overlaps the
    desired extents of the subset array.

    :param x_sub: (np.array) [x] cell center x coords, hdf5 subset
    :param y_sub:(np.array) [y] cell center y coords, hdf5 subset
    :param x_tile: (np.array) [x] cell center x coords, MODIS tile
    :param y_tile: (np.array) [y] cell center x coords, MODIS tile
    :return:(int,int,int,int) x_start,x_end,y_start,y_end
            x start & end, y start & end indices
            if, for example:
            a = [y,x] modis tile,
            a_sub = [y,x] portion of subset array overlapping that tile,
            then:
            a_sub = a[y_start:y_end,x_start:x_end]
    '''
    xmin = x_tile[0];xmax = x_tile[-1]
    pad = 0.2
    try:
        x_s = int(np.where((x_sub < xmin+pad) & (x_sub > xmin-pad))[0][0])
    except:
        x_s = None
    try:
        x_e = np.where((x_sub < xmax+pad) & (x_sub > xmax-pad))[0][0] 
        x_e = int(x_e+ 1)
    except:
        x_e = None
    
    ymin = y_tile[-1];ymax = y_tile[0]
    try:
        y_s = int(np.where((y_sub < ymax+pad) & (y_sub > ymax-pad))[0][0])
    except:
        y_s = None
    try:
        y_e = np.where((y_sub < ymin+pad) & (y_sub > ymin-pad))[0][0] 
        y_e = int(y_e + 1)
    except:
        y_e = None
    
    return x_s,x_e,y_s,y_e
    
### BUILDER FUNCTIONS
def Start_modis_hdf(project, hdfp):
    '''Start aggregating MODIS dates and tiles to one hdf5 array file.

    Defines subset array extents, and passes parameters on to Mk_hdf

    :param project: (dict) py-stint project parameters
    :param hdfp: (dict) conversion parameters
    :return: None
    '''
    x_var, y_var = Get_mosaic_extent( project, hdfp['dset'] )
    x_s,x_e,y_s,y_e = Get_spatial_indexes(x_var, y_var, 
                                          project['aoi'])
    x_sub = x_var[x_s:x_e]
    y_sub = y_var[y_s:y_e]

    modis_md = Get_modis_md( project['modis_dir'],
                             hdfp['dset'],hdfp['dnum'],
                             project['modis_tiles'][0] )
    Mk_hdf( hdfp, x_sub, y_sub, modis_md)


def Mk_hdf( hdfp, x_var, y_var, modis_md ):
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
        # arr_out[:] = modis_md['fill_value'] # todo go back in later, find these locations, and replace with fill_value
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


### AOI functions
def Get_sin_aoi( aoi ):
    '''Return sinusoidal bounding coordinates of AOI.

    :param aoi: (dict) aoi extents
    :return: (float,float,float,float) xmin,ymin,xmax,ymax
             aoi extents in sinusoidal coordinates
    '''
    # Load aoi extents
    aoi_bbox = Mk_bbox(aoi['xmin'],aoi['ymin'],aoi['xmax'],aoi['ymax'])

    # transform aoi from native prj to MODIS sinusoidal
    sin_transform = osr.CoordinateTransformation(aoi['srs'],sin_srs) #src,dst
    # if CoordinateTransformation fails, it will return null:
    if sin_transform == None:
        print '[ERROR] Could not reproject AOI box to MODIS sinusoidal'
        sys.exit( 1 )

    aoi_bbox.Transform(sin_transform)  
    
    xmin,xmax,ymin,ymax = aoi_bbox.GetEnvelope() 
    return xmin,ymin,xmax,ymax


def Get_spatial_indexes(x_var, y_var, aoi):
    '''Return start and end indices for x and y coordinates.

    indices identify the subset of x_var and y_var intersecting aoi extents

    index interpretation:
    if x_s==None:
        print 'aoi extends west of grid'
    if x_e==None:
        print 'aoi extends east of grid'
    if y_s==None:
        print 'aoi extends north of grid'
    if y_e==None:
        print 'aoi extends south of grid'


    :param x_var: (np.array) [x] cell center x coordinates, MODIS tile
    :param y_var: (np.array) [y] cell center y coordinates, MODIS tile
    :param aoi: (dict) aoi extents
    :return:(int,int,int,int) x_start,x_end,y_start,y_end
            x start & end, y start & end indices
            if, for example:
            a = [y,x] modis tile,
            a_sub = [y,x] portion of tile overlapping the aoi bounding box,
            then:
            a_sub = a[y_start:y_end,x_start:x_end]
            while the cell center corodinates of the tile subset are:
            x_sub = x_var[x_s:x_e]
            y_sub = y_var[y_s:y_e]
    '''
    aoi_xmin,aoi_ymin,aoi_xmax,aoi_ymax = Get_sin_aoi(aoi)
    
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


def Mod2hdf( project, dset, dnum ):
    '''Coordinates the aggregation of a MODIS dataset to hdf5 array file.

    If the hdf5 has been created, Mod2hdf will pick up where it left off.

    :param project: (dict) py-stint project parameters
    :param dset: (string) MODIS product, ie 'MCD43A3'
    :param dnum: (integer) subdataset number, ie 7 -> BSA_vis
    :return: None
    '''
    dname = project['modis']
    modis_days = project['modis_days']

    hdfp = {} # hdf parameters
    hdfp['appendnum']   = 5
    hdfp['dset']        = dset
    hdfp['dnum']        = dnum
    hdfp['sds']         = dname[dset][dnum]
    hdfp['modis_dir']   = project['modis_dir']
    hdfp['modis_days']  = project['modis_days']
    hdfp['modis_tiles'] = project['modis_tiles']
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

    modis_md = Get_modis_md(hdfp['modis_dir'],dset,dnum,hdfp['modis_tiles'][0])
    hdfp['fill_value'] = modis_md['fill_value']
    if not os.path.isfile(hdfp['h5f']):
        Start_modis_hdf(project, hdfp)
    
    Continue_modis_hdf( project, hdfp )
    print 'Bing!!'
    
