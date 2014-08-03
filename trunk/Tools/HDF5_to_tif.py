#  HDF5 to GEOTIFF conversion for climate and MODIS collections 
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

import os,sys,ogr,osr,gdal,h5py
import numpy as np
from ORG_tools import Daynum2date, Parse_input

input_fn = '../INPUT_dev.txt'

def Ld_src_array(h5f,i,sds):
    '''
    mod_scene = Ld_src_array(h5f,i)
    mod_scene.keys: 'a'   : array
                    'x'   : x axis coords
                    'y'   : y axis coords
                    'day' : modis start date
                    'prj' : projection
                    'fill': NODATA value
    '''
    with h5py.File(h5f,'r') as hdf:
        out = {}
        out['a']   = hdf[sds][i,:,:]
        out['x']   = hdf['x'][:]
        out['y']   = hdf['y'][:]
        out['prj'] = hdf[sds].attrs['projection']
        out['scale'] = float(hdf[sds].attrs['scale_factor'])
        out['offset'] = float(hdf[sds].attrs['add_offset'])
        out['fill']= int(hdf[sds].attrs['fill_value'])

        basedatestr = hdf['time'].attrs['basedate']
        daynum   = int(hdf['time'][i])
        out['day'] = Daynum2date(daynum,basedatestr)
        
    return out


def Write_tif(scene,outf):
    format = "GTiff"
    driver = gdal.GetDriverByName( format )
    # Apply Scale and Offset
    a = np.where( scene['a'] == scene['fill'], np.nan, scene['a'] )
    a = a * scene['scale'] + scene['offset']
    a = np.where( a == np.nan, scene['fill'], a )
    # Get array parameters
    numRows,numCols = scene['a'].shape
    minx            = np.min(scene['x'])
    maxy            = np.max(scene['y'])
    dx =  scene['x'][1] - scene['x'][0]
    dy =  scene['y'][1] - scene['y'][0]

    # Make outfile
    dst_ds = driver.Create( outf, numCols, numRows, 1, gdal.GDT_Int16 )
    #! Change if MCD43A2 or other byte word dataset?
    #! dst_ds = driver.Create( outf, numCols, numRows, 1, gdal.GDT_Byte )

    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    dst_ds.SetGeoTransform( [ minx, dx, 0, maxy, 0, dy ] )

    dst_ds.SetProjection( scene['prj'] )

    # write the band
    band = dst_ds.GetRasterBand(1)
    band.SetNoDataValue(int(scene['fill']))
    band.WriteArray( a )
    del dst_ds


if __name__ == '__main__':
    '''
    CLI arguments: 
    1) path/to/sds.hdf5
    2) timestep_index
    $ python dset2polygrid.py BSA_vis.hdf5 0
    '''
    project = Parse_input(input_fn)
    
    if len( sys.argv ) != 3: # first argument is 'ras2polygrid.py'
        print "[ ERROR ] you must supply two arguments: <dataset> <timestep>"
        sys.exit( 1 )
    
    processing_dir = os.path.join(project['prj_directory'],'Processing/')
    hdf_dir = os.path.join(processing_dir,'HDF/')
    tif_dir = os.path.join(processing_dir,'TIF/')
    
    sds     = sys.argv[1]
    h5f     = os.path.join(hdf_dir,project['prj_name']+'_'+sds+'.hdf5')
    i       = int(sys.argv[2])

    scene   = Ld_src_array(h5f,i,sds)
    outfn   = sds+scene['day']+'.tif'
    outf    = os.path.join(tif_dir,outfn)
    
    Write_tif(scene,outf)
