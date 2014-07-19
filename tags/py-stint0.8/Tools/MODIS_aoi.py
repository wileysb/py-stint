#  Tilefinder for MODIS datasets within a bounding box
#   part of MODIS/climate/landcover intersection toolset
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

import osr,ogr,os,sys
import numpy as np

from SPATIAL_tools import Parse_extents,Mk_bbox
# Use this tool to find modis tiles within a project's bounding box AOI
# load bounding box - from args, from shp, or from rasterS 

### DEFINE modis shapefile
modis_tile_fn = '../Data/MODIS_tiles/modis_sinusoidal_grid_world.shp'


def Check_mod_tiles(xmin,ymin,xmax,ymax,dx,dy,srs):
    mod_grid  = ogr.Open(modis_tile_fn)
    mod_layer = mod_grid.GetLayer()
    sin_srs   = mod_layer.GetSpatialRef()

    aoi_box = Mk_bbox(xmin,ymin,xmax,ymax)

    # transform aoi from native prj to MODIS sinusoidal
    sin_transform = osr.CoordinateTransformation(srs,sin_srs) #src_srs,dst_srs
    # if CoordinateTransformation fails, it will return null:
    if sin_transform == None:
        print '[ERROR] Could not reproject AOI box to MODIS sinusoidal'
        sys.exit( 1 )

    aoi_box.Transform(sin_transform)  

    # loop through modis tiles to see which intersect aoi 
    tiles_in_aoi = []
    for tile in mod_layer:
        if aoi_box.Intersects(tile.GetGeometryRef()):
            h = tile.GetFieldAsInteger("h")
            v = tile.GetFieldAsInteger("v")
            tile_str = "h%02dv%02d" % (h,v)
            tiles_in_aoi.append(tile_str)

    return tiles_in_aoi
    

if __name__ == '__main__':
    '''
    $ python MODIS_aoi.py aoi.shp/.tif/.nc/.hdf
    OR
    $ python MODIS_aoi.py xmin ymin xmax ymax prj_EPSG_code
    '''
    if len( sys.argv ) == 2:
        src_fn = sys.argv[1]
        aoi = Parse_extents(src_fn)
    elif len (sys.argv ) == 6:
        aoi = {}
        getsrs = osr.SpatialReference()
        getsrs.ImportFromEPSG( sys.argv[5] )
        aoi['srs']  = getsrs
        #! exit with error if these are not int or float
        #! and xmax>xmin, ymax>ymin:
        aoi['xmin'] = sys.argv[1]
        aoi['ymin'] = sys.argv[2]
        aoi['xmax'] = sys.argv[3]
        aoi['ymax'] = sys.argv[4]
        aoi['dx']   = 'user input'
        aoi['dy']   = 'user input'
    else:
        print '[ ERROR ] Input extent file, or "xmin ymin xmax ymax epsg_number"'
        sys.exit( 1 )
    tiles_in_aoi = Check_mod_tiles(**aoi)
    if len(tiles_in_aoi)==1:
        print "1 MODIS tiles encompassing area of interest:"
    elif len(tiles_in_aoi)==0:
        print '[ ERROR ] No intersections. Check AOI geometry and projections'
    else:
        print "%s MODIS tiles encompassing area of interest:" % len(tiles_in_aoi)
    for tile in tiles_in_aoi:
        print tile
