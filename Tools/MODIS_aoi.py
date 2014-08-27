#  Tilefinder for MODIS datasets within a bounding box
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
'''Module or CLI script for finding MODIS tiles intersecting a bounding box.

Usage:
$ python MODIS_aoi [raster|shp]
 OR
$ python MODIS_aoi <xmin> <ymin> <xmax> <ymax> <epsg_number>
'''
import osr,ogr,sys,os

from SPATIAL_tools import Parse_extents,Mk_bbox


def Check_mod_tiles(xmin,ymin,xmax,ymax,dx,dy,srs,modis_tile_fn=None):
    '''Returns list of MODIS tiles (strings, 'h__v__') intersecting bounding
    box or raster/vector data source supplied.

    :param xmin: (float or int) bounding box minimum x coordinate
    :param ymin: (float or int) bounding box minimum y coordinate
    :param xmax: (float or int) bounding box maximum x coordinate
    :param ymax: (float or int) bounding box maximum y coordinate
    :param dx: not used (dict unpacking placeholder)
    :param dy: not used (dict unpacking placeholder)
    :param srs: osr.SpatialReference().ImportFromEPSG('EPSG_num')
    :param modis_tile_fn: (string) path to shapefile with modis tile boundaries
    :return: (list) e.g. ['h19v02','h19v03',...]
    '''
    if modis_tile_fn==None:
        tools_dir = os.path.split(__file__)[0]
        modis_tile_fn = os.path.join(tools_dir,'Data','MODIS_tiles.shp')


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


def Args_are_floats(*args):
    '''Return True if all args can be converted to float.
    Else, return False.
    '''
    for arg in args:
        try:
            float(arg)
        except ValueError:
            return False
    return True


def Valid_bounds(xmin,ymin,xmax,ymax):
    '''Return False unless xmax>xmin & ymax>ymin.'''
    if not float(xmax)>float(xmin):
        return False
    if not float(ymax)>float(ymin):
        return False
    return True


if __name__ == '__main__':
    '''
    $ python MODIS_aoi.py aoi.shp/.tif/.nc/.hdf
    OR
    $ python MODIS_aoi.py xmin ymin xmax ymax prj_EPSG_code
    '''
    if len( sys.argv ) == 2:
        src_fn = sys.argv[1]
        if not os.path.isfile(src_fn):
            print '[ERROR] File not found:',src_fn
            sys.exit(1)
        try:
            aoi = Parse_extents(src_fn)
        except:
            print '[ERROR] Could not parse bounding box from file:',src_fn
            sys.exit(1)
    elif len (sys.argv ) == 6:
        aoi = {}
        try:
            getsrs = osr.SpatialReference()
            getsrs.ImportFromEPSG( int(sys.argv[5]) )
            aoi['srs']  = getsrs
        except:
            print '[ERROR] EPSG not found: ',sys.argv[5]
            sys.exit(1)
        if not Args_are_floats(*sys.argv[1:5]):
            print '[ERROR] bounding box values are not all numbers'
            sys.exit(1)
        if not Valid_bounds(*sys.argv[1:5]):
            print '[ERROR] bounding box maximums should be greater than minimums'
            sys.exit(1)
        aoi['xmin'] = float(sys.argv[1])
        aoi['ymin'] = float(sys.argv[2])
        aoi['xmax'] = float(sys.argv[3])
        aoi['ymax'] = float(sys.argv[4])
        aoi['dx']   = 'user input'
        aoi['dy']   = 'user input'

    else:
        print "[ ERROR ] bad arguments"
        print 'Input extent file, or "<xmin> <ymin> <xmax> <ymax> <epsg_number>"'
        sys.exit( 1 )

    # If no errors in parsing bbox extents, continue:
    tiles_in_aoi = Check_mod_tiles(**aoi)
    if len(tiles_in_aoi)==1:
        print "1 MODIS tiles encompassing area of interest:"
    elif len(tiles_in_aoi)==0:
        print '[ ERROR ] No intersections. Check AOI geometry and projections'
    else:
        print "%s MODIS tiles encompassing area of interest:" % len(tiles_in_aoi)
    for tile in tiles_in_aoi:
        print tile
