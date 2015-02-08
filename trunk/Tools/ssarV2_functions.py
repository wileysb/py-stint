__author__ = 'wiley'

import os
import ogr
import osr
import math
import csv
from rtree import index
from SPATIAL_tools import FastRtree, Ogr_open, Mk_proj
from ORG_tools import Countdown
from MODIS_aoi import Mk_bbox, Parse_extents

sin_prj = 'PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
sin_srs = osr.SpatialReference()
sin_srs.ImportFromWkt(sin_prj)

utm33n = osr.SpatialReference()
utm33n.ImportFromEPSG(32633)
utm33n_string = utm33n.ExportToWkt()

sin2utm33n = osr.CoordinateTransformation(sin_srs, utm33n) #src,dst

def Isect_mod_clim_ssar(project):
    '''Isect_mod_clim_ssar(project, project_paths)

    project_paths = {'climate_fn'   :climate_fn,
                     'modis_fn'     :modis_fn,
                     'ssarV1_dir'   : '/space/wib_data/LANDCOVER/ss_ar_shp/',
                     'ssarV1_tiles' : tiles_dsn
    }

    # Define blocks of 30x30 MODIS cells (cells MUST NOT repeat!!)
    # for given unique, unrepeated block (30x30) of MODIS cells:
        # which regions ('tile_id') in tile_bounds.shp intersect block
        # for those tile_id which intersect block:
            # intersect modis features, climate features, lc features
            # export lc isect to csv
            # export applicable modis and climate cells to csv'''

    # load climate bbox
    climate_params = Parse_extents(project['paths']['climate_fn'])
    climate_bbox = Mk_bbox(climate_params['xmin'], climate_params['ymin'],climate_params['xmax'], climate_params['ymax'])


    # Define blocks of 30x30 MODIS cells (cells MUST NOT repeat!!)
    mod_params = Parse_extents(project['paths']['modis_fn'])

    mod_xmin = float(mod_params['xmin'])
    mod_xmax = float(mod_params['xmax'])
    mod_ymin = float(mod_params['ymin'])
    mod_ymax = float(mod_params['ymax'])
    mod_dx   = mod_params['dx']
    mod_dy   = mod_params['dy']
    mod_prj  = mod_params['srs']

    tile_dx = mod_dx * 30
    tile_dy = mod_dy * 30


    # load ssarV1 tile bounds
    ssarV1_ds, ssarV1_lyr  = Ogr_open(project['paths']['ssarV1_tiles'])
    ssarV1_r = FastRtree(project['paths']['ssarV1_tiles'])

    # prepare tile_bounds out
    tiles_out = os.path.join(project['shp_dir'],'ssarV2_tile_bounds')

     # Define shapefile path for tile_bounds.shp, with feature type polygon
    driver = ogr.GetDriverByName('Esri Shapefile')
    tiles_out_ds = driver.CreateDataSource(tiles_out+'.shp')
    tiles_out_layer = tiles_out_ds.CreateLayer('',None,ogr.wkbPolygon)
    tiles_out_layer.CreateField(ogr.FieldDefn('id',ogr.OFTInteger))

    # Define text field
    field_name = ogr.FieldDefn("tile_name", ogr.OFTString)
    field_name.SetWidth(24)
    tiles_out_layer.CreateField(ogr.FieldDefn("tile_name", ogr.OFTString))
    # tiles_out_layer.CreateField(ogr.FieldDefn(field_name)) # yind_xind

    defn = tiles_out_layer.GetLayerDefn()

    idVar   = 0
    modis_idVar = 0
    count        = 0
    count_max    = int(math.ceil((mod_ymax-mod_ymin) / tile_dy))
    # count_update = count_max * 0.05 # print progress every 5%!
    progress_bar = Countdown(count_max)

    tile_uly = mod_ymax
    tile_y_ind = 0
    while round(tile_uly,3) > round(mod_ymin,3):
        tile_ulx = mod_xmin
        tile_x_ind = 0
        while round(tile_ulx,3) < round(mod_xmax,3):
            new_tile = True

            tile_xmin = tile_ulx
            tile_ymin = tile_uly - tile_dy
            tile_xmax = tile_ulx + tile_dx
            tile_ymax = tile_uly

            # if creating from scratch:
            # tile_y_ind = 83; tile_x_ind = 34
            #tile_ymin = tile_uly - tile_y_ind*tile_dy
            #tile_xmax = tile_ulx + tile_x_ind*tile_dx

            tile_bbox_utm33 = Mk_bbox(tile_xmin, tile_ymin, tile_xmax, tile_ymax)
            tile_bbox_utm33.Transform(sin2utm33n)

            txmin,txmax,tymin,tymax = tile_bbox_utm33.GetEnvelope()

            if tile_bbox_utm33.Intersects(climate_bbox):
                hits = ssarV1_r.intersection((txmin,tymin,txmax,tymax)) # (gxmin,gymin,gxmax,gymax)
                for hit_fid in hits:
                    if new_tile == True:
                        tile_id = '{0}_{1}'.format(tile_y_ind,tile_x_ind)

                        # check for intersection with climate
                        feat = ogr.Feature(defn)
                        feat.SetField('id',idVar)
                        feat.SetField('tile_name',tile_id)
                        feat.SetGeometry(tile_bbox_utm33)

                        tiles_out_layer.CreateFeature(feat)
                        feat =  None
                        idVar += 1

                        tile_out_fmt = os.path.join(project['csv_dir'], 'ssarV2_{0}_'+tile_id+'.csv') # .format(sds)
                        modis_rows_to_write = set()
                        climate_rows_to_write = set()

                        # generate modis features in tile
                        # - transform to utm33
                        # - get areas
                        # generate climate features intersecting tile
                        modis = [(tile_xmin, tile_ymin,tile_xmax,tile_ymax), mod_params, modis_idVar]
                        climate = project['paths']['climate_shp']# [(txmin, tymin, txmax, tymax), climate_params]
                        mod_clim_isect, modis_idVar   = Mk_mod_clim_tile(modis, climate, tile_x_ind, tile_y_ind, tile_dx, tile_dy)
                        mod_clim_isect_r = mod_clim_isect['idx']

                        new_tile = False


                    ssarV1_tile = ssarV1_lyr.GetFeature(hit_fid)
                    geom2 = ssarV1_tile.GetGeometryRef()
                    if tile_bbox_utm33.Intersects(geom2):
                        ssarV1_tile_id = ssarV1_tile.GetField('tile_id')

                        # load tile'
                        ssarV1_tile_dsn = os.path.join(project['paths']['ssarV1_dir'], 'ss_ar_'+ssarV1_tile_id)
                        ssarV1_tile_ds, ssarV1_tile_lyr = Ogr_open(ssarV1_tile_dsn)
                        for fid1 in range(0,ssarV1_tile_lyr.GetFeatureCount()):
                            ssarV1_feat = ssarV1_tile_lyr.GetFeature(fid1)
                            ssar_geom = ssarV1_feat.GetGeometryRef()
                            fxmin,fxmax,fymin,fymax = ssar_geom.GetEnvelope()

                            final_hits = mod_clim_isect_r.intersection((fxmin,fymin,fxmax,fymax))
                            for mod_clim_id in final_hits:
                                mod_clim_feat = mod_clim_isect['geom'][mod_clim_id]
                                if ssar_geom.Intersects(mod_clim_feat):
                                    isect = ssar_geom.Intersection()
                                    isect_area = isect.GetArea()

                                    climate_id = mod_clim_isect['climate_id'][mod_clim_id]
                                    climate_x_ind = mod_clim_isect['climate_x_ind'][mod_clim_id]
                                    climate_y_ind = mod_clim_isect['climate_y_ind'][mod_clim_id]
                                    modis_id = mod_clim_isect['modis_id'][mod_clim_id]
                                    modis_area = mod_clim_isect['modis_area'][mod_clim_id]
                                    modis_ctr_x = mod_clim_isect['modis_ctr_x'][mod_clim_id]
                                    modis_ctr_y = mod_clim_isect['modis_ctr_y'][mod_clim_id]
                                    modis_x_ind = mod_clim_isect['modis_x_ind'][mod_clim_id]
                                    modis_y_ind = mod_clim_isect['modis_y_ind'][mod_clim_id]

                                    climate_rows_to_write.add((climate_id, climate_x_ind, climate_y_ind))
                                    modis_rows_to_write.add((modis_id, modis_area, modis_ctr_x, modis_ctr_y, modis_x_ind, modis_y_ind))

                                    # todo
                                    # write ssarV1_attributes (plus area, modis_area, mod_id, climate_id) to csv


            # todo:
            Write_modis_tile(project, modis_rows_to_write, tile_id)
            Write_climate_tile(project, climate_rows_to_write, tile_id)
            # Write modis and climate datasets to CSVs, for all cells which had hits
            tile_ulx += tile_dx
            tile_x_ind+=1
        tile_uly -= tile_dy
        tile_y_ind+=1
        count += 1
        progress_bar.check(count)

    # Save and close everything
    ds = layer = feat = perim = polygon = None
    progress_bar.flush()
    prj = project['srs'].ExportToWkt()
    Mk_proj( utm33n_string,tiles_out )



def Mk_polygrid_memory(params, tile_ulx_ind, tile_uly_ind, mk_idx=True, record_ctr_coords=True, record_area=True, transform=None, idVar=None):

    xmin = float(params['xmin'])
    xmax = float(params['xmax'])
    ymin = float(params['ymin'])
    ymax = float(params['ymax'])
    dx   = params['dx']
    dy   = params['dy']

    boxes      = []
    x_ind_list = []
    y_ind_list = []

    if mk_idx==True:
        idx = index.Index(interleaved=True)
    if record_ctr_coords==True:
        ctr_xy = []
    if record_area==True:
        arealist = []
    if idVar != None:
        idVar_list = []

    grid_idVar   = 0

    y = ymax
    y_ind = tile_uly_ind
    while round(y,3) > round(ymin,3):
        x = xmin
        x_ind = tile_ulx_ind
        while round(x,3) < round(xmax,3):
            feat = Mk_bbox(x, y - dy, x + dx, y)

            x_ind_list.append(x_ind)
            y_ind_list.append(y_ind)

            if transform!=None:
                feat.Transform(transform)
            boxes.append(feat)
            if record_area==True:
                arealist.append(feat.GetArea())
            if mk_idx==True:
                 idx.insert(grid_idVar,(x, y - dy, x + dx, y))
            if record_ctr_coords==True:
                ctr_xy.append((feat.GetX(), feat.GetY()))
            if idVar != None:
                idVar_list.append(idVar)

            feat = geom = None
            grid_idVar += 1
            if idVar != None:
                idVar+=1
            x += dx
            x_ind+=1
        y -= dy
        y_ind+=1

    # Save and close everything
    ds = layer = feat = perim = polygon = None
    out = {}
    out['boxes'] = boxes
    out['x_ind'] = x_ind_list
    out['y_ind'] = y_ind_list
    if mk_idx==True:
        out['idx']    = idx
    if record_ctr_coords==True:
        out['ctr_xy'] = ctr_xy
    if record_area==True:
        out['area'] = arealist
    if idVar != None:
        out['idVar'] = idVar_list


    return out


def Mk_mod_clim_tile(modis, climate_dsn, tile_x_ind, tile_y_ind, tile_dx, tile_dy ):
    # corner full-grid indices (modis)
    mod_tile_ulx_ind = tile_x_ind*tile_dx
    mod_tile_uly_ind = tile_y_ind*tile_dy

    # unpack parameters
    tile_xmin, tile_ymin, tile_xmax, tile_ymax = modis[0]
    mod_params = modis[1]
    modis_idVar = modis[2]

    modis_tile_params = {
    'xmin' : float(tile_xmin),
    'xmax' : float(tile_xmax),
    'ymin' : float(tile_ymin),
    'ymax' : float(tile_ymax),
    'dx'   : mod_params['dx'],
    'dy'   : mod_params['dy']
    }

    modis_tile  = Mk_polygrid_memory(modis_tile_params, mod_tile_ulx_ind, mod_tile_uly_ind, mk_idx=False, record_ctr_coords=True, record_area=True, transform=sin2utm33n, idVar=modis_idVar)
    modis_idVar = modis_tile['idVar'][-1]

    ## CLIMATE PARAMETERS - JUST USE CLIMATE_DSN WITH FASTRTREE IDX??? ###
    climate_ds, climate_lyr = Ogr_open(climate_dsn)
    climate_idx             = FastRtree(climate_dsn)

    # initiate lists for everything which will go into isect_out{}
    idx = index.Index(interleaved=True)
    isect_geom    = []
    climate_id    = []
    climate_x_ind = []
    climate_y_ind = []
    modis_id      = []
    modis_area    = []
    modis_x_ind   = []
    modis_y_ind   = []
    modis_ctr_x   = []
    modis_ctr_y   = []

    tile_idVar = 0

    for mod_i, mod_cell in enumerate(modis_tile['boxes']):
        mod_id = modis_tile['idVar'][mod_i]
        mod_area = mod_cell.GetArea()
        mod_x_ind = modis_tile['x_ind'][mod_i]
        mod_y_ind = modis_tile['y_ind'][mod_i]
        mod_ctr_x = modis_tile['ctr_xy'][mod_i][0]
        mod_ctr_y = modis_tile['ctr_xy'][mod_i][1]
        mxmin,mxmax,mymin,mymax = mod_cell.GetEnvelope()
        climate_in_modis = climate_idx.intersection((mxmin,mymin,mxmax,mymax))
        for climate_fid in climate_in_modis:
            climate_feat = climate_lyr.GetFeature(climate_fid)
            climate_cell = climate_feat.GetGeometryRef()
            if mod_cell.Intersects(climate_cell):
                isect = mod_cell.Intersection(climate_cell)
                # Insert to index
                xmin,xmax,ymin,ymax = isect.GetEnvelope()
                idx.insert(tile_idVar,(xmin,ymin,xmax,ymax))

                climate_id.append(climate_feat.GetField('id'))
                climate_x_ind.append(climate_feat.GetField('x_ind'))
                climate_y_ind.append(climate_feat.GetField('y_ind'))

                modis_id.append(mod_id)
                modis_area.append(mod_area)
                modis_x_ind.append(mod_x_ind)
                modis_y_ind.append(mod_y_ind)
                modis_ctr_x.append(mod_ctr_x)
                modis_ctr_y.append(mod_ctr_y)

                isect_geom.append(mod_cell.Intersection)
                tile_idVar+=1


    isect_out = {
    'idx'           : idx,
    'geom'          : isect_geom,
    'climate_id'    : climate_id,
    'climate_x_ind' : climate_x_ind,
    'climate_y_ind' : climate_y_ind,
    'modis_id'      : modis_id,
    'modis_area'    : modis_area,
    'modis_x_ind'   : modis_x_ind,
    'modis_y_ind'   : modis_y_ind,
    'modis_ctr_x'   : modis_ctr_x,
    'modis_ctr_y'   : modis_ctr_y}

    return isect_out, modis_idVar


def Write_modis_tile(project, modis_rows_to_write, tile_id):
    # todo
    todo = 'fill in the structure with code'
    # Sort rows by id, ascending

    for dset in project['modis'].keys():
        for dnum in project['modis'][dset].keys():
            # Create csv

            # Write header

            # Open hdf5

            # write rows:
            for row in range(len(modis_rows_to_write)):
                modis_series = hdf[:,y,x]
                # write modis_rows_elements, followed by modis_series


def Write_climate_tile(project, climate_rows_to_write, tile_id):
    # todo
    todo = 'fill in the structure with code'
    # Sort rows by id, ascending
    climate_ds_list = ['fsw','rr','sd','swe','tam']
    for climate_ds in climate_ds_list:
        # Create csv

        # Write header

        # Open hdf5

        # write rows:
        for row in range(len(climate_rows_to_write)):
            climate_series = hdf[:,y,x]
            # write climate_rows_elements, followed by climate_series