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
    tiles_out = os.path.join(project['shp_dir'],'ssarV2_tiles')

     # Define outfile as shapefile, with feature type polygon
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
                        modis = [(tile_xmin, tile_ymin,tile_xmax,tile_ymax), mod_params]
                        climate = [(txmin, tymin, txmax, tymax), climate_params]
                        mod_clim_isect   = Mk_mod_clim_tile(modis, climate) # todo
                        mod_clim_isect_r = 'rtree idx for mod_clim_isect'
                        # each feature in mod_clim_isect should have id, geom, modis_area, modis_ & climate_id & x_ind & y_ind
                        # mod_clim_isect_idx??

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
                            geom1 = ssarV1_feat.GetGeometryRef()
                            fxmin,fxmax,fymin,fymax = geom1.GetEnvelope() # todo

                            final_hits = mod_clim_isect_r.intersection((fxmin,fymin,fxmax,fymax))
                            for mod_clim_id in final_hits:
                                TODO = 'THIS:' # todo
                                # get intersection area
                                # write ssarV1_attributes to csv
                                # add modis_ & climate_ id, x&y_ind, modis_area to 'out_rows' lists

#! todo in writing CSV datasets, climate and modis_ids need to continue to increase
#! over the whole output, not just within each tile then resetting to 0

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



def Mk_polygrid_memory(params, mk_idx=True, record_ctr_coords=True, record_area=True, transform=None):

    xmin = float(params['xmin'])
    xmax = float(params['xmax'])
    ymin = float(params['ymin'])
    ymax = float(params['ymax'])
    dx   = params['dx']
    dy   = params['dy']
    outf = params['outf']

    boxes      = []
    x_ind_list = []
    y_ind_list = []

    if mk_idx==True:
        idx = index.Index(interleaved=True)
    if record_ctr_coords==True:
        ctr_xy = []
    if record_area==True:
        arealist = []

    id_list = 'implicit; idVar only needed for rtree'

    idVar   = 0

    y = ymax
    y_ind = 0
    while round(y,3) > round(ymin,3):
        x = xmin
        x_ind = 0
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
                 idx.insert(idVar,(x, y - dy, x + dx, y))
            if record_ctr_coords==True:
                # todo get ctr-x ctr-y out of Centroid
                ctr_xy.append((feat.GetX(), feat.GetY()))

            feat = geom = None
            idVar += 1
            x += dx
            x_ind+=1
        y -= dy
        y_ind+=1

    # Save and close everything
    ds = layer = feat = perim = polygon = None
    # sys.stdout.write("\n")
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


    out['other params'] = 'other_lists'



def Mk_mod_clim_tile(modis, climate ):

    # unpack parameters
    tile_xmin, tile_ymin, tile_xmax, tile_ymax = modis[0]
    mod_params = modis[1]

    txmin, tymin, txmax, tymax = climate[0] # These are utm33 projected coordinates of tile_xmin, ...
    climate_params = climate[1]

    climate_xmin = float(climate_params['xmin'])
    climate_xmax = float(climate_params['xmax'])
    climate_ymin = float(climate_params['ymin'])
    climate_ymax = float(climate_params['ymax'])
    climate_dx   = climate_params['dx']
    climate_dy   = climate_params['dy']

    # todo assign climate xy minmax such that climate tile encloses tile min max
    '''params interpreted like this:
    xmin = float(params['xmin'])
    xmax = float(params['xmax'])
    ymin = float(params['ymin'])
    ymax = float(params['ymax'])
    dx   = params['dx']
    dy   = params['dy']
    outf = params['outf']'''

    climate_tile_ = Mk_polygrid_memory(climate_tile_params, mk_idx=True, record_ctr_coords=True, record_area=False, transform=None)


    mod_xmin = float(mod_params['xmin'])
    mod_xmax = float(mod_params['xmax'])
    mod_ymin = float(mod_params['ymin'])
    mod_ymax = float(mod_params['ymax'])
    mod_dx   = mod_params['dx']
    mod_dy   = mod_params['dy']

    modis_tile = Mk_polygrid_memory(modis_tile_params, mk_idx=False, record_ctr_coords=True, record_area=True, transform=sin2utm33n)

    for mod_idx, mod_cell in enumerate(modis_tile['boxes'])