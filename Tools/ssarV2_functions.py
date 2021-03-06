__author__ = 'wiley'

import os
import ogr
import osr
import math
import csv
import h5py
import glob
import numpy as np
import datetime as dt
from rtree import index
from SPATIAL_tools import FastRtree, Ogr_open, Mk_proj, Mk_bbox, Parse_extents
from ORG_tools import Countdown
import zipfile
try:
    import zlib
    compression = zipfile.ZIP_DEFLATED
except:
    compression = zipfile.ZIP_STORED

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

    # get ssarV1 attributes
    ssar_attribs = Get_lc_attribs(project)
    ssar_hdr = ['area','modis_area','modis_id','climate_id']
    for attrib in ssar_attribs:
        ssar_hdr.append(attrib)

    # Make climate csv header
    start_date = dt.datetime.strptime(str(project['modis_days'][0]), '%Y%j').date()
    end_date   = dt.datetime.strptime(str(project['modis_days'][-1]), '%Y%j').date()
    numdays = (end_date-start_date).days + 1
    daterange = [(start_date + dt.timedelta(days=x)).strftime('%Y%j') for x in range(0, numdays)]
    project['climate_hdr'] = ['climate_id', 'x_ind', 'y_ind'] + daterange

    # load climate bbox
    climate_params = Parse_extents(project['paths']['climate_fn'])
    climate_bbox = Mk_bbox(float(climate_params['xmin']), float(climate_params['ymin']), float(climate_params['xmax']), float(climate_params['ymax']))


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

     # Define shapefile path for tile_bounds.shp, with feature type polygon
    if 'restart' in project.keys():
        modis_idVar = project['restart']['modis_idVar']
        modis_nanmask = project['restart']['modis_nanmask']

    else:
        # Make a mask for MODIS cells which are worth processing (contain any non-nan values from a cloud-free image)
        modis_nanmask = Mk_modis_nanmask(project)
        modis_idVar = 0


    idVar = 0
    count_max    = modis_nanmask.sum()
    progress_bar = Countdown(count_max, update_interval=.01)

    tile_uly = mod_ymax
    tile_y_ind = 0
    while round(tile_uly,3) > round(mod_ymin,3):
        tile_ulx = mod_xmin
        tile_x_ind = 0
        while round(tile_ulx,3) < round(mod_xmax,3):
            if modis_nanmask[tile_y_ind,tile_x_ind]:
                idVar += 1
                new_tile = True
                new_ssar_csv = True

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
                modis_rows_to_write = [] # no list if no intersections
                climate_rows_to_write = [] # no list if no intersections
                if tile_bbox_utm33.Intersects(climate_bbox):
                    hits = ssarV1_r.intersection((txmin,tymin,txmax,tymax)) # (gxmin,gymin,gxmax,gymax)
                    for hit_fid in hits:
                        if new_tile == True:
                            tile_id = '{0}_{1}'.format(tile_y_ind,tile_x_ind)

                            tile_out_fmt = os.path.join(project['csv_dir'], '{0}/ssarV2_{0}_'+tile_id+'.csv') # .format(sds)
                            modis_rows_to_write = set()
                            climate_rows_to_write = set()

                            # generate modis features in tile
                            # - transform to utm33
                            # - get areas
                            # generate climate features intersecting tile
                            modis = [(tile_xmin, tile_ymin,tile_xmax,tile_ymax), mod_params, modis_idVar]
                            climate = project['paths']['climate_dsn']# [(txmin, tymin, txmax, tymax), climate_params]
                            mod_clim_isect, modis_idVar   = Mk_mod_clim_tile(modis, climate, tile_x_ind, tile_y_ind)
                            mod_clim_isect_r = mod_clim_isect['idx']



                            new_tile = False


                        ssarV1_tile = ssarV1_lyr.GetFeature(hit_fid)
                        geom2 = ssarV1_tile.GetGeometryRef()
                        if tile_bbox_utm33.Intersects(geom2):
                            ssarV1_tile_id = ssarV1_tile.GetField('tile_id')

                            # load tile'
                            ssarV1_tile_dsn = os.path.join(project['paths']['ssarV1_dir'], 'ss_ar_'+str(ssarV1_tile_id))
                            if os.path.isfile(ssarV1_tile_dsn+'.shp'):
                                ssarV1_tile_ds, ssarV1_tile_lyr = Ogr_open(ssarV1_tile_dsn)
                                for fid1 in range(0,ssarV1_tile_lyr.GetFeatureCount()):
                                    ssarV1_feat = ssarV1_tile_lyr.GetFeature(fid1)
                                    ssar_geom = ssarV1_feat.GetGeometryRef()
                                    fxmin,fxmax,fymin,fymax = ssar_geom.GetEnvelope()

                                    # get attribute fields
                                    ssar_row_proto = [ssarV1_feat.GetField(attrib) for attrib in ssar_attribs]

                                    final_hits = mod_clim_isect_r.intersection((fxmin,fymin,fxmax,fymax))
                                    for mod_clim_id in final_hits:
                                        mod_clim_feat = mod_clim_isect['geom'][mod_clim_id]
                                        if ssar_geom.Intersects(mod_clim_feat):
                                            isect = ssar_geom.Intersection(mod_clim_feat)
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

                                            if new_ssar_csv==True:
                                                # open ssar csv, write header
                                                ssar_fn  = tile_out_fmt.format('lc')
                                                ssar_f   = open(ssar_fn,'wt')
                                                ssar_csv = csv.writer(ssar_f)
                                                ssar_csv.writerow(ssar_hdr)
                                                new_ssar_csv = False

                                            ssar_row = [isect_area,modis_area,modis_id,climate_id] + ssar_row_proto
                                            ssar_csv.writerow(ssar_row)
                            else:
                                print ssarV1_tile_dsn, 'absent; moving on'

                if new_ssar_csv==False:
                    ssar_f.close()
                # Write modis and climate datasets to CSVs, for all cells which had hits
                if len(modis_rows_to_write)>0:
                    Write_modis_tile(project, modis_rows_to_write, tile_out_fmt)
                if len(climate_rows_to_write)>0:
                    Write_climate_tile(project, climate_rows_to_write, tile_out_fmt)

            tile_ulx += tile_dx
            tile_x_ind+=1
        tile_uly -= tile_dy
        tile_y_ind+=1
        progress_bar.check(idVar)

    progress_bar.flush()


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
                ctr = feat.Centroid()
                ctr_xy.append((ctr.GetX(), ctr.GetY()))
                del ctr
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


def Mk_mod_clim_tile(modis, climate_dsn, tile_x_ind, tile_y_ind):
    # corner full-grid indices (modis)
    mod_tile_ulx_ind = tile_x_ind*30
    mod_tile_uly_ind = tile_y_ind*30

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

                isect_geom.append(isect)
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

    del climate_ds, climate_lyr, climate_feat, climate_cell

    return isect_out, modis_idVar


def Write_modis_tile(project, modis_rows_to_write, tile_out_fmt):
    # Sort rows by id, ascending
    modis_rows = sorted(modis_rows_to_write, key=lambda rows: rows[0])

    for dset in project['modis'].keys():
        for dnum in project['modis'][dset].keys():
            sds = project['modis'][dset][dnum]
            hdf_fn = os.path.join(project['hdf_dir'],project['prj_name']+ \
                               '_'+sds+'.hdf5')

            # Create csv
            csv_fn = tile_out_fmt.format(sds)
            csv_f   = open(csv_fn,'wt')
            mod_csv = csv.writer(csv_f)

            # Write header
            hdr = ['modis_id', 'x_ind', 'y_ind', 'ctr_x', 'ctr_y', 'area'] + project['modis_days']
            mod_csv.writerow(hdr)

            # Open hdf5
            hdf = h5py.File(hdf_fn, 'r')

            # Get offset, scale, nan
            scale  = hdf[sds].attrs['scale_factor']
            offset = hdf[sds].attrs['add_offset']
            ds_nan = hdf[sds].attrs['fill_value']

            # write rows:
            for row in range(len(modis_rows)):
                modis_id, area, ctr_x, ctr_y, x_ind, y_ind = modis_rows[row]
                modis_series = hdf[sds][:,y_ind,x_ind]
                modis_row =  [modis_id, x_ind, y_ind, ctr_x, ctr_y, area]  + list(np.where(modis_series==ds_nan,np.nan,modis_series)*scale+offset)
                mod_csv.writerow(modis_row)

            hdf.close()
            csv_f.close()


def Write_climate_tile(project, climate_rows_to_write, tile_out_fmt):
    # Sort rows by id, ascending
    climate_rows = sorted(climate_rows_to_write, key=lambda rows: rows[0])

    climate_ds_list = ['fsw','rr','sd','swe','tam']
    for sds in climate_ds_list:
        # Create csv
        csv_fn = tile_out_fmt.format(sds)
        csv_f   = open(csv_fn,'wt')
        climate_csv = csv.writer(csv_f)

        # Write header
        hdr = project['climate_hdr']
        climate_csv.writerow(hdr)

        # Open hdf5
        hdf_fn = os.path.join(project['climate_dir'],project['prj_name']+ \
                               '_'+sds+'.hdf5')
        hdf = h5py.File(hdf_fn, 'r')

        # Get offset, scale, nan
        scale  = hdf[sds].attrs['scale_factor']
        offset = hdf[sds].attrs['add_offset']
        ds_nan = hdf[sds].attrs['fill_value']

        # write rows:
        for row in range(len(climate_rows)):
            climate_id, x_ind, y_ind = climate_rows[row]
            climate_series = hdf[sds][:,y_ind,x_ind]
            climate_row = [climate_id, x_ind, y_ind]  + list(np.where(climate_series==ds_nan,np.nan,climate_series)*scale+offset)
            climate_csv.writerow(climate_row)

        hdf.close()
        csv_f.close()


def Get_lc_attribs(project):
    lc_dir = project['paths']['ssarV1_dir']

    proto_dsn = glob.glob(os.path.join(lc_dir, '*.shp'))[0].split('.')[0]

    proto_ds, proto_lyr = Ogr_open(proto_dsn)

    layerDefinition = proto_lyr.GetLayerDefn()

    attribs = []
    for i in range(layerDefinition.GetFieldCount()):
        attribs.append(layerDefinition.GetFieldDefn(i).GetName())

    del layerDefinition, proto_lyr, proto_ds, proto_dsn

    attribs_to_omit = ['cell_id', 'feature_id', 'ctr_x_cell', 'ctr_y_cell', 'ctr_x_poly', 'ctr_y_poly', 'area_m2']
    for attr in attribs_to_omit:
        attribs.remove(attr)
    return attribs


def Gen_ssarV2_tiles(project, out_dsn):

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

    # prepare tile_bounds out
    tiles_out = os.path.join(project['shp_dir'],'ssarV2_tile_bounds')
    Mk_proj( utm33n_string,out_dsn )

     # Define shapefile path for tile_bounds.shp, with feature type polygon
    driver = ogr.GetDriverByName('Esri Shapefile')
    tiles_out_ds = driver.CreateDataSource(out_dsn+'.shp')
    tiles_out_layer = tiles_out_ds.CreateLayer('',None,ogr.wkbPolygon)
    tiles_out_layer.CreateField(ogr.FieldDefn('id',ogr.OFTInteger))


    # Define text field
    field_name = ogr.FieldDefn("tile_name", ogr.OFTString)
    field_name.SetWidth(24)
    tiles_out_layer.CreateField(ogr.FieldDefn("tile_name", ogr.OFTString))

    defn = tiles_out_layer.GetLayerDefn()

    idVar   = 0
    count        = 0
    count_max    = int(math.ceil((mod_ymax-mod_ymin) / tile_dy))
    progress_bar = Countdown(count_max, update_interval=.01)

    tile_uly = mod_ymax
    tile_y_ind = 0
    while round(tile_uly,3) > round(mod_ymin,3):
        tile_ulx = mod_xmin
        tile_x_ind = 0
        while round(tile_ulx,3) < round(mod_xmax,3):
            tile_xmin = tile_ulx
            tile_ymin = tile_uly - tile_dy
            tile_xmax = tile_ulx + tile_dx
            tile_ymax = tile_uly

            tile_bbox_utm33 = Mk_bbox(tile_xmin, tile_ymin, tile_xmax, tile_ymax)
            tile_bbox_utm33.Transform(sin2utm33n)

            tile_id = '{0}_{1}'.format(tile_y_ind,tile_x_ind)

            feat = ogr.Feature(defn)
            feat.SetField('id',idVar)
            feat.SetField('tile_name',tile_id)
            feat.SetGeometry(tile_bbox_utm33)

            tiles_out_layer.CreateFeature(feat)
            feat =  None
            idVar += 1


            tile_ulx += tile_dx
            tile_x_ind+=1
        tile_uly -= tile_dy
        tile_y_ind+=1
        count += 1
        progress_bar.check(count)

    # Save and close everything
    ds = layer = feat = perim = polygon = None
    progress_bar.flush()


def Mk_modis_nanmask(project):
    sds = 'BSA_nir'
    fn = hdf_fn = os.path.join(project['hdf_dir'], project['prj_name']+'_'+sds+'.hdf5')
    hdf = h5py.File(fn,'r')
    nanval = hdf[sds].attrs['fill_value']
    i = 18

    scene = hdf[sds][i,:,:]
    scene = np.where(scene==nanval,np.nan,scene)

    ylen, xlen = scene.shape

    mask = np.ones((110, 109))

    tx = ty = 30

    for yi in range(110):
        for xi in range(109):
            sx = tx * xi
            sy = ty * yi

            ex = tx * (xi+1)
            ey = ty * (yi+1)

            if ex>xlen:
                ex=xlen
            if ey>ylen:
                ey=ylen

            num_not_nan = np.count_nonzero(~np.isnan(scene[sy:ey,sx:ex]))
            if num_not_nan==0:
                mask[yi,xi]=0

    return mask


def Check_tile(tile_id, csv_fmt='/home/wiley/wrk/ntnu/ssarV2/CSV/{0}/ssarV2_{0}_{1}.csv', tolerance=1, return_details=False):

    lc_fn = csv_fmt.format('lc',tile_id)

    dtypes = {'names': ('feat_area', 'mod_area', 'mod_id', 'climate_id'),
             'formats': ('float',    'float',    'int',    'int')}

    lc_tile = np.loadtxt(lc_fn, delimiter=',', skiprows=1, usecols=[0,1,2,3], dtype=dtypes)

    mod_ids = np.unique(lc_tile['mod_id'])
    climate_ids = np.unique(lc_tile['climate_id'])

    ds_diff   = {}
    area_diff = []
    if return_details==True:
        details = {}
    for mod_id in mod_ids:
        where_id = np.where(lc_tile['mod_id']==mod_id)
        matches = lc_tile[where_id]
        feat_sum = np.sum(matches['feat_area'])
        mod_area = matches['mod_area'][0]
        area_diff.append(mod_area - feat_sum)
        if return_details==True:
            if np.abs(mod_area - feat_sum)>tolerance:
                details[mod_id] = 100*np.abs(mod_area - feat_sum)/mod_area # %!


    ds_diff['area'] = np.max(np.abs(np.array(area_diff)))<tolerance
    if not ds_diff['area']:
        ds_diff['count'] = (np.abs(np.array(area_diff))>tolerance).sum() # count vals outside tolerance



    modis_datasets   =  ['BSA_ancill', 'BSA_nir', 'BSA_sw', 'BSA_band', 'BSA_quality', 'BSA_vis']
    climate_datasets =  ['fsw', 'sd', 'tam', 'rr', 'swe',]

    for ds in climate_datasets:
        ds_csv = csv_fmt.format(ds, tile_id)
        climate_id = np.loadtxt(ds_csv, delimiter=',', skiprows=1, usecols=[0,])
        ds_diff[ds] = (climate_id==climate_ids).all()

    for ds in modis_datasets:
        ds_csv = csv_fmt.format(ds, tile_id)
        modis_id = np.loadtxt(ds_csv, delimiter=',', skiprows=1, usecols=[0,])
        ds_diff[ds] = (modis_id==mod_ids).all()

    if return_details==True:
        return details
    else:
        return ds_diff, len(mod_ids)


def Check_output(csv_dir, tolerance=1):
    '''missing_files, data_mismatch = Check_output(csv_dir)

    '''
    missing_files = []
    data_mismatch = {}

    csv_format = os.path.join(csv_dir, '{0}/ssarV2_{0}_{1}.csv')

    tiles = glob.glob(os.path.join(csv_dir, 'lc/*.csv'))
    tile_ids = ['_'.join(tile.split('.')[0].split('_')[-2:]) for tile in tiles]
    rows_cols = sorted([tile_id.split('_') for tile_id in tile_ids], key=lambda rows: rows[0])
    tile_ids = ['_'.join(row_col) for row_col in rows_cols]

    dsets = ['BSA_ancill', 'BSA_nir', 'BSA_sw', 'BSA_band', 'BSA_quality', 'BSA_vis', 'fsw', 'sd', 'tam', 'rr', 'swe']

    cells_missing_features = 0
    total_modis_cells      = 0

    for tile_id in tile_ids:
        all_files = True
        for dset in dsets:
            all_files = all_files * os.path.isfile(csv_format.format(dset, tile_id))
        if all_files:
            ds_diff, num_mod_cells = Check_tile(tile_id, csv_fmt=csv_format, tolerance=tolerance)
            total_modis_cells+=num_mod_cells
            if False in ds_diff.values():
                data_mismatch[tile_id] = ds_diff
                if 'count' in ds_diff.keys():
                    cells_missing_features+=ds_diff['count']
        else:
            missing_files.append(tile_id)

    count_fmt = '{0}/{1} modis cells (in {2} tiles) not completely filled by landcover'
    print 'missing files in {} tiles'.format(len(missing_files))
    print count_fmt.format(cells_missing_features, total_modis_cells, len(data_mismatch))

    return missing_files, data_mismatch


def Restart_isect(project):
    project['restart'] = {}
    csv_format = os.path.join(project['csv_dir'], '{0}/ssarV2_{0}_{1}.csv')
    dsets = ['BSA_ancill', 'BSA_nir', 'BSA_sw', 'BSA_band', 'BSA_quality', 'BSA_vis', 'fsw', 'sd', 'tam', 'rr', 'swe', 'lc']


    # Find last tile completed
    tiles = glob.glob(os.path.join(project['csv_dir'], 'lc/*.csv'))
    tile_ids = ['_'.join(tile.split('.')[0].split('_')[-2:]) for tile in tiles]
    rows_cols = np.array([[int(val) for val in tile_id.split('_')] for tile_id in tile_ids])

    last_row = rows_cols[np.where(rows_cols[:,0]==np.max(rows_cols[:,0]))]
    last_tile = last_row[np.where(last_row[:,1]==np.max(last_row[:,1]))][0]

    last_tile_id = '_'.join((str(last_tile[0]), str(last_tile[1])))

    # Find first mod_idVar in last tile
    lc_csv = csv_format.format('lc',last_tile_id)
    with open(lc_csv,'r') as lc:
        hdr = lc.readline()
        start = lc.readline()
    project['restart']['modis_idVar'] = int(start.split(',')[2])


    # Clean up last tile
    print 'Removing output:',last_tile_id
    for ds in dsets:
        if os.path.isfile(csv_format.format(ds,last_tile_id)):
            os.remove(csv_format.format(ds,last_tile_id))

    # Remake rows_cols
    tiles = glob.glob(os.path.join(project['csv_dir'], 'lc/*.csv'))
    tile_ids = ['_'.join(tile.split('.')[0].split('_')[-2:]) for tile in tiles]
    rows_cols = np.array([[int(val) for val in tile_id.split('_')] for tile_id in tile_ids])

    # Update nanmask to prevent duplication of already finished tiles
    project['restart']['modis_nanmask'] = Mk_modis_nanmask(project)
    project['restart']['modis_nanmask'][rows_cols[:,0], rows_cols[:,1]] = 0

    Isect_mod_clim_ssar(project)
    # todo find which ogr op is throwing empty intersection warnings and logic it


def Pack_it_up(project):
    csv_dir = project['csv_dir']
    csv_format = os.path.join(project['csv_dir'], '{0}/ssarV2_{0}_{1}.csv')# .format(sds,tile_id)
    dsets = ['BSA_ancill', 'BSA_nir', 'BSA_sw', 'BSA_band', 'BSA_quality', 'BSA_vis', 'fsw', 'sd', 'tam', 'rr', 'swe', 'lc']

    # Get list of tiles
    tiles = glob.glob(os.path.join(project['csv_dir'], 'lc/*.csv'))
    tile_ids = ['_'.join(tile.split('.')[0].split('_')[-2:]) for tile in tiles]

    # make / clean directory for unfilled_modis_cells
    # make / clean directory for zip archives

    print 'Number of tiles: ', len(tile_ids)
    for tile_id in tile_ids:
        unfilled = False
        unfilled_modis = Check_tile(tile_id, csv_fmt=csv_format, tolerance=1, return_details=True)
        if len(unfilled_modis.keys())>0:
            # add a csv listing modis_id and % filled
            unfilled = True
            unfilled_csv = csv_format.format('unfilled_mod_cells',tile_id)
            Write_out_unfilled(unfilled_csv, unfilled_modis)
        file_list = [csv_format.format(ds, tile_id) for ds in dsets]
        if unfilled==True:
            file_list.append(unfilled_csv)
        archive_fn = os.path.join(project['csv_dir'], 'ZIP', 'ssarV2_'+tile_id+'.zip')
        Zip_files(archive_fn, file_list)


def Write_out_unfilled(csv_fn, details):
    csv_f   = open(csv_fn,'wt')
    details_csv = csv.writer(csv_f)
    hdr = ['modis_id','modis_area - Sum(feat_areas) (m2)']
    details_csv.writerow(hdr)
    for i in range(len(details.keys())):
        details_csv.writerow([details.keys()[i], details.values()[i]])
    csv_f.close()


def Zip_files(archive_fn, file_list):

    zf = zipfile.ZipFile(archive_fn, mode='w')
    for fn in file_list:
        zf.write(fn, compress_type=compression)
    zf.close()
