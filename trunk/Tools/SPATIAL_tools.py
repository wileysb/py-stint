#  Intersection, reprojection, indexing, and other spatial functions
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
'''Module containing functions for basic shapefile and raster manipulation.

'''
import os
import cPickle
import sys
import rtree
import ogr
import osr
import gdalconst
import numpy as np
import sqlite3
from ORG_tools import Countdown

class FastRtree(rtree.Rtree):
    '''Accelerate Rtree implementation by using cPickle.'''
    def dumps(self, obj):
        return cPickle.dumps(obj, -1)


def Reprj_and_idx( src_dsn, dst_dsn, dst_srs, fields, area=False ):
    '''
    Reprj_and_idx(src_dsn,dst_dsn,dst_srs,fields_to_preserve)
    in_shp  = 'path/to/in' # no extension 
    dst_dsn = 'path/to/out' # no extension
    dst_srs = osr.SpatialReference(), initialized to output ref system
    :param area: (bool) True if 'area' attribute should be output
    '''
    # Load source dataset
    src_ds  = ogr.Open(src_dsn+'.shp',gdalconst.GA_ReadOnly)
    src_lyr = src_ds.GetLayer(0)
    src_srs = src_lyr.GetSpatialRef()
    
    # Open index and prepare coordinate transformation
    r = FastRtree(dst_dsn,interleaved=True)
    t = osr.CoordinateTransformation(src_srs,dst_srs) #src_srs,dst_srs
    # if CoordinateTransformation fails, it will return null:
    if t == None:
        print '[ERROR] Could not reproject between given reference systems'
        sys.exit( 1 )
    
    # Prepare output dataset
    driver = ogr.GetDriverByName('Esri Shapefile')
    dst_ds = driver.CreateDataSource(dst_dsn+'.shp')
    dst_lyr = dst_ds.CreateLayer('',None,ogr.wkbPolygon)
    dst_lyr.CreateField(ogr.FieldDefn('id',ogr.OFTInteger))
    
    # Define extra fields linking shapefile to raster&dataset
    src_lyr_defn = src_lyr.GetLayerDefn()
    for ind,field_name in enumerate(fields):
        try:
            field_i    = src_lyr_defn.GetFieldIndex(field_name)
            field_defn = src_lyr_defn.GetFieldDefn(field_i)
            dst_lyr.CreateField(field_defn)
        except:
            print 'Source dataset has no field named',field_name
            jnk = fields.pop(ind)
    if area == True:
        dst_lyr.CreateField(ogr.FieldDefn('area',ogr.OFTReal))

    defn = dst_lyr.GetLayerDefn()

    # Initialize progress updater
    count_max    = float(src_lyr.GetFeatureCount())
    # count_update = count_max * 0.05 # print progress every 5%!
    progress_bar = Countdown(count_max)

    for fid in range(0,src_lyr.GetFeatureCount()):
        src_feat = src_lyr.GetFeature(fid)
        geom = src_feat.GetGeometryRef()
        
        geom.Transform(t)
        
        # Insert to index
        xmin,xmax,ymin,ymax = geom.GetEnvelope()
        r.insert(fid,(xmin,ymin,xmax,ymax))
        
        # Write to output shapefile
        dst_feat = ogr.Feature(defn)
        dst_feat.SetField('id',fid)
        for field_name in fields:
            field_val = src_feat.GetField(field_name)
            dst_feat.SetField(field_name,field_val)
        if area == True:
                    dst_feat.SetField('area',geom.GetArea())

        dst_feat.SetGeometry(geom)

        dst_lyr.CreateFeature(dst_feat)
        dst_feat = geom = None
        # Print progress:
        progress_bar.check(fid)
        #if int( math.fmod( fid, count_update ) ) == 0:
        #    prog = int( fid / count_max * 100 )
        #    report = '%s%% . . ' % prog
        #   sys.stdout.write( report )
        #    sys.stdout.flush()

    progress_bar.flush()
    #sys.stdout.write("\n")
    dst_prj = dst_srs.ExportToWkt()
    Mk_proj(dst_prj,dst_dsn)
    # Close, flush, save output files
    r.close()
    dst_ds = dst_lyr = dst_feat = geom = defn = None


def Mk_proj(prj,outf):
    '''Write an ESRI style .prj file defining
    spatial reference for a shapefile dataset.

    :param prj: (str) wkt representation of dataset projection
    :param outf: (str) path to shapefile dataset, no extensions
    :return: None
    '''
    with open("%s.prj" % outf, "w") as proj:
        prj_out = osr.SpatialReference()
        val = prj_out.ImportFromWkt(prj)
        if val==0:
            prj_out.MorphToESRI()
            prj = prj_out.ExportToPrettyWkt()
            if len(prj)>0:
                proj.write(prj)
            else:
                print '[ERROR] osr projection string is empty:',outf
                sys.exit(1)
        else:
            print '[ERROR] Could not determine projection for',outf
            sys.exit(1)


def Ogr_open(dsn):
    '''Error handling for opening shapefiles.

    :param dsn: (str) path to shapefile, no extensions
    :return: (osgeo.ogr.DataSource) src_ds, (osgeo.ogr.Layer) src_lyr
    '''
    try:
        src_ds  = ogr.Open(dsn+'.shp',gdalconst.GA_ReadOnly)
        src_lyr = src_ds.GetLayer(0)
    # AttributeError indicates src_ds==None
    except AttributeError:
        print '[ERROR] could not access shapefile',dsn
        sys.exit(1)

    return src_ds,src_lyr


def Isect_poly_idx( src1_dsn, src1_pre, src1_id, src1_fields, area,
                    src2_dsn, src2_pre, src2_id, src2_fields, dst_dsn ):
    '''Intersect two shapefiles, assuming both datasets:
    * contain only polygon geometries
    * consist of a complete set of shapefile and rtree idx extensions
    * are already in the same projection or spatial reference system

    :param src1_dsn: (str) path to first input shapefile, no extensions
    :param src1_pre: (str) characters to prepend to output fields preserved from 1st src dsn
    :param src1_id: (str) characters to prepend to id field in output dataset
    :param src1_fields: (list) list of field names to preserve from 1st input dataset
    :param area: (bool) True if 'area' attribute should be output
    :param src2_dsn: (str) path to second input shapefile, no extensions
    :param src2_pre: (str) characters to prepend to output fields preserved from 2nd src dsn
    :param src2_id: (str) characters to prepend to id field in output dataset
    :param src2_fields: (list) list of field names to preserve from 2nd input dataset
    :param dst_dsn: (str) path to output shapefile, no extensions
    :return: None
    '''

    # Load first dataset, and projection
    src_ds1, src_lyr1  = Ogr_open(src1_dsn)
    srs      = src_lyr1.GetSpatialRef()

    # Load second source dataset, plus index
    src_ds2, src_lyr2  = Ogr_open(src2_dsn)
    src_r2   = FastRtree(src2_dsn)
    
    # Prepare output dataset
    driver = ogr.GetDriverByName('Esri Shapefile')
    dst_ds = driver.CreateDataSource(dst_dsn+'.shp')
    dst_lyr = dst_ds.CreateLayer('',None,ogr.wkbPolygon)
    
    # Basic id fields
    dst_lyr.CreateField(ogr.FieldDefn('id',ogr.OFTInteger))
    dst_lyr.CreateField(ogr.FieldDefn(src1_id+'id',ogr.OFTInteger))
    dst_lyr.CreateField(ogr.FieldDefn(src2_id+'id',ogr.OFTInteger))
    
    # Additional fields from src1
    src_lyr1_defn = src_lyr1.GetLayerDefn()
    for ind,field_name in enumerate(src1_fields):
        try:
            field_i    = src_lyr1_defn.GetFieldIndex(field_name)
            field_defn = src_lyr1_defn.GetFieldDefn(field_i)
            field_defn.SetName(src1_pre+field_name)
            dst_lyr.CreateField(field_defn)
        except:
            print 'First source dataset has no field named',field_name
            jnk = src1_fields.pop(ind)
    
    # Additional fields from src2
    src_lyr2_defn = src_lyr2.GetLayerDefn()
    for ind,field_name in enumerate(src2_fields):
        try:
            field_i    = src_lyr2_defn.GetFieldIndex(field_name)
            field_defn = src_lyr2_defn.GetFieldDefn(field_i)
            field_defn.SetName(src2_pre+field_name)
            dst_lyr.CreateField(field_defn)
        except:
            print 'First source dataset has no field named',field_name
            jnk = src2_fields.pop(ind)
    if area == True:
        dst_lyr.CreateField(ogr.FieldDefn('area',ogr.OFTReal))
    defn = dst_lyr.GetLayerDefn()
    
    # Open output index
    dst_r  = FastRtree(dst_dsn,interleaved=True)

    # ID counter for output shapefile
    idVar = 0
    
    # Initialize progress updater
    count_max    = float(src_lyr1.GetFeatureCount())
    # count_update = count_max * 0.05 # print progress every 5%!
    progress_bar = Countdown(count_max)

    # Loop through features in src1
    for fid1 in range(0,src_lyr1.GetFeatureCount()):
        src_feat1 = src_lyr1.GetFeature(fid1)
        geom1 = src_feat1.GetGeometryRef()
        
        gxmin,gxmax,gymin,gymax = geom1.GetEnvelope()
        
        # use src2 index to find intersections
        hits = src_r2.intersection((gxmin,gymin,gxmax,gymax))
        for hit_fid in hits:
            src_feat2 = src_lyr2.GetFeature(hit_fid)
            geom2 = src_feat2.GetGeometryRef()
            if geom1.Intersects(geom2):
                isect = geom1.Intersection(geom2)
                # Insert to index
                xmin,xmax,ymin,ymax = isect.GetEnvelope()
                dst_r.insert(idVar,(xmin,ymin,xmax,ymax))
                
                # Write to output shapefile
                dst_feat = ogr.Feature(defn)
                dst_feat.SetField('id',idVar)
                dst_feat.SetField(src1_id+'id',fid1)
                dst_feat.SetField(src2_id+'id',hit_fid)
                for field_name in src1_fields:
                    field_val = src_feat1.GetField(src1_pre+field_name)
                    dst_feat.SetField(src1_pre+field_name,field_val)
                for field_name in src2_fields:
                    field_val = src_feat2.GetField(src2_pre+field_name)
                    dst_feat.SetField(src2_pre+field_name,field_val)
                if area == True:
                    dst_feat.SetField('area',isect.GetArea())
                dst_feat.SetGeometry(isect)

                dst_lyr.CreateFeature(dst_feat)
                dst_feat = isect = None
                idVar+=1
        progress_bar.check(fid1)
    progress_bar.flush()
    dst_r.close()
    dst_r = dst_ds = dst_lyr = dst_feat = isect = defn = None
    dst_prj = srs.ExportToWkt()
    Mk_proj(dst_prj,dst_dsn)


def Get_spatial_indexes(x_var,y_var,bounds):
    '''x_s,x_e,y_s,y_e = Get_spatial_indexes(x_var,y_var)
    x_sub = x_var[x_s:x_e]
    y_sub = y_var[y_s:y_e]
    
    # index interpretation:
    if x_s==None:
        print 'aoi extends west of grid'
    if x_e==None:
        print 'aoi extends east of grid'
    if y_s==None:
        print 'aoi extends north of grid'
    if y_e==None:
        print 'aoi extends south of grid'
    '''
    xmin,xmax,ymin,ymax = bounds
    
    try:
        x_s = np.where(x_var<xmin)[0][-1]
    except:
        x_s = 0 # slicing [0:end] is equivalent to slicing[None:end]
    try:
        x_e = np.where(x_var>xmax)[0][0]
    except:
        x_e = None
    try:
        y_s = np.where(y_var>ymax)[0][-1]
    except:
        y_s = 0 # slicing [0:end] is equivalent to slicing[None:end]
    try:
        y_e = np.where(y_var<ymin)[0][0]
    except:
        y_e = None

    return x_s,x_e,y_s,y_e


def Isect_ras_poly(ras_fn,poly_dsn,dst_fn):
    '''
    Isect_ras_poly(ras,poly_dsn,dst_dsn)
    raster should have extension
    poly_ and dst_dsn should have no extension
    raster and poly should already be in the same projection
    Function loops through features in poly, finds ras cells
    wholly within the feature, and intersect ras cells on the border

    if dst_fn.split('.')=='p', output to pickle
    ## use a pickle output if the dataset is smallish
    ## (less than millions of output lcm features)
    ## The output file must load completely to memory
    ## while writing csv output
    if dst_fn.split('.')=='db', output to sqlite
    ## Use sqlite db output for larger regions or finer resolutions.
    ## THe output data gets loaded step by step, so this puts less
    ## strain on the system's memory during csv output
    ## and scales better to larger output datsets
    '''
    out = {}
    ras      = Parse_extents(ras_fn)
    poly_ds  = ogr.Open(poly_dsn+'.shp',gdalconst.GA_ReadOnly)
    poly_lyr = poly_ds.GetLayer(0)

    ras_x    = np.arange(ras['xmin'],ras['xmax'],ras['dx'])
    ras_y    = np.arange(ras['ymax'],ras['ymin'],-1*ras['dy'])
    cell_size = np.abs(ras['dy']) * np.abs(ras['dx'])

    ras_box  = Mk_bbox(min(ras_x),min(ras_y),max(ras_x),max(ras_y))

    # Initialize progress updater
    count_max    = float(poly_lyr.GetFeatureCount())
    # count_update = count_max * 0.05 # print progress every 5%!
    progress_bar = Countdown(count_max)

    if dst_fn.split('.')[-1]=='db':
        dst = 'db'
        conn = sqlite3.connect(dst_fn)
        c = conn.cursor()
        c.execute('''CREATE TABLE isect
                    (fid integer, px integer, py integer, area real)''')
        conn.commit()
    else:
        dst='p' # default to pickle

    for fid in range(0,poly_lyr.GetFeatureCount()):
        # out[fid] = [[within],[intersecting]]
        # within   = [fid,x,y,area], x,y in raster coordinates (ncol,nrow)
        # isecting = [fid,x,y,area]
        feat = poly_lyr.GetFeature(fid)
        geom = feat.GetGeometryRef()
        if geom.Intersects(ras_box):
            out[fid] = [[],[]]
            bounds = geom.GetEnvelope()
            x_s,x_e,y_s,y_e = Get_spatial_indexes(ras_x,ras_y,bounds)

            cx_range   = ras_x[x_s:x_e]
            cy_range   = ras_y[y_s:y_e]

            for i,cx in enumerate(cx_range):
                px = x_s+i
                for j,cy in enumerate(cy_range):
                    py = y_s+j
                    # xmin,ymin,xmax,ymax:
                    ras_cell = Mk_bbox(cx,cy-ras['dy'],cx+ras['dx'],cy)
                    if ras_cell.Intersects(geom):
                        if ras_cell.Within(geom):
                            if dst=='p':
                                coords=(px,py)
                            elif dst=='db':
                                coords=(fid,px,py,cell_size)
                            out[fid][0].append(coords)
                        else:
                            isect = ras_cell.Intersection(geom)
                            if dst=='p':
                                coords=(px,py,isect.Area())
                            elif dst=='db':
                                coords=(fid,px,py,isect.Area())
                            out[fid][1].append(coords)

            # Tuple to Array if there were any intersecting cells
            # Remove this item if there were not

            # This conversion slows down the code, but should shrink
            # the pickle size at the end?
            keep = False
            if len(out[fid][0]) > 0:
                if dst=='p':
                    dtype = [('x',np.int64),('y',np.int64)]
                    out[fid][0] = np.array(out[fid][0],dtype)
                    keep = True
                elif dst=='db':
                    c.executemany('INSERT INTO isect VALUES (?,?,?,?)',out[fid][0])
            if len(out[fid][1]) > 0:
                if dst=='p':
                    dtype=[('x',np.int64),('y',np.int64),('area',np.float64)]
                    out[fid][1] = np.array(out[fid][1],dtype)
                    keep = True
                elif dst=='db':
                    c.executemany('INSERT INTO isect VALUES (?,?,?,?)',out[fid][1])
            if not keep:
                jnk = out.pop(fid)
        conn.commit()
        progress_bar.check(fid)

    progress_bar.flush()

    if dst=='p':
        cPickle.dump(out,open(dst_fn,'w'))
    elif dst=='db':
        ### CREATE INDEX
        idx_sql = '''CREATE INDEX IF NOT EXISTS %s ON %s (%s)'''
        idx_name = 'fid_idx'
        ds_name = 'isect'
        idx_col = 'fid'
        sql = idx_sql % (idx_name, ds_name, idx_col)
        c.execute(sql)
        ### CLOSE CONNECTION
        conn.close()
    del out
    

def Parse_extents(src_fn):
    out = {}
    src_ext = os.path.splitext(src_fn)[1]
    if src_ext == None:
        print '[ ERROR ] input should be MODIS hdf, MODIS or ERA hdf5, ERA netcdf (CF-1.0), shapefile, or geoTIFF' 
        sys.exit( 1 )
    elif src_ext == '.shp':
        toprint =  'reading shapefile extents'
        src         = ogr.Open(src_fn)
        layer       = src.GetLayer()
        extent      = layer.GetExtent()
        out['dx']   = 'shapefile'
        out['dy']   = 'shapefile'
        out['srs']  = layer.GetSpatialRef()
        out['xmin'] = extent[0]
        out['xmax'] = extent[1]
        out['ymin'] = extent[2]
        out['ymax'] = extent[3]
    elif src_ext in ('.tif', '.hdf'):
        import gdal
        ras = gdal.Open(src_fn)
        if src_ext=='.hdf':
            toprint =  'reading MODIS hdf grid extents'
            # get subdataset (all have the same params in modis)
            dset0 = ras.GetSubDatasets()[0][0] # first listed subdataset
            del ras
            ras = gdal.Open(dset0)
        elif src_ext=='.tif':
            toprint =  'reading GeoTIFF extents'
        srs = osr.SpatialReference()
        tmp = srs.ImportFromWkt(ras.GetProjection())
        if tmp == 1:
            print '[ERROR] FAILED to import projection from raster file'
            sys.exit( 1 )
        gt          = ras.GetGeoTransform()
        out['srs']  = srs
        out['dx']   = gt[1]
        out['dy']   = abs(gt[-1]) # abs if this should be negative
        out['xmin'] = gt[0]
        out['ymax'] = gt[3]
        out['xmax'] = out['xmin'] + ras.RasterXSize*out['dx']
        out['ymin'] = out['ymax'] - ras.RasterYSize*out['dy']
        del ras
    elif src_ext == '.hdf5':
        toprint = 'reading hdf5 stack'
        import h5py
        fn = os.path.split(src_fn)[-1]
        prj_name = fn.split('_')[0]
        sds = os.path.splitext(fn)[0].split(prj_name+'_')[1]
        hdf = h5py.File(src_fn,'r')
        out['dx']   = hdf[sds].attrs['dx']
        out['dy']   = np.abs(hdf[sds].attrs['dy'])
        out['srs']  = osr.SpatialReference()
        out['srs'].ImportFromWkt(hdf[sds].attrs['projection'])
        out['xmin'] = np.min(hdf['x'])
        out['xmax'] = np.max(hdf['x'])+out['dx']
        out['ymin'] = np.min(hdf['y'])-out['dy']
        out['ymax'] = np.max(hdf['y'])
    elif src_ext == '.nc':
        #from Scientific.IO.NetCDF import NetCDFFile
        from scipy.io.netcdf import netcdf_file as NetCDFFile
        #Can this simply drop in like this?
        toprint =  'converting NetCDF Grid'
        toprint =  'assuming Geographic coordinates, WGS84'
        get4326 = osr.SpatialReference();get4326.ImportFromEPSG(4326)
        epsg4326 = get4326
        out['srs'] = epsg4326
        ras =  NetCDFFile(src_fn,'r')
        lat = ras.variables['latitude'][:]
        lon = ras.variables['longitude'][:]
        ras.close()
        out['dx']   = float(lat[0] - lat[1])
        out['dy']   = float(np.abs(lon[0] - lon[1]))
        out['xmin'] = float(np.min(lon)) - 0.5 * out['dx']
        out['ymax'] = float(np.max(lat)) + 0.5 * out['dy']
        out['xmax'] = float(np.max(lon)) + 0.5 * out['dx']
        out['ymin'] = float(np.min(lat)) - 0.5 * out['dy']
    else:
        print '[ ERROR ] input file not of type tif|hdf|nc' 
        sys.exit( 1 )

    return out


def Mk_bbox(xmin,ymin,xmax,ymax):
    perim = ogr.Geometry(ogr.wkbLinearRing)
    perim.AddPoint(xmin,ymin)
    perim.AddPoint(xmax,ymin)
    perim.AddPoint(xmax,ymax)
    perim.AddPoint(xmin,ymax)
    perim.AddPoint(xmin,ymin)

    bbox = ogr.Geometry(ogr.wkbPolygon)
    bbox.AddGeometry(perim)
    
    return bbox


def Mk_polygrid(params):
    '''Mk_polygrid(**params) 
    after params = Get_grid_params(rasterfile)
    
    QGIS fTools approach (VectorGrid tool); 
    modified to replace QGIS stuff with OGR
    and my little extent grabber functions'''
    xmin = float(params['xmin'])
    xmax = float(params['xmax'])
    ymin = float(params['ymin'])
    ymax = float(params['ymax'])
    dx   = params['dx']
    dy   = params['dy']
    prj  = params['srs']
    outf = params['outf']
    
    # Define outfile as shapefile, with feature type polygon
    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(outf+'.shp')
    layer = ds.CreateLayer('',None,ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn('id',ogr.OFTInteger))

    # Define extra fields linking shapefile to raster&dataset
    layer.CreateField(ogr.FieldDefn('ctr_x',ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('ctr_y',ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('x_ind',ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('y_ind',ogr.OFTInteger))

    # Define text field
    # field_name = ogr.FieldDefn("Name", ogr.OFTString)
    # field_name.SetWidth(24)
    # layer.CreateField(field_name)

    defn = layer.GetLayerDefn()

    idVar   = 0
    count        = 0
    count_max    = (ymax-ymin) / dy
    # count_update = count_max * 0.05 # print progress every 5%!
    progress_bar = Countdown(count_max)

    y = ymax
    y_ind = 0
    while round(y,3) > round(ymin,3):
        x = xmin
        x_ind = 0
        while round(x,3) < round(xmax,3):
            perim = ogr.Geometry(ogr.wkbLinearRing)
            perim.AddPoint(x, y)
            perim.AddPoint(x + dx, y)
            perim.AddPoint(x + dx, y - dy)
            perim.AddPoint(x, y - dy)
            perim.AddPoint(x, y)

            polygon = ogr.Geometry(ogr.wkbPolygon)
            polygon.AddGeometry(perim)

            feat = ogr.Feature(defn)
            feat.SetField('id',idVar)
            feat.SetField('ctr_x',x+0.5*dx)
            feat.SetField('ctr_y',y-0.5*dy)
            feat.SetField('x_ind',x_ind)
            feat.SetField('y_ind',y_ind)
            feat.SetGeometry(polygon)

            layer.CreateFeature(feat)
            feat = geom = None
            idVar += 1
            x += dx
            x_ind+=1
        y -= dy
        y_ind+=1
        count += 1
        progress_bar.check(count)
        #if int( math.fmod( count, count_update ) ) == 0:
        #    prog = int( count / count_max * 100 )
        #   report = '%s%% . . ' % prog
        #   sys.stdout.write( report )
        #    sys.stdout.flush()
    
    # Save and close everything
    ds = layer = feat = perim = polygon = None
    # sys.stdout.write("\n")
    progress_bar.flush()
    prj = params['srs'].ExportToWkt()
    Mk_proj( prj, params['outf'] )


def Dbf2db(poly_dsn, fields='all', idxs=None):
    '''
    poly_dsn = os.path.join( project['shp_dir'],'k3tif_mc' )
    fields = ['mod_id','era_id','mod_area']
    idxs = [['fid_idx','fid'],['mid_idx','mod_id']]

    :param src_dsn:
    :param dst_fn:
    :param fields:
    :return:
    '''
    # Derive some paths
    ds_path =  os.path.splitext(poly_dsn)[0]
    dst_db = ds_path + '.db'
    ds_name = os.path.basename(ds_path)

    # create connection objects
    poly_ds, poly_lyr = Ogr_open(poly_dsn)
    poly_lyr_defn = poly_lyr.GetLayerDefn()

    conn = sqlite3.connect(dst_db)
    c    = conn.cursor()

    # Get attribute data types:
    if fields=='all':
        attribs = {}
        for i in range(poly_lyr_defn.GetFieldCount()):
            field_defn = poly_lyr_defn.GetFieldDefn(i)
            aname = field_defn.GetName()
            fieldTypeCode = field_defn.GetType()
            attribs[aname] = {}
            attribs[aname]['type'] = field_defn.GetFieldTypeName(fieldTypeCode)
            attribs[aname]['width'] = field_defn.GetWidth()
            attribs[aname]['precision'] = field_defn.GetPrecision()
    else:
        attribs = {}
        for ind,aname in enumerate(fields):
            try:
                field_i    = poly_lyr_defn.GetFieldIndex(aname)
                field_defn = poly_lyr_defn.GetFieldDefn(field_i)
                fieldTypeCode = field_defn.GetType()
                attribs[aname] = {}
                attribs[aname]['type'] = field_defn.GetFieldTypeName(fieldTypeCode)
                attribs[aname]['width'] = field_defn.GetWidth()
                attribs[aname]['precision'] = field_defn.GetPrecision()
            except:
                print 'First source dataset has no field named',aname
                jnk = fields.pop(ind)

    # Make Database and Table
    sql_fmt = 'CREATE TABLE ' + ds_name + ' (%s)'
    col_fmt = 'id integer'
    for i in range(len(attribs.keys())):
        aname = attribs.keys()[i]
        col_fmt = col_fmt + ', %s %s' %  (aname, attribs[aname]['type'])

    sql = sql_fmt % col_fmt
    c.execute(sql)

    # LOOP THROUGH FEATURES AND POPULATE DATABASE
    out = []
    for fid in range(0,poly_lyr.GetFeatureCount()):
        feat = poly_lyr.GetFeature(fid)
        loop_out = [fid,]
        for aname in attribs.keys():
            loop_out.append(feat.GetField(aname))
        out.append(loop_out)
    insert_fmt = 'INSERT INTO ' + ds_name + ' VALUES (%s)'
    insert_sql = insert_fmt % ','.join('?'*(1+len(attribs.keys())))
    c.executemany(insert_sql,out)
    conn.commit()
    # make index on mod_id

    if idxs!=None:
        for idx in idxs:
            idx_name = idx[0]
            idx_col  = idx[1]
            idx_sql = '''CREATE INDEX IF NOT EXISTS %s ON %s (%s)'''
            sql = idx_sql % (idx_name, ds_name, idx_col)
            c.execute(sql)
        conn.commit()
    conn.close()
    del conn,c,poly_lyr,poly_ds
