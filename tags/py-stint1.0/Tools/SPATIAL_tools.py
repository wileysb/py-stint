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

import os,sys,math,cPickle
import rtree,ogr,osr,gdalconst
import numpy as np

class FastRtree(rtree.Rtree):
    def dumps(self, obj):
        return cPickle.dumps(obj, -1)


def Reprj_and_idx( src_dsn, dst_dsn, dst_srs, fields ):
    '''
    Reprj_and_idx(src_dsn,dst_dsn,dst_srs,fields_to_preserve)
    in_shp  = 'path/to/in' # no extension 
    dst_dsn = 'path/to/out' # no extension
    dst_srs = osr.SpatialReference(), initialized to output ref system
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
    
    defn = dst_lyr.GetLayerDefn()

    # Initialize progress updater
    count_max    = float(src_lyr.GetFeatureCount())
    count_update = count_max * 0.05 # print progress every 5%!

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
        dst_feat.SetGeometry(geom)

        dst_lyr.CreateFeature(dst_feat)
        dst_feat = geom = None
        # Print progress:
        if int( math.fmod( fid, count_update ) ) == 0:
            prog = int( fid / count_max * 100 )
            report = '%s%% . . ' % prog
            sys.stdout.write( report )
            sys.stdout.flush()
    
    sys.stdout.write("\n")
    dst_prj = dst_srs.ExportToWkt()
    Mk_proj(dst_prj,dst_dsn)
    # Close, flush, save output files
    r.close()
    dst_ds = dst_lyr = dst_feat = geom = defn = None


def Mk_proj(prj,outf):
    with open("%s.prj" % outf, "w") as proj:
        prj_out = osr.SpatialReference()
        prj_out.ImportFromWkt(prj)
        prj_out.MorphToESRI()
        prj = prj_out.ExportToPrettyWkt()
        proj.write(prj)


def Isect_poly_idx( src1_dsn, src1_pre, src1_id, src1_fields, area,
                    src2_dsn, src2_pre, src2_id, src2_fields, dst_dsn ):
    '''
    Isect(src1,src2,dst)
    assume src1 and src2:
    * are polygon shapefiles
    * have shapefile and index extensions
    * have same projection
    '''
    # Load first dataset, and projection
    src_ds1  = ogr.Open(src1_dsn+'.shp',gdalconst.GA_ReadOnly)
    src_lyr1 = src_ds1.GetLayer(0)
    srs      = src_lyr1.GetSpatialRef()
    
    # Load second source dataset, plus index
    src_ds2  = ogr.Open(src2_dsn+'.shp',gdalconst.GA_ReadOnly)
    src_lyr2 = src_ds2.GetLayer(0)
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
    count_update = count_max * 0.05 # print progress every 5%!
    
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
        # Print progress:
        if int( math.fmod( fid1, count_update ) ) == 0:
            prog = int( fid1 / count_max * 100 )
            report = '%s%% . ' % prog
            sys.stdout.write( report )
            sys.stdout.flush()

    sys.stdout.write("\n")
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


def Isect_ras_poly(ras_fn,poly_dsn,dst_p):
    '''
    Isect_ras_poly(ras,poly_dsn,dst_dsn)
    raster should have extension
    poly_ and dst_dsn should have no extension
    raster and poly should already be in the same projection
    Function loops through features in poly, finds ras cells
    wholly within the feature, and intersect ras cells on the border
    '''
    out = {}
    ras      = Parse_extents(ras_fn)
    poly_ds  = ogr.Open(poly_dsn+'.shp',gdalconst.GA_ReadOnly)
    poly_lyr = poly_ds.GetLayer(0)
    
    ras_x    = np.arange(ras['xmin'],ras['xmax'],ras['dx'])
    ras_y    = np.arange(ras['ymax'],ras['ymin'],-1*ras['dy'])
    
    ras_box  = Mk_bbox(min(ras_x),min(ras_y),max(ras_x),max(ras_y))
    
    # Initialize progress updater
    count_max    = float(poly_lyr.GetFeatureCount())
    count_update = count_max * 0.05 # print progress every 5%!
    
    for fid in range(0,poly_lyr.GetFeatureCount()):
        # out[fid] = [[within],[intersecting]]
        # within   = [x,y], x,y in raster coordinates (ncol,nrow)
        # isecting = [x,y,area]
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
                            out[fid][0].append((px,py))
                        else:
                            isect = ras_cell.Intersection(geom)
                            out[fid][1].append((px,py,isect.Area()))
            
            # Tuple to Array if there were any intersecting cells
            # Remove this item if there were not
            
            # This conversion slows down the code, but should shrink
            # the pickle size at the end?
            keep = False
            if len(out[fid][0]) > 0:
                dtype = [('x',np.int64),('y',np.int64)]
                out[fid][0] = np.array(out[fid][0],dtype)
                keep = True
            if len(out[fid][1]) > 0:
                dtype=[('x',np.int64),('y',np.int64),('area',np.float64)]
                out[fid][1] = np.array(out[fid][1],dtype)
                keep = True
            if not keep:
                jnk = out.pop(fid)
            
        # Print progress:
        if int( math.fmod( fid, count_update ) ) == 0:
            prog = int( fid / count_max * 100 )
            toprint = '%s%% . . ' % prog
            print toprint
    
    cPickle.dump(out,open(dst_p,'w'))
    

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
    progVar = 0
    count        = 0
    count_max    = (ymax-ymin) / dy
    count_update = count_max * 0.05 # print progress every 5%!

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
        if int( math.fmod( count, count_update ) ) == 0:
            prog = int( count / count_max * 100 )
            report = '%s%% . . ' % prog
            sys.stdout.write( report )
            sys.stdout.flush()
    
    # Save and close everything
    ds = layer = feat = perim = polygon = None
    sys.stdout.write("\n")
    prj = params['srs'].ExportToWkt()
    Mk_proj( prj, params['outf'] )


