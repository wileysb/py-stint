#  Organizational tools for ERA climate archive 
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

import osr
import datetime as dt
import numpy as np
#from Scientific.IO.NetCDF import NetCDFFile
from scipy.io.netcdf import netcdf_file as NetCDFFile
#Can this simply drop in like this?

from MODIS_aoi import Parse_extents, Mk_bbox
from ORG_tools import Yearday2hrnum
from ERA_subsetter import Build_src_nc, wgs84


def Val_era( project ):
    flaws = {}
    flaws['extents']    = False
    flaws['aoi']        = False
    flaws['timeseries'] = False
    flaws['modis_days'] = False

    ## Verify that all era extents are the same
    era = []
    for era_dset in project['era'].keys():
        for era_fn in project['era'][era_dset]:
            era.append( Parse_extents(era_fn) )
            jnk = era[-1].pop('srs') # srs == srs returns false
    
    for i in range(1,len(era)):
        if era[i] != era[0]:
            flaws['extents'] = True

    ## Verify that aoi extents are within era extents
    era_bbox = Mk_bbox( era[0]['xmin'], era[0]['ymin'],
                era[0]['xmax'], era[0]['ymax'] )
    
    aoi_bbox = Get_wgs84_aoi( project['aoi'] )
    if not aoi_bbox.Within( era_bbox ):
        flaws['aoi'] = True
    
    ## Verify timeseries do not overlap or jump between files
    era_t = {}
    for era_dset in project['era'].keys():
        for era_fn in project['era'][era_dset]:
            era_nc = NetCDFFile(era_fn,'r')
            try:
                era_t[era_dset] = np.hstack((era_t[era_dset], 
                                   era_nc.variables['time'][:]))
            except:
                era_t[era_dset] = era_nc.variables['time'][:]
            dtime     = np.diff(era_t[era_dset])
            errors = np.where( dtime != dtime[0] )
            
            # check that d_time is consistent
            if len(errors[0]) > 0:
                flaws['timeseries'] = True
    
    ## Verify timeseries contain modis_days
    basedate= dt.datetime(1900,01,01,00)
    mstart = project['modis_days'][0]
    m_end  = project['modis_days'][-1]

    estart = Yearday2hrnum(basedate,mstart)
    e_end  = Yearday2hrnum(basedate,m_end)
    for era_dset in era_t.keys():
        t_var = era_t[era_dset]
        if (np.min(t_var) > estart) or (np.max(t_var) < e_end):
            flaws['modis_days'] = True
    
    return flaws
    

def Get_wgs84_aoi( aoi ):
    aoi_bbox = Mk_bbox(aoi['xmin'], aoi['ymin'], 
                       aoi['xmax'], aoi['ymax'])

    # transform aoi from native prj to MODIS sinusoidal
    wgs84_transform = osr.CoordinateTransformation( aoi['srs'], 
                                                    wgs84) #src,dst
    # if CoordinateTransformation fails, it will return null:
    if wgs84_transform == None:
        print '[ERROR] Could not reproject AOI box to MODIS sinusoidal'
        sys.exit( 1 )

    aoi_bbox.Transform(wgs84_transform) 
    
    return aoi_bbox 

