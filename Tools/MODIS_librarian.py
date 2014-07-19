#  Organizational tools for MODIS archive 
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

import os,sys,glob
import gdal,gdalconst
import numpy as np


def Retrieve_archive(project):
    skip = 'code this next'


def Get_modis_fn(modis_dir, dset, modis_day, mtile):
    '''
    fn       = '*A%s.%s.*.hdf' % (modis_day,mtile)
    dset_dir = os.path.join(modis_dir,dset)
    src_fn   = glob.glob(os.path.join(dset_dir,fn))[0]
    '''
    fn       = '*A%s.%s.*.hdf' % (modis_day,mtile)
    dset_dir = os.path.join(modis_dir,dset)
    tile_dir = os.path.join(dset_dir,mtile)
    src_fn   = glob.glob(os.path.join(tile_dir,fn))[0]
    
    return src_fn


def Val_mod(modis_dir, dset, mtile, subdsets, mdays):
    missing_files = []
    corrupt_files = []
    corrupt_dsets = []

    for modis_day in mdays:
        try:
            fn = Get_modis_fn(modis_dir, dset, modis_day, mtile)
            try:
                ds = gdal.Open(fn)
                ds_sds = ds.GetSubDatasets()
                del ds
                for dnum in subdsets.keys():
                    try:
                        sds = gdal.Open(ds_sds[dnum][0])
                        del sds
                    except:
                        corrupt_dsets.append((fn,dnum))
            except:
                corrupt_files.append(fn)
        except:
            missing_files.append((dset,mtile,modis_day))
    
    out = {
           'missing':missing_files,
           'corrupt':corrupt_files,
           'partial':corrupt_dsets
           }
    return out


def Gather_mod_flaws(project):
    flaws = {}
    for mtile in project['modis_tiles']:
        flaws[mtile] = {}
        for mdset in project['modis'].keys():
            subdsets = project['modis'][mdset]
            flaws[mtile][mdset] = Val_mod( project['modis_dir'],
                                           mdset, mtile, subdsets,
                                           project['modis_days'] )
    return flaws


def Record_flaws(flaws,out_dsn):
    if len(flaws['missing'])>0:
        with open(out_dsn+'.missing','w') as outf:
            for i in range(len(flaws['missing'])):
                outf.write('%s %s %s \n' % flaws['missing'][i])
    if len(flaws['corrupt'])>0:
        with open(out_dsn+'.corrupt','w') as outf:
            for fn in flaws['corrupt']:
                outf.write(fn+'\n')
    if len(flaws['partial'])>0:
        with open(out_dsn+'.partial','w') as outf:
            for i in range(len(flaws['partial'])):
                outf.write('%s %s \n' % flaws['partial'][i])


if __name__ == '__main__':
    flaw_sum = '%s %s had %s missing, %s corrupt, and %s partial files'
    for dset in MOD_dsets.keys():
        for mtile in ['h18v01','h18v02','h18v03','h19v01','h19v02','h19v03']:
            flaws = Val_mod(dset, mtile, mdays)
            flaw_args = (dset,mtile,len(flaws['missing']),
                    len(flaws['corrupt']),len(flaws['partial']))
            print flaw_sum % flaw_args
            Record_flaws(flaws,dset+'_'+mtile)

