#  Control tool for Spatial/Temporal Intersection Toolset
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

def _cmdhelp():
    '''Documentation for py-stint run.py command line options.

    Usage:
    $ python run.py -help
    prints this documentation

    $ python run.py
    no arguments: executes stages 1-7 in order

    $ python run.py <integer>
    Integer from 1 to 7
    Executes specified stage independently


    Stages
    -------
    These are the stages in the py-stint workflow
    +------+---------------------+----------------------+--------------------+
    |Stage | Summary             | Requires             | Produces           |
    +======+=====================+======================+====================+
    |      |                     | MODIS archive,       |                    |
    |      |                     | ERA archive          |                    |
    | 1    |Parse input,         | INPUT.txt,           |  directories       |
    |      |check archives       | aoi/landcover,       |                    |
    |      |                     | MODIS tiles shp      |                    |
    +------+---------------------+----------------------+--------------------+
    | 2    | src -> hdf5         | As stage 1           | datasets.hdf5      |
    +------+---------------------+----------------------+--------------------+
    | 3    | hdf5 -> shp         | ERA & MODIS          | ERA & MODIS shp    |
    |      |                     | hdf5 array files     | -native projections|
    +------+---------------------+----------------------+--------------------+
    | 4    | reproject+index     | ERA & MODIS shp:     | ERA & MODIS shp+idx|
    |      | ERA & MODIS shp     | -native projections  | -aoi/landcover prj |
    +------+---------------------+----------------------+--------------------+
    | 5    | lc + climate -> lcc |landcover/aoi shp     |lcc.shp             |
    |      |                     |era-reproj shp+idx    |                    |
    +------+---------------------+----------------------+--------------------+
    | 6    | lcc + modis -> lcm  | lcc.shp,             |lcm.shp             |
    |      |                     | modis-reproj shp+idx |                    |
    +------+---------------------+----------------------+--------------------+
    | 7    |  Output hdf5+lcm    | lcm.shp,             | datasets.csv       |
    |      |  -> CSV             | datasets.hdf5        |                    |
    +------+---------------------+----------------------+--------------------+



      Stage 1
        Loads project parameters from INPUT.txt
        Loads landcover/AOI shapefile
        Checks MODIS archive for corrupt or missing files
          within in the region, timeframe, and datasets specified
        Checks ERA archive for gaps in timeseries
          within the timeframe and datasets specified

      Stage 2
        Constructs array file (hdf5) with dimensions [d,y,x] for each
        MODIS and ERA dataset
        * x and y are spatial coordinates in dataset's native projection
        * d represents index number of modis interval
        * d=0 is the first modis interval in the timeframe, eg 2000049
        * hdf5 arrays are subset to the project timeframe and minimum
          rectangle necessary to completely contain the landcover region (AOI)

      Stage 3
        Produce MODIS & ERA shapefiles from hdf5 subset arrays
        * Shapefiles are in native projection (Sinusoidal & WGS84)
        * consist of a polygon grid, each raster cell represented by
          rectangle with attributes id,x_ctr,y_ctr,x_ind,y_ind

      Stage 4
        Reproject MODIS & ERA shapefiles to landcover(aoi) projection
        Construct spatial rtree index (idx) for each

      Stage 5
        Intersect reprojected climate (ERA) + landcover (aoi) shapefiles
        to produce lcc.shp
          * lcc = land cover climate
          * attributes specified by lc_ in INPUT are preserved from aoi
          * era idx to speed up feature matching

      Stage 6
        Intersect lcc + reprojected modis shapefiles to produce lcm.shp
        * lcm = landcover climate modis
        * modis idx to speed up feature matching
        * lcm.shp is effectively the central product of this workflow, linking
          breaking each landcover feature into polygons wholly within and
          linked to individual MODIS and ERA cells

      Stage 7
        lcm.shp used as guide to export modis and era timeseries to CSV
        * lc.csv: each row is a landcover feature with MODIS and ERA links
        * MODIS and ERA dataset csvs: each row is a full timeseries for an
          individual cell
        * datasets with more than 20k landcover features are broken up into
          regions made up of 50 MODIS cells

    '''
    pass

import sys,os
import numpy as np
import Tools.SPATIAL_tools
from Tools.ORG_tools import Parse_input, Check_input,Get_modis_days
from Tools.MODIS_aoi import Parse_extents, Check_mod_tiles
from Tools.MODIS_librarian import Gather_mod_flaws
from Tools.ERA_librarian import Val_era
from Tools.MODIS_subsetter import Mod2hdf
from Tools.ERA_subsetter import Era2hdf
from Tools.STINT_outputs import Veclc2csv
from Tools.pyModis_downmodis import downModis

# Load INPUT variables to 'project' dictionary
project = Parse_input('INPUT.txt')

# Check for errors in input
err = Check_input(project)
if len(err):
    print 'FIX ERRORS IN INPUT.txt:'
    for problem in err:
        print problem
    sys.exit( 1 )

### DEFINE AND RESOLVE HDF DIRECTORY
processing_dir = os.path.join(project['prj_directory'],'Processing/')
hdf_dir = os.path.join(processing_dir,'HDF/')
if not os.path.isdir(processing_dir):
    os.mkdir(processing_dir)

if not os.path.isdir(hdf_dir):
    os.mkdir(hdf_dir)

## DEFINE AND RESOLVE SHP DIRECTORY
shp_dir = os.path.join(processing_dir,'SHP/')
if not os.path.isdir(shp_dir):
    os.mkdir(shp_dir)

## DEFINE AND RESOLVE TIF DIRECTORY
tif_dir = os.path.join(processing_dir,'TIF/')
if not os.path.isdir(tif_dir):
    os.mkdir(tif_dir)

## DEFINE AND RESOLVE OUTPUT DIRECTORY
out_dir = os.path.join(project['prj_directory'],'Output/')
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

## DEFINE AND RESOLVE CSV OUTPUT DIRECTORY
csv_dir = os.path.join(out_dir,'CSV/')
if not os.path.isdir(csv_dir):
    os.mkdir(csv_dir)

# Expand 'project' dictionary
project['NODATA']  = np.nan
project['aoi']     = Parse_extents(project['lc_src'])
project['hdf_dir'] = hdf_dir
project['shp_dir'] = shp_dir
project['tif_dir'] = tif_dir
project['out_dir'] = out_dir
project['csv_dir'] = csv_dir
project['modis_tiles'] = Check_mod_tiles(project['modis_tile_fn'],**project['aoi'])
project['modis_days']  = Get_modis_days( project['start_year'],
                                         project['end_year'])

def Proceed(): # not coded yet!
    #num_stages = 7
    #last_stage = Ld_stage_file() # Read stage file
    # Check products up through last_stage
    # Clean up anything beyond last_stage and check directory structure
    #start_stage = last_stage
    #for stage_num in range(start_stage,num_stages+1):
    #    Run_stage(stage_num)
    TODO = 'code this'


def Get_modis(tiles,dates,modis_dir,dset='MCD32A2'):
    dest_dir = os.path.join(modis_dir,dset)
    dataset  = dset+'.005'
    lday = dates[-1]  # last date in range, datetime
    sday = dates[0]   # first date in range, datetime
    dm = downModis(dest_dir,product=dataset,
                   today=lday,enddate=sday,delta=None)
    dm.tiles = tiles
    dm.connect()
    dm.downloadsAllDay(clean=True)
    

def Run_stage(stage_num):
    if stage_num == 1:
        # Validate MODIS datasets for ID'd tiles
        mflaws = Gather_mod_flaws(project)
          # Write modis flaws to file so they can be used to fix archive
          # report flaw summary
        # Validate ERA extents, date range, and file overlap
        eflaws = Val_era(project)
        if sum(eflaws.values()) == True:
            for k in eflaws.keys():
                if eflaws[k] == True:
                    print 'Some issues with era',k
    elif stage_num == 2:
        # Collect MODIS data to 3d subset grids
        # save to .hdf5 in Processing/HDF/
        for dset in project['modis'].keys():
            for dnum in project['modis'][dset].keys():
                Mod2hdf( project, dset, dnum )

        # Collect ERA data to 3d subset grids
        # save to .hdf5 in Processing/HDF/
        for sds in project['era'].keys():
            Era2hdf( project, sds )
        # tools for examining hdf data: hdf2tif,dset2grid
    elif stage_num == 3:
        pre  = project['prj_name']+'_' # start to all the project files

        # Create MODIS shapefile from sample hdf5 stack
        msds = project['modis'].values()[0].values()[0]
        modis_fn           = os.path.join( project['hdf_dir'],
                                           pre + msds + '.hdf5' )
        modis_shp          = os.path.join( project['shp_dir'],
                                           pre + 'modis' )
        mod_params         = Parse_extents( modis_fn )
        mod_params['outf'] = modis_shp
        print 'Creating '+pre+'modis.shp from '+modis_fn
        Tools.SPATIAL_tools.Mk_polygrid( mod_params )
                        
        # Create ERA shapefile from sample hdf5 stack
        esds = project['era'].keys()[0]
        era_fn             = os.path.join( project['hdf_dir'],
                                           pre + esds + '.hdf5' )
        era_shp            = os.path.join( project['shp_dir'],
                                           pre + 'era' )
        era_params         = Parse_extents( era_fn )
        era_params['outf'] = era_shp
        print 'Creating '+pre+'era.shp from '+era_fn
        Tools.SPATIAL_tools.Mk_polygrid( era_params )
        
    elif stage_num == 4:
        pre  = project['prj_name']+'_' # start to all the project files
        # Reproject and index MODIS shapefile (to land cover projection)
        modis_shp          = os.path.join( project['shp_dir'],
                                           pre + 'modis' )                                   
        mod_tx             = os.path.join( project['shp_dir'],
                                           pre+'modis_reprj' )
        mod_reprj = {'dst_dsn':mod_tx, 'src_dsn':modis_shp,
                     'dst_srs':project['srs'], 
                     'fields':['ctr_x','ctr_y','x_ind','y_ind']}
        print 'Reprojecting '+pre+'modis.shp to project reference system and creating rtree index'
        Tools.SPATIAL_tools.Reprj_and_idx( **mod_reprj)
        
        # Reproject and index ERA shapefile (to land cover projection)
        era_shp            = os.path.join( project['shp_dir'],
                                           pre + 'era' )
        era_tx             = os.path.join( project['shp_dir'],
                                           pre+'era_reprj' )
        era_reprj = {'dst_dsn':era_tx, 'src_dsn':era_shp,
                     'dst_srs':project['srs'], 
                     'fields':['ctr_x','ctr_y','x_ind','y_ind']}
        print 'Reprojecting '+pre+'era.shp to project reference system and creating rtree index'
        Tools.SPATIAL_tools.Reprj_and_idx( **era_reprj)
    elif stage_num == 5:
        pre  = project['prj_name']+'_' # start to all the project files
        # Create landcover-climate intersection (lcc)
        lcc_dsn = os.path.join( project['shp_dir'],
                                pre+'lcc' )
        era_dsn             = os.path.join( project['shp_dir'],
                                           pre+'era_reprj' )
        lcc_params = {'src1_dsn':os.path.splitext(project['lc_src'])[0],
                      'src1_pre': 'lc_',
                      'src1_id' : 'lc_',
                      'src1_fields':project['lc'].values(),
                      'src2_dsn': era_dsn,
                      'src2_pre': 'era_',
                      'src2_id' : 'era_',
                      'src2_fields': ['ctr_x','ctr_y','x_ind','y_ind'],
                      'dst_dsn' : lcc_dsn,
                      'area'    : False  }
        print 'Intersecting LandCover & Climate grids -> lcc.shp'
        Tools.SPATIAL_tools.Isect_poly_idx( **lcc_params )
        
    elif stage_num == 6:
        pre  = project['prj_name']+'_' # start to all the project files
        # Create lcc-MODIS intersection (final product, shapefile form!)
        lcc_dsn             = os.path.join( project['shp_dir'],
                                            pre+'lcc' )
        lcc_fields = [ 'era_ctr_x', 'era_ctr_y', 'era_id',
                       'era_x_ind', 'era_y_ind', 'lc_id']
        for val in project['lc'].values():
            lcc_fields.append('lc_'+val)
        mod_dsn             = os.path.join( project['shp_dir'],
                                            pre+'modis_reprj' )
        lcm_dsn           = os.path.join( project['shp_dir'],
                                            pre+'lcm' )
        lcm_params = {'src1_dsn':lcc_dsn,
                      'src1_pre': '',
                      'src1_id' : 'lcc_',
                      'src1_fields':lcc_fields,
                      'src2_dsn': mod_dsn,
                      'src2_pre': 'mod_',
                      'src2_id' : 'mod_',
                      'src2_fields':['ctr_x','ctr_y','x_ind','y_ind'],
                      'dst_dsn' : lcm_dsn,  
                      'area'    : True      }
        print 'Intersecting lcc.shp with MODIS grid for final Landcover-Climate-MODIS (lcm) shapefile'
        Tools.SPATIAL_tools.Isect_poly_idx( **lcm_params )

    elif stage_num == 7:
        # Export era, modis, and intersected landcover to csv!
        Veclc2csv(project)

    elif stage_num > 7:
        print "OK, OK, you're done already!"
        print "only 7 stages..."


if __name__ == '__main__':
    '''
    CLI arguments: 
    <none>: proceed from last completed stage
    1 : run stage 1
    2 : run stage 2
    3 ...
    '''
    if len( sys.argv ) > 2: # first argument is 'run.py'
        print "[ ERROR ] bad arguments"
        print "$ python run.py to proceed automatically; "
        print "eg $ python run.py <stage_num> to run specified stage"
        sys.exit( 1 )
    
    if len( sys.argv ) == 1:
        Proceed()
    elif len( sys.argv ) == 2:
        try:
            stage_num = int(sys.argv[1])
            Run_stage(stage_num)
        except:
            # Encourage the user to provide an integer
            print '[ ERROR ] executing run %s; is it a normal integer?' % sys.argv[1]
