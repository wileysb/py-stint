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

import sys,os
import Tools.SPATIAL_tools
from Tools.ORG_tools import Parse_input, Check_input,Get_modis_days
from Tools.MODIS_aoi import Parse_extents, Check_mod_tiles
from Tools.MODIS_librarian import Val_mod, Gather_mod_flaws
from Tools.ERA_librarian import Val_era
from Tools.MODIS_subsetter import Mod2hdf
from Tools.ERA_subsetter import Era2hdf

# Load INPUT variables to 'project' dictionary
project = Parse_input('INPUT_dev.txt')

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
# Expand 'project' dictionary
project['aoi']         = Parse_extents(project['lc_src'])
project['hdf_dir'] = hdf_dir
project['shp_dir'] = shp_dir
project['tif_dir'] = tif_dir
project['modis_tiles'] = Check_mod_tiles(**project['aoi'])
project['modis_days']  = Get_modis_days(project['start_year'], 
                                       project['end_year'])

def Proceed(): # not coded yet!
    num_stages = 7
    last_stage = Ld_stage_file() # Read stage file
    # Check products up through last_stage
    # Clean up anything beyond last_stage and check directory structure
    start_stage = last_stage
    for stage_num in range(start_stage,numstages+1):
        Run_stage(stage_num)    


def Get_modis(tiles,dates,modis_dir,mflaws):
    TODO = 'Download MODIS dsets'
    TODO+='Try to get Missing datasets'
    TODO+='Delete corrupt datasets, try to replace them'
    

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
            print '[ ERROR ] argument, if any, should be normal integer'
