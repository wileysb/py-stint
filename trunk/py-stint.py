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

'''Commandline control for Spatial/Temporal Intersection Toolset

    Usage:
    python run.py --help
      prints this documentation

    python run.py <integer>
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
    | 1    | Parse input,        | INPUT.txt,           | directories        |
    |      | check archives      | aoi/landcover,       |                    |
    |      |                     | MODIS tiles shp      |                    |
    +------+---------------------+----------------------+--------------------+
    | 2    | hdf/netcdf -> hdf5  | As stage 1           | datasets.hdf5      |
    +------+---------------------+----------------------+--------------------+
    | 3    | hdf5 -> shp         | ERA & MODIS          | ERA & MODIS shp    |
    |      |                     | hdf5 array files     | -native projections|
    +------+---------------------+----------------------+--------------------+
    | 4    | reproject+index     | ERA & MODIS shp:     | ERA & MODIS shp+idx|
    |      | ERA & MODIS shp     | -native projections  | -aoi/landcover prj |
    +------+---------------------+----------------------+--------------------+
    | 5    | lc + climate -> lcc | landcover/aoi shp    | lcc.shp            |
    |      |                     | era-reproj shp+idx   |                    |
    +------+---------------------+----------------------+--------------------+
    | 6    | lcc + modis -> lcm  | lcc.shp,             | lcm.shp            |
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

def prj_mkdir(dir_path):
    '''Create a directory, with a bit of extra syntax.

    Check if directory exists, print error and exit if directory cannot be created.
    This error probably indicates a bad path, pointing to a directory structure
    which doesn't exist yet.

    :param dir_path: path to directory to create
    :return: Exit(1) on failure, None on success
    '''
    if not os.path.isdir(dir_path):
        try:
            os.mkdir(dir_path)
        except:
            print '[ERROR] Problem creating ',dir_path
            sys.exit(1)


def Load_params(input_fn):
    '''Populate project parameters dictionary.

    :param input_fn: path to INPUT.txt
    :return: (dict)
    '''
    in_fn = input_fn
    try:
        project = Parse_input(in_fn)
    except:
        print '[ERROR] Problem loading parameters from', in_fn
        sys.exit(1)

    # Check for errors in input
    err = Check_input(project)
    if len(err):
        print 'FIX ERRORS IN %s:' % in_fn
        for problem in err:
            print problem
        sys.exit( 1 )

    ### DEFINE AND RESOLVE HDF DIRECTORY
    hdf_dir = os.path.join(project['prj_directory'],'HDF/')
    prj_mkdir(hdf_dir)

    ## DEFINE AND RESOLVE SHP DIRECTORY
    shp_dir = os.path.join(project['prj_directory'],'SHP/')
    prj_mkdir(shp_dir)

    ## DEFINE AND RESOLVE TIF DIRECTORY
    tif_dir = os.path.join(project['prj_directory'],'TIF/')
    prj_mkdir(tif_dir)

    ## DEFINE AND RESOLVE CSV OUTPUT DIRECTORY
    csv_dir = os.path.join(project['prj_directory'],'CSV/')
    prj_mkdir(csv_dir)

    # Expand 'project' dictionary
    project['NODATA']  = np.nan
    project['hdf_dir'] = hdf_dir
    project['shp_dir'] = shp_dir
    project['tif_dir'] = tif_dir
    project['csv_dir'] = csv_dir

    try:
        project['aoi']     = Parse_extents(project['lc_src'])
    except:
        print '[ERROR] Problems loading bounding box for AOI'
        print project['lc_src']
        sys.exit(1)

    try:
        project['modis_tiles'] = Check_mod_tiles(modis_tile_fn=project['modis_tile_fn'],**project['aoi'])
    except:
        print '[ERROR] Problems identifying MODIS tiles intersecting project AOI'

    try:
        project['modis_days']  = Get_modis_days( project['start_year'],
                                                 project['end_year'])
    except:
        print '[ERROR] Problems deriving MODIS interval days from specified timeframe'

    return project


def Run_stage(stage_num):
    '''Main workflow function.

    :param stage_num: (INT 1-7) workflow stage to execute
    :return: None
    '''
    if stage_num == 1:
        # Validate MODIS and ERA datasets
        stage_1()
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
        # Create MODIS and ERA shapefiles from sample hdf5 stack
        stage_3()
        
    elif stage_num == 4:
        # Reproject/index MODIS and ERA shapefiles --> land cover projection
        stage_4()

    elif stage_num == 5:
        if project['lc_type']=='shp':
            # landcover + climate -> lcc
            shp_stage_5()
        elif project['lc_type']=='tif_dir':
            # modis + climate -> mc
            tif_stage_5()
        
    elif stage_num == 6:
        if project['lc_type']=='shp':
            # lcc + modis -> lcm
            shp_stage_6()
        elif project['lc_type']=='tif_dir':
            # LandCover + mc -> lcmc
            tif_stage_6()

    elif stage_num == 7:
        if project['lc_type']=='shp':
            # Export era, modis, and intersected landcover to csv!
            Veclc2csv(project)
        elif project['lc_type']=='tif_dir':
            codethis = 'raslc2csv(project)'

    elif stage_num > 7:
        print "OK, OK, you're done already!"
        print "only 7 stages..."


def stage_1():
    print 'Looking for flaws in MODIS and ERA archives'
    # Validate MODIS datasets for ID'd tiles
    mflaws = Gather_mod_flaws(project)
    flawnum={'missing':0,'corrupt':0,'partial':0}
    for mtile in mflaws.keys():
        for mprod in mflaws[mtile].keys():
            for flaw in flawnum.keys():
                flawnum[flaw] += len(mflaws[mtile][mprod][flaw])

    for flaw in flawnum.keys():
        if flawnum[flaw]>0:
            print flawnum[flaw],flaw,'files in MODIS archive'
      # Write modis flaws to file so they can be used to fix archive
      # report flaw summary
    # Validate ERA extents, date range, and file overlap
    eflaws = Val_era(project)
    if sum(eflaws.values()) == True:
        for k in eflaws.keys():
            if eflaws[k] == True:
                print 'Some issues with era',k


def stage_3():
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


def stage_4():
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


def shp_stage_5():
    pre  = project['prj_name']+'_' # start to all the project files
    # Create landcover-climate intersection (lcc)
    lcc_dsn = os.path.join( project['shp_dir'],pre+'lcc' )
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


def tif_stage_5():
    pre  = project['prj_name']+'_'
    era_dsn             = os.path.join( project['shp_dir'],
                                       pre+'era_reprj' )
    mod_dsn             = os.path.join( project['shp_dir'],
                                       pre+'modis_reprj' )
    mc_dsn           = os.path.join( project['shp_dir'],
                                       pre+'mc' )
    lcc_params = {'src1_dsn': mod_dsn,
                  'src1_pre': 'mod_',
                  'src1_id' : 'mod_',
                  'src1_fields':['ctr_x','ctr_y','x_ind','y_ind'],
                  'src2_dsn': era_dsn,
                  'src2_pre': 'era_',
                  'src2_id' : 'era_',
                  'src2_fields': ['ctr_x','ctr_y','x_ind','y_ind'],
                  'dst_dsn' : mc_dsn,
                  'area'    : True  }

    print 'Intersecting Modis & Climate grids -> mc.shp'
    Tools.SPATIAL_tools.Isect_poly_idx( **lcc_params )


def shp_stage_6():
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


def tif_stage_6():
    pre  = project['prj_name']+'_' # start to all the project files
    ras_fn = project['lc_src']
    poly_dsn = os.path.join( project['shp_dir'],pre+'mc' )
    dst_p = os.path.join(project['prj_directory'],'lcmc.p')
    Tools.SPATIAL_tools.Isect_ras_poly(ras_fn,poly_dsn,dst_p)


def is_stagenum(s,numstages=7):
    '''Check if commandline argument is an acceptable stage number.

    :param s: (string) arg to check against numstages
    :param numstages: (int) Maximum stage_num for Run_stage(); Default 7
    :return: True or False
    '''
    try:
        argnum = int(s)
        if argnum in range(1,numstages+1):
            return True
        else:
            return False
    except ValueError:
        return False


def bad_arg_exit():
    '''Print standard error for bad CLI arguments and exit.'''
    print "[ ERROR ] bad arguments"
    print "to run a specified stage, try $ python run.py <stage_num>"
    print "for documentation, type $ python run.py --help"
    sys.exit( 1 )


if __name__ == '__main__':
    '''
    CLI arguments: 
    --help: print help __doc__
    1 : run stage 1
    2 : run stage 2
    3 ...
    '''

    if len(sys.argv)==2:
        if sys.argv[1] in ['-h','--help']:
            print __doc__
        elif is_stagenum(sys.argv[1]):
            project = Load_params('INPUT.txt')
            stage_num = int(sys.argv[1])
            Run_stage(stage_num)
        else:
            bad_arg_exit()
    elif len(sys.argv)==3:
        if os.path.isfile(sys.argv[1]):
            input_fn = sys.argv[1]
            if is_stagenum(sys.argv[2]):
                project = Load_params(input_fn)
                stage_num = int(sys.argv[2])
                Run_stage(stage_num)
            else:
                bad_arg_exit()
        else:
            bad_arg_exit()
    else:
        bad_arg_exit()
