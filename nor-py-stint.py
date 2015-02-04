import os
import sys
import osr
import Tools.SPATIAL_tools
from Tools.ORG_tools import Parse_input,Get_modis_days
from Tools.MODIS_aoi import Parse_extents, Check_mod_tiles
from Tools.MODIS_subsetter import Mod2hdf

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


if __name__ == '__main__':

    input = '/space/wib_data/ssarV2/ssarV2_INPUT.txt'
    project = Parse_input(input)

    project['srs']  = osr.SpatialReference()
    project['srs'].ImportFromEPSG(32633) # utm33n

    project['modis_days']  = Get_modis_days( project['start_year'],
                                             project['end_year'])

    project['climate_dir'] = '/space/wib_data/CLIMATE'

    # project['prj_directory'] = '/space/wib_data/ssarV2' # is this already loaded in parse_input?

    # MAKE HDF, SHP, CSV
    hdf_dir = os.path.join(project['prj_directory'],'HDF/')
    prj_mkdir(hdf_dir)

    shp_dir = os.path.join(project['prj_directory'],'SHP/')
    prj_mkdir(shp_dir)

    csv_dir = os.path.join(project['prj_directory'],'CSV/')
    prj_mkdir(csv_dir)

    project['hdf_dir'] = hdf_dir
    project['shp_dir'] = shp_dir
    project['csv_dir'] = csv_dir

    ### AGGREGATE NVE AND METNO CLIMATE DIRECTORIES TO HDF5, WHOLE NORWAY
    from Tools.CLIMATE_aggregator import Aggregate_climate_grids
    # todo done Aggregate_climate_grids(project)

    ### CLIMATE shapefile grid + idx
    pre = project['prj_name']
    cds = 'sd'
    climate_fn          = '/space/wib_data/CLIMATE/ssarV2_sd.hdf5'
    climate_shp          = os.path.join( project['shp_dir'],
                                       pre + 'climate' )
    climate_params         = Tools.SPATIAL_tools.Parse_extents( climate_fn )
    climate_params['outf'] = climate_shp
    climate_params['idx'] = 'yes, please'
    print 'Creating '+pre+'climate.shp from '+os.path.split(climate_fn)[1]
    # todo done Tools.SPATIAL_tools.Mk_polygrid(climate_params )

    # define CLIMATE shapefile as project['aoi']
    project['aoi']     = Parse_extents(climate_shp+'.shp')

    ### MODIS hdf5 subset (climate NORWAY extents)
    # project['modis_tiles'] = Check_mod_tiles(modis_tile_fn=project['modis_tile_fn'],**project['aoi'])

    project['modis_tiles'] = ['h18v01', 'h18v02', 'h18v03', 'h19v01', 'h19v02', 'h19v03']

    for dset in project['modis'].keys():
        for dnum in project['modis'][dset].keys():
            Mod2hdf( project, dset, dnum )

    ### MODIS shapefile grid
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

    ### REPROJECT + IDX MODIS shapefile grid
    pre  = project['prj_name']+'_' # start to all the project files
    # Reproject and index MODIS shapefile (to land cover projection)
    modis_shp          = os.path.join( project['shp_dir'],
                                       pre + 'modis' )
    mod_tx             = os.path.join( project['shp_dir'],
                                       pre+'modis_reprj' )
    mod_reprj = {'dst_dsn':mod_tx, 'src_dsn':modis_shp,
                 'dst_srs':project['srs'],
                 'fields':['ctr_x','ctr_y','x_ind','y_ind']}
    mod_reprj['area'] = True
    print 'Reprojecting '+pre+'modis.shp to project reference system and creating rtree index'
    Tools.SPATIAL_tools.Reprj_and_idx( **mod_reprj)

    ### intersections:
    # todo major
    # for given unique, unrepeated block (25x25) of MODIS cells:
        # for those cells which intersect landcover data:
            # intersect modis features, climate features, lc features
            # export to csv
