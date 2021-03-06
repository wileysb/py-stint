# '#' denotes comments: remainder of the line is not parsed
# ',' is used to identify ERA and MODIS dataset lines
# ',' should teherefore be used in every ERA and MODIS dataset line
# ',' should NOT be used in any other line

# spaces between words and symbols are meaningless, and will be stripped
# same goes for empty lines; they'll be ignored!
# however, spaces must not appear within variables.
# use '_' instead!


# PROJECT INFO #########################################################
start_year      = YY   # starts 1 January, 20YY (late Feb if '00)
end_year        = YY   # runs through 31 December, 20YY
# mod_per_region  = 500
# lcm_thresh      = 10000

# PATHS  ###############################################################
aoi           = /path/to/aoi.shp # shp or raster directory
modis_tile_fn = /path/to/MODIS_tiles.shp
prj_name      = project_name # no spaces
prj_directory = /path/to/project_output_directory/ # this directory should exist before beginning processing
modis_dir     = /path/to/MODIS_ARCHIVE/
era_dir       = /path/to/CLIMATE/ERA/


# ERA datasets #########################################################
# short_name  = file1,[file2,file3] ( in order of dates contained )    #
example       = example_jan2000_Dec2013.nc,
sd            = snowdepth_Jan2000_Dec2006.nc, snowdepth_Jan2007_Dec2013.nc
t2m           = T2m_Feb2000_Dec2004.nc, T2m_Jan2005_Dec2008.nc, T2m_Jan2009_Dec2013.nc


# MODIS datasets #######################################################
# dset   =  dnum:short_name,  [dnum:short_name,  dnum:short_name]      #
example  =     2:example_subdataset,
MCD43A2  =     0:BSA_quality,   2:BSA_ancill,    3:BSA_band
MCD43A3  =     7:BSA_vis,       8:BSA_nir,       9:BSA_sw
MOD10A2  =     0:mse,           1:dsc


# Landcover datasets ###################################################
# lc_short_name  =  filename.tif OR shapefile field name               #
# short_names should all start with 'lc_' for 'land cover'             #
lc_category = code_06
# lc_age      = shp_age_field