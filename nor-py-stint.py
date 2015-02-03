from Tools.ORG_tools import Parse_input,Get_modis_days

if __name__ == '__main__':

    nor_py_cmds = '''running below!'''
    input = '/space/wib_data/ssarV2/ssarV2_INPUT.txt'
    project = Parse_input(input)
    project['modis_days']  = Get_modis_days( project['start_year'],
                                             project['end_year'])
    ### AGGREGATE NVE AND METNO CLIMATE DIRECTORIES TO HDF5, WHOLE NORWAY
    from Tools.CLIMATE_aggregator import Aggregate_climate_grids
    Aggregate_climate_grids(project)
