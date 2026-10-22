#!/usr/bin/env python

"""
nemo_to_regional_bisicles.py

For each BISICLES region to be used, cylc will run a copy of this task,
regridding the global ocean/iceshelf outputs from ocean_average_for_ice into the local
domain. It creates a BISICLES hdf5 boundary condition file.

nemo3.6 version uses UKESM NEMO output variable names and sign conventions

Called by suite as ocean_to_regional_ice_sheet
"""

import os
import sys
from netCDF4 import Dataset
import numpy as np
import cf
#import mapping_class as mp
from common_arg_to_file_exist import arg_to_file_exist
from common_params_and_constants import days_in_year, secs_in_day


def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input, name of nemo input nc file")
    parser.add_argument("--output",
                        help="input, name of bisicles output hdf file")
    parser.add_argument("--region", help="input, name of ice sheet")
    args = parser.parse_args()

    err = 0
    input, err = arg_to_file_exist(args.input, mandatory=True, err=err)
    if err != 0:
        sys.exit(err)
    output, err = arg_to_file_exist(args.output, io="out", err=err)

    region = "GrIS"
    if args.region is not None:
        region = args.region

    if err > 0:
        parser.print_help()
        sys.exit(err)

    return input, output, region

if __name__ == "__main__":

    # Read the command line arguments
    nctgridfile, regrid_hdf5file, region = parse_commandline()

    regrid_ncfile = os.path.splitext(os.path.basename(regrid_hdf5file))[0] + \
                    ".nc"

    bike_ncgridfile = "cf_bikegridfile.nc"
    nemo_ncgridfile = "cf_gridfile.nc"

    #NEMO is kg/s, BISICLES wants m/yr
    water_unit_factor = secs_in_day * days_in_year / 1e3
    #NEMO is J/s, BISICLES wants J/yr
    heat_unit_factor = secs_in_day * days_in_year

    h = cf.read(nctgridfile)

    melt_water = h.select_field('ice_shelf_melting').array.squeeze()
    melt_heat = h.select_field('ice_shelf_heat_flux').array.squeeze()

    # NEMO output just has zero where there is no ocean (eg grounded ice)
    # We don't try to mask this out cf regridding ends up putting missing_data
    # flags into the real cavity even where there is mostly good information.
    # For now, we just let cf treat the no ocean 0s as if they were valid info
    # to regrid

    # Check for array size- we used to optionally do time averaging, but
    # that's now in process_ocean
    if melt_water.ndim == 3:
        print("3D meltwater array found. This routine shouldn't have to do time averaging",
              "- have I picked up a 3d volume field by accident?!")
        sys.exit(1)

    g = cf.read(bike_ncgridfile)[0]

    print(" regridding melt field")
    f = cf.read(nemo_ncgridfile)[0]
    f.set_data(cf.Data(melt_water, units='kg/m2/s'), axes=('Y', 'X'))
    w_regrid = f.regrids(g, src_cyclic=True, method='linear')

    print(" regridding heat field")
    f = cf.read(nemo_ncgridfile)[0]
    f.set_data(cf.Data(melt_heat, units='J/m2/s'), axes=('Y', 'X'))
    h_regrid = f.regrids(g, src_cyclic=True, method='linear')

    # Like BIKEtoNEMO, needs a rollaxis due to the array ordering when the lon,lat
    # information for BISICLES was generated. We need to switch x and y here to 
    # get it the way BISICLES wants it, and should check this step when new BISICLES
    # domain information is generated
    print("SWAPPING AXES OF REGRIDDED FIELD FOR BISICLES. THIS WAS NEEDED FOR THE ",
          "ORIGINAL UKESM1 AIS BISICLES GRID DEFINITION, IF YOU HAVE MOVED ON FROM ",
          "THAT CONFIG YOU SHOULD CHECK THAT YOUR REGRIDDED ICE SHELF INPUTS TO ",
          "BISICLES MATCH THE GEOGRAPHY YOU EXPECT") 
    heat = np.swapaxes(h_regrid.array, 0, 1)
    melt = np.swapaxes(w_regrid.array, 0, 1)

    # Do the unit conversion, get rid of any missing data flags in the regridded array
    heat.fill_value = 0.
    heat = heat.filled() * heat_unit_factor
    melt.fill_value = 0.
    melt = melt.filled() * water_unit_factor

    # Antony Siahaan wrote a clause here that adjusted collections of cells in the 
    # regridded melt field to conserve with respect to the original NEMO cell their
    # centres map to, using the precomputed bikegridmap_file.map2d file.  
    # Removed, since Robin found places where it made very large +ve and -ve
    # adjustments, preserving the mean but distorting the patterns. This might have
    # come from the calculation or errors in the mapping file. The idea should be
    # revisited, for code see earlier revisions in the repository.

    # write out a netCDF version of the regridded melt fields

    x_bike = np.load("bike_xcoords.dump", allow_pickle=True, fix_imports=True,
                     encoding='latin1')
    y_bike = np.load("bike_ycoords.dump", allow_pickle=True, fix_imports=True,
                     encoding='latin1')

    ncfile_out = Dataset(regrid_ncfile, 'w', format='NETCDF3_CLASSIC')
    ncfile_out.createDimension('x', x_bike.shape[0])
    ncfile_out.createDimension('y', x_bike.shape[1])

    x_nc = ncfile_out.createVariable('x', np.dtype('float64').char, ('x'))
    y_nc = ncfile_out.createVariable('y', np.dtype('float64').char, ('y'))

    x_nc[:] = x_bike[0, :]
    y_nc[:] = y_bike[:, 0]

    heat_nc = ncfile_out.createVariable('melt_heat',
                                        np.dtype('float64').char, ('y', 'x'))
    melt_nc = ncfile_out.createVariable('melt_water',
                                        np.dtype('float64').char, ('y', 'x'))
    heat_nc[:] = heat
    melt_nc[:] = melt

    ncfile_out.close()

    # Create the HDF file
    print(" Calling nctoamr")
    status = os.system(
        './nctoamr2d.ex %(regrid_ncfile)s %(regrid_hdf5file)s melt_water melt_heat' %
        locals())
    if status != 0:
        print("ERROR: failed to create HDF file")
        sys.exit(1)
