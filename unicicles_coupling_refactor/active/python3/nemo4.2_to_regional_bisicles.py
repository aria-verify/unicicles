#!/usr/bin/env python

"""
nemo_to_regional_bisicles.py

For each BISICLES region to be used, cylc will run a copy of this task,
regridding the global output from ocean_average_for_ice into ther local
domain. It creates a BISICLES hdf5 boundary condition file.

nemo4.2 version uses BISICLES-NEMO ASE output variable names and sign
conventions

Called by suite as ocean_to_regional_ice_sheet
"""

import os
import sys
from netCDF4 import Dataset
import numpy as np
import cf
from common_arg_to_file_exist import arg_to_file_exist
from common_params_and_constants import days_in_year, secs_in_day


def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input, name of nemo input nc file")
    parser.add_argument("--output",
                        help="input, name of bisicles output hdf file")
    parser.add_argument("--region", help="input, name of bisicles region")
    args = parser.parse_args()

    err = 0
    input, err = arg_to_file_exist(args.input, mandatory=True, err=err)
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

    # NEMO melt rate is positive for melting and negative for freezing,
    # BISICLES has the opposite sign convention, so multiply NEMO melt by -1
    # NEMO is kg/s, BISICLES wants m/yr
    water_unit_factor = -1 * secs_in_day * days_in_year / 1e3
    # NEMO is J/s, BISICLES wants J/yr
    heat_unit_factor = secs_in_day * days_in_year

    h = cf.read(nctgridfile)

    melt_water = h.select_field(
        'long_name=Ice shelf fresh water flux ( from isf to oce )')\
                  .array.squeeze()
    melt_heat = h.select_field(
        'long_name=Ice shelf ocean  heat flux ( from isf to oce )')\
                 .array.squeeze()

    # NEMO output just has zero where there is no ocean (eg grounded ice)
    # We don't try to mask this out cf regridding ends up putting missing_data
    # flags into the real cavity even where there is mostly good information.
    # For now, we just let cf treat the no ocean 0s as if they were valid info
    # to regrid

    # Check for array size - we used to optionally do time averaging, but
    # that's now in process_ocean
    if melt_water.ndim == 3:
        print("Doing time average", np.shape(melt_water))
        print("I shouldn't have to do this - have I picked up a 3d volume",
              "field by accident?!")
        sys.exit(1)

    g = cf.read(bike_ncgridfile)[0]

    print(" regridding melt field")
    f = cf.read(nemo_ncgridfile)[0]
    f.set_data(cf.Data(melt_water, units='kg/m2/s'))
    w_regrid = f.regrids(g, src_cyclic=True, method='linear')

    print(" regridding heat field")
    f = cf.read(nemo_ncgridfile)[0]
    f.set_data(cf.Data(melt_heat, units='kg/m2/s'), axes=('Y', 'X'))
    h_regrid = f.regrids(g, src_cyclic=True, method='linear')

    heat = h_regrid.array
    melt = w_regrid.array

    # Do the unit conversion, get rid of any missing data flags in the regridded array
    try: 
      heat.fill_value = 0.
      heat = heat.filled() * heat_unit_factor
      melt.fill_value = 0.
      melt = melt.filled() * water_unit_factor
    except AttributeError:
      heat = heat * heat_unit_factor
      melt = melt * water_unit_factor

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
    ncfile_out.createDimension('x', x_bike.shape[1])
    ncfile_out.createDimension('y', x_bike.shape[0])

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
