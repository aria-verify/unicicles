"""
regional_icecouple_to_global

Subroutine for merging fields from the Glint icecouple files
from the regional ice sheet domains into unified global ones.

Used by unicicles_cap_global_to_um
"""

import sys
import os
import numpy as np
from common_arg_to_file_exist import arg_to_file_exist


def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_GrIS",
                        help="input, name of regional Greenland ancil")
    parser.add_argument("--input_AIS",
                        help="input, name of regional Antarctic ancil")
    parser.add_argument("--output_global",
                        help="output, name of modified global ancil")
    args = parser.parse_args()

    err = 0

    GrIS_filename, err = arg_to_file_exist(args.input_GrIS, mandatory=False,
                                           err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)
    AIS_filename, err = arg_to_file_exist(args.input_AIS, mandatory=False,
                                          err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)
    global_filename_out, err = arg_to_file_exist(args.output_global, io="out",
                                                 err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)

    return GrIS_filename, AIS_filename, global_filename_out


def splice_icecouple(GrIS_filename, AIS_filename, global_filename_out):
    from netCDF4 import Dataset

    found_file = 0
    fields_to_splice = ["surface_elevation", "cell_calving_flux",
                        "tile_land_fraction", "tile_ice_fraction",
                        "nonice_snowdepth", "ice_stemp"]

    # we need info from either GrIS and/or AIS, but at least one
    if GrIS_filename is not None:
        gris_file = Dataset(GrIS_filename)
        found_file += 1

    if AIS_filename is not None:
        ais_file = Dataset(AIS_filename)
        found_file += 2

    if found_file == 0:
        print("Must specify at least one of GrIS, AIS i.icecouple files",
              GrIS_filename, AIS_filename)
        sys.exit(1)

    # if only info from one ice sheet, then the single regional icecouple
    # file we have can be used as the global one
    if found_file == 1:
        os.system("cp " + GrIS_filename + " " + global_filename_out)

    if found_file == 2:
        os.system("cp " + AIS_filename + "  " + global_filename_out)

    # if info from both ice sheets, then we need to combine the regional files
    # into a global one. 
    if found_file == 3:

        out_file = Dataset(global_filename_out, 'w')

        # start setting up the global file by copying structure from the 
        # pre-existing regional ones. 

        # Echo global dimensions
        for dimname in gris_file.dimensions:
            dim = gris_file.dimensions[dimname]
            out_file.createDimension(dimname, len(dim))

        # Determine where in the 2D lat-lon field each ice sheet can provide valid data
        # for each field. The total coverage is 1 for grid cells entirely in the ISM domain
        # but may be less for those that overlap the boundary
        coverage_gris = np.sum(gris_file.variables['tile_land_fraction'][:],
                               axis=0)
        coverage_ais = np.sum(ais_file.variables['tile_land_fraction'][:],
                              axis=0)
        valid = np.where((coverage_gris > 0) | (coverage_ais > 0))

        # Echo individual variables, with attributes
        for varname in gris_file.variables:
            print(varname)
            var = gris_file.variables[varname]
            var_ais = ais_file.variables[varname]

            # construct a global ice volume and make new variables
            # for each separate ice sheet
            if varname == "total_ice_volume":
                outfield = out_file.createVariable(varname, var.dtype,
                                                   var.dimensions)
                outfield[:] = var[:] + var_ais[:]
                outfield = out_file.createVariable('total_ice_volume_GrIS',
                                                   var.dtype, var.dimensions)
                outfield[:] = var[:]
                outfield = out_file.createVariable('total_ice_volume_AIS',
                                                   var.dtype, var.dimensions)
                outfield[:] = var_ais[:]

            # construct a global change in ice volume and make new variables
            # for each separate ice sheet
            elif varname == "change_in_ice_volume":
                outfield = out_file.createVariable(varname, var.dtype,
                                                   var.dimensions)
                outfield[:] = var[:] + var_ais[:]
                outfield = out_file.createVariable('change_in_ice_volume_GrIS',
                                                   var.dtype, var.dimensions)
                outfield[:] = var[:]
                outfield = out_file.createVariable('change_in_ice_volume_AIS',
                                                   var.dtype, var.dimensions)
                outfield[:] = var_ais[:]

            # loop over each of the fields we just want to add together
            elif varname in fields_to_splice:
                outfield = out_file.createVariable(varname, var.dtype,
                                                   var.dimensions)
                merge = np.zeros_like(var[:])

                # check to see if we're dealing with a field on elevation tiles
                l_tile = len(np.shape(var[:])) == 3

                if l_tile:
                    # tiled field, take each elevation level separately and combine the contributions
                    # from each ice sheet in their valid areas, weighted by the coverage fraction 
                    for tile in range(np.shape(var[:])[0]):
                        var_2d = var[:][tile, :, :]
                        var_ais_2d = var_ais[:][tile, :, :]

                        merge_2d = np.zeros_like(var_2d)
                        merge_2d[valid] = ((var_2d[valid] *
                                            coverage_gris[valid]) +
                                           (var_ais_2d[valid] *
                                            coverage_ais[valid])) / \
                            (coverage_gris[valid] + coverage_ais[valid])
                        merge[tile, :, :] = merge_2d

                else:
                    # non-tiled field, just combine the contributions from each ice sheet in their 
                    # valid areas, weighted by the coverage fraction 
                    merge[valid] = ((var[:][valid] * coverage_gris[valid]) +
                                    (var_ais[:][valid] *
                                     coverage_ais[valid])) / \
                        (coverage_gris[valid] + coverage_ais[valid])

                outfield[:] = merge
            else:
                # fields not dealt with by the above cases contain generic global grid info
                # just copy acros the GrIS version
                outfield = out_file.createVariable(varname, var.dtype,
                                                   var.dimensions)
                outfield[:] = var[:]
                outfield.setncatts({k: var.getncattr(k)
                                    for k in var.ncattrs()

        # this effectively writes the file out
        out_file.close()


if __name__ == "__main__":

    # Read the command line arguments
    GrIS_filename, AIS_filename, global_filename_out = parse_commandline()

    # Create a global field from regional icecouple files
    splice_icecouple(GrIS_filename, AIS_filename, global_filename_out)
