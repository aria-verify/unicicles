"""
um_landmask_to_regional

EXPERIMENTAL: part of being able to change landsea mask in the ice
coupling workflow. Cuts global UM mask file down to a regional version
that corresponds to the specified ice sheet domain

Called as a task directly from the suite, lets us feed the (possibly changed) 
mask generated from the NEMO coastlines to the regional CAP instances
"""

import sys
import numpy as np
import mule
import common_mule_rss as mule_rss
from common_arg_to_file_exist import arg_to_file_exist
from common_params_and_constants import ncols_AIS_n96, ncols_GrIS_n96, \
    nrows_AIS_n96, nrows_GrIS_n96, xoff_AIS_n96, xoff_GrIS_n96, \
    yoff_AIS_n96, yoff_GrIS_n96

def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--region", help="input, which ice sheet")
    parser.add_argument("--flip_array",
                        help="input, flip array N-S (not meta)",
                        action="store_true")
    parser.add_argument("--binary",
                        help="input, round fractional mask to binary",
                        action="store_true")
    parser.add_argument("--input_global", help="input, name of global file")
    parser.add_argument("--output_regional",
                        help="output, name of regional file")
    args = parser.parse_args()

    region = "GrIS"
    if args.region is not None:
        region = args.region

    flip_array = False
    if args.flip_array:
      #flip the array rows of the mask to run from N->S or vice versa
      flip_array = True

    binary = False
    if args.binary:
      #transform a fractional mask into binary one
      binary = True

    err = 0
    input_global, err = arg_to_file_exist(args.input_global, err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)
    output_regional, err = arg_to_file_exist(args.output_regional, err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)

    return input_global, output_regional, region, flip_array, binary

if __name__ == "__main__":

    # Read the command line arguments
    infile, outfile, region, flip, binary = parse_commandline()

    fieldg = mule.load_umfile(infile)

    #add this header, empty, to make a valid ancil file
    fieldg.level_dependent_constants \
        = mule.ff.FF_LevelDependentConstants(
            np.zeros([fieldg.integer_constants.num_levels, 4]))

    # get grid resolution from the input file
    nx_orig = fieldg.integer_constants.raw[6]
    ny_orig = fieldg.integer_constants.raw[7]
    dx = fieldg.fields[0].raw[62]
    dy = fieldg.fields[0].raw[60]

    # load  up the parameters for the ice sheet domain we want to cut out
    if region == "AIS":
        # Cutout Ant from N96e
        nx = ncols_AIS_n96
        ny = nrows_AIS_n96
        x_offset = xoff_AIS_n96
        y_offset = yoff_AIS_n96
    elif region == "GrIS":
        # Cutout GrIS from N96e (-1 at each edge to keep CAP domain from
        # failing/wrapping?)
        nx = ncols_GrIS_n96
        ny = nrows_GrIS_n96
        x_offset = xoff_GrIS_n96
        y_offset = yoff_GrIS_n96
    else:
        print("Don't know parameters for region", region, " exiting")
        sys.exit(1)


    # Report what we've got
    #print("X", nx, dx, fieldg.real_constants.raw[4] + (x_offset * dx),
    #      fieldg.real_constants.raw[4] + (x_offset * dx) + nx * dx)
    #print("Y", ny, dy, fieldg.real_constants.raw[3] + (y_offset * dy),
    #      fieldg.real_constants.raw[3] + (y_offset * dy) + ny * dy)

    # Change from global to non-global latlon grid without wrap
    if (nx < fieldg.integer_constants.raw[6]) or \
       (ny < fieldg.integer_constants.raw[7]):
        print("Setting to non-global grid type")
        fieldg.fixed_length_header.raw[4] = 3

    # update the headers with the new grid dimensions
    # Switch to newer version of coding the UM version
    fieldg.fixed_length_header.raw[12] = 901
    # 1 should always be MDI
    fieldg.fixed_length_header.raw[1] = -32768

    fieldg.integer_constants.raw[6] = nx
    fieldg.integer_constants.raw[7] = ny

    fieldg.real_constants.raw[3] = fieldg.real_constants.raw[3] + \
                                   (y_offset * dy)
    fieldg.real_constants.raw[4] = fieldg.real_constants.raw[4] + \
                                   (x_offset * dx)


    # cut down the fields themselves, modify the lookup headers to match the new grid at the same time
    for i in range(len(fieldg.fields)):
        fieldg.fields[i].raw[17] = 3
        fieldg.fields[i].raw[18] = ny
        fieldg.fields[i].raw[19] = nx
        fieldg.fields[i].raw[59] = fieldg.fields[i].raw[59] + (y_offset * dy)
        fieldg.fields[i].raw[61] = fieldg.fields[i].raw[61] + (x_offset * dx)

        datag = fieldg.fields[i].get_data()
        if flip:
            print("Flipping input array NS, apparently it was generated wrong")
            datag = datag[-1:0:-1, :]
        if binary:
            print("I think I've been given a fraction landmask, make it binary")
            # By eyeballing what I currently use the binary mask is around
            # this cf fractional!?
            fieldg.fields[0].raw[23] = 38
            THRESHOLD = 0.5
            datag[np.where(datag >= THRESHOLD)] = 1.
            datag[np.where(datag < THRESHOLD)] = 0.
            datag = datag.astype(int)
            fieldg.fields[i].raw[42] = 30

        # Write as logical data
        fieldg.fields[i].raw[38] = 9011111
        fieldg.fields[0].raw[39] = 3

        # This is the line that actually cuts down the global data to the regional domain
        fieldg.fields[i] \
            = mule_rss.overwrite_data(fieldg.fields[i],
                                      datag[y_offset:y_offset + ny,
                                            x_offset:x_offset + nx])

    # write out the output
    fieldg.to_file(outfile)
