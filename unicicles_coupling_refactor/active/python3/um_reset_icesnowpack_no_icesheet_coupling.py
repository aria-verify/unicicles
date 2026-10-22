#!/usr/bin/env python
"""
um_reset_icesnowpack_no_icesheet_coupling, based on unicicles_cap_global_to_um

Modifies the UM restart file that will be used by the next cycle of the
atmosphere/land model. Without information from the evolving ice sheet, can
only reset the snow mass on icesheet tile surfaces (using the same routine as 
unicicles_cap_global_to_um does) so that SMB can continue to be diagnosed.

Called by suite as ice_sheets_to_atmos_RESETSNOW
"""

import sys
from mule import load_umfile
from netCDF4 import Dataset
from common_arg_to_file_exist import arg_to_file_exist
from um_reset_icesnowpack import reset_icetile_snowpack

# the usual headers to downgrade MULE's strict internal checks to warnings
# so we can work with a range of ancils and dumps from different eras
def validate_warn(self, filename=None, warn=True):
    self.validate_errors(filename=filename, warn=True)
    print("(no more than one warning is noted)")
    return


def disable_mule_validators(umf):
    umf.validate_errors = umf.validate
    print("OVERWRITING MULE'S VALIDATION FUNCTION")
    print("(the original is saved as self.validate_errors)")
    funcType = type(umf.validate)
    print("VALIDATION ERRORS DOWNGRADED TO WARNINGS")
    umf.validate = funcType(validate_warn, umf)
    return umf


def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--um_input",
                        help="input, name of UM dump to be modified")
    parser.add_argument("--toice_input",
                        help="input, name of file going to Unicicles")
    parser.add_argument("--um_output", help="output, name of modified UM dump")
    parser.add_argument("--coupling_period",
                        help="climate-ice coupling period in secs")
    args = parser.parse_args()

    err = 0
    um_input, err = arg_to_file_exist(args.um_input, mandatory=True, err=err)
    toice_input, err = arg_to_file_exist(args.toice_input, mandatory=True,
                                         err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)
    um_output, err = arg_to_file_exist(args.um_output, io="out", err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)

    if args.coupling_period is not None:
        coupling_period = float(args.coupling_period)
    else:
        err = 4
        print("ERROR: specify climate-ice coupling period")

    if err > 0:
        parser.print_help()
        sys.exit(err)

    return um_input, toice_input, um_output, coupling_period


# Process filenames etc from the command line
if __name__ == "__main__":

    """
    add last year's SMB back into the snowpacks on the icesheet tiles in the dump
    now that that mass change has been applied to the ice sheet, keeping a standard
    depth of snowpack on top of the ice sheet surface
    """

    um_input, toice_input, um_output, coupling_period = parse_commandline()

    # Read in the UM dump to be modified
    print("Reading UM input ", um_input)

    umf = load_umfile(um_input)
    um_dump = umf.copy()
    #19-11-24: I totally don't recall why I made it remove all the fields from
    #the copy and then add them back in manually!
    um_dump.fields = []
    um_dump = disable_mule_validators(um_dump)
    for i in range(umf.fixed_length_header.raw[152]):
        um_dump.fields.append(umf.fields[i])

    # Read in the wrapper-processed output from unicicles
    print("Reading toice input ", toice_input)
    toice_file = Dataset(toice_input, 'r')

    # 4. elev_ice snowpacks
    print("Resetting ice tile snow mass")
    um_dump = reset_icetile_snowpack(um_dump, toice_file, coupling_period,
                                     real_tile_only=True)

    # Output the modified dump
    print("Writing UM output")
    um_dump.to_file(um_output)
