#!/usr/bin/env python

"""
unicicles_cap_global_to_um

Reads in regional CAP orography files and Glint icecouple outputs from
previous stages, unifies them into global fields and inserts them into
the UM restart file that will be used by the next cycle of the 
atmosphere/land model

Called by suite as ice_sheets_ancils_to_atmos
"""

import sys
from mule import load_umfile
import numpy as np
from netCDF4 import Dataset
from cap_orog_to_um import insert_cap_ice_orog
from tile_fracs_to_um import fix_polarrow_tilefrac, insert_tile_fracs
from um_reset_icesnowpack import reset_icetile_snowpack
from tsurf_elev_to_um import insert_tsurf_elev
from interface_fortran_jules_snowpack import adjust_snow
from um_relandmask_restart import new_lsm
from regional_cap_ancil_to_global import splice_cap_ancil
from regional_icecouple_to_global import splice_icecouple
from common_arg_to_file_exist import arg_to_file_exist
from common_params_and_constants import stashcode_frac_surf_type

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


def save_originals(um_dump):
    from common_um_to_np import um_to_np_3d
    frac_elev0, _ = um_to_np_3d(um_dump, stashcode_frac_surf_type, elev=True)

    return frac_elev0


def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--um_input",
                        help="input, name of UM dump to be modified")
    parser.add_argument("--um_ref_input",
                        help="input, name of source of UM reference fields")
    parser.add_argument("--toice_input",
                        help="input, name of file going to Unicicles")
    parser.add_argument("--um_output", help="output, name of modified UM dump")
    parser.add_argument("--coupling_period",
                        help="climate-ice coupling period in secs")
    parser.add_argument("--do_snow_shuffling",
                        help="(re)enable non-ice snow shuffling, in "
                        "coupling, experimental",
                        action="store_true")
    parser.add_argument("--do_temp_boundary",
                        help="(re)enable passing BISICLES surface "
                        "temperature to the botton of ice tile snowpacks, "
                        "experimental",
                        action="store_true")
    parser.add_argument("--input_new_landseamask",
                        help="alter UM land area if CAP has given us a new "
                        "ice shelf land sea mask, experimental")
    parser.add_argument("--input_orog_GrIS",
                        help="input, name of regional Greenland ancil")
    parser.add_argument("--input_orog_AIS",
                        help="input, name of regional Antarctic ancil")
    parser.add_argument("--output_orog_global",
                        help="output, name of modified global ancil")
    parser.add_argument("--input_icecouple_GrIS",
                        help="input, name of regional Greenland icecouple file")
    parser.add_argument("--input_icecouple_AIS",
                        help="input, name of regional Antarctic icecouple file")
    parser.add_argument("--output_icecouple_global_out",
                        help="output, name of modified global icecouple file")
    parser.add_argument("--fix_polarrow_landfrac",
                        help="experimental, fix weird AIS polar row land "
                        "frac artefacts",
                        action="store_true")

    args = parser.parse_args()

    err = 0

    um_input, err = arg_to_file_exist(args.um_input, mandatory=True, err=err)
    if err != 0:
        sys.exit(err)
    um_ref_input, err = arg_to_file_exist(args.um_ref_input,
                                          mandatory=True, err=err)
    if err != 0:
        sys.exit(err)
    toice_input, err = arg_to_file_exist(args.toice_input,
                                         mandatory=True, err=err)
    if err != 0:
        sys.exit(err)
    um_output, err = arg_to_file_exist(args.um_output, io="out", err=err)
    if err != 0:
        sys.exit(err)
    GrIS_orog_filename, err = arg_to_file_exist(args.input_orog_GrIS,
                                                mandatory=False, err=err)
    if err != 0:
        sys.exit(err)
    AIS_orog_filename, err = arg_to_file_exist(args.input_orog_AIS,
                                               mandatory=False, err=err)
    if err != 0:
        sys.exit(err)
    global_orog_filename_out, err = arg_to_file_exist(args.output_orog_global,
                                                      io="out", err=err)
    if err != 0:
        sys.exit(err)
    GrIS_icecouple_filename, err = arg_to_file_exist(args.input_icecouple_GrIS,
                                                     mandatory=False)
    AIS_icecouple_filename, err = arg_to_file_exist(args.input_icecouple_AIS,
                                                    mandatory=False)
    global_icecouple_filename_out, err \
        = arg_to_file_exist(args.output_icecouple_global_out, io="out")

    if args.coupling_period is not None:
        coupling_period = float(args.coupling_period)
    else:
        err = 4

    do_snow_shuffling = False
    if args.do_snow_shuffling:
        do_snow_shuffling = True

    do_temp_boundary = False
    if args.do_temp_boundary:
        do_temp_boundary = True

    fix_polarrow_landfrac = False
    if args.fix_polarrow_landfrac:
        fix_polarrow_landfrac = True

    new_landseamask = None
    if args.input_new_landseamask is not None:
      new_landseamask, err = arg_to_file_exist(args.input_new_landseamask,
                                               mandatory=True, err=err)

    if err > 0:
        parser.print_help()
        sys.exit(err)

    return um_input, um_ref_input, toice_input, um_output, \
        coupling_period, do_snow_shuffling, do_temp_boundary, \
        new_landseamask, GrIS_orog_filename, AIS_orog_filename, \
        global_orog_filename_out, GrIS_icecouple_filename, \
        AIS_icecouple_filename, global_icecouple_filename_out, \
        fix_polarrow_landfrac


# Process filenames etc from the command line
if __name__ == "__main__":

    # Read the command line arguments
    um_input, um_ref_input, toice_input, um_output, \
        coupling_period, do_snow_shuffling, do_temp_boundary, \
        new_landseamask, GrIS_orog_filename, AIS_orog_filename, \
        global_orog_filename_out, GrIS_icecouple_filename, \
        AIS_icecouple_filename, global_icecouple_filename_out, \
        fix_polarrow_landfrac = parse_commandline()

    final_orog = splice_cap_ancil(GrIS_orog_filename, AIS_orog_filename,
                                  um_ref_input)

    final_orog = disable_mule_validators(final_orog)

    final_orog.to_file(global_orog_filename_out)

    splice_icecouple(GrIS_icecouple_filename, AIS_icecouple_filename,
                     global_icecouple_filename_out)

    # Read in the UM dump to be modified.
    print("Reading UM input ", um_input)
    umf = load_umfile(um_input)
    um_dump = umf.copy()
    #19-11-24: I totally don't recall why I made it remove all the fields from
    #the copy and then add them back in manually!
    um_dump.fields = []
    um_dump = disable_mule_validators(um_dump)
    for i in range(umf.fixed_length_header.raw[152]):
        um_dump.fields.append(umf.fields[i])

    # Read in the UM reference fields
    print("Reading UM reference input ", um_ref_input)
    umf = load_umfile(um_ref_input)
    um_ref_dump = umf.copy()
    um_ref_dump.fields = []
    for i in range(umf.fixed_length_header.raw[152]):
        um_ref_dump.fields.append(umf.fields[i])

    # Read in the *icecouple fields passed to/from unicicles this cycle
    print("Reading toice input ", toice_input)
    toice_file = Dataset(toice_input, 'r')

    print("Reading fromice input ", global_icecouple_filename_out)
    fromice_file = Dataset(global_icecouple_filename_out, 'r+')

    if fix_polarrow_landfrac:
        print("(and fixing the AIS polar row land fractions)")
        fromice_file = fix_polarrow_tilefrac(fromice_file)

    # Read in the CAP-processed orography fields
    print("Reading CAP input ", global_orog_filename_out)
    fromcap_file = load_umfile(global_orog_filename_out)

    # Save the originals of some fields for use later.
    # frac_elev0 is needed to do snow_shuffling modifications, we have to
    # update the original ice area fractions before we can shuffle
    frac_elev0 = save_originals(um_dump)

    # Apply the new landsea mask from the file specified on the command line.
    if new_landseamask is not None:
        print("EXPERIMENTAL WARNING!-----------")
        print("Putting in new land sea mask")
        print("EXPERIMENTAL WARNING!-----------")
        print("Reading mask", new_landseamask)
        lsm_file = Dataset(new_landseamask, 'r')
        um_dump = new_lsm(um_dump, lsm_file)

    # Insert the orography fields from the CAP/ANTS output
    print("Putting in CAP orog fields")
    um_dump = insert_cap_ice_orog(um_dump, um_ref_dump,
                                  fromcap_file, fromice_file)

    # Update tile fractions based on area fractions from UniCiCles
    print("Update tile fractions")
    um_dump = insert_tile_fracs(um_dump, um_ref_dump, fromice_file)

    # Reset snowpack mass on ice elevation tiles to the default
    print("Resetting ice tile snow mass")
    um_dump = reset_icetile_snowpack(um_dump, toice_file, coupling_period)

    # Adjust snowpacks on non-ice elevated tiles to account for changes in
    # area fraction or where UniCiCles has initialised new ice from deep
    # snowpacks, or deglaciating ice areas that have thinned to near 0

    if do_snow_shuffling:
        print("EXPERIMENTAL WARNING!-----------")
        print("SNOW SHUFFLING ON")
        print("Adjusting the rock elev snow on the new fracs")
        print("This ALWAYS goes wrong eventually! FIX ME!")
        print("Currently this option JUST SHOWS WHAT THE CORRECTION WOULD BE!")
        print("We haven't written code to modify the existing SMB field ",
              "with the change yet!")
        print("EXPERIMENTAL WARNING!-----------")
        um_dump, smb_from_nisnow = adjust_snow(
            um_dump, fromice_file, frac_elev0)

        if np.sum(smb_from_nisnow) > 0.0:
            print("Modify the SMB file with this overflow before we run ",
                  "BISICLES again")
            import cf
            f = cf.read(toice_input)
            g = f[0].copy()
            g.set_property('long_name', 'smb_from_nisnow')
            g.id = 'smb_from_nisnow'
            g.nc_set_variable('smb_from_nisnow')
            for prop in ['stash_code', 'um_stash_source']:
                if g.has_property(prop):
                    g.del_property(prop)

            g.set_data(cf.Data(smb_from_nisnow / coupling_period,
                               units='kg/m2/s'))
            f.append(g)
            cf.write(f, toice_input)

    else:
        print("SNOW SHUFFLING OFF")
        print("Not adjusting elev_rock snowpacks")
        print("This implies snow mass conservation problems when the ice ",
              "sheet moves")
        print("but the alternative is complex and currently broken - see ",
              "code under do_snow_shuffling")

        # elev_tile subsurface temperatures
    if do_temp_boundary:
        print("EXPERIMENTAL WARNING!-----------")
        print("Updating sub-icesnowpack boundary conditions")
        print("we haven't tried configurations of UKESM with temperature ",
              "updating")
        print("EXPERIMENTAL WARNING!-----------")
        um_dump = insert_tsurf_elev(um_dump, um_ref_dump, fromice_file)

    else:
        # BISICLES not using temps at the moment - just getting 0s back
        print("NOT updating sub-icesnowpack boundary conditions")
        print("we haven't tried configurations of UKESM with temperature ",
              "updating")

    # Output the modified dump
    print("Writing UM output")
    um_dump.to_file(um_output)
