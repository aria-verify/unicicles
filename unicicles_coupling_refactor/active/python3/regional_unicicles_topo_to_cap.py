#!/usr/bin/env python

"""
regional_unicicles_topo_to_cap

The plot.*hdf5 output directly from BISICLES is read in and used to create
a high-resolution UM-format regional ancillary file containing the current
icesheet topography

Called by suite as regional_ice_sheet_topo_to_ancil
"""

import sys
import os
import time
import numpy as np
from common_arg_to_file_exist import arg_to_file_exist
from common_params_and_constants import dx_regrid_coarse, dy_regrid_coarse, \
    dx_regrid_fine, dy_regrid_fine, \
    epsg_ASE, epsg_GrIS, epsg_global, x0_AIS, x0_GrIS, \
    xcyclic_AIS, xcyclic_GrIS, xmax_AIS, xmax_GrIS, xmin_AIS, xmin_GrIS, \
    y0_AIS, y0_GrIS, ymax_AIS, ymax_GrIS, ymin_AIS, ymin_GrIS


# some versions of UM mule library have over zealous internal consistency checks
# that fail on valid, older format files. The next two functions are needed to 
# downgrade mule validation errors to warnings for those versions

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
                        help="input, name of UM ancil to be modified")
    parser.add_argument("--hdf5_input", help="input, name of unicicles hdf5")
    parser.add_argument("--um_output",
                        help="output, name of modified UM ancil")
    parser.add_argument("--region", help="GrIS or AIS?")
    args = parser.parse_args()

    err = 0
    um_input, err = arg_to_file_exist(args.um_input, err=err)
    hdf5_input, err = arg_to_file_exist(args.hdf5_input, err=err)
    um_output, err = arg_to_file_exist(args.um_output, io="out", err=err)

    region = "GrIS"
    if args.region is not None:
        region = args.region

    if err > 0:
        parser.print_help()
        sys.exit(err)

    return um_input, hdf5_input, um_output, region


def extract_bisicles_topography(hdf5_input, region="GrIS"):
    from amrfile import io as amrio
    from pyproj import Transformer

    level = 0  # level of grid refinement, 0 = lowest
    order = 0  # interpolation order, 0 for piecewise constant, 1 for linear

    # Load BISICLES hdf5
    amrID = amrio.load(hdf5_input)
    lo, hi = amrio.queryDomainCorners(amrID, level)
    nx = hi[0] - lo[0] + 1
    ny = hi[1] - lo[1] + 1

    # Extract cartesian grid and variables. x and y are the box centres
    # distance from origin
    x_1d, y_1d, z_surf = amrio.readBox2D(amrID, level, lo, hi, "Z_surface",
                                         order)

    # Make 2D meshes of x,y as coordinate arrays to match the data
    dx = x_1d[1] - x_1d[0]
    dy = y_1d[1] - y_1d[0]
    y_2d, x_2d = np.mgrid[y_1d[0]:y_1d[-1] + dy:dy, x_1d[0]:x_1d[-1] + dx:dx]

    # select projection defintions for pyproj to use
    outProj = epsg_global

    # If you set up a new BISICLES domain, check that coastlines etc look
    # right either side of the transformation - these offsets have changed
    # in the past, odd reflections have been applied etc.
    if region == "GrIS":
        inProj = epsg_GrIS
        x0 = x0_GrIS
        y0 = y0_GrIS
    elif region == "AIS":
        inProj = epsg_ASE
        x0 = x0_AIS
        y0 = y0_AIS
    else:
        print("extract_bisicles_topography: region unknown ", region)

    # Transform from the BISICLES grid info to lons and lats
   
    # Original AIS definition needed longitude reflection
    # if region == "AIS": mfact = -1
    mfact = 1.

    t = Transformer.from_crs(inProj, outProj, always_xy=True)
    xl_2d, yl_2d = t.transform((x_2d + x0) * mfact, (y_2d + y0))

    # Original AIS definition needed *additional* longitude offset 
    # if region == "AIS": steph_offset=90.
    steph_offset = 0.

    xl_2d = xl_2d + steph_offset

    return z_surf, xl_2d, yl_2d, x_1d, y_1d


def construct_CF(field, lon, lat, x_1d, y_1d):

    import cf

    """
    Construct a minimal CF Field from a numpy array and some lon-lat data
    so that cf-python can then use it in regridding routines
    """

    field_CF = cf.Field()

    dimx = cf.DimensionCoordinate(data=cf.Data(x_1d, 'm'))
    field_CF.set_construct(cf.DomainAxis(size=len(x_1d)), key="X")
    field_CF.set_construct(dimx, axes="X")

    dimy = cf.DimensionCoordinate(data=cf.Data(y_1d, 'm'))
    field_CF.set_construct(cf.DomainAxis(size=len(y_1d)), key="Y")
    field_CF.set_construct(dimy, axes="Y")
    lats = cf.AuxiliaryCoordinate(data=cf.Data(lat, 'degrees_north'),
                                  properties={'standard_name': 'latitude'})
    field_CF.set_construct(lats, axes=('Y', 'X'))
    lons = cf.AuxiliaryCoordinate(data=cf.Data(lon, 'degrees_east'),
                                  properties={'standard_name': 'longitude'})
    field_CF.set_construct(lons, axes=('Y', 'X'))

    field_ma = np.ma.masked_equal(field, 0.)
    field_CF.set_data(cf.Data(field_ma, units='m'), axes=('Y', 'X'))
    field_CF.set_property('_FillValue', -99.)

    return field_CF


def regrid_to_GLOBE30(field_CF, region="GrIS"):

    import cf

    """
    Regrid the CF Field of BISICLES topography to the UM CAP ancil
    high res lat-lon grid. 

    Originally this procedure used two separate regriddings - cf-python
    ESMF regridding is very quick for the basic transform onto a regular 
    lat-lon grid but memory limitations meant it couldn't handle going 
    to the 30min resolution needed. We did a coarser regrid in cf-python, 
    then call cdo (which is very slow for the first regrid) to refine it.

    New releases of cf-python *should* be able to cope with the memory issue
    so this function was modified to work with a single stage regrid. In
    practice we've still not had this work in anger on the AIS domain on the
    platforms available to us, so the "one stage" lines are here but commented 
    out and the two stage is procedure is still in force. KEEP UNDER REVIEW
    """

    # Construct the regular grid in CF

    #set region bounds
    if region == "GrIS":
        ymin = ymin_GrIS
        ymax = ymax_GrIS
        xmin = xmin_GrIS
        xmax = xmax_GrIS
        xcyclic = xcyclic_GrIS
    elif region == "AIS":
        ymin = ymin_AIS
        ymax = ymax_AIS
        xmin = xmin_AIS
        xmax = xmax_AIS
        xcyclic = xcyclic_AIS

    # set target resolution for cf-python regrid
    # One stage regrid /still/ running out of memory for AIS
    # dy = dy_regrid_fine
    # dx = dx_regrid_fine
    # Coarser initial regrid for 2-stage still required right now
    dy = dy_regrid_coarse
    dx = dx_regrid_coarse

    x_1d = np.arange(xmin + dx / 2., xmax - dx / 2., dx)
    y_1d = np.arange(ymin + dx / 2., ymax - dx / 2., dy)

    g = cf.Field()

    dimx = cf.DimensionCoordinate(data=cf.Data(x_1d, 'degrees_east'),
                                  properties={'axis': 'X',
                                              'standard_name': 'longitude'})
    dimy = cf.DimensionCoordinate(data=cf.Data(y_1d, 'degrees_north'),
                                  properties={'axis': 'Y',
                                              'standard_name': 'latitude'})
    g.set_construct(cf.DomainAxis(size=len(x_1d)), key="X")
    g.set_construct(dimx, axes='X')
    g.set_construct(cf.DomainAxis(size=len(y_1d)), key="Y")
    g.set_construct(dimy, axes='Y')
    g.set_data(cf.Data(np.zeros([dimy.shape[0], dimx.shape[0]])),
               axes=('Y', 'X'))

    # Do the cf-python regrid from BISICLES to the regular grid
    t0 = time.time()
    field_REGRID = field_CF.regrids(g, dst_cyclic=xcyclic, method='linear',
                                    src_axes={'X': 'key%X', 'Y': 'key%Y'})
    print("coarse cf-python regrid: %.1fs" % (time.time() - t0))

    # For the one-stage regrid, return the array now
    # field_REGRID_array = field_REGRID.array

    #######################################
    # Two_stage extension - REMOVE WHEN NEWER CF AVAILABLE AND USE FINER
    # DX, DY ABOVE
     
    # write out the coarsely regridded field
    cf.write(field_REGRID, 'cf_0.083deglatlonfield.nc', fmt="NETCDF3_CLASSIC")

    # Write out cdo grid definitions for coarse and fine grids
    gf = open('gridfile-Reg0.083', 'w')
    gf.write('#0.083 lon lat grid file \n')
    gf.write('gridtype  =  lonlat\n')
    gf.write('xsize =' + str(len(x_1d)) + '\n')
    gf.write('ysize =' + str(len(y_1d)) + '\n')
    gf.write('xfirst =' + str(xmin + dx / 2.) + '\n')
    gf.write('xinc =' + str(dx) + '\n')
    gf.write('yfirst =' + str(ymin + dy / 2.) + '\n')
    gf.write('yinc =' + str(dy) + '\n')
    gf.close()

    resfact = dy_regrid_coarse/dy_regrid_fine
    gf = open('gridfile-Reg0.0083', 'w')
    gf.write('#0.0083 res lon lat grid file \n')
    gf.write('gridtype = lonlat\n')
    gf.write('xsize =' + str(int(resfact * len(x_1d))) + '\n')
    gf.write('ysize =' + str(int(resfact * len(y_1d))) + '\n')
    gf.write('xfirst =' + str(xmin + dx / (2. * resfact)) + '\n')
    gf.write('xinc =' + str(dx / resfact) + '\n')
    gf.write('yfirst =' + str(ymin + dy / (2. * resfact)) + '\n')
    gf.write('yinc =' + str(dy / resfact) + '\n')
    gf.close()

    # Get cdo to refine to the desired resolution
    # Newer cdo's use *much* less memory using no_remap_weights option
    cmd = "cdo --no_remap_weights remapbil,gridfile-Reg0.0083 " + \
          "-setgrid,gridfile-Reg0.083 cf_0.083deglatlonfield.nc " + \
          "cdo_0.0083deglatlonfield.nc"
    t0 = time.time()
    status = os.system(cmd)
    print("cdo refine: %.1fs" % (time.time() - t0))
    if status != 0:
        print("ERROR: failed to run cdo command")
        sys.exit(1)

    # Read the cdo results back in so we can return the array from this function
    t0 = time.time()
    field_REGRID_array = cf.read("cdo_0.0083deglatlonfield.nc")[0].array
    print("read back cdo output: %.1fs" % (time.time() - t0))
    cmd = "rm cdo_0.0083deglatlonfield.nc gridfile-Reg0.0083 " + \
          "gridfile-Reg0.083 cf_0.083deglatlonfield.nc"
    status = os.system(cmd)
    if status != 0:
        print("ERROR: failed to remove intermediate cdo files")
        sys.exit(1)
    # End two_stage extension
    #######################################

    # centre of the first (south-west) cell of the fine regional grid, needed
    # later to work out where this regional data sits within a global ancil
    lon0 = xmin + dx_regrid_fine / 2.
    lat0 = ymin + dy_regrid_fine / 2.

    return field_REGRID_array, lon0, lat0


def feather_weights(topog_data, feather_deg):
    """Weight 0 at the edge of valid BISICLES data rising linearly to 1
    feather_deg degrees inside it (array edges count as edges)."""
    from scipy.ndimage import distance_transform_edt

    valid = ~np.ma.getmaskarray(topog_data)
    if feather_deg <= 0:
        return valid.astype(np.float32)
    padded = np.pad(valid, 1, constant_values=False)
    dist = distance_transform_edt(padded)[1:-1, 1:-1]
    nfeather = feather_deg / dy_regrid_fine
    return np.clip(dist / nfeather, 0., 1.).astype(np.float32)


def splice_into_ancil_template(topog_data, lon0, lat0, ancil_template,
                               ancil_output, feather_deg=4):
    import mule
    import common_mule_rss as mule_rss

    """
    Update the topography file CAP will use as input only where we have
    data from BISICLES.

    ancil_template may either be cropped to exactly the regional domain
    (in which case topog_data lines up 1:1 with it), or it may be a full
    global field, in which case we work out where the regional data
    belongs using the grid definitions in the ancil's own header and the
    known origin (lon0, lat0) of the regional grid, rather than assuming
    the two arrays share an origin.
    """

    t0 = time.time()
    fieldg = mule.load_umfile(ancil_template)
    field0 = fieldg.fields[0]
    globe30 = field0.get_data()
    print("load ancil template + get_data: %.1fs" % (time.time() - t0))

    # Latitude axes in the usual CAP input field ordered in the opposite 
    # direction, so flip our data array N-S
    topog_data = topog_data[::-1, :]
    ny, nx = topog_data.shape

    # lat/lon of the centre of each row/column of the regional data, in the
    # same (north to south) order as topog_data after the flip above
    lat_reg = lat0 + dy_regrid_fine * np.arange(ny)
    lat_reg = lat_reg[::-1]
    lon_reg = lon0 + dx_regrid_fine * np.arange(nx)

    nrow, ncol = globe30.shape

    # blend weights so the BISICLES edge fades into the existing GLOBE30 data
    weights = feather_weights(topog_data, feather_deg)

    t0 = time.time()
    if (nrow, ncol) == (ny, nx):
        # template is already cropped to the regional domain, splice directly
        update = np.where(topog_data.mask == False)
        w = weights[update]
        globe30[update] = w * topog_data[update] + (1. - w) * globe30[update]
    else:
        # template is a full global field - use the file-level real_constants
        # for the grid definition. Some ancils (e.g. qrparm.orog) carry a
        # bogus, rescaled column spacing in the per-field lookup headers
        # (bzx/bdx), so those can't be trusted for this calculation.
        # Latitude never wraps, so the regional rows map onto one contiguous
        # block of global rows. Longitude can wrap around the date line (and,
        # since the regional grid deliberately overlaps a bit for safety, may
        # even wrap more than once), so that's handled a contiguous segment
        # at a time. Plain slicing is used throughout (rather than fancy
        # indexing) as it is vastly faster for arrays this size.
        rc = fieldg.real_constants

        row0 = int(round((rc.start_lat - lat_reg[0]) / rc.row_spacing))
        row1 = row0 + ny
        r_lo, r_hi = max(row0, 0), min(row1, nrow)
        t_row_lo, t_row_hi = r_lo - row0, ny - (row1 - r_hi)

        if r_lo < r_hi:
            col_start = int(round(np.mod(lon_reg[0] - rc.start_lon, 360.)
                                  / rc.col_spacing)) % ncol
            src_col = 0
            remaining = nx
            while remaining > 0:
                take = min(remaining, ncol - col_start)

                sub_globe = globe30[r_lo:r_hi, col_start:col_start + take]
                sub_topog = topog_data[t_row_lo:t_row_hi,
                                       src_col:src_col + take]
                update = np.where(sub_topog.mask == False)
                w = weights[t_row_lo:t_row_hi,
                            src_col:src_col + take][update]
                sub_globe[update] = w * sub_topog[update] \
                    + (1. - w) * sub_globe[update]
                globe30[r_lo:r_hi, col_start:col_start + take] = sub_globe

                src_col += take
                remaining -= take
                col_start = 0

    print("splice loop: %.1fs" % (time.time() - t0))

    # poke the updated array back into the CAP ancil
    t0 = time.time()
    fieldg.fields[0] = mule_rss.overwrite_data(fieldg.fields[0], globe30)
    fieldg = disable_mule_validators(fieldg)
    fieldg.to_file(ancil_output)
    print("write output ancil: %.1fs" % (time.time() - t0))

    return


if __name__ == "__main__":

    # Read the command line arguments
    ancilin_file, hdf5in_file, ancilout_file, region = parse_commandline()

    # Get what we need from the BISICLES plot file
    print("extract fields from BISICLES hdf5")
    surface_height, lon, lat, x_dist, y_dist \
        = extract_bisicles_topography(hdf5in_file, region=region)

    # Make topography ancil for CAP
    print("regrid, splice CAP fields")
    surface_height_CF = construct_CF(surface_height, lon, lat, x_dist, y_dist)
    surface_height_REGRID, lon0, lat0 = regrid_to_GLOBE30(surface_height_CF,
                                                          region=region)
    splice_into_ancil_template(surface_height_REGRID, lon0, lat0,
                               ancilin_file, ancilout_file)
