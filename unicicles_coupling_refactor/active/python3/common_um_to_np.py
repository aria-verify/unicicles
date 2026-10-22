"""
common_um_to_np

utilities to read/write numpy arrays into / out of UM dumps and coupling
fields.
"""

import sys
import numpy as np
import common_mule_rss as mule_rss
from common_params_and_constants import start_tile_id_elev_ice, \
    start_tile_id_elev_ice_snowlayers, start_tile_id_snowlayers

SMALL_THRESHOLD = 1e-4


def um_to_np_2d(um_dump, stash):
    """
    extract a single 2d field from a UM dump by stash code and return as a numpy array.
    """
    try:
        #find the location of the field by searching for the stash code
        #assume user has checked there will only be one occurrence of this
        #code (so only one level), since this will only return one 2d field
        #from the first location 
        index = mule_rss.findindex(um_dump, stash)[0]
    except IndexError:
        print("")
        print("ERROR in um_to_np_2d. Field not found for STASH: ", stash)
        print("")
        sys.exit(2)

    return um_dump.fields[index].get_data()


def find_pslevs(um_dump, indices):
    """
    find the pseudolevel ID codes associated with the fields at a set of indices
    return them as a numpy array.
    """
    pslev = []
    for i in indices:
        pslev.append(um_dump.fields[i].lbuser5)

    pslev = np.asarray(pslev)
    return pslev


def select_elev_indices(um_dump, indices):
    """
    return the sub-set of locations (and their pseudolevel ID codes) of just the 
    elevated tiles levels in a multi-level set of UM fields. Just used internally
    to this library by the 3d and 4d um_to_np and np_to_um routines.
    
    Needs to cope with ID codes for 3d fields [lon, lat, tile] and 4d ones [lon, lat, tile, snow_level]
    """
    pslev = find_pslevs(um_dump, indices)

    indices = indices[np.where( ((pslev >= start_tile_id_elev_ice) &
                                 (pslev < start_tile_id_snowlayers)) |
                                (pslev >= start_tile_id_elev_ice_snowlayers) )]

    pslev = pslev[np.where( ((pslev >= start_tile_id_elev_ice) &
                             (pslev < start_tile_id_snowlayers)) |
                            (pslev >= start_tile_id_elev_ice_snowlayers) )]

    return indices, pslev


def um_to_np_3d(um_dump, stash, elev=False):
    """
    extract a 3d (ie stack of 2d levels) field from a UM dump by stash code and 
    return as a numpy array. Optionally, return just the elevated tile levels, which
    is often all the ice sheet coupling cares about.
    """
    try:
        indices = mule_rss.findindex(um_dump, stash)
    except IndexError:
        print("")
        print("ERROR in um_to_np_3d. Field not found for STASH: ", stash)
        print("")
        sys.exit(2)

    #having got all the levels, get the ID codes that identify which pseudolevels, if
    #that's what they are
    pslev = find_pslevs(um_dump, indices)
    if elev:
        #sub-select the locations of just the elevated ice tiles, if desired
        indices, pslev = select_elev_indices(um_dump, indices)

    return mule_rss.get_data3d(um_dump, indices), pslev


def um_to_np_4d(um_dump, stash, elev=False):
    """
    extract a 4d field from a UM dump by stash code and return as a numpy array. 
    Usually these are internal snowpack fields with possible dimensions of 
    [lon, lat, tile_id, snow_level] which the UM dump stores as a big stack of 
    2d [lon, lat] fields in which format it's difficult to find the field for the individual 
    layer on the individual tile you want. Optionally, return just the elevated tile levels, which
    is often all the ice sheet coupling cares about.
    """
    try:
        indices = mule_rss.findindex(um_dump, stash)
    except IndexError:
        print("")
        print("ERROR in um_to_np_4d. Field not found for STASH: ", stash)
        print("")
        sys.exit(2)

    #having got all the levels and tiles, get the ID codes that identify which pseudolevels, if
    #that's what they are
    pslev = find_pslevs(um_dump, indices)
    if elev:
        #sub-select the locations of snow levels on just the elevated ice tiles, if desired
        indices, pslev = select_elev_indices(um_dump, indices)

    #pull the field out of the dump as 3d stack, before discrimintaing along the 3rd dimension 
    #into snow levels and tiles
    field3d = mule_rss.get_data3d(um_dump, indices)

    #work out the size of each of the 4 dimensions we'll need. JULES snow level id codes are
    #constructed as [(1000*tile id code) + snow level index]
    nsl = int(max(pslev - (pslev // 1000) * 1000))
    nt = int(len(pslev) / nsl)
    ny = int(np.shape(field3d)[0])
    nx = int(np.shape(field3d)[1])

    field4d = np.zeros([ny, nx, nt, nsl])

    #populate the 4d template by splitting the 3rd dimension of field3d
    #by snow_level for the different tiles
    for i in range(len(pslev)):
        tile_id = pslev[i] // 1000
        snowlayer_index = pslev[i] - (tile_id * 1000) - 1
        tile_index = np.mod(i, nt)
        field4d[:, :, tile_index, snowlayer_index] = field3d[:, :, i]

    return field4d, pslev


def np_to_um_2d(um_dump, stash, array):
    """
    modify a single 2d field in a UM dump by stash code from a numpy array.
    """
    try:
        #find the existing location of the field in the dump by searching for the stash code
        #assume user has checked there will only be one occurrence of this
        #code (so only one level), since this will only modify one 2d field
        #in the first location
        index = mule_rss.findindex(um_dump, stash)[0]
    except IndexError:
        print("")
        print("ERROR in np_to_um_2d. Field not found for STASH: ", stash)
        print("")
        sys.exit(2)

    um_dump.fields[index] = mule_rss.overwrite_data(um_dump.fields[index],
                                                    array)

    return um_dump


def np_to_um_3d(um_dump, stash, array, elev=False):
    """
    modify a 3d (ie stack of 2d levels) field in a UM dump by stash code 
    by inserting a numpy array. Optionally, modify just the elevated tile levels, which
    is often all the ice sheet coupling cares about.
    """
    try:
        indices = mule_rss.findindex(um_dump, stash)
    except IndexError:
        print("")
        print("ERROR in np_to_um_3d. Field not found for STASH: ", stash)
        print("")
        sys.exit(2)

    if elev:
        #sub-select the locations of just the elevated ice tiles, if desired
        indices, pslev = select_elev_indices(um_dump, indices)

    if len(indices) != np.shape(array)[2]:
        print("ERROR in np_to_um_3d. Array doesn't match 3d axis length")
        print("")
        sys.exit(2)
    else:
        um_dump = mule_rss.put_data3d(um_dump, indices, array)

    return um_dump


def np_to_um_4d(um_dump, stash, field4d, elev=False):
    """
    modify a 4d (ie stack of 2d levels by tile id and snow level) field in a UM dump by stash code 
    by inserting a numpy array. Optionally, modify just the elevated tile levels, which
    is often all the ice sheet coupling cares about.
    """
    try:
        indices = mule_rss.findindex(um_dump, stash)
    except IndexError:
        print("")
        print("ERROR in np_to_um_4d. Field not found for STASH: ", stash)
        print("")
        sys.exit(2)

    #having got all the levels, get the ID codes that identify which pseudolevels
    pslev = find_pslevs(um_dump, indices)
    if elev:
        #sub-select the locations of just the elevated ice tiles, if desired
        indices, pslev = select_elev_indices(um_dump, indices)

    ny = np.shape(field4d)[0]
    nx = np.shape(field4d)[1]
    nt = np.shape(field4d)[2]
    nsl = np.shape(field4d)[3]

    #UM stack of 2d fields needs to combine the 3rd and 4th dimension (tile id and snowpack level)
    #of the numpy array by interleaving them into the stack in the right order. UM stack order
    #varies snow level faster than tile id. We use the pseudolevel id to identify the layers we're
    #looking at, JULES snow level id codes are constructed as [(1000*tile id code) + snow level index]
    array = np.zeros([ny, nx, nt * nsl])

    for i in range(len(pslev)):
        tile_id = pslev[i] // 1000
        snowlayer_index = pslev[i] - (tile_id * 1000) - 1
        tile_index = np.mod(i, nt)
        array[:, :, i] = field4d[:, :, tile_index, snowlayer_index]

    if len(indices) != np.shape(array)[2]:
        print("ERROR in np_to_um_4d. Array doesn't match 3d axis length")
        print("")
        sys.exit(2)
    else:
        #put the final 3d stack of fields into the dump
        um_dump = mule_rss.put_data3d(um_dump, indices, array)

    return um_dump


def ice_to_np(ice_file, ncvar_name, ns_flip=True):
    """
    extract a field from a Glint format netCDF file by name, and reorder it
    to match the axes order etc we have in the fields extracted from UM dumps
    """
    try:
        field = ice_file.variables[ncvar_name][:].squeeze()
    except KeyError:
        print("")
        print("ERROR in ice_to_np: field ", ncvar_name, " not found in file")
        print("")
        sys.exit(2)

    if field.ndim > 2:
        # Use UM convention: ice_file arrays are [(nsnow),ntile,ny,nx], um
        # dumps are [ny,nx,ntile]
        field = shift_tile_dim(field)
    if field.ndim > 3:
        # The above has done the snow. Now move the tile fraction back past the
        # 2D
        field = shift_tile_dim(field, end=0)

    # The glint-produced ice_file arrays are N-S flipped wrt UM
    # The atmos.iceceouple file fed *to* glint is the right way up though
    if ns_flip:
        field = field[::-1, ...]

    return field


def shift_tile_dim(field, end=1):
    """
    Simple interface to allow some common array rearranging between the two conventions being used in 
    the models in the ice coupling. Fields extracted from UM dumps are ordered 
    [ny,nx,ntile,(nsnow)], ice_file arrays are [(nsnow),ntile,ny,nx].
    """

    # Shuffle a leading-dimension (ie tiles) to the back, eg
    # [ntile,ny,nx] -> [ny,nx,ntile] (this is my mule convention)
    # this swaps position 1 -> 3
    field = np.swapaxes(field, 0, 1)
    field = np.swapaxes(field, 1, 2)

    # A 4D array needs one more swap to go right to the back
    if (field.ndim == 4) & (end == 1):
        field = np.swapaxes(field, 2, 3)

    return field


def get_ice_um_grid_overlap(ice_file):
    """
    Pull 3d tile land fraction from ice model and collapse it down to say what fraction
    of each UM [lon,lat] gridcell exists in the ice model regional domain.
    We need this for the possibly fractional overlap of grid coverage at the edges of the
    ice domain.
    """

    land_frac_ice = ice_to_np(ice_file, 'tile_land_fraction')
    fmask_ice = np.sum(land_frac_ice, axis=2)

    return fmask_ice


def merge_ice_um_arrays(um, ice, overlap, points=None, noFrac=False):
    """
    merge an ice-model regional field into an existing global UM one
    weighted by the partial grid coverage of the ice grid at the edges

    Optionally only do this on a specified subset of (ie ice cell) points
    in the global field. Optionally, for cells with only partial coverage of the 
    global domain cell at the edges of the regional ice model domain, weight 
    the information being merged in by the coverage fraction, or ignore the
    weighting.
    """

    if noFrac:
        overlap[np.where(np.abs(overlap - 1.0) < SMALL_THRESHOLD)] = 1.
        overlap[np.where(np.abs(overlap - 1.0) >= SMALL_THRESHOLD)] = 0.

    # Splice in the new field from the ice, weighted by the UM/ice grid overlap
    if points is None:
        um = ice * overlap + um * (1. - overlap)
    else:
        um[points] = ice[points] * overlap[points] + \
            um[points] * (1. - overlap[points])

    return um
