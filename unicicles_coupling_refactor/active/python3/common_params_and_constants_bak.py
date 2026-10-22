"""
common_params_and_constants

central place to store commonly reused parameters and magic numbers
"""

import numpy as np

# This matches the value in Glimmer
rho_ice_glimmer = 910.

# These match the value in the UM/JULES run
nsmax = 10
number_elevs = 10
dzsnow = np.asarray([0.04, 0.12, 0.34, 0.5, 0.5, 0.5, 1., 1., 2., 2.])
tile_elevations = np.asarray([31.96, 301.57, 550.11, 850.85, 1156.89,
                              1453.16, 1806.74, 2268.25, 2753.70, 3341.4])
# Standard constants
planet_radius = 6371229.0
pi = 3.14159265358979323846

degC_to_Kelvin = -273.15
deg_to_rad = pi / 180.
# UKESM1 still uses the 360 day calendar. UKESM2 will move to
# Gregorian, with leap years, and we'll have an accuracy issue
# however we set this. iirc BISICLES actually has its own internal
# hardcoded factor for seconds->years conversions for some
# variables, so we whatever we do it'll be inconsistent somewhere
days_in_year = 360.
secs_in_day = 24. * 60. * 60.


# These match jules/src/control/shared/max_dimensions.F90
max_number_elevs = 25

# These match jules/src/params/um/land_tile_ids.F90
start_tile_id_elev_ice = 901
start_tile_id_elev_rock = 926
start_tile_id_snowlayers = 1 * 1000
start_tile_id_elev_ice_snowlayers = start_tile_id_elev_ice * 1000

# These match um/src/control/stash/um_stashcode_mod.F90
stashcode_orog_x_grad = 5
stashcode_orog_y_grad = 6
stashcode_sil_orog_rough = 17
stashcode_hlf_pk_to_trf_ht = 18
stashcode_snow_amount = 23
stashcode_lsm = 30
stashcode_orog = 33
stashcode_orog_var = 34
stashcode_orog_gdxx = 35
stashcode_orog_gdxy = 36
stashcode_orog_gdyy = 37
stashcode_vol_smc_sat = 43
stashcode_frac_surf_type = 216
stashcode_rgrain = 231
stashcode_snow_tile = 240
stashcode_snowdep_grd_tile = 376
stashcode_snowpack_bk_dens = 377
stashcode_nsnow_layrs_tiles = 380
stashcode_snow_laythk_tiles = 381
stashcode_snow_ice_tile = 382
stashcode_snow_liq_tile = 383
stashcode_snow_t_tile = 384
stashcode_snow_laydns_tiles = 385
stashcode_snow_grnsiz_tiles = 386
stashcode_land_frac = 505
stashcode_tsurf_elev_surft = 576
stashcode_tile_fractions = 3317
stashcode_land_frac_diag = 3395
stashcode_tile_snow_mass = 8236
stashcode_tile_snow_depth = 8376
stashcode_tile_snow_temp = 8576
stashcode_tile_snow_htflux = 8577
stashcode_tile_smb = 8578

# these are standard coordinate systems for different projections
# used by pyproj
epsg_global = 'epsg:4326'

# Ask CPOM BISICLES people for offsets for different domains
# when they create them
epsg_GrIS = 'epsg:3413'
x0_GrIS = -654650.0
y0_GrIS = -3385950.0

epsg_AIS = 'epsg:3031'
x0_AIS = -3072000
y0_AIS = -3072000

# earlier (BedMap2?) versions of the BISICLES Amundsen Sea Embayment
# used something a little different, check what you're coupling to
epsg_ASE = 'epsg:3031'
x0_ASE = -1831000
y0_ASE = -904000

# Offsets of regional domains within the UM N96 global grid
xmin_GrIS = 250.
xmax_GrIS = 360.
ymin_GrIS = 45.
ymax_GrIS = 90.
xcyclic_GrIS = False

xmin_AIS = 0.
xmax_AIS = 360.
ymin_AIS = -90.
ymax_AIS = -55.
xcyclic_AIS = True

# resolutions of finer grids for feeding orography to CAP
dx_regrid_coarse = 0.083333
dy_regrid_coarse = 0.083333

dx_regrid_fine = 0.0083333
dy_regrid_fine = 0.0083333

# Minimum depth tests for bisicles_global_to_nemo
bgtn_thresh = 1e-3
bgtn_mindep_isf = 1e-3
bgtn_mindep_bath = 1e-3
bgtn_mincav = 1e-3
bgtn_topopen = 12
bgtn_max_jindex = 300
bgtn_pflag = 99999.
bgtn_nflag = -99999.
bgtn_bathy_threshold = -5000.

# CAP regional, resolution params
ncols_GrIS_n216 = 122
ncols_AIS_n216 = 432
assig_n216 = 1e-4
h_n216 = 0.3
xoff_GrIS_n216 = 309
yoff_GrIS_n216 = 266
xoff_AIS_n216 = 0
yoff_AIS_n216 = 1

ncols_GrIS_n96 = 54
nrows_GrIS_n96 = 25
ncols_AIS_n96 = 192
nrows_AIS_n96 = 25
assig_n96 = 5e-5
h_n96 = 0.15
xoff_GrIS_n96 = 137
yoff_GrIS_n96 = 118
xoff_AIS_n96 = 0
yoff_AIS_n96 = 1
las_k = 0.4
las_cd = 0.3
las_a = 5e-4
las_bigsig = 100.
