import os

from lcbuilder.star.starinfo import StarInfo

from tkmatrix.tkmatrix_class import MATRIX

matrix_user_properties = {}
matrix_user_properties['TARGET'] = "KIC 11615498"
star_properties = {}
star_properties['MASS'] = 1.977
star_properties['MASS_LOWER_ERROR'] = 0.062
star_properties['MASS_UPPER_ERROR'] = 0.062
star_properties['RADIUS'] = 9.863
star_properties['RADIUS_LOWER_ERROR'] = 0.015
star_properties['RADIUS_UPPER_ERROR'] = 0.015
star_properties['TEFF'] = 4978.6
star_properties['LD_COEFFICIENTS'] = [0.4356, 0.2367]
star_info = StarInfo(matrix_user_properties['TARGET'],
         tuple(star_properties["LD_COEFFICIENTS"]) if "LD_COEFFICIENTS" in star_properties else None,
         star_properties["TEFF"] if "TEFF" in star_properties else None,
         star_properties["LUM"] if "LUM" in star_properties else None,
         star_properties["LOGG"] if "LOGG" in star_properties else None,
         star_properties["LOGG_ERR"] if "LOGG_ERR" in star_properties else None,
         star_properties["RADIUS"] if "RADIUS" in star_properties else None,
         star_properties["RADIUS_LOWER_ERROR"] if "RADIUS_LOWER_ERROR" in star_properties else None,
         star_properties["RADIUS_UPPER_ERROR"] if "RADIUS_UPPER_ERROR" in star_properties else None,
         star_properties["MASS"] if "MASS" in star_properties else None,
         star_properties["MASS_LOWER_ERROR"] if "MASS_LOWER_ERROR" in star_properties else None,
         star_properties["MASS_UPPER_ERROR"] if "MASS_UPPER_ERROR" in star_properties else None,
         star_properties["RA"] if "RA" in star_properties else None,
         star_properties["DEC"] if "DEC" in star_properties else None)
matrix_user_properties['SECTORS'] = 'all'
matrix_user_properties['PERIOD_GRID'] = None
matrix_user_properties['RADIUS_GRID'] = None
matrix_user_properties['PERIOD_GRID_GEOM'] = 'lin'
matrix_user_properties['RADIUS_GRID_GEOM'] = 'lin'
matrix_user_properties['PHASES'] = 2
matrix_user_properties["CPUS"] = os.cpu_count() - 1
matrix_user_properties['MIN_PERIOD'] = 4
matrix_user_properties['MAX_PERIOD'] = 20
matrix_user_properties['MIN_PERIOD_SEARCH'] = 0.5
matrix_user_properties['MAX_PERIOD_SEARCH'] = 40
matrix_user_properties['STEPS_PERIOD'] = 2
matrix_user_properties['MIN_RADIUS'] = 4
matrix_user_properties['MAX_RADIUS'] = 22
matrix_user_properties['STEPS_RADIUS'] = 2
matrix_user_properties['PHASES'] = 2
matrix_user_properties['DETREND_WS'] = 3.5
matrix_user_properties['OVERSAMPLING'] = 1
matrix_user_properties['SIGNAL_SELECTION_MODE'] = 'period'
matrix_user_properties['FIT_METHOD'] = 'bls-periodogram'
matrix_user_properties['SNR_THRESHOLD'] = 5
matrix_user_properties['USE_SEARCH_CACHE'] = True
matrix_user_properties['DETREND_METHOD'] = 'biweight'
matrix_user_properties['RUN_LIMIT'] = 2
target = matrix_user_properties["TARGET"]
file = None
author = None
custom_search = None
prepare_algorithm = None
initial_mask = None
initial_smooth_enabled = False
initial_transit_mask = None
auto_detrend_period = None
auto_detrend_ratio = 1
auto_detrend_method = 'biweight'
auto_detrend_enabled = False
oscillation_reduction = False
oscillation_min_snr = 3
oscillation_amplitude_threshold = 0.1
oscillation_ws_percent = 60
oscillation_min_period = 0.0006875
oscillation_max_period = 0.25
high_rms_bin_hours = 4
high_rms_threshold = 3
high_rms_enabled = False
outliers_sigma = 3
exptime = 1800
eleanor_corr_flux = "pca_flux"
cache_dir = None
if cache_dir is None:
    cache_dir = os.path.expanduser('~') + "/"
search_engine = 'cpu'
run_dir = './'
ir = MATRIX(target, matrix_user_properties["SECTORS"], author, run_dir, False, star_info, file, exptime,
                initial_mask, initial_transit_mask, eleanor_corr_flux, outliers_sigma, high_rms_enabled,
                high_rms_threshold, high_rms_bin_hours, initial_smooth_enabled, auto_detrend_enabled,
                auto_detrend_method, auto_detrend_ratio, auto_detrend_period, prepare_algorithm, cache_dir,
                oscillation_reduction, oscillation_min_snr, oscillation_amplitude_threshold, oscillation_ws_percent,
                oscillation_min_period, oscillation_max_period, matrix_user_properties["CPUS"],
                search_engine=search_engine)
inject_dir = None
inject_dir, period_grid, radius_grid = ir.inject(matrix_user_properties["PHASES"],
                           matrix_user_properties["MIN_PERIOD"], matrix_user_properties["MAX_PERIOD"],
                           matrix_user_properties["STEPS_PERIOD"],
                           matrix_user_properties["MIN_RADIUS"], matrix_user_properties["MAX_RADIUS"],
                           matrix_user_properties["STEPS_RADIUS"],
                           period_grid=matrix_user_properties['PERIOD_GRID'],
                           radius_grid=matrix_user_properties['RADIUS_GRID'],
                           period_grid_geom=matrix_user_properties["PERIOD_GRID_GEOM"],
                           radius_grid_geom=matrix_user_properties["RADIUS_GRID_GEOM"],
                           inject_dir=inject_dir)
ir.recovery(inject_dir, matrix_user_properties["SNR_THRESHOLD"],
            matrix_user_properties["DETREND_METHOD"],
            matrix_user_properties["DETREND_WS"], matrix_user_properties["FIT_METHOD"],
            matrix_user_properties["RUN_LIMIT"],
            None, matrix_user_properties["MAX_PERIOD_SEARCH"], matrix_user_properties["OVERSAMPLING"],
            matrix_user_properties["SIGNAL_SELECTION_MODE"],
            use_search_cache=matrix_user_properties["USE_SEARCH_CACHE"])