import os
import sys
import pickle
import yaml
from argparse import ArgumentParser
from os import path
from lcbuilder.star.starinfo import StarInfo
import importlib.util
from tkmatrix.tkmatrix_class import MATRIX
import datetime
from pathlib import Path


def load_module(module_path):
    spec = importlib.util.spec_from_file_location("customs", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module

def extract_custom_class(module_path):
    class_module = None
    if module_path is not None:
        class_module = load_module(module_path)
        class_name = Path(module_path.replace(".py", "")).name
        class_module = getattr(class_module, class_name)
        globals()[class_name] = class_module
        pickle.dumps(class_module)
        class_module = class_module()
    return class_module

# Check the variable CPUS from user-properties.
# If the value is greater than the real cpus number, we replace it.
def get_cpus():
    user_properties_cpus = matrix_user_properties["CPUS"]
    cpu_count = os.cpu_count()
    if user_properties_cpus > cpu_count:
        final_cpus = cpu_count
        print("Property change [CPUS]: You only have " + str(cpu_count) + " CPUs.")
    else:
        final_cpus = matrix_user_properties["CPUS"]

    return final_cpus

def get_star_info(properties, id):
    input_star_info = None
    if properties["STAR"] is not None and properties["STAR"][id] is not None:
        star_properties = properties["STAR"][id]
        input_star_info = StarInfo(id, tuple(star_properties["LD_COEFFICIENTS"]) if "LD_COEFFICIENTS" in star_properties else None,
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
    return input_star_info


if __name__ == '__main__':
    # We save the start time:
    start_time = datetime.datetime.now()

    current_path = os.path.dirname(os.path.realpath(__file__))
    ap = ArgumentParser(description='Sherlock Inject&Recovery tool')
    ap.add_argument('--dir', default="./", help="Working directory (if empty your current dir will be assumed)",
                    required=False)
    ap.add_argument('--properties', help="Configuration file", required=True)
    ap.add_argument('--preserve', help="Preserve the inject file. By default they should be removed and only kept if the flag is enabled.",
                    action="store_true", required=False)
    args = ap.parse_args()
    resources_dir = path.join(path.dirname(__file__))
    file_dir = resources_dir + "/" + 'properties.yaml' if resources_dir != "" and resources_dir is not None \
        else 'properties.yaml'
    print("The resource dir is: " + str(resources_dir))
    matrix_user_properties = yaml.load(open(file_dir), yaml.SafeLoader)
    user_properties = yaml.load(open(args.properties), yaml.SafeLoader)
    preserve = args.preserve
    matrix_user_properties.update(user_properties)
    matrix_user_properties["CPUS"] = get_cpus()
    target = matrix_user_properties["TARGET"]
    file = matrix_user_properties["FILE"]
    star_info = get_star_info(matrix_user_properties, target)
    custom_search = extract_custom_class(matrix_user_properties["CUSTOM_SEARCH_ALGORITM"])
    custom_clean = extract_custom_class(matrix_user_properties["CUSTOM_CLEAN_ALGORITM"])

    ir = MATRIX(target, matrix_user_properties["SECTORS"], args.dir, args.preserve, star_info, file,
                matrix_user_properties["EXPOSURE_TIME"])
    inject_dir = ir.inject(matrix_user_properties["PHASES"], matrix_user_properties["MIN_PERIOD"],
                           matrix_user_properties["MAX_PERIOD"], matrix_user_properties["STEPS_PERIOD"],
                           matrix_user_properties["MIN_RADIUS"], matrix_user_properties["MAX_RADIUS"],
                           matrix_user_properties["STEPS_RADIUS"], matrix_user_properties["PERIOD_GRID_GEOM"],
                           matrix_user_properties["RADIUS_GRID_GEOM"])
    ir.recovery(matrix_user_properties["CPUS"], inject_dir, matrix_user_properties["SNR_THRESHOLD"],
                matrix_user_properties["SHERLOCK_DEEPNESS"], matrix_user_properties["KNOWN_TRANSITS"],
                matrix_user_properties["DETREND_WS"], matrix_user_properties["FIT_METHOD"],
                matrix_user_properties["RUN_LIMIT"],
                matrix_user_properties["DETREND_PERIOD"], matrix_user_properties["DETREND_PERIOD_METHOD"],
                custom_clean, custom_search)
    ir.plot_results(target, inject_dir, period_grid_geom=matrix_user_properties["PERIOD_GRID_GEOM"],
                    radius_grid_geom=matrix_user_properties["RADIUS_GRID_GEOM"])
    # print the execution time:
    print("Execution time: " + str(datetime.datetime.now() - start_time))
