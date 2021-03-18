import os

import yaml
from argparse import ArgumentParser
from os import path
from tkmatrix.tkmatrix_class import MATRIX
import datetime


# Check the variable CPUS from user-properties.
# If the value is greater than the real cpus number, we replace it.
def get_cpus():
    user_properties_cpus = tirma_user_properties["CPUS"]
    cpu_count = os.cpu_count()
    if user_properties_cpus > cpu_count:
        final_cpus = cpu_count
        print("Property change [CPUS]: You only have " + str(cpu_count) + " CPUs.")
    else:
        final_cpus = tirma_user_properties["CPUS"]

    return final_cpus


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
    tirma_user_properties = yaml.load(open(file_dir), yaml.SafeLoader)
    user_properties = yaml.load(open(args.properties), yaml.SafeLoader)
    preserve = args.preserve
    tirma_user_properties.update(user_properties)

    tirma_user_properties["CPUS"] = get_cpus()

    ir = MATRIX(tirma_user_properties["TARGET"], tirma_user_properties["SECTORS"], args.dir, args.preserve) \
        .inject(tirma_user_properties["PHASES"], tirma_user_properties["MIN_PERIOD"],
                tirma_user_properties["MAX_PERIOD"], tirma_user_properties["STEP_PERIOD"],
                tirma_user_properties["MIN_RADIUS"], tirma_user_properties["MAX_RADIUS"],
                tirma_user_properties["STEP_RADIUS"], tirma_user_properties["EXPOSURE_TIME"]) \
        .recovery(tirma_user_properties["CPUS"], tirma_user_properties["SHERLOCK_DEEPNESS"],
                  tirma_user_properties["KNOWN_TRANSITS"], tirma_user_properties["DETREND_WS"],
                  tirma_user_properties["FIT_METHOD"])

    # print the execution time:
    print("Execution time: " + str(datetime.datetime.now() - start_time))
