import os

import yaml
from argparse import ArgumentParser
from os import path
from tkmatrix.tkmatrix_class import MATRIX


if __name__ == '__main__':
    current_path = os.path.dirname(os.path.realpath(__file__))
    ap = ArgumentParser(description='Sherlock Inject&Recovery tool')
    ap.add_argument('--dir', default="./", help="Working directory (if empty your current dir will be assumed)",
                    required=False)
    ap.add_argument('--properties', help="Configuration file", required=True)
    args = ap.parse_args()
    resources_dir = path.join(path.dirname(__file__))
    file_dir = resources_dir + "/" + 'properties.yaml' if resources_dir != "" and resources_dir is not None \
        else 'properties.yaml'
    print(resources_dir)
    tirma_user_properties = yaml.load(open(file_dir), yaml.SafeLoader)
    user_properties = yaml.load(open(args.properties), yaml.SafeLoader)
    tirma_user_properties.update(user_properties)
    ir = MATRIX(tirma_user_properties["TARGET"], tirma_user_properties["SECTORS"], args.dir) \
        .inject(tirma_user_properties["PHASES"], tirma_user_properties["MIN_PERIOD"],
                tirma_user_properties["MAX_PERIOD"], tirma_user_properties["STEP_PERIOD"],
                tirma_user_properties["MIN_RADIUS"], tirma_user_properties["MAX_RADIUS"],
                tirma_user_properties["STEP_RADIUS"]) \
        .recovery(tirma_user_properties["CPUS"], tirma_user_properties["SHERLOCK_DEEPNESS"],
                  tirma_user_properties["KNOWN_TRANSITS"], tirma_user_properties["DETREND_WS"])
