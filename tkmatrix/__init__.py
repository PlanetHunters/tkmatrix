"""Root module for the MATRIX ToolKit classes"""
import os
import shutil
import subprocess
import sys

# Determine the path to libellc.so
ellc_path = os.path.join(os.path.dirname(__file__), 'ellc')
lib_path = os.path.join(ellc_path, 'ellc','libellc.so')
# Check if it exists
if not os.path.exists(lib_path):
    print("[ellc] libellc.so not found, running make...")
    try:
        subprocess.check_call(['make', '-B'], cwd=ellc_path)
        shutil.copy(ellc_path + '/libellc.so', os.path.join(ellc_path, 'ellc') + '/libellc.so')
        print("[ellc] libellc.so built successfully.")
    except Exception as e:
        print(f"Could not build libellc.so. Please ensure make and dependencies are available: {e}")

#Patching ellc with submodule
import tkmatrix.ellc.ellc as _mypackage_ellc

# Override the 'ellc' name in sys.modules to point to your internal package module
sys.modules['ellc'] = _mypackage_ellc