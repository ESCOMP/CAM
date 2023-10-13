#! /usr/bin/env python3

"""
Script run from CIME when calling case.setup
Expects 3 arguments:
   (1) case root path
   (2) cam root path
   (3) cam configuration options
"""

import sys, os, shutil

case_root = sys.argv[1]
cam_root = sys.argv[2]
cam_options = sys.argv[3]

# If using GEOS-Chem chemistry then copy GEOS-Chem configuration files from source code to case
if '-chem geoschem' in cam_options:
    geoschem_src = os.path.join(cam_root,'src','chemistry','geoschem','geoschem_src')
    if not os.path.isdir(geoschem_src):
        raise SystemExit("ERROR: Did not find path to GEOS-Chem source code at {:s}".format(geoschem_src))
    for fname in ['species_database.yml', 'geoschem_config.yml',
                  'HISTORY.rc', 'HEMCO_Config.rc', 'HEMCO_Diagn.rc']:
        file1 = os.path.join(geoschem_src, 'run','CESM', fname)
        if not os.path.exists(file1):
            raise SystemExit("ERROR: GEOS-Chem configuration file does not exist: {}".format(file1))
        file2 = os.path.join(case_root, fname)
        shutil.copy(file1,file2)
