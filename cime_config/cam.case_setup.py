#! /usr/bin/env python3

"""Copy GEOS-Chem configuration files from source to the case directory.
This script is run from CIME when calling case.setup"""

import logging
import os
import shutil
import sys

_CIMEROOT = os.environ.get("CIMEROOT")
if _CIMEROOT is None:
    raise SystemExit("ERROR: must set CIMEROOT environment variable")
# end if
_LIBDIR = os.path.join(_CIMEROOT, "CIME", "Tools")
sys.path.append(_LIBDIR)
sys.path.insert(0, _CIMEROOT)

#pylint: disable=wrong-import-position
from CIME.case import Case

logger = logging.getLogger(__name__)

if len(sys.argv) != 3:
    raise SystemExit(f"Incorrect call to {sys.argv[0]}, need CAM root and case root")
# end if
cam_root = sys.argv[1]
case_root = sys.argv[2]

with Case(case_root) as case:
    cam_config = case.get_value('CAM_CONFIG_OPTS')
    # Gather case information (from _build_usernl_files in case_setup.py)
    comp_interface = case.get_value("COMP_INTERFACE")

    if comp_interface == "nuopc":
        ninst = case.get_value("NINST")
    elif ninst == 1:
        ninst = case.get_value("NINST_CAM")
    # end if
# end with

# GEOS-Chem only: copy config files to case
if '-chem geoschem' in cam_config:
    geoschem_config_src = os.path.join(cam_root, 'src', 'chemistry',
                                       'geoschem', 'geoschem_src', 'run', 'CESM')
    if not os.path.isdir(geoschem_config_src):
        raise SystemExit(f"ERROR: Did not find path to GEOS-Chem config files at {geoschem_config_src}")
    # end if
    for fileName in ['species_database.yml', 'geoschem_config.yml', 'HISTORY.rc',
                     'HEMCO_Config.rc', 'HEMCO_Diagn.rc']:
        source_file = os.path.join(cam_root, geoschem_config_src, fileName)
        if not os.path.exists(source_file):
            raise SystemExit(f"ERROR: Did not find source file, {source_file}")
        # end if
        spaths = os.path.splitext(source_file)
        for inst_num in range(ninst):
            if ninst > 1:
                target_file = f"{spaths[0]}_{inst_num+1:04d}{spaths[1]}"
            else:
                target_file = os.path.join(case_root, fileName)
            # end if
            if not os.path.exists(target_file):
                logger.info("CAM namelist one-time copy of GEOS-Chem run directory files: source_file %s target_file %s ",
                            source_file, target_file)
                shutil.copy(source_file, target_file)
            # end if
        # end for
    # end for
# end if
