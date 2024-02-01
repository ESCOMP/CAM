"""
CIME MGP test.  This class inherits from SystemTestsCompareTwo

This is a changing config options test to compare between MG3 and 
PUMAS in camdev. The use of MG3 or PUMAS should be bfb.
This is just like an ERC test and it's meant for CAM only 
as it only does a single build.
 
(1) Do an initial run with microphys setup as MG3 (suffix MG3)
(2) Do an initial run with microphys setup as PUMAS (suffix PUMAS)
"""

import sys
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)

class MGP(SystemTestsCompareTwo):

    def __init__(self, case,
        separate_builds=True,
        run_one_suffix="mg3",
        run_two_suffix="pumas",
        run_one_description="MG3 microphysics",
        run_two_description="PUMAS microphysics",
        multisubmit=False,
        **kwargs
    ):
        SystemTestsCompareTwo.__init__(self, case,
                             separate_builds=separate_builds,
                             run_one_suffix=run_one_suffix,
                             run_two_suffix=run_two_suffix,
                             run_one_description=run_one_description,
                             run_two_description=run_two_description,
                             multisubmit=multisubmit,
                             **kwargs
                        )
    def _case_one_setup(self):
        stop_n = self._case1.get_value("STOP_N")
        expect(stop_n >= 3, "STOP_N must be at least 3, STOP_N = {}".format(stop_n))
        self._case.set_value("CAM_CONFIG_OPTS","-phys cam_dev -microphys mg3")

    def _case_two_setup(self):
        self._case.set_value("CAM_CONFIG_OPTS","-phys cam_dev -microphys pumas")


