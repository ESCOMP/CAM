"""
Implementation of the CAM physics load balancing test.

This is a CAM specific test:
Verifies that changeing pyhysics load balancing doesn't change answers
(1) do a run with the default physics load balancing 
(2) Do a run with the new physics load balancing from the namelist change 

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
from CIME.utils import append_testlog

logger = logging.getLogger(__name__)

class PLB(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       ignore_fieldlist_diffs = True,
                                       run_two_suffix = 'default',
                                       run_one_description = 'Default phys_loadbalance',
                                       run_two_description = 'Changed phys_loadbalance')

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "phys_loadbalance = 2")
        comments = "Setting phys_loadbalance to default 2."
        append_testlog(comments, self._orig_caseroot)

    def _case_two_setup(self):
        CAM_CONFIG_OPTS = self._case1.get_value("CAM_CONFIG_OPTS")
        self._case.set_value("CAM_CONFIG_OPTS",CAM_CONFIG_OPTS.replace(' -nadv_tt 5 -cppdefs -DTRACER_CHECK',''))
        comments = "Setting phys_loadbalance to value that's set in the testmod."
        append_testlog(comments, self._orig_caseroot)


