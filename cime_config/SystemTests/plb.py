"""
Implementation of the CAM physics load balancing test.

This is a CAM specific test:
Verifies that changing physics load balancing doesn't change answers
(1) do a run with the default physics load balancing
(2) Do a run with the new physics load balancing from the namelist change

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
from CIME.status import append_testlog

logger = logging.getLogger(__name__)

class PLB(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       ignore_fieldlist_diffs = False,
                                       run_two_suffix = 'default',
                                       run_one_description = 'Default phys_loadbalance',
                                       run_two_description = 'Changed phys_loadbalance')

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "phys_loadbalance = 2")
        comments = "Overriding phys_loadbalance to 2 (usual default)."
        append_testlog(comments, self._orig_caseroot)

    def _case_two_setup(self):
        comments = "Leave phys_loadbalance at value that's set in the testmod."
        append_testlog(comments, self._orig_caseroot)
