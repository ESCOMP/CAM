"""
Implementation of the CAM sub column test.

This is a CAM specific test:
Verifies that turning sub columns on and off doesn't change answers
(1) do a run with sub columns on
(2) Do a run with sub columns off

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
from CIME.utils import append_testlog


logger = logging.getLogger(__name__)

class SUB(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'default',
                                       run_one_description = 'Sub columns on',
                                       run_two_description = 'Sub columns off')

    def _case_one_setup(self):
        cam_config_opts = self._case.get_value("CAM_CONFIG_OPTS")
        cam_config_opts = self._case.set_value("CAM_CONFIG_OPTS",cam_config_opts+" -psubcols 3")
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "pbuf_global_allocate = .false.")

        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "microp_uniform = .false.")

        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "use_subcol_microp = .true.")

        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "subcol_scheme='tstcp'")

        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "subcol_tstcp_noAvg= .true.")
        comments = "Sub columns on."
        append_testlog(comments, self._orig_caseroot)

    def _case_two_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "cam",
                                contents = "pbuf_global_allocate=.false.")
        comments = "Sub columns off."
        append_testlog(comments, self._orig_caseroot)

