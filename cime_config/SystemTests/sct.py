"""
Implementation of the CAM single column test.

This is a CAM specific test:
Verifies that single column is working by checking that T and Q
only have round off differences.

(1) 3D
(2) scam

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
from CIME.test_status import *
from CIME.utils import append_testlog


logger = logging.getLogger(__name__)

class SCT(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = 'default',
                                       run_one_description = '3D CAM run',
                                       run_two_description = 'SCAM run')

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "inithist = 'CAMIOP'")

        CAM_CONFIG_OPTS = self._case1.get_value("CAM_CONFIG_OPTS")



    def _case_two_setup(self):
        case_name = self._case.get_value("CASE")
        RUN_STARTDATE = self._case1.get_value("RUN_STARTDATE")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "ncdata = '../"+case_name+".cam.i."+RUN_STARTDATE+"-00000.nc'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "NDENS    = 1,1,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "MFILT    = 1,7,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "nhtfrq   = 1,1,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "fincl2='T','Q','TDIFF','QDIFF','LANDFRAC'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "iopfile = '../"+case_name+".cam.h1."+RUN_STARTDATE+"-00000.nc'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "inithist = 'YEARLY'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scm_cambfb_mode                = .true.")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scm_use_obs_uv         = .true.")
        for comp in self._case.get_values("COMP_CLASSES"):
            self._case.set_value("NTASKS_{}".format(comp), 1)
            self._case.set_value("NTHRDS_{}".format(comp), 1)
            self._case.set_value("ROOTPE_{}".format(comp), 0)

        self._case.set_value("PTS_MODE","TRUE")
        self._case.set_value("PTS_LAT",-20.0)
        self._case.set_value("PTS_LON",140.0)

        CAM_CONFIG_OPTS = self._case1.get_value("CAM_CONFIG_OPTS")
        self._case.set_value("CAM_CONFIG_OPTS","{} -scam ".format(CAM_CONFIG_OPTS))
        self._case.case_setup(test_mode=True, reset=True)

    def _component_compare_test(self, suffix1, suffix2,
                                success_change=False,
                                ignore_fieldlist_diffs=False):
        with self._test_status:
            stat,netcdf_filename,err=run_cmd('ls ./run/case2run/*h1*8400.nc ')
            stat,DIFFs,err=run_cmd('ncdump -ff -p 9,17 -v QDIFF,TDIFF '+netcdf_filename+' | egrep //\.\*DIFF | sed s/^\ \*// | sed s/\[,\;\].\*// | sed s/^0/0.0/  | uniq')
            array_of_DIFFs=DIFFs.split("\n")
            answer=max([abs(float(x)) for x in array_of_DIFFs])
            comments="Checking QDIFF,TDIFF in SCAM run."
            append_testlog(comments, self._orig_caseroot)
            # Test for greater that round off changes.
            if answer < 1e-10:
                self._test_status.set_status("{}_{}_{}".format(COMPARE_PHASE, self._run_one_suffix, self._run_two_suffix), TEST_PASS_STATUS)
                comments="QDIFF,TDIFF: PASS"
            else:
                self._test_status.set_status("{}_{}_{}".format(COMPARE_PHASE, self._run_one_suffix, self._run_two_suffix), TEST_FAIL_STATUS)
                comments="QDIFF,TDIFF: Difference greater than round off."
            append_testlog(comments, self._orig_caseroot)
