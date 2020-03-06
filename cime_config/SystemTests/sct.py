"""
Implementation of the CAM single column test.
not working yet

This is a CAM specific test:
Verifies that single column is working
(1) 3d
(2) scam 

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
from CIME.XML.machines import Machines


logger = logging.getLogger(__name__)

class SCT(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = True,
                                       run_two_suffix = 'default',
                                       run_one_description = 'Default phys_loadbalance',
                                       run_two_description = 'Changed phys_loadbalance')

    def _case_one_setup(self):
#        self._case.set_value("phys_loadbalance".cam,2)
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "NDENS = 1,1phys_loadbalance = 2")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "MFILT = 1,10phys_loadbalance = 2")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "nhtfrq = 0,1phys_loadbalance = 2")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "inithist = 'CAMIOP'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "inithist_all = .true.")

        CAM_CONFIG_OPTS = self._case1.get_value("CAM_CONFIG_OPTS")
        self._case.set_value("CAM_CONFIG_OPTS","{} -camiop".format(CAM_CONFIG_OPTS))



    def _case_two_setup(self):
        case_name = self._case.get_value("CASE")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "ncdata         = '../"+case_name+".cam.i.0000-09-01-00000.nc'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "NDENS    = 1,1,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "MFILT    = 1,7,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "nhtfrq   = 1,1,1,1,1,1")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "fincl2='TDIFF','QDIFF','LANDFRAC'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scmlon= 140.")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scmlat= -20.")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "iopfile = '../"+case_name+".cam.h1.0000-09-01-00000.nc'")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "inithist_all = .true.")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scm_cambfb_mode                = .true.")
        append_to_user_nl_files(caseroot = self._get_caseroot(), component = "cam", contents = "scm_use_obs_uv         = .true.")
        for comp in self._case.get_values("COMP_CLASSES"):
            self._case.set_value("NTASKS_{}".format(comp), 1)
            self._case.set_value("NTHRDS_{}".format(comp), 1)
            self._case.set_value("ROOTPE_{}".format(comp), 0)

        mach_name = self._case.get_value("MACH")
        mach_obj = Machines(machine=mach_name)
#        if mach_obj.is_valid_MPIlib("mpi-serial"):
#            self._case.set_value("MPILIB","mpi-serial")

        self._case.set_value("PTS_MODE","TRUE")
        self._case.set_value("PTS_LAT",-20.0)
        self._case.set_value("PTS_LON",140.0)

        CAM_CONFIG_OPTS = self._case1.get_value("CAM_CONFIG_OPTS")
        self._case.set_value("CAM_CONFIG_OPTS","{} -scam".format(CAM_CONFIG_OPTS))
        self._case.case_setup(test_mode=True, reset=True)



