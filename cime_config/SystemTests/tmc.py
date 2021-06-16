"""
CAM mass conservation test  This class inherits from SystemTestsCommon
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.test_status import *
from CIME.utils import append_testlog
import glob, gzip


class TMC(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the TMC system test
        """
        SystemTestsCommon.__init__(self, case)

    def run_phase(self):

        with self._test_status:
            self._test_status.set_status("COMPARE_MASS", TEST_PEND_STATUS)
        self.run_indv()
        cpllog = ''.join(self._get_latest_cpl_logs())
        atmlog  = cpllog.replace("cpl.log","atm.log")
        if '.gz' == atmlog[-3:]:
            fopen = gzip.open
        else:
            fopen = open

        f = fopen(atmlog,'r')
        lines = f.readlines()
        first_val = -9.0
        with self._test_status:
            self._test_status.set_status("COMPARE_MASS", TEST_PASS_STATUS)
        use_this_tt_un = False 
        for line in lines:
            if re.search('vvvvv gmean_mass: before tphysbc DRY',line.decode('utf-8')): 
                use_this_tt_un = True
            if re.search('TT_UN ',line.decode('utf-8')) and use_this_tt_un:
                tt_un_flt=re.findall("\d+\.\d+",line.decode('utf-8'))
                if len(tt_un_flt) > 0:
                    if first_val == -9.0:
                        first_val = tt_un_flt[0]
                    if first_val != tt_un_flt[0]:
                        with self._test_status:
                            self._test_status.set_status("COMPARE_MASS", TEST_FAIL_STATUS, comments="Mass Not Conserved")
                        comments = "CAM mass conservation test FAILED."
                        append_testlog(comments, self._orig_caseroot)
                use_this_tt_un = False 
        if first_val == -9.0:
                with self._test_status:
                    self._test_status.set_status("COMPARE_MASS", TEST_FAIL_STATUS, comments="Failed to find TT_UN in atm.log")
                comments = "CAM mass conservation test FAILED."
                append_testlog(comments, self._orig_caseroot)

