"""
CAM mass conservation test  This class inherits from SystemTestsCommon
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.test_status import *
from CIME.utils import append_testlog
import glob, gzip


logger = logging.getLogger(__name__)

class TMC(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the TMC system test
        """
        SystemTestsCommon.__init__(self, case)

    def run_phase(self):

        self.run_indv()
        cesmlog = ''.join(self._get_latest_cpl_logs())
        if '.gz' == cesmlog[-3:]:
            fopen = gzip.open
        else:
            fopen = open

        with fopen(cesmlog, "rb") as f:
            lines = [line.rstrip('\n') for line in f]
        first_val = 0.0
        with self._test_status:
            self._test_status.set_status(RUN_PHASE, TEST_PASS_STATUS)
        for line in lines:
            if re.search('between DRY m=30 name=TT_UN gavg dry, wet, min, max',line):
                if first_val == 0.0:
                    first_val = re.findall('\s*[\d]+ name=TT_UN [^0-9]+([\S]+)',line)
                if first_val != re.findall('\s*[\d]+ name=TT_UN [^0-9]+([\S]+)',line):
                    with self._test_status:
                        self._test_status.set_status(RUN_PHASE, TEST_FAIL_STATUS)
                    comments = "CAM mass conservation test FAILED."
                    CIME.utils.append_testlog(comments, self._orig_caseroot)
            first_var = 1.1

    def _get_latest_cesm_logs(self):
        """
        find and return the latest cesm log file in the run directory
        """
        coupler_log_path = self._case.get_value("RUNDIR")
        cesmlogs = glob.glob(os.path.join(coupler_log_path, 'cesm*.log.*'))
        lastcesmlogs = []
        if cesmlogs:
            lastcesmlogs.append(max(cesmlogs, key=os.path.getctime))
            basename = os.path.basename(lastcesmlogs[0])
            suffix = basename.split('.',1)[1]
            for log in cesmlogs:
                if log in lastcesmlogs:
                    continue

                if log.endswith(suffix):
                    lastcesmlogs.append(log)

        return lastcesmlogs

