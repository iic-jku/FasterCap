"""Small end-to-end checks that the external safety guard actually stops work."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

GUARD = Path(__file__).resolve().parents[1] / 'tools' / 'run_guarded.py'


class ResourceGuardTests(unittest.TestCase):
    def run_guard(self, code, *options):
        with tempfile.TemporaryDirectory() as root:
            result = subprocess.run([sys.executable, str(GUARD), '--root', root,
                                     '--seconds', '5', *options, '--',
                                     sys.executable, '-c', code], capture_output=True,
                                    text=True, timeout=12)
            report = json.loads(result.stdout.splitlines()[0])
            self.assertFalse([p for p in Path(root).glob('run-*') if p.is_dir()])
            return result, report

    def test_success_and_child_status(self):
        result, report = self.run_guard('print("bounded output")')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIsNone(report['guard'])
        result, _ = self.run_guard('raise SystemExit(7)')
        self.assertEqual(result.returncode, 7)

    def test_timeout(self):
        result, report = self.run_guard('import time; time.sleep(10)', '--seconds', '0.3')
        self.assertEqual(result.returncode, 125)
        self.assertEqual(report['guard'], 'elapsed time limit')

    def test_memory_limit(self):
        result, report = self.run_guard('import time; a=bytearray(32*1024*1024); time.sleep(3)',
                                        '--rss-mib', '24')
        self.assertEqual(result.returncode, 125)
        self.assertEqual(report['guard'], 'process-group RSS limit')

    def test_log_limit(self):
        result, report = self.run_guard('print("x"*100000)', '--log-mib', '0.01')
        self.assertEqual(result.returncode, 125)
        self.assertEqual(report['guard'], 'log output limit')

    def test_disk_preflight(self):
        result, report = self.run_guard('raise SystemExit("must not start")',
                                        '--reserve-mib', '1000000000')
        self.assertEqual(result.returncode, 125)
        self.assertEqual(report['guard'], 'free disk below reserve')

    def test_temp_limit(self):
        code = ('import os,time; from pathlib import Path; '
                'Path(os.environ["TMPDIR"],"data").write_bytes(b"x"*2000000); time.sleep(3)')
        result, report = self.run_guard(code, '--temp-mib', '1')
        self.assertEqual(result.returncode, 125)
        self.assertEqual(report['guard'], 'temporary-file limit')


if __name__ == '__main__':
    unittest.main()
