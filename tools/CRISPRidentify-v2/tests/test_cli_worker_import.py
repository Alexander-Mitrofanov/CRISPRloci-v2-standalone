"""A spawned matcher must not re-enter CLI setup or alter its working directory."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class CliWorkerImportTests(unittest.TestCase):
    def test_worker_import_does_not_parse_reexec_import_pipeline_or_create_temp(self):
        probe = '''
import json, os, pathlib, runpy, sys
root = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(root))
sys.argv = ['foreign-worker', '--invalid-for-crispridentify', 'missing-relative.fa']
def forbidden(*args, **kwargs):
    raise AssertionError('worker import tried to re-exec')
os.execve = forbidden
before = os.getcwd()
names = runpy.run_path(str(root/'CRISPRidentify.py'), run_name='__mp_main__')
assert os.getcwd() == before
assert 'args' not in names and 'pid_work_directory' not in names
assert 'components.pipeline' not in sys.modules
assert not list(pathlib.Path(os.environ['TMPDIR']).iterdir())
print(json.dumps({'imported': True}))
'''
        for hash_seed in (None, '19', '0'):
            with self.subTest(hash_seed=hash_seed), tempfile.TemporaryDirectory() as directory:
                environment = dict(os.environ, TMPDIR=directory, PYTHONDONTWRITEBYTECODE='1')
                environment.pop('_CRISPRIDENTIFY_HASH_REEXEC', None)
                if hash_seed is None:
                    environment.pop('PYTHONHASHSEED', None)
                else:
                    environment['PYTHONHASHSEED'] = hash_seed
                result = subprocess.run([sys.executable, '-c', probe, str(ROOT)],
                                        cwd=directory, env=environment,
                                        text=True, capture_output=True, timeout=30)
                self.assertEqual(0, result.returncode, result.stderr)
                self.assertEqual({'imported': True}, json.loads(result.stdout))


if __name__ == '__main__':
    unittest.main()
