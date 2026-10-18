import importlib.util
import json
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch

source = Path('/home/alex/CodexProjects/CRISPRtracrRNA/new_version/scripts/dependencies.py')
spec = importlib.util.spec_from_file_location('dependencies_probe', source)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def git(path, *args):
    return subprocess.check_output(['git', '-C', str(path), *args], text=True, stderr=subprocess.DEVNULL).strip()

def commit(path, text):
    (path / 'runtime.py').write_text(text)
    git(path, 'add', 'runtime.py')
    git(path, '-c', 'user.name=Test', '-c', 'user.email=test@example.invalid', 'commit', '-m', text)
    return git(path, 'rev-parse', 'HEAD')

with tempfile.TemporaryDirectory(prefix='tracrrna-tooling-review-') as temporary:
    work = Path(temporary)
    upstream = work / 'upstream'
    upstream.mkdir()
    git(upstream, 'init', '-b', 'main')
    first = commit(upstream, 'v1')
    root = work / 'project'
    root.mkdir()
    m.ROOT = root
    lock = {'CasAndra': {'url': str(upstream), 'commit': first}}
    lock_path = root / 'upstream-lock.json'
    original_bytes = (json.dumps(lock) + '\n').encode()
    lock_path.write_bytes(original_bytes)
    m.sync('CasAndra', lock['CasAndra'])
    target = root / 'external/CasAndra'
    second = commit(upstream, 'v2')
    original_git = m.git
    original_unlink = Path.unlink
    for rollback_fails in (False, True):
        original_git('checkout', '--detach', first, cwd=target)
        def fail_replace(self, destination):
            raise OSError('PRIMARY: lock publication failed')
        def fail_unlink(self, *args, **kwargs):
            if self.name.startswith('.upstream-lock-'):
                raise PermissionError('CLEANUP: cannot remove staged lock')
            return original_unlink(self, *args, **kwargs)
        def conditional_git(*args, **kwargs):
            if rollback_fails and args == ('checkout', '--detach', first):
                raise OSError('RESTORATION: original HEAD cannot be restored')
            return original_git(*args, **kwargs)
        with patch.object(Path, 'replace', fail_replace), patch.object(Path, 'unlink', fail_unlink), patch.object(m, 'git', conditional_git):
            try:
                m.update('CasAndra', 'main', lock)
            except Exception as exc:
                print(json.dumps({
                    'rollback_fails': rollback_fails,
                    'surfaced_type': type(exc).__name__,
                    'surfaced_message': str(exc),
                    'context_type': type(exc.__context__).__name__,
                    'context_message': str(exc.__context__),
                    'primary_lock_unchanged': lock_path.read_bytes() == original_bytes,
                    'mapping_unchanged': lock['CasAndra']['commit'] == first,
                    'head_is_original': git(target, 'rev-parse', 'HEAD') == first,
                    'head_is_candidate': git(target, 'rev-parse', 'HEAD') == second,
                    'staged_lock_count': len(list(root.glob('.upstream-lock-*')))
                }, sort_keys=True))
        for staged in root.glob('.upstream-lock-*'):
            staged.unlink()
