"""Apply the predeclared RNG-only patch and preserve an exact build manifest."""
from pathlib import Path
import difflib
import hashlib
import json
import subprocess

HERE = Path(__file__).resolve().parent
TREE = HERE / 'fasta36-fasta-v36.3.8g'
archive = HERE / 'fasta-v36.3.8g.tar.gz'
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert digest(archive) == 'fa5318b6f8d6a3cfdef0d29de530eb005bfd3ca05835faa6ad63663f8dce7b2e'
header = '''/* CRISPRtracrRNA council amendment, 2026-09-27.
 * Explicit, untuned RNG seeds for the tested Linux pthread, -T1 profile.
 * The RNG recurrence, shuffling and statistical thresholds are unchanged.
 */
#ifndef CRISPRTRACRRNA_COUNCIL_SEED_H
#define CRISPRTRACRRNA_COUNCIL_SEED_H
#ifdef PCOMPLIB
#error "Council seed profile does not support MPI builds"
#endif
#define FASTA_COUNCIL_MASTER_SEED 20260927
void *my_srand(int set);
#endif
'''
changes = {}
for name in ('comp_lib9.c', 'work_thr2.c'):
    path = TREE / 'src' / name
    before = path.read_text()
    if name == 'comp_lib9.c':
        assert before.count('void *my_srand();') == 1
        after = before.replace('void *my_srand();', '#include "council_seed.h"')
        assert after.count('rand_state = my_srand();') == 1
        after = after.replace('rand_state = my_srand();', 'rand_state = my_srand(FASTA_COUNCIL_MASTER_SEED);')
    else:
        assert before.count('extern void *my_srand();') == 1
        after = before.replace('extern void *my_srand();', '#include "council_seed.h"')
        assert after.count('  my_rand_state=my_srand();\n') == 1
        after = after.replace('  my_rand_state=my_srand();\n', '')
        assert after.count('  my_worker = work_info->worker;') == 1
        after = after.replace('  my_worker = work_info->worker;', '''  my_worker = work_info->worker;
  if (my_worker < 0 || my_worker > 2147483645 - FASTA_COUNCIL_MASTER_SEED) {
    fprintf(stderr, "Council FASTA worker seed exceeds RNG range\\n");
    exit(1);
  }
  my_rand_state = my_srand(FASTA_COUNCIL_MASTER_SEED + my_worker + 1);''')
    changes[name] = {'before_sha256': hashlib.sha256(before.encode()).hexdigest(), 'after_sha256': hashlib.sha256(after.encode()).hexdigest()}
    path.write_text(after)
    with (HERE / 'explicit-seed.patch').open('a') as out:
        out.writelines(difflib.unified_diff(before.splitlines(True), after.splitlines(True), fromfile='a/src/' + name, tofile='b/src/' + name))
(TREE / 'src/council_seed.h').write_text(header)
with (HERE / 'explicit-seed.patch').open('a') as out:
    out.writelines(difflib.unified_diff([], header.splitlines(True), fromfile='/dev/null', tofile='b/src/council_seed.h'))
command = ['make', '-f', '../make/Makefile.linux_sse2', 'LIB_M=-lm /usr/lib/x86_64-linux-gnu/libz.so.1', 'fasta36']
manifest = {'source_url': 'https://github.com/wrpearson/fasta36/archive/refs/tags/fasta-v36.3.8g.tar.gz',
            'source_archive_sha256': digest(archive), 'master_seed': 20260927,
            'worker_seed': '20260927 + zero_based_worker_index + 1', 'supported_native_threads': 1,
            'scope': 'Linux SSE2 pthread build only; no MPI or arbitrary-thread determinism claim',
            'patch_sha256': digest(HERE / 'explicit-seed.patch'), 'changes': changes,
            'source_files_sha256': {str(p.relative_to(TREE)): digest(p) for p in sorted(TREE.rglob('*')) if p.is_file()},
            'compiler': subprocess.check_output(['gcc', '--version'], text=True),
            'make_version': subprocess.check_output(['make', '--version'], text=True),
            'build_command': command, 'build_cwd': str(TREE / 'src')}
(HERE / 'build_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
with (HERE / 'build.log').open('w') as handle:
    result = subprocess.run(command, cwd=TREE / 'src', stdout=handle, stderr=subprocess.STDOUT)
manifest['build_exit_code'] = result.returncode
if result.returncode == 0:
    binary = HERE / 'fasta36-council-seed20260927'
    binary.write_bytes((TREE / 'bin/fasta36').read_bytes())
    binary.chmod(0o755)
    manifest['binary_sha256'] = digest(binary)
    manifest['binary'] = str(binary)
manifest['build_log_sha256'] = digest(HERE / 'build.log')
(HERE / 'build_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
raise SystemExit(result.returncode)
