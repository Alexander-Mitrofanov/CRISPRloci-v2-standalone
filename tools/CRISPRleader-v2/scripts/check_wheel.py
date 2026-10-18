"""Install a built wheel offline and check the CLI outside every source checkout."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import venv


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    wheels = sorted((root / 'dist').glob('crisprleader-*.whl'))
    if len(wheels) != 1:
        raise SystemExit('Expected exactly one wheel in dist; remove stale builds')
    with tempfile.TemporaryDirectory(prefix='crisprleader-wheel-') as temporary:
        cwd = Path(temporary)
        environment = cwd / 'venv'
        venv.EnvBuilder(with_pip=True).create(environment)
        python = environment / 'bin' / 'python'
        subprocess.run([str(python), '-m', 'pip', 'install', '--no-index', '--no-deps', str(wheels[0])], check=True, cwd=cwd)
        for name in ('synthetic.fa', 'synthetic.arrays.json', 'benchmark.synthetic.dataset.json',
                     'benchmark.synthetic.predictions.json'):
            shutil.copyfile(root / 'examples' / name, cwd / name)
        base = [str(python), '-I', '-m', 'crisprleader']
        subprocess.run(base + ['--version'], check=True, cwd=cwd)
        subprocess.run(base + ['extract', '--fasta', 'synthetic.fa', '--arrays', 'synthetic.arrays.json',
                              '--output', 'result', '--flank-length', '12'], check=True, cwd=cwd)
        subprocess.run(base + ['validate', '--fasta', 'synthetic.fa', '--result', 'result/result.json'], check=True, cwd=cwd)
        manifest = json.loads((cwd / 'result/manifest.json').read_text())
        for name, receipt in manifest['files'].items():
            payload = (cwd / 'result' / name).read_bytes()
            assert len(payload) == receipt['size_bytes']
            assert hashlib.sha256(payload).hexdigest() == receipt['sha256']
        failure = subprocess.run(base + ['predict', '--fasta', 'synthetic.fa', '--arrays', 'synthetic.arrays.json',
                                         '--output', 'prediction'], cwd=cwd, capture_output=True, text=True)
        assert failure.returncode == 2 and not (cwd / 'prediction').exists()
        assert 'unavailable' in json.loads(failure.stderr)['message']
        subprocess.run(base + ['benchmark', '--dataset', 'benchmark.synthetic.dataset.json',
                              '--predictions', 'benchmark.synthetic.predictions.json',
                              '--output', 'benchmark.json', '--bootstrap-replicates', '100'],
                       check=True, cwd=cwd)
        assert json.loads((cwd / 'benchmark.json').read_text())['scientific_superiority_established'] is False
        subprocess.run([str(python), '-I', '-c',
                        'import importlib.resources as r; import crisprleader.inference.motifs; '
                        'assert all(r.files("crisprleader").joinpath(name).is_file() for name in '
                        '["data/REFERENCE_SCHEMA.md", "schemas/annotations.schema.json", '
                        '"schemas/context-result.schema.json"])'],
                       check=True, cwd=cwd)
    print('Offline installed-wheel CLI checks passed; biological gates remain unevaluated.')


if __name__ == '__main__':
    main()
