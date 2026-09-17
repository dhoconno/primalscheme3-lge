#!/usr/bin/env python3
"""Build a source-bound LGE wheel from a clean git commit, without network access.

Prepare tooling with: python -m pip install -r requirements-build.txt
Run twice into separate directories and compare SHA256SUMS.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tarfile
import tempfile
import tomllib

UPSTREAM_COMMIT = 'dd13ec5cb1cf375f052640355c73101c0c4bf839'


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args], text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if git(root, 'status', '--porcelain'):
        raise SystemExit('Build requires a clean committed source tree')
    tools = {}
    for line in (root / 'requirements-build.txt').read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        name, expected = line.split('==')
        actual = importlib.metadata.version(name)
        if actual != expected:
            raise SystemExit(f'Build tooling mismatch: {name} expected {expected}, found {actual}')
        tools[name] = actual
    commit = git(root, 'rev-parse', 'HEAD')
    epoch = int(os.environ.get('SOURCE_DATE_EPOCH', git(root, 'show', '-s', '--format=%ct', 'HEAD')))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise SystemExit('Output directory must be empty')
    project = tomllib.loads((root / 'pyproject.toml').read_text())['project']
    version = project['version']
    metadata = {
        'schemaVersion': 1,
        'distribution': 'primalscheme3',
        'version': version,
        'sourceRepository': 'https://github.com/dhoconno/primalscheme3-lge',
        'sourceCommit': commit,
        'upstreamRepository': 'https://github.com/artic-network/primalscheme3',
        'upstreamCommit': UPSTREAM_COMMIT,
        'sourceDateEpoch': epoch,
        'buildPython': platform.python_version(),
        'buildTools': tools,
    }
    with tempfile.TemporaryDirectory(prefix='primalscheme-lge-build-') as temporary:
        stage = Path(temporary)
        archive = stage / 'source.tar'
        with archive.open('wb') as handle:
            subprocess.run(['git', '-C', str(root), 'archive', '--format=tar', 'HEAD'], stdout=handle, check=True)
        source = stage / 'source'
        source.mkdir()
        with tarfile.open(archive) as handle:
            handle.extractall(source, filter='data')
        (source / 'primalscheme3' / 'lge-build.json').write_text(json.dumps(metadata, indent=2, sort_keys=True) + '\n')
        environment = os.environ.copy()
        environment.update(SOURCE_DATE_EPOCH=str(epoch), PYTHONHASHSEED='0')
        subprocess.run([sys.executable, '-m', 'hatchling', 'build', '-t', 'wheel', '-d', str(output)], cwd=source, env=environment, check=True)
    wheels = list(output.glob('*.whl'))
    if len(wheels) != 1:
        raise SystemExit('Expected exactly one wheel')
    wheel = wheels[0]
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    metadata.update(wheelFilename=wheel.name, wheelSHA256=digest, wheelSize=wheel.stat().st_size)
    (output / 'build-evidence.json').write_text(json.dumps(metadata, indent=2, sort_keys=True) + '\n')
    (output / 'SHA256SUMS').write_text(f'{digest}  {wheel.name}\n')
    print(json.dumps(metadata, sort_keys=True))


if __name__ == '__main__':
    main()
