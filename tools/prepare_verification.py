"""Restore losslessly compressed source files for the computational programs."""
from pathlib import Path
import gzip
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifest = json.loads((ROOT / 'distribution.json').read_text())
    for item in manifest['compressed_sources']:
        archive = (ROOT / item['archive']).read_bytes()
        assert hashlib.sha256(archive).hexdigest() == item['archive_sha256']
        raw = gzip.decompress(archive)
        assert len(raw) == item['bytes']
        assert hashlib.sha256(raw).hexdigest() == item['sha256']
        target = ROOT / item['path']
        if target.exists():
            if target.read_bytes() != raw:
                raise RuntimeError(f'Refusing to overwrite changed data: {item["path"]}')
        else:
            target.write_bytes(raw)
        print('Verified and restored', item['path'], flush=True)
    # Control verifiers read these archived task directories directly.
    archives = list((ROOT / 'research/symmetrizable').glob('rank*/constant_tasks.zip'))
    archives += [ROOT / 'research/symmetrizable/control-outputs.zip']
    for path in archives:
        destination = path.parent / path.stem if path.name == 'constant_tasks.zip' else path.parent
        with zipfile.ZipFile(path) as archive:
            for name in archive.namelist():
                if name == 'manifest.json':
                    continue
                target = (destination / name).resolve()
                if not target.is_relative_to(destination.resolve()):
                    raise ValueError(f'Unsafe archive member: {name}')
                raw = archive.read(name)
                if target.exists():
                    if target.read_bytes() != raw:
                        raise RuntimeError(f'Refusing to overwrite changed task: {name}')
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(raw)
        print('Restored verification inputs from', path.relative_to(ROOT), flush=True)


if __name__ == '__main__':
    main()
