"""Distribute compressed SMT queries in independently readable ZIP archives."""
import argparse
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import zipfile

HERE = Path(__file__).resolve().parent
LIMIT = 48 * 1024 * 1024


@lru_cache(maxsize=None)
def member_archives(directory):
    manifest = directory / 'smt_queries_manifest.json'
    if not manifest.exists():
        return None
    data = json.loads(manifest.read_text())
    result = {}
    for archive in data['archives']:
        for member in archive['members']:
            if member in result:
                raise ValueError(f'Duplicate query: {member}')
            result[member] = directory / archive['file']
    return result


def archive_path(directory, member):
    mapping = member_archives(directory)
    return mapping[member] if mapping is not None else directory / 'smt_queries.zip'


def package(directory, names, read_member, original_sha256=None):
    """Preserve each .gz member exactly; verify each resulting ZIP after writing."""
    groups, group, size = [], [], 22
    for name in sorted(names):
        if Path(name).name != name or not name.endswith('.smt2.gz'):
            raise ValueError(f'Invalid query filename: {name}')
        raw = read_member(name)
        cost = len(raw) + 76 + 2 * len(name.encode('utf-8'))
        if cost + 22 > LIMIT:
            raise ValueError(f'Query exceeds archive size limit: {name}')
        if group and size + cost > LIMIT:
            groups.append(group)
            group, size = [], 22
        group.append((name, raw))
        size += cost
    if group:
        groups.append(group)
    if not groups:
        raise ValueError('No queries to package')
    archives = []
    for number, members in enumerate(groups, 1):
        name = 'smt_queries.zip' if len(groups) == 1 else f'smt_queries-{number:02d}.zip'
        target = directory / name
        temporary = target.with_suffix('.zip.tmp')
        with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_STORED) as archive:
            for member, raw in members:
                info = zipfile.ZipInfo(member, date_time=(1980, 1, 1, 0, 0, 0))
                archive.writestr(info, raw)
        with zipfile.ZipFile(temporary) as archive:
            assert archive.namelist() == [member for member, _ in members]
            for member, raw in members:
                assert archive.read(member) == raw
        assert temporary.stat().st_size <= LIMIT
        temporary.replace(target)
        archives.append({'file': name, 'bytes': target.stat().st_size,
                         'sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
                         'members': [member for member, _ in members]})
    manifest = {'schema_version': 1, 'query_count': len(names), 'archives': archives}
    if original_sha256:
        manifest['original_archive_sha256'] = original_sha256
    (directory / 'smt_queries_manifest.json').write_text(
        json.dumps(manifest, indent=2) + '\n', encoding='utf-8', newline='\n')
    member_archives.cache_clear()
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rank', type=int)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--source', type=Path)
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    directory = HERE / f'rank{args.rank}'
    if args.verify:
        manifest = json.loads((directory / 'smt_queries_manifest.json').read_text())
        queries = [json.loads(line) for line in (directory / 'verification.jsonl').read_text().splitlines()]
        expected = {query['query_file']: query for query in queries}
        assert len(expected) == len(queries) == manifest['query_count']
        seen = set()
        for part in manifest['archives']:
            path = directory / part['file']
            assert path.stat().st_size == part['bytes'] <= LIMIT
            assert hashlib.sha256(path.read_bytes()).hexdigest() == part['sha256']
            with zipfile.ZipFile(path) as archive:
                assert archive.namelist() == part['members']
                for member in archive.namelist():
                    assert member not in seen
                    seen.add(member)
                    query = expected[member]
                    assert query['result'] == 'unsat'
                    raw = gzip.decompress(archive.read(member))
                    assert hashlib.sha256(raw).hexdigest() == query['sha256']
            print(f'Verified {path.name}: {len(part["members"])} queries.', flush=True)
        assert seen == set(expected)
        print(f'PASS: {len(seen)} exact query hashes and complete archive coverage.')
        return
    with zipfile.ZipFile(args.source) as source:
        report = package(directory, source.namelist(), source.read,
                         hashlib.sha256(args.source.read_bytes()).hexdigest())
    print(json.dumps({'queries': report['query_count'],
                      'archives': [{k: v for k, v in a.items() if k != 'members'}
                                   for a in report['archives']]}, indent=2))


if __name__ == '__main__':
    main()
