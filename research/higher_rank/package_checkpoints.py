"""Deterministic, verified archives of completed checkpoint files."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
from query_archives import package

HERE=Path(__file__).resolve().parent


def archive(directory,target,pattern,selected=None):
    files=sorted(directory.glob(pattern),key=lambda p:p.name)
    if selected is not None:
        files=[p for p in files if p.name in selected]
        assert {p.name for p in files}==set(selected)
    assert files,(directory,pattern)
    temporary=target.with_suffix(target.suffix+'.tmp')
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in files:
            info=zipfile.ZipInfo(p.name,date_time=(1980,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(info,p.read_bytes())
    with zipfile.ZipFile(temporary) as z:
        assert z.testzip() is None and z.namelist()==[p.name for p in files]
        for p in files:assert z.read(p.name)==p.read_bytes()
    temporary.replace(target)
    return {'file':target.name,'members':len(files),'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}


def main():
    p=argparse.ArgumentParser();p.add_argument('rank',type=int);p.add_argument('stage',choices=['constants','queries','controls']);a=p.parse_args()
    directory=HERE/f'rank{a.rank}'
    if a.stage=='constants':
        assert json.loads((directory/'constant-progress.json').read_text())['complete']
        reports=[archive(directory/'constant_tasks',directory/'constant_tasks.zip','upper-*.json')]
    elif a.stage=='queries':
        queries=[json.loads(line) for line in (directory/'verification.jsonl').read_text().splitlines()]
        count=json.loads((directory/'constant_candidates.json').read_text())['count']
        assert len(queries)==len({q['id'] for q in queries})==count and all(q['result']=='unsat' for q in queries)
        import gzip
        for q in queries:
            raw=gzip.decompress((directory/'smt_queries'/q['query_file']).read_bytes())
            assert hashlib.sha256(raw).hexdigest()==q['sha256']
        distribution=package(directory,[q['query_file'] for q in queries],lambda name:(directory/'smt_queries'/name).read_bytes())
        reports=[{k:v for k,v in a.items() if k!='members'} for a in distribution['archives']]
    else:
        reports=[archive(HERE/name,HERE/(name+'.zip'),'upper-*.json')
                 for n in (4,5) for name in (f'control-rank{n}-columns',f'control-rank{n}-hereditary')]
    print(json.dumps(reports),flush=True)


if __name__=='__main__':main()
