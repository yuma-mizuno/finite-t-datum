"""Check every stored query hash and archive the complete proof inputs."""
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED

HERE=Path(__file__).resolve().parent
report=json.loads((HERE/'verification.json').read_text())
assert report['replay_counts']=={'unsat':4865}
assert report['lift_status_counts']=={'sat':37,'unsat':4828}
assert len(report['queries'])==4865
with ZipFile(HERE/'smt_queries.zip','w',compression=ZIP_DEFLATED,compresslevel=9) as archive:
    for item in report['queries']:
        raw=(HERE/'smt_queries'/item['file']).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==item['sha256']
        info=ZipInfo(item['file'],date_time=(1980,1,1,0,0,0));info.compress_type=ZIP_DEFLATED
        archive.writestr(info,raw,compresslevel=9)
with ZipFile(HERE/'smt_queries.zip') as archive:
    assert len(archive.namelist())==4865
    for item in report['queries']:
        assert hashlib.sha256(archive.read(item['file'])).hexdigest()==item['sha256']
print('Verified and archived all 4865 query hashes;', (HERE/'smt_queries.zip').stat().st_size,'bytes.')
