"""Check source links, hashes, generated record consistency and HTML packaging."""
import hashlib
import json
from pathlib import Path
import re
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
data=json.loads((HERE/'catalogue.json').read_text(encoding='utf-8'))
schema=json.loads((HERE/'record.schema.json').read_text(encoding='utf-8'))
try:
    import jsonschema
except ImportError:
    jsonschema=None
if jsonschema is not None:
    # Older jsonschema releases cannot resolve nested fragment references under
    # a URN base. A file URI gives the same local schema an offline-resolvable base.
    local_schema=dict(schema, **{'$id':(HERE/'record.schema.json').as_uri()})
    validator=jsonschema.Draft202012Validator(local_schema)
    validator.check_schema(local_schema)
    for record in data['records']:validator.validate(record)
    print('PASS: JSON Schema draft 2020-12 validation for all 61 records.')
html=(HERE/'index.html').read_text(encoding='utf-8')
embedded=re.search(r'<script id="catalogue-data" type="application/json">(.*?)</script>',html,re.S).group(1)
assert json.loads(embedded)==data
assert not re.search(r'<script[^>]+src=|<link[^>]+href=',html)
assert not re.search(r'/\*__[A-Z]+__\*/',html)
files=list((HERE/'records').glob('*.json'))
assert len(files)==len(data['records'])==61
hashes={}
with zipfile.ZipFile(ROOT/'research/rank4/smt_queries.zip') as archive:
    archive_names=archive.namelist()
    for r in data['records']:
        assert set(schema['required'])<=set(r)
        assert json.loads((HERE/'records'/f'{r["id"]}.json').read_text(encoding='utf-8'))==r
        for source in r['provenance']['sources']+r['provenance']['enrichment_sources']:
            p=ROOT/source['path']
            if p not in hashes: hashes[p]=hashlib.sha256(p.read_bytes()).hexdigest()
            assert hashes[p]==source['sha256']
        assert (ROOT/r['exponents']['plot_path']).is_file()
        assert r['verification']['jacobian_spectrum']=='independent interval-certified match'
        assert sum(x['multiplicity'] for x in r['exponents']['multiplicities'])==r['slice']['vertices']
        for key in ('manuscript','pdf','query_path'):
            if r['provenance'][key]:assert (ROOT/r['provenance'][key]).is_file(),r['provenance'][key]
        if r['rank']<3:continue
        if r['rank']==4:
            member=r['provenance']['query_member']
            names=[x for x in archive_names if x==member or x.endswith('/'+member)]
            assert len(names)==1,(member,names)
            query=archive.read(names[0])
        else:
            query=(ROOT/r['provenance']['query_path']).read_bytes()
        assert hashlib.sha256(query).hexdigest()==r['family']['coverage']['sha256']
print('PASS: 61 standalone records, embedded catalogue, source links and SHA-256 provenance, including every coverage query.')
