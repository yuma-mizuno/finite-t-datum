"""Refresh the distributed source archive, retaining historical enumerators."""
from pathlib import Path
import hashlib,json,zipfile
P=Path(__file__).resolve().parents[1]
old=P/'research/symmetrizable/certificate-sources.zip'
with zipfile.ZipFile(old) as z:
    old_manifest=json.loads(z.read('manifest.json'))
    historical={name:z.read(name) for name in z.namelist() if name.startswith('historical/')}
sources=[]
for folder in ('research','docs/catalogue','tools'):
    for p in (P/folder).rglob('*'):
        if p.is_file() and '.qa' not in p.parts and p.suffix in ('.py','.cpp','.inc','.js','.css') and not any(x in p.parts for x in ('__pycache__','_deps','node_modules')):sources.append(p.relative_to(P).as_posix())
sources += ['docs/catalogue/index.template.html','docs/catalogue/proofs.json','research/symmetrizable/methods.template.html','research/symmetrizable/proofs.json']
sources += [p.relative_to(P).as_posix() for p in (P/'research/symmetrizable').glob('principal-constants*.txt*')]
manifest={'source_commit':old_manifest['source_commit'],'scope':'Distributed verification sources and original historical constant-enumerator sources. The imported revision identifies the research origin; current source bytes are identified by their hashes. External SageMath, SymPy and Z3 runtimes are separate dependencies.','files':[{'path':name,'sha256':hashlib.sha256((P/name).read_bytes()).hexdigest()} for name in sorted(set(sources))],'historical_enumerators':old_manifest['historical_enumerators'],'rank_reports':old_manifest['rank_reports']}
with zipfile.ZipFile(old,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for item in manifest['files']:z.write(P/item['path'],item['path'])
    for name,raw in historical.items():z.writestr(name,raw)
    z.writestr('manifest.json',json.dumps(manifest,indent=2)+'\n')
manifest['archive_sha256']=hashlib.sha256(old.read_bytes()).hexdigest()
(P/'research/symmetrizable/source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')
print('Archived',len(sources),'distributed verification sources')
