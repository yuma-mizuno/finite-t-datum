"""Replay valued witnesses and assemble verified symmetrizable reader records."""
import contextlib,gzip,hashlib,io,json,subprocess,sys,zipfile
from pathlib import Path
from sage.all import matrix,QQ,identity_matrix,diagonal_matrix
from quiver_types import verify
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'catalogue'));import plot_exponents
sys.path.insert(0,str(ROOT/'docs/catalogue'));from build_catalogue import source
def main(rank):
    directory=HERE/f'rank{rank}';records=json.loads((directory/'base-records.json').read_text());report=json.loads((directory/'classification.json').read_text());assert report['complete']
    spectra,notes,quivers,verified=[json.loads((directory/p).read_text()) for p in ('spectral-data.json','family-notes.json','quiver-data.json','spectral-verification.json')]
    assert all(set(d)=={r['id'] for r in records} for d in (spectra,notes,quivers,verified))
    from folded_families import verify_relation
    from quiver_context import enrich as enrich_quiver
    parents={r['id']:r for r in json.loads((ROOT/'docs/catalogue/catalogue.json').read_text())['records'] if r['scope']['symmetrizer']=='identity'}
    for r in records:
        s=spectra[r['id']];a,b=(matrix(QQ,r['datum'][key]) for key in ('A_plus_1','A_minus_1'))
        k,ki=(matrix(QQ,s['matrix_ratios'][key]) for key in ('A_plus_inverse_A_minus','A_minus_inverse_A_plus'))
        assert a*k==b and b*ki==a and k*ki==identity_matrix(QQ,rank)
        D=diagonal_matrix(QQ,r['datum']['symmetrizer']);assert matrix(QQ,s['exponents']['fixed_point_matrix'])==D.inverse()*k*D
        assert verify(r,quivers[r['id']]);verified[r['id']]['mutation_class_path']='exact replay'
        enrich_quiver(r,quivers[r['id']])
        folds=[note for note in notes[r['id']]['identifications'] if note.get('category')=='Fold']
        assert all(verify_relation(r,note,parents) for note in folds)
        verified[r['id']]['exact_equitable_fold_replays']=len(folds)
    (directory/'quiver-data.json').write_text(json.dumps(quivers,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')
    (directory/'enrichment-verification.json').write_text(json.dumps(verified,indent=2)+'\n',encoding='utf8',newline='\n')
    plot_exponents.HERE=directory
    with contextlib.redirect_stdout(io.StringIO()):plot_exponents.main()
    for folder in ('constant_tasks','constant_parts','smt_queries','slice_tasks'):
        path=directory/folder
        if path.exists():
            with zipfile.ZipFile(directory/(folder+'.zip'),'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
                for file in sorted(path.iterdir()):
                    if file.is_file():archive.write(file,file.name)
    commit=subprocess.check_output(['git','log','-1','--format=%H','--','research/symmetrizable'],cwd=ROOT,text=True).strip()
    for r in records:
        r.update(spectra[r['id']]);r['notes']={'family':notes[r['id']],'quiver':quivers[r['id']]};r['verification']=verified[r['id']]
        r['exponents']['plot_path']=f'research/symmetrizable/rank{rank}/plots/{r["id"]}.svg'
        r['provenance']['source_commit']=commit
        r['provenance']['sources']=[source(s['path'],s['role']) for s in r['provenance']['sources']]+[source(f'research/symmetrizable/rank{rank}/slice_signatures.json','complete weighted slice distinctness')]
        r['provenance']['enrichment_sources']=[source(f'research/symmetrizable/rank{rank}/'+p,role) for p,role in
            [('spectral-data.json','certified thesis spectra and rational matrix ratios'),('family-notes.json','weighted named-family matches'),
             ('quiver-data.json','valued mutation-class witnesses'),('enrichment-verification.json','independent Jacobian spectrum and mutation replay')]]
    (directory/'catalogue-records.json').write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')
    report['enrichment_complete']=True;(directory/'classification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print(rank,len(records),'new weighted records fully enriched and archived',flush=True)
if __name__=='__main__':
    for rank in map(int,sys.argv[1:]):main(rank)
