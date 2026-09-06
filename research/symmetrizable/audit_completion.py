"""Audit closure of every weighted rank and agreement with the preserved identity catalogue."""
import gzip,hashlib,itertools,json,sys,zipfile
from functools import lru_cache
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"higher_rank"))
from query_archives import archive_path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def readlines(path):return {r['id']:r for r in map(json.loads,path.read_text().splitlines())}
@lru_cache(None)
def file_hash(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def constant_key(p,m):
    n=len(p)
    return min(tuple(a[i][j] for a in pair for i in order for j in order) for order in itertools.permutations(range(n)) for pair in ((p,m),(m,p)))
def main():
    catalogue=json.loads((ROOT/'docs/catalogue/catalogue.json').read_text());out=[]
    for n in range(2,7):
        directory=HERE/f'rank{n}';report=json.loads((directory/'classification.json').read_text())
        assert report['complete'] and report['enrichment_complete']
        constants=json.loads((directory/'constant_candidates.json').read_text());cs={c['id']:c for c in constants['candidates']}
        assert constants['enumeration_complete'] and set(cs)==set(range(1,constants['count']+1))
        w,f,v=(readlines(directory/name) for name in ('lift_feasibility.jsonl','families.jsonl','verification.jsonl'))
        assert set(w)==set(v)==set(cs) and all(r['status'] in ('sat','unsat') for r in w.values())
        assert set(f)=={cid for cid,r in w.items() if r['status']=='sat'}
        assert all(r['result']=='unsat' for r in v.values())
        for r in f.values():
            assert r['coverage_status']=='unsat' and len(r['spaces'])==1
            s=r['spaces'][0];c=s['certificate']
            assert s['dimension']==n and s['scaling_and_shifts'] and c['positive'] and c['negative'] and c['labelled_tropical_seed_period']
        old=[r for r in catalogue['records'] if r['rank']==n and r['scope']['symmetrizer']=='identity']
        expected={constant_key(*[[[2*(i==j)-r['datum'][s][i][j] for j in range(n)] for i in range(n)] for s in ('A_plus_1','A_minus_1')]) for r in old}
        identity={constant_key(cs[cid]['N_plus_1'],cs[cid]['N_minus_1']) for cid in f if cs[cid]['symmetrizer']==[1]*n}
        assert identity==expected and len(identity)==len(old)
        new=[r for r in catalogue['records'] if r['rank']==n and r['scope']['symmetrizer']=='positive_diagonal']
        assert {r['constant_id'] for r in new}=={cid for cid in f if cs[cid]['symmetrizer']!=[1]*n}
        assert len(new)==report['new_nonidentity_families'] and len(f)==report['polynomial_families']
        from contextlib import ExitStack
        with ExitStack() as stack:
            archives={}
            for row in v.values():
                path=archive_path(directory,row['query_file'])
                if path not in archives: archives[path]=stack.enter_context(zipfile.ZipFile(path))
                query=gzip.decompress(archives[path].read(row['query_file']))
                assert hashlib.sha256(query).hexdigest()==row['sha256']
        for r in new:
            assert sum(x['multiplicity'] for x in r['exponents']['multiplicities'])==r['slice']['vertices']
            assert r['notes']['quiver']['certificate_replayed']
            for group in ('sources','enrichment_sources'):
                for item in r['provenance'][group]:assert file_hash(item['path'])==item['sha256']
        out.append({'rank':n,'constant_triples':len(cs),'complete_families':len(f),'identity_families_exactly_recovered':len(identity),'new_primitive_nonidentity_families':len(new),'unresolved_lifts_or_replays':0,'archived_query_hashes_checked':len(v)})
        print(n,'completion audit passed',out[-1],flush=True)
    (HERE/'completion-audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':main()
