"""Compare permutation-certified orbit closure with the previously complete labelled orbits."""
import gzip,json
from pathlib import Path
from unlabelled_quiver_orbits import canonical,close_orbit,verify_orbit
HERE=Path(__file__).resolve().parent
def main():
    directory=HERE/'control-quiver-orbits';directory.mkdir(exist_ok=True);rows=[]
    records=json.loads((HERE/'rank4/base-records.json').read_text())
    data=json.loads((HERE/'rank4/quiver-data.json').read_text())
    for r in records:
        q=data[r['id']]
        if q['status']!='certified-finite-orbit':continue
        d=r['slice']['symmetrizer'];archive=HERE.parents[1]/q['certificate']['orbit_archive']
        old=json.loads(gzip.decompress(archive.read_bytes()))
        expected={canonical(tuple(map(tuple,b)),d)[0] for b in old['states']}
        target=directory/(r['id']+'.json.gz');new=close_orbit(r['slice']['B'],d,target)
        assert verify_orbit(r,new)
        actual={tuple(map(tuple,b)) for b in json.loads(gzip.decompress(target.read_bytes()))['states']}
        assert actual==expected
        rows.append({'record':r['id'],'labelled_matrices':len(old['states']),'matrices_up_to_relabelling':len(actual),'equal_after_independent_canonicalization':True,'every_integer_transition_replayed':True})
    path=HERE/'quiver-orbit-controls.json';path.write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n');print(rows)
if __name__=='__main__':main()
