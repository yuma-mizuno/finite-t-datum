"""Compare exact C++ minima and orbit sizes with independently saved results."""
import json
from pathlib import Path
from slice_cpp import orbit

HERE=Path(__file__).resolve().parent


def main():
    items=[s for s in json.loads((HERE/'rank5/slice_signatures.json').read_text()) if s['canonical_signature']['orbit_minimum'] is not None]
    items += [json.loads(p.read_text())['slice'] for p in sorted((HERE/'rank6/slice_tasks').glob('*.json'))]
    reports=[]
    for i,s in enumerate(items):
        result=orbit(s['B'],s['mutation_word'],s['relabel_old_to_new'],HERE/'slice_cpp_controls',str(i))
        assert result['minimum']==s['canonical_signature']['orbit_minimum'],s['id']
        assert result['states']==s['equivalent_canonical_states_including_sign_exchange'],s['id']
        reports.append({'case':s['id'],'vertices':len(s['B']),'states':result['states'],'exact_minimum_and_count_match':True})
        print(reports[-1],flush=True)
    (HERE/'slice-cpp-controls.json').write_text(json.dumps(reports,indent=2)+'\n')


if __name__=='__main__':main()
