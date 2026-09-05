"""All labelled lower-rank constants for exact hereditary branch pruning."""
import hashlib
import itertools as it
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    source=HERE.parents[1]/'docs/catalogue/catalogue.json'
    records=json.loads(source.read_text(encoding='utf8'))['records']
    lines=[];counts={}
    for n in range(1,6):
        labelled=set()
        selected=[r for r in records if r['rank']==n]
        assert len(selected)==[2,6,16,37,55][n-1]
        for r in selected:
            matrices=[[[2*(i==j)-r['datum'][a][i][j] for j in range(n)] for i in range(n)] for a in ('A_plus_1','A_minus_1')]
            for perm in it.permutations(range(n)):
                for pair in (matrices,matrices[::-1]):
                    labelled.add(tuple(b[i][j] for b in pair for i in perm for j in perm))
        counts[n]=len(labelled)
        lines.extend(' '.join(map(str,(n,*k))) for k in sorted(labelled))
    target=HERE/'principal-constants.txt'
    target.write_text('\n'.join(lines)+'\n',encoding='utf8',newline='\n')
    report={'source':'docs/catalogue/catalogue.json','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'rank_counts':counts,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
            'dependency':'Complete simultaneous-positive polynomial-lift classifications in ranks 1 through 5.'}
    (HERE/'principal-constants-provenance.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print(report)


if __name__=='__main__':main()
