"""Replay every pulse exclusion with an independent exhaustive C++ algorithm."""
import hashlib
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent


def main():
    directory=HERE/'rank6';target=directory/'pulse_replays';target.mkdir(exist_ok=True)
    records={r['id']:r for r in json.loads((HERE.parents[1]/'docs/catalogue/catalogue.json').read_text())['records']}
    queries=[json.loads(x) for x in (directory/'verification.jsonl').read_text().splitlines()]
    results=[]
    for query in queries:
        c=query.get('analytic_exclusion',{})
        if c.get('lemma')!='finite pulse enumeration for a single attachment to a classified principal core':continue
        datum=records[c['core_record']]['datum'];n=len(datum['delays']);h=max(datum['delays']);a=c['core_attachment_species']
        matrices=[]
        for name in ('N_plus','N_minus'):
            terms=[[[0]*n for _ in range(n)] for _ in range(h+1)]
            for i,r in enumerate(datum['delays']):terms[0][i][i]=1;terms[r][i][i]=1
            for i,row in enumerate(datum[name]):
                for j,entry in enumerate(row):
                    for coefficient,degree in entry:terms[degree][i][j]-=coefficient
            matrices.append(terms)
        if c['core_sign_exchange']:matrices.reverse()
        values=[n,h,a,c['diagonal_coefficient_sum'],datum['delays'][a]]
        values.extend(x for row in c['target_coefficient_sums'] for x in row)
        values.extend(matrix[d][i][j] for d in range(1,h+1) for i in range(n) for j in range(n) for matrix in matrices)
        source=target/f'{query["id"]}.txt';output=target/f'{query["id"]}.json';source.write_text(' '.join(map(str,values))+'\n')
        subprocess.run([str(HERE/'bin/verify_pulses'),str(source),str(output)],check=True)
        result=json.loads(output.read_text());assert result['solutions']==0 and result['delay_bound']==c['proved_delay_bound'],query['id']
        result.update({'constant_id':query['id'],'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
        results.append(result);print(query['id'],result['pulse_positions_checked'],'positions; independently excluded',flush=True)
    report={'source_sha256':hashlib.sha256((HERE/'verify_pulses.cpp').read_bytes()).hexdigest(),
            'binary_sha256':hashlib.sha256((HERE/'bin/verify_pulses').read_bytes()).hexdigest(),'replays':results}
    (directory/'pulse-verification.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
