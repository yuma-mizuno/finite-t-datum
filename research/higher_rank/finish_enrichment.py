"""Verify, render and assemble higher-rank catalogue records."""
import contextlib
import io
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'catalogue'))
import verify_enrichment
import plot_exponents
sys.path.insert(0,str(ROOT/'docs/catalogue'))
from build_catalogue import source

def main(rank):
    directory=HERE/f'rank{rank}';records=json.loads((directory/'base-records.json').read_text())
    verify_enrichment.HERE=directory;verify_enrichment.records=lambda:records
    with contextlib.redirect_stdout(io.StringIO()):verify_enrichment.main()
    print(f'{len(records)} matrix ratios, quiver witnesses and independent Jacobian spectra verified.',flush=True)
    plot_exponents.HERE=directory
    with contextlib.redirect_stdout(io.StringIO()):plot_exponents.main()
    print(f'{len(records)} standalone SVG exponent plots rendered.',flush=True)
    spectra=json.loads((directory/'spectral-data.json').read_text());notes=json.loads((directory/'family-notes.json').read_text())
    quivers=json.loads((directory/'quiver-data.json').read_text());verified=json.loads((directory/'enrichment-verification.json').read_text())
    for r in records:
        r.update(spectra[r['id']]);r['notes']={'family':notes[r['id']],'quiver':quivers[r['id']]};r['verification']=verified[r['id']]
        assert r['notes']['quiver']['status'].startswith('certified-')
        r['exponents']['plot_path']=f'research/higher_rank/rank{rank}/plots/{r["id"]}.svg'
        r['provenance']['enrichment_sources']=[source(f'research/higher_rank/rank{rank}/'+p,role) for p,role in
             [('spectral-data.json','certified multiplicities and exact matrix ratios'),('family-notes.json','exact named-family identifications'),
              ('quiver-data.json','exact quiver mutation-class witnesses'),('enrichment-verification.json','independent quiver and Jacobian replay')]]
    (directory/'catalogue-records.json').write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')
    print(f'{len(records)} enriched rank-{rank} records ready.',flush=True)
if __name__=='__main__':main(int(sys.argv[1]))
