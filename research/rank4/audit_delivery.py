"""Check independent controls, the chain algebra, and document references."""
import itertools as it
import json
from pathlib import Path
import re

import sympy as sp

HERE=Path(__file__).resolve().parent


def key(p,m):
    return min(tuple(a[q[i]][q[j]] for a in pair for i in range(4) for j in range(4))
               for q in it.permutations(range(4)) for pair in ((p,m),(m,p)))


def main():
    classes=json.loads((HERE/'classification.json').read_text())
    degree_two=json.loads((HERE/'degree_two.json').read_text())
    lookup={key(*[[[2*(i==j)-c[field][i][j] for j in range(4)] for i in range(4)]
                  for field in ('A_plus_1','A_minus_1')]):c['class'] for c in classes}
    mapping=sorted(lookup[key(x['N_plus_1'],x['N_minus_1'])] for x in degree_two['classes'])
    assert mapping==[1,2,3,28,29,30,31,34,35]
    signatures=json.loads((HERE/'slice_signatures.json').read_text())
    assert len(signatures)==37 and len({tuple(x['canonical_signature']) for x in signatures})==37
    q,u,v,b,s,X,Y,Z,W=sp.symbols('q u v b s X Y Z W')
    F=1+q;H=b*(1+1/s)
    plus=sp.Matrix([[F,0,0,0],[0,F,0,0],[0,0,F,0]])
    minus=sp.Matrix([[F,-u,0,-b],[-q/u,F,-v,0],[0,-q/v,F,0]])
    apstar=sp.Matrix([-X,-Y/u,-Z/(u*v),1+1/s])
    amstar=sp.Matrix([0,0,-W/(u*v),1+1/s])
    raw=plus*amstar-minus*apstar
    expected=sp.Matrix([F*X-Y+H,(F*Y-q*X-Z)/u,(F*(Z-W)-q*Y)/(u*v)])
    assert all(sp.expand(x)==0 for x in raw-expected)
    assert sp.expand(F*(F*X+H)-q*X-((1+q+q*q)*X+F*H))==0
    tex=(HERE/'rank4-classification.tex').read_text()
    labels=re.findall(r'\\label\{([^}]+)\}',tex)
    refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',tex)
    assert len(labels)==len(set(labels)) and set(refs)<=set(labels)
    report={'degree_two_classes':mapping,'distinct_slice_classes':37,
            'chain_equation_derivation':'symbolically verified','labels':len(labels),
            'unresolved_labels':[],'constant_matrices':74}
    (HERE/'delivery_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(report)


if __name__=='__main__':main()
