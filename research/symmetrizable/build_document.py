"""Build the human-readable proof and the reader's weighted proof summaries."""
import json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main():
    reports={}
    for path in HERE.glob('rank*/classification.json'):
        r=json.loads(path.read_text())
        if r['complete'] and r['enrichment_complete']:reports[r['rank']]=r
    common={'scope':'Positive diagonal symmetrizers, recorded primitively, and diagonal N₀. The rank is the number of species.',
            'theorem':'Every finite datum in the stated scope belongs to one recorded admissible scale-and-shift family, up to species permutations, sign exchange and the recorded slice equivalences. Scalar multiples of primitive D are recovered separately.',
            'definition':'<p>A± = diag(1+zʳⁱ)−N± have nonnegative integral N±, strict row support and disjoint opposite signs. D is positive diagonal; D⁻¹N±D is integral and A₊DA₋† = A₋DA₊†. Primitive normalization means gcd(d₁,…,dₙ)=1. General leading permutations and automatic Langlands-dual identification are outside this scope.</p>',
            'status':'Computer-assisted classification with finite coefficient and symmetrizer bounds, exact lift coverage, periodicity and weighted slice certificates.',
            'references':[{'title':'Mizuno: Difference equations arising from cluster algebras','url':'https://arxiv.org/abs/1912.05710','use':'General symmetrizable T-data, reconstruction, positivity and periodicity.'},
                          {'title':'Mizuno: periodic Y-systems of rank two','url':'https://arxiv.org/abs/2301.13239','use':'Published identity-symmetrizer rank-two subcatalogue.'},
                          {'title':'Mizuno thesis','url':'https://yuma-mizuno.github.io/thesis.pdf','use':'Dual Nahm fixed point and the determinant formula for exponents.'}]}
    proofs={'common':common};rows=['<tr><td>1</td><td>2</td><td>0</td><td>2</td><td>Complete</td></tr>']
    for n in range(2,7):
        if n not in reports:
            rows.append(f'<tr><td>{n}</td><td>{[0,2,6,16,37,55,108][n]}</td><td>Pending</td><td>Pending</td><td>Search running</td></tr>');continue
        r=reports[n];count=r['polynomial_families'];constants=r['constant_candidates'];new=r['new_nonidentity_families']
        rows.append(f'<tr><td><a href="../../docs/catalogue/index.html#s{n}-c01/matrices">{n}</a></td><td>{r["identity_families"]}</td><td>{new}</td><td>{count}</td><td>Complete</td></tr>')
        proofs[str(n)]={'count':count,'status':'computer-assisted-complete','statistics':[[constants,'retained necessary constant triples'],[constants-count,'no polynomial lift'],[count,'complete periodic families'],[new,'with nonidentity primitive D']],
          'decomposable':'Take direct sums of indecomposable blocks with independent time scales, sign choices and positive integer scalar factors on block symmetrizers.',
          'steps':[
            {'title':'1. Bound both coefficients and the symmetrizer','kind':'Necessary positivity and integrality','html':f'<p>Original and dual off-diagonal coefficient sums are at most {2**n-1}. A directed spanning tree bounds every primitive weight by {(2**n-1)**(n-1)}. There are {math.factorial(n+1)} upper-triangle possibilities per sign. The exact enumeration leaves {constants} necessary triples, retaining D in every equivalence comparison.</p>'},
            {'title':'2. Exhaust all polynomial lifts','kind':'Unbounded linear arithmetic','html':f'<p>Entry (i,j) uses atoms of coefficient dᵢ/gcd(dᵢ,dⱼ). Strict support, dual integrality and weighted Laurent identities are encoded exactly. Every exclusion and outside-space query is resolved and replayed. Each of the {count} surviving spaces has dimension {n}, spanned by scaling and relative shifts.</p><p>Eligible lower-rank principal data supply necessary equations and exclusions; the original identity-only obstructions are used only for D = I.</p>'},
            {'title':'3. Verify periodicity and slice distinctness','kind':'Exact weighted mutation certificates','html':'<p>The reconstructed exchange matrix satisfies BD = −DBᵀ with the stored vertex weights. Every representative has a labelled tropical return and two-sided negative permutation C-matrices. Weighted terminal cycles distinguish singleton groups; complete cyclic-loop orbits resolve every collision group.</p>'},
            {'title':'4. Record names, valued quivers and exponents','kind':'Independent enrichment certificates','html':'<p>Family notes retain exact constructor transformations. Quiver notes contain mutation paths or complete finite orbits. The thesis fixed-point equation uses K∨ = D⁻¹A₊(1)⁻¹A₋(1)D, while the determinant uses the original A matrices. Interval-certified multiplicities agree with an independent mutation-Jacobian calculation.</p>'}]}
    table='<table><thead><tr><th>Rank</th><th>D = I</th><th>New D ≠ I</th><th>Total families</th><th>Status</th></tr></thead><tbody>'+''.join(rows)+'</tbody></table>'
    complete=max(reports,default=1)
    table='<p class="status">The symmetrizable classification is complete through rank '+str(complete)+('. Rank six remains under computation; no complete weighted rank-six count is asserted.' if complete<6 else ', including rank six.')+'</p>'+table
    (HERE/'methods.html').write_text((HERE/'methods.template.html').read_text(encoding='utf8').replace('<!--RESULTS-->',table),encoding='utf8',newline='\n')
    (HERE/'proofs.json').write_text(json.dumps(proofs,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')
    print('Built weighted proof document through rank',complete)
if __name__=='__main__':main()
