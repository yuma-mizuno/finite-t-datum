"""Thesis determinant spectrum using the dual Nahm equation for nonidentity D."""
import sys,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'catalogue'))
from exponents import *
from sage.all import diagonal_matrix
HERE=Path(__file__).resolve().parent
def compute(r):
 n=r['rank'];ap=matrix(QQ,r['datum']['A_plus_1']);am=matrix(QQ,r['datum']['A_minus_1']);K=ap.inverse()*am
 D=diagonal_matrix(QQ,r["datum"]["symmetrizer"]);KD=K*D
 assert KD==KD.transpose() and all(KD[:i,:i].det()>0 for i in range(1,n+1))
 dual=D.inverse()*K*D
 shift=centered_gauge(r);f,box=solution_box(dual)
 bar=interval_determinant(polynomial_matrix(r,f));t=r['slice']['components'];omega=r['periodicity']['labelled_period']
 assert omega%t==0;N=omega//t;degree=sum(r['datum']['delays'])//t
 assert degree==r['slice']['vertices'] and bar.degree()==degree*t
 assert all(bar[i].contains_zero() for i in range(bar.degree()+1) if i%t)
 coeff=[bar[t*i] for i in range(degree+1)];assert coeff[-1]==1 and coeff[0]==1
 powers=[RI(degree)]
 for k in range(1,N):
  s=sum(coeff[degree-j]*powers[k-j] for j in range(1,min(k,degree+1)))
  s+=k*coeff[degree-k] if k<=degree else 0
  powers.append(-s)
 root=CI.zeta(N);multiplicities=[];bounds=[]
 for m in range(N):
  value=sum(CI(powers[k])*root**(-m*k) for k in range(N))/N
  lo=value.real().lower().ceil();hi=value.real().upper().floor()
  assert lo==hi and lo>=0 and value.imag().contains_zero(),(r['id'],m,value)
  multiplicities.append(int(lo));bounds.append({'m':m,'multiplicity':int(lo),'real_lower':str(value.real().lower()),'real_upper':str(value.real().upper()),
                                             'imag_lower':str(value.imag().lower()),'imag_upper':str(value.imag().upper())})
 assert sum(multiplicities)==degree
 CF=CyclotomicField(N);w=CF.gen();PP=PolynomialRing(CF,'x');x=PP.gen();tau=PP(1)
 for m,mult in enumerate(multiplicities):tau*=(x-w**m)**mult
 # The exact product and independently enclosed determinant coefficients agree.
 # Evaluate the canonical power-basis expression at the certified enclosure of
 # exp(2*pi*i/N). Interval equality cannot serve Sage's exact homomorphism check.
 embedding=lambda v: sum(CI(QQ(a))*root**j for j,a in enumerate(v.list()))
 assert all(coeff[i].overlaps(embedding(tau[i]).real()) and embedding(tau[i]).imag().contains_zero() for i in range(degree+1))
 return {'matrix_ratios':{'A_plus_inverse_A_minus':[[str(v) for v in row] for row in K.rows()],
                         'A_minus_inverse_A_plus':[[str(v) for v in row] for row in K.inverse().rows()]},
  'exponents':{'status':'certified using interval arithmetic and recorded periodicity', 'source':'https://yuma-mizuno.github.io/thesis.pdf#page=56',
   'formula':'bar_tau(z) = det(A_plus(z) diag(1-f) + A_minus(z) diag(f)); tau(z) = bar_tau(z^(1/t))',
   'fixed_point_equation':'log(f) = (D^(-1) A_plus(1)^(-1) A_minus(1) D) log(1-f)',
   'fixed_point_matrix':[[str(v) for v in row] for row in dual.rows()],
   'cartan_like_centering_shifts':[str(s) for s in shift], 'centering_integer_dilation':int(lcm(s.denominator() for s in shift)),
   'components':t,'representative_period':omega,'root_order':N,'degree':degree,
   'values':[m for m,mult in enumerate(multiplicities) for _ in range(mult)],
   'multiplicities':[{'m':m,'multiplicity':mult} for m,mult in enumerate(multiplicities) if mult],
   'fixed_point_decimal':[str(v.center()) for v in f],
   'tau_coefficients_decimal':[str(v.center()) for v in coeff],
   'tau_exact_cyclotomic_coefficients':[str(v) for v in tau.list()],
   'cyclotomic_generator':f'zeta{N} = exp(2*pi*i/{N})',
   'certificate':{'fixed_point_box':box,'multiplicity_bounds':bounds,'method':'Interval Newton sums and discrete Fourier inversion on the exact period grid',
                  'dependency':'Thesis Theorem 3.4.5 and the record’s universal periodicity certificate guarantee integral multiplicities.'}}}

def main(rank):
 directory=HERE/f'rank{rank}';rs=json.loads((directory/'base-records.json').read_text());results={};start=time.monotonic()
 for r in rs:
  results[r['id']]=compute(r);print(r['id'],'weighted thesis spectrum certified',flush=True)
 (directory/'spectral-data.json').write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')
 with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:out.write(json.dumps({'rank':rank,'stage':'weighted-thesis-spectra','wall_seconds':time.monotonic()-start,'records':len(results)})+'\n')
if __name__=='__main__':main(int(sys.argv[1]))
