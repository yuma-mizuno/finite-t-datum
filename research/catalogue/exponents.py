"""Mizuno thesis (3.4.1), Definition 3.4.4: certified exponent multiplicities.

Run under SageMath. A Krawczyk box encloses the unique positive Nahm solution.
Interval polynomial arithmetic + Newton sums + the discrete Fourier inversion
enclose each (a priori integral) root-of-unity multiplicity in one integer.
The root-of-unity premise is the recorded periodicity and thesis Theorem 3.4.5.
"""
from pathlib import Path
import itertools as it
import json
import mpmath as mp
import sys
from fractions import Fraction
from sage.all import (QQ, ZZ, matrix, vector, identity_matrix, RealField,
                      RealIntervalField, ComplexIntervalField, PolynomialRing,
                      CyclotomicField, lcm)

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
BITS=512;RF=RealField(BITS);RI=RealIntervalField(BITS);CI=ComplexIntervalField(BITS)
POLY=PolynomialRing(RI,'z');z=POLY.gen()

def records():
 return json.loads((HERE/'lower-rank-records.json').read_text())+[r for r in json.loads((ROOT/'docs/catalogue/catalogue.json').read_text(encoding='utf-8'))['records'] if r['rank']>=3]

def centered_gauge(record):
 n=record['rank'];r=record['datum']['delays'];rows=[];rhs=[]
 for key in ['N_plus','N_minus']:
  for i,row in enumerate(record['datum'][key]):
   for j,entry in enumerate(row):
    if entry:
     target=QQ(r[i]-entry[0][1]-entry[-1][1])/2
     rows.append([int(k==i)-int(k==j) for k in range(n)]);rhs.append(target)
 rows.append([int(i==n-1) for i in range(n)]);rhs.append(0)
 shift=list(matrix(QQ,rows).solve_right(vector(QQ,rhs)))
 for key in ['N_plus','N_minus']:
  for i,row in enumerate(record['datum'][key]):
   for j,entry in enumerate(row):
    transformed={QQ(p)+shift[i]-shift[j]:c for c,p in entry}
    assert all(0<p<r[i] and transformed.get(r[i]-p)==c for p,c in transformed.items())
 return shift

def solution_box(K):
 n=K.nrows();mp.mp.dps=185;k=mp.matrix([[mp.mpf(str(K[i,j].numerator()))/mp.mpf(str(K[i,j].denominator())) for j in range(n)] for i in range(n)])
 eye=mp.eye(n);x=mp.matrix([0]*n)
 def fun(a):return a+(k-eye)*mp.matrix([mp.log(1+mp.exp(t)) for t in a])
 for _ in range(120):
  f=mp.matrix([1/(1+mp.exp(-t)) for t in x]);g=fun(x)
  if max(abs(v) for v in g)<mp.mpf('1e-170'):break
  J=eye+(k-eye)*mp.diag(f);v=mp.lu_solve(J,g);rate=mp.mpf(1)
  while max(abs(v) for v in fun(x-rate*v))>=max(abs(v) for v in g):rate/=2
  x-=rate*v
 else:raise RuntimeError('Nahm iteration did not converge')
 centers=[mp.nstr(v,180) for v in f]
 fractions=[Fraction(v) for v in centers]
 exact=vector(QQ,[QQ(v.numerator)/v.denominator for v in fractions]);radius=QQ(1)/2**350
 X=vector(RI,[RI(v-radius,v+radius) for v in exact]);point=vector(RI,exact)
 assert all(0<v.lower() and v.upper()<1 for v in X)
 G=vector(RI,[v.log() for v in point])-matrix(RI,K)*vector(RI,[(1-v).log() for v in point])
 def jac(values):return matrix(RI,[[RI(int(i==j))/values[i]+RI(K[i,j])/(1-values[j]) for j in range(n)] for i in range(n)])
 J=jac(X);J0=matrix(RF,[[RF(jac(point)[i,j].center()) for j in range(n)] for i in range(n)])
 C=matrix(RI,matrix(QQ,J0.inverse()))
 box=point-C*G+(identity_matrix(RI,n)-C*J)*(X-point)
 assert all(X[i].lower()<box[i].lower() and box[i].upper()<X[i].upper() for i in range(n))
 return X,{'method':'Krawczyk strict interior inclusion','center_decimal':centers,'radius':'2^-350','arithmetic_bits':BITS,
           'maximum_image_radius':str(max((box[i]-point[i]).absolute_diameter()+abs((box[i]-point[i]).center()) for i in range(n)))}

def polynomial_matrix(r,f):
 d=r['datum'];n=r['rank']
 return [[(1+z**d['delays'][i] if i==j else POLY(0))
          -sum(RI(c)*(1-f[j])*z**p for c,p in d['N_plus'][i][j])
          -sum(RI(c)*f[j]*z**p for c,p in d['N_minus'][i][j]) for j in range(n)] for i in range(n)]

def interval_determinant(a):
 n=len(a);out=POLY(0)
 for perm in it.permutations(range(n)):
  sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n));term=POLY(sign)
  for i in range(n):term*=a[i][perm[i]]
  out+=term
 return out

def compute(r):
 n=r['rank'];ap=matrix(QQ,r['datum']['A_plus_1']);am=matrix(QQ,r['datum']['A_minus_1']);K=ap.inverse()*am
 assert K==K.transpose() and all(K[:i,:i].det()>0 for i in range(1,n+1))
 shift=centered_gauge(r);f,box=solution_box(K)
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
   'fixed_point_equation':'log(f) = (A_plus(1)^(-1) A_minus(1)) log(1-f)',
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

def main():
 results={}
 for r in records():
  results[r['id']]=compute(r)
  print(r['id'],results[r['id']]['exponents']['root_order'],results[r['id']]['exponents']['values'],flush=True)
  (HERE/'spectral-data.json').write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
 print('All 61 fixed-point boxes and exponent multiplicities certified.')

if __name__=='__main__':main()
