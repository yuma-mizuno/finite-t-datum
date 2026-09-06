/* Exact arithmetic and mutation functions; shared by the reader and tests. */
(function (root) {
  'use strict';
  const gcd = (a,b) => { a=a<0n?-a:a; b=b<0n?-b:b; while(b){const t=a%b;a=b;b=t;} return a; };
  class Rational {
    constructor(n,d=1n) { n=BigInt(n);d=BigInt(d);if(!d)throw Error('A denominator cannot be zero.');if(d<0n){n=-n;d=-d;}const g=gcd(n,d);this.n=n/g;this.d=d/g; }
    static parse(x) {
      if(x instanceof Rational)return x;
      const s=String(x).trim();if(s.length>60)throw Error('Use at most 60 characters per number.');
      if(/^[+-]?\d+(\/\d+)?$/.test(s)){const [n,d='1']=s.split('/');return new Rational(n,d);}
      if(/^[+-]?\d+\.\d+$/.test(s)){const [a,b]=s.split('.');return new Rational(BigInt(a+b),10n**BigInt(b.length));}
      throw Error('Enter an integer or a fraction, such as 3/2.');
    }
    add(x){x=Rational.parse(x);return new Rational(this.n*x.d+x.n*this.d,this.d*x.d);}
    sub(x){x=Rational.parse(x);return new Rational(this.n*x.d-x.n*this.d,this.d*x.d);}
    mul(x){x=Rational.parse(x);return new Rational(this.n*x.n,this.d*x.d);}
    cmp(x){x=Rational.parse(x);const d=this.n*x.d-x.n*this.d;return d<0n?-1:d>0n?1:0;}
    toString(){return this.d===1n?String(this.n):`${this.n}/${this.d}`;}
    toJSON(){return this.toString();}
  }
  function transform(record,scale,shifts) {
    const errors=[];let lambda,s;
    try{lambda=Rational.parse(scale);s=shifts.map(Rational.parse);if(s.length!==record.rank-1)throw Error('Wrong number of species shifts.');s.push(new Rational(0));}
    catch(e){return {valid:false,errors:[e.message]};}
    if(lambda.cmp(0)<=0)errors.push('The time scale λ must be positive.');
    const delays=record.datum.delays.map(r=>lambda.mul(r));
    delays.forEach((r,i)=>{if(r.d!==1n||r.cmp(0)<=0)errors.push(`Delay r${i+1} = ${r} must be a positive integer.`);});
    const convert=(a,sign)=>a.map((row,i)=>row.map((entry,j)=>entry.map(([c,p])=>{
      const t=lambda.mul(p).add(s[i]).sub(s[j]);
      if(t.d!==1n)errors.push(`N${sign}[${i+1},${j+1}] has nonintegral exponent ${t}.`);
      else if(t.cmp(0)<=0||t.cmp(delays[i])>=0)errors.push(`N${sign}[${i+1},${j+1}]: exponent ${t} must lie strictly between 0 and ${delays[i]}.`);
      return [c,t];
    })));
    const plus=convert(record.datum.N_plus,'+'),minus=convert(record.datum.N_minus,'−');
    return {valid:errors.length===0,errors:[...new Set(errors)],lambda,shifts:s,delays,symmetrizer:record.datum.symmetrizer||Array(record.rank).fill(1),N_plus:plus,N_minus:minus};
  }
  function normalize(terms) {
    const map=new Map();for(const [c,p] of terms){const k=Rational.parse(p).toString();map.set(k,(map.get(k)||0)+c);}
    return [...map].filter(([,c])=>c).map(([p,c])=>[c,Rational.parse(p)]).sort((a,b)=>a[1].cmp(b[1]));
  }
  function polynomialMatrix(datum,sign,kind='A') {
    const a=datum[sign==='plus'?'N_plus':'N_minus'];
    if(kind==='N')return a.map(row=>row.map(normalize));
    return a.map((row,i)=>row.map((entry,j)=>normalize([
      ...(i===j?[[1,0],[1,datum.delays[i]]]:[]),...entry.map(([c,p])=>[-c,p])])));
  }
  function atOne(datum,sign){return polynomialMatrix(datum,sign).map(row=>row.map(terms=>terms.reduce((s,[c])=>s+c,0)));}
  function mutate(b,k){return b.map((row,i)=>row.map((x,j)=>i===k||j===k?(x===0?0:-x):x+
    Math.max(b[i][k],0)*Math.max(b[k][j],0)-Math.max(-b[i][k],0)*Math.max(-b[k][j],0)));}
  function relabel(b,p){const inv=p.map((_,i)=>p.indexOf(i));return inv.map(i=>inv.map(j=>b[i][j]));}
  function sliceFrames(slice){let b=slice.B;const frames=[b];for(const k of slice.mutation_word){b=mutate(b,k);frames.push(b);}frames.push(relabel(b,slice.relabel_old_to_new));return frames;}
  function polyText(terms) {
    terms=normalize(terms);if(!terms.length)return '0';
    return [...terms].reverse().map(([c,p],i)=>{
      const a=Math.abs(c),power=p.cmp(0)===0?'':p.cmp(1)===0?'z':`z^(${p})`;
      const coefficient=power&&a===1?'':String(a);
      return (c<0?(i?' − ':'−'):(i?' + ':''))+coefficient+power;
    }).join('');
  }
  function latexPolynomial(terms){return polyText(terms).replace(/−/g,'-').replace(/z\^\(([^)]+)\)/g,(_,p)=>'z^{'+(p.includes('/')?'\\frac{'+p.split('/')[0]+'}{'+p.split('/')[1]+'}':p)+'}');}
  function latexMatrix(a){return '\\begin{pmatrix}'+a.map(row=>row.map(latexPolynomial).join(' & ')).join('\\\\')+'\\end{pmatrix}';}
  const api={Rational,transform,normalize,polynomialMatrix,atOne,mutate,relabel,sliceFrames,polyText,latexPolynomial,latexMatrix};
  if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TDataCore=api;
})(typeof globalThis!=='undefined'?globalThis:this);
