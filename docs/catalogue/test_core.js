'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const C=require('./core.js'),data=JSON.parse(fs.readFileSync(path.join(__dirname,'catalogue.json'),'utf8'));
const frozen=JSON.stringify(data);
const strMatrix=a=>a.map(row=>row.map(terms=>terms.map(([c,p])=>[c,String(p)])));
const isPermutation=p=>assert.deepEqual([...p].sort((a,b)=>a-b),p.map((_,i)=>i));
assert.equal(data.records.length,116);
assert.equal(new Set(data.records.map(r=>r.id)).size,116);
assert.deepEqual([1,2,3,4].map(rank=>data.records.filter(r=>r.rank===rank).length),[2,6,16,37]);
assert.equal(data.records.filter(r=>r.rank===5).length,55);
for(const r of data.records){
  const d=r.datum;
  assert.equal(d.delays.length,r.rank);
  for(const sign of ['plus','minus']){
    assert.equal(d['N_'+sign].length,r.rank);
    assert.deepEqual(C.atOne(d,sign),d['A_'+sign+'_1'],r.id+' specialization');
    d['N_'+sign].forEach((row,i)=>{
      assert.equal(row.length,r.rank);
      row.forEach(entry=>{let last=0;for(const [c,p] of entry){assert(Number.isSafeInteger(c)&&c>0);assert(Number.isSafeInteger(p)&&last<p&&p<d.delays[i]);last=p;}});
    });
  }
  const base=C.transform(r,'1',Array(r.rank-1).fill('0'));
  assert(base.valid,r.id);assert.deepEqual(base.delays.map(String),d.delays.map(String));
  for(const sign of ['plus','minus'])assert.deepEqual(strMatrix(base['N_'+sign]),strMatrix(d['N_'+sign]));
  const doubled=C.transform(r,'2',Array(r.rank-1).fill('0'));
  assert(doubled.valid);assert.deepEqual(doubled.delays.map(String),d.delays.map(x=>String(2*x)));
  for(const sign of ['plus','minus'])assert.deepEqual(C.atOne(doubled,sign),d['A_'+sign+'_1']);
  const frames=C.sliceFrames(r.slice);
  assert.equal(frames.length,r.slice.mutation_word.length+2);
  assert.deepEqual(frames.at(-1),r.slice.B,r.id+' slice loop');
  assert.equal(r.slice.B.length,r.slice.vertices);
  for(const b of frames){
    b.forEach((row,i)=>{assert.equal(row.length,b.length);row.forEach((x,j)=>{assert(Number.isSafeInteger(x));assert(x===-b[j][i]);});});
    for(let k=0;k<b.length;k++)assert.deepEqual(C.mutate(C.mutate(b,k),k),b);
  }
  for(const p of [r.slice.relabel_old_to_new,r.exchange.relabel_old_to_new,r.periodicity.positive_negative_permutation,r.periodicity.negative_negative_permutation])isPermutation(p);
  for(const v of r.slice.mutation_word)assert(Number.isInteger(v)&&v>=0&&v<r.slice.vertices);
  for(const v of r.exchange.mutation_vertices)assert(Number.isInteger(v)&&v>=0&&v<r.exchange.vertices.length);
  assert(r.family.rref.every(row=>row.length===r.family.variable_names.length));
  assert.equal(r.family.representative_values.length,r.family.variable_names.length);
  for(const row of r.family.rref)assert.equal(row.reduce((s,a,i)=>s.add(C.Rational.parse(a).mul(r.family.representative_values[i])),C.Rational.parse(0)).toString(),'0');
  if(r.rank>=3)assert.equal(r.family.coverage.result,'unsat');
}
const r=data.records.find(r=>r.id==='r4-c01');
const shifted=C.transform(r,'3/2',['1/2','1/2','1/2']);
assert(shifted.valid);assert.deepEqual(shifted.delays.map(String),['3','3','3','3']);
assert.deepEqual(strMatrix(shifted.N_minus)[0][3],[[1,'2']]);
assert.deepEqual(strMatrix(shifted.N_minus)[3][0],[[1,'1']]);
assert.equal(C.transform(r,'1/2',['0','0','0']).valid,false);
assert.equal(C.transform(r,'1',['2','0','0']).valid,false);
for(const value of ['0','-1','1/0','1e2','abc','', '1'.repeat(61)])assert.equal(C.transform(r,value,['0','0','0']).valid,false);
assert.equal(C.Rational.parse('-0.125').toString(),'-1/8');
assert.equal(C.Rational.parse('9007199254740993').add(1).toString(),'9007199254740994');
assert(C.transform(r,'9007199254740993',['0','0','0']).valid);
assert.equal(C.latexPolynomial([[1,2],[-1,1],[1,0]]),'z^{2} - z + 1');
assert.deepEqual(data.records.filter(r=>r.rank===4&&r.datum.delays.every(x=>x===2)).map(r=>r.class_number),[1,2,3,28,29,30,31,34,35]);
assert.equal(JSON.stringify(data),frozen,'The reader must not mutate the source records.');
console.log(`PASS: ${data.records.length} exact records; specializations, admissible lifts, RREF witnesses, slice loops, permutations and arithmetic edge cases.`);
