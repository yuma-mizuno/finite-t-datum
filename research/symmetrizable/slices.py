"""Complete slice comparisons retaining primitive vertex symmetrizers."""
from collections import defaultdict,deque
import hashlib,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
def invariant(s):
    p=s['relabel_old_to_new'];d=s['symmetrizer'];seen=set();out=[]
    for root in range(len(p)):
        if root in seen:continue
        at=root;length=0
        while at not in seen:
            assert d[at]==d[root];seen.add(at);length+=1;at=p[at]
        out.append((length,d[root]))
    return tuple(sorted(out))
def cpp(s,directory,identifier):
    directory.mkdir(exist_ok=True);inp=directory/(identifier+'.txt');out=directory/(identifier+'.json')
    values=[len(s['mutation_word']),s['vertices']]+[x for row in s['B'] for x in row]+s['mutation_word']+s['relabel_old_to_new']+s['symmetrizer']
    inp.write_text(' '.join(map(str,values))+'\n',encoding='utf8',newline='\n')
    subprocess.run([str(HERE/'bin/slice_orbit'),str(inp),str(out)],check=True)
    result=json.loads(out.read_text());assert result['complete'] and result['states']==result['processed']
    result.update({'source_sha256':hashlib.sha256((HERE/'slice_orbit.cpp').read_bytes()).hexdigest(),
                   'binary_sha256':hashlib.sha256((HERE/'bin/slice_orbit').read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest()})
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');return result
def mutation(b,k):
    return tuple(tuple(-b[i][j] if i==k or j==k else b[i][j]+max(b[i][k],0)*max(b[k][j],0)-max(-b[i][k],0)*max(-b[k][j],0) for j in range(len(b))) for i in range(len(b)))
def python_orbit(s):
    word=tuple(s['mutation_word']);p=tuple(s['relabel_old_to_new']);d=tuple(s['symmetrizer']);b=tuple(map(tuple,s['B']))
    pending=deque([(b,word,p,d),(tuple(tuple(-x for x in row) for row in b),word,p,d)]);seen=set()
    while pending:
        b,word,p,d=pending.popleft();order=[];lengths=[]
        for root in word:
            at=root;length=0
            while True:
                order.append(at);length+=1;at=p[at]
                if at==root:break
            lengths.append(length)
        assert sorted(order)==list(range(len(b)))
        signature=tuple(lengths)+tuple(d[root] for root in word)+tuple(b[i][j] for i in order for j in order)
        if signature in seen:continue
        seen.add(signature);inverse={v:i for i,v in enumerate(order)}
        b=tuple(tuple(b[i][j] for j in order) for i in order);word=tuple(inverse[v] for v in word);p=tuple(inverse[p[v]] for v in order);d=tuple(d[v] for v in order)
        inv=tuple(p.index(i) for i in range(len(p)))
        pending.append((mutation(b,word[0]),word[1:]+(inv[word[0]],),p,d));prefix=b
        for i in range(len(word)-1):
            if not prefix[word[i]][word[i+1]]:
                other=list(word);other[i],other[i+1]=other[i+1],other[i];pending.append((b,tuple(other),p,d))
            prefix=mutation(prefix,word[i])
    return {'states':len(seen),'minimum':list(min(seen))}
def main(rank):
    start=time.monotonic();directory=HERE/f'rank{rank}';path=directory/'base-records.json';rs=json.loads(path.read_text());groups=defaultdict(list)
    for r in rs:groups[invariant(r['slice'])].append(r['id'])
    minima=set();results={};controls={}
    for index,r in enumerate(rs):
        s=r['slice'];key=invariant(s);result=None
        if len(groups[key])>1 or (rank<=4 and index<3):
            result=cpp(s,directory/'slice_tasks',r['id'])
            if rank<=4 and index<3:
                check=python_orbit(s);assert all(result[k]==check[k] for k in ('states','minimum'))
                controls[r['id']]=check
        signature={'weighted_terminal_cycle_lengths':key,'orbit_minimum':result['minimum'] if len(groups[key])>1 else None}
        encoded=json.dumps(signature,sort_keys=True,separators=(',',':'));assert encoded not in minima,r['id'];minima.add(encoded)
        s.update({'canonical_signature':signature,'signature_sha256':hashlib.sha256(encoded.encode()).hexdigest(),
                  'signature_convention':'Weighted terminal cycle lengths and complete cyclic-loop minima within collision groups.',
                  'distinctness':{'weighted_terminal_cycle_lengths':key,'records_with_same_invariant':groups[key],
                                  'method':'exhaustive weighted cyclic-loop orbit' if len(groups[key])>1 else 'unique weighted terminal cycle lengths'},
                  'equivalent_canonical_states_including_sign_exchange':result['states'] if result else None})
        results[r['id']]=s;print(r['id'],'slice distinct;',result['states'] if result else 'unique invariant',flush=True)
    path.write_text(json.dumps(rs,indent=2)+'\n',encoding='utf8',newline='\n')
    (directory/'slice_signatures.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
    (directory/'slice-controls.json').write_text(json.dumps(controls,indent=2)+'\n',encoding='utf8',newline='\n')
    report=json.loads((directory/'classification.json').read_text());report['slice_distinctness_complete']=True;report['complete']=True
    (directory/'classification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:out.write(json.dumps({'rank':rank,'stage':'weighted-slice-distinctness','wall_seconds':time.monotonic()-start,'families':len(rs),'complete':True})+'\n')
if __name__=='__main__':main(int(sys.argv[1]))
