"""Finite pulse enumeration for a single attachment to a classified core.

Only a necessary subsystem is tested. An empty result is an exact exclusion;
a surviving pulse pattern does not assert a complete T-datum.
"""
from functools import lru_cache
from obstructions import attachments
from principal_relations import labelled_records


def core(candidate,s,ell,a,top):
    counts=[candidate['N_plus_1'],candidate['N_minus_1']]
    key=tuple(b[i][j] for b in counts for i in top for j in top)
    found=labelled_records().get(len(top),{}).get(key)
    if not found:return None
    record,exchange,positions=found;permutation=[top[i] for i in positions]
    target=tuple(tuple(counts[t][ell][i] for i in permutation) for t in (s,1-s))
    return record,s^exchange,permutation.index(a),target,counts[s][ell][ell],permutation


@lru_cache(maxsize=4096)
def patterns(record_id,flip,a,target,diagonal):
    records={r[0]['id']:r[0] for table in labelled_records().values() for r in table.values()}
    datum=records[record_id]['datum'];n=len(target[0]);height=max(datum['delays']);mass=sum(map(sum,target))
    assert diagonal in (0,1) and mass>0
    # Each h consecutive zero coefficients would terminate the response.
    # Before the last pulse that would leave an uncancelled single pulse.
    bound=(mass+1)*height
    ops=[]
    for i,d in enumerate(datum['delays']):ops.append((d,i,i,0,1))
    names=('N_plus','N_minus') if not flip else ('N_minus','N_plus')
    for sign,name in enumerate(names):
        for i,row in enumerate(datum[name]):
            for j,terms in enumerate(row):
                for coefficient,d in terms:ops.append((d,i,j,1 if sign else -1,coefficient))
    def response(history,t,pulse=0):
        w=[0]*n;w[a]=pulse
        for d,i,j,kind,c in ops:
            if t>=d:
                value=history[t-d][j] if t-d<len(history) else 0
                w[i]+=(-value if kind==0 else max(value,0) if kind==1 else min(value,0))*c
        return w
    def add(counts,w):
        for j,v in enumerate(w):counts[0 if v>=0 else 1][j]+=abs(v)
        return any(counts[s][j]>target[s][j] for s in range(2) for j in range(n))
    initial=[];counts=[[0]*n for _ in range(2)];first_excess=None
    for t in range(bound+1):
        w=response(initial,t,-1 if t==0 else 0);initial.append(w)
        if add(counts,w):first_excess=t;break
    assert first_excess is not None,'Single-pulse mass bound failed'
    choices=range(1,first_excess+1) if diagonal else (None,)
    solutions=[]
    for middle in choices:
        history=[];counts=[[0]*n for _ in range(2)];zeros=0;stopped=False
        for t in range(bound+1):
            pulse=-1 if t==0 else 1 if t==middle else 0
            w=response(history,t,pulse)
            if t>(middle or 0) and counts==[list(x) for x in target] and w==[int(i==a) for i in range(n)]:
                # The last negative pulse cancels w. A full zero memory state
                # after h further coefficients proves the whole tail is zero.
                if all(not any(response(history,u)) for u in range(t+1,t+height+1)):
                    solutions.append({'delay':t,'middle_pulse':middle,
                                      'diagonal_exponent':t-middle if middle is not None else None,
                                      'row_exponents':[[[t-degree for degree in range(len(history)-1,-1,-1)
                                                       for _ in range(max((1 if s==0 else -1)*history[degree][j],0))]
                                                      for j in range(n)] for s in range(2)]})
            history.append(w)
            if add(counts,w):stopped=True;break
            zeros=0 if any(w) else zeros+1
            if zeros>=height and (middle is None or t>=middle):stopped=True;break
        assert stopped,'Pulse search reached its proved bound without stopping'
    return {'core_record':record_id,'core_sign_exchange':bool(flip),'core_attachment_species':a,
            'target_coefficient_sums':[list(x) for x in target],'diagonal_coefficient_sum':diagonal,
            'core_height':height,'row_mass':mass,'proved_delay_bound':bound,
            'middle_pulses_tested':len(choices),'patterns':solutions}


def exclusion(candidate):
    for s,ell,a,top in attachments(candidate):
        found=core(candidate,s,ell,a,top)
        if found is None:continue
        record,flip,local_a,target,diagonal,permutation=found
        result=patterns(record['id'],flip,local_a,target,diagonal)
        if not result['patterns']:
            return {'lemma':'finite pulse enumeration for a single attachment to a classified principal core',
                    'extra_species':ell,'attachment_species':a,'sign_exchange':bool(s),
                    'core_to_original_species':permutation,**result}
    return None
