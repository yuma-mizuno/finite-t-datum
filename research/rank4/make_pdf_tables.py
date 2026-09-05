"""Typeset all representatives and their constant matrices from exact data."""
import json
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
z=sp.Symbol('z')


def matrix(a):
    return r'\begin{pmatrix}'+r'\\'.join(' & '.join(sp.latex(sp.expand(x)) for x in row) for row in a.tolist())+r'\end{pmatrix}'


def main():
    verification=json.loads((HERE/'verification.json').read_text())
    assert verification['replay_counts']=={'unsat':4865}
    fs=[json.loads(line) for line in (HERE/'families.jsonl').read_text().splitlines()]
    fs=sorted(fs,key=lambda x:x['id'])
    assert len(fs)==37 and all(x['coverage_status']=='unsat' and len(x['spaces'])==1 for x in fs)
    constants={x['id']:x for x in json.loads((HERE/'constant_candidates.json').read_text())['candidates']}
    slices={x['id']:x for x in json.loads((HERE/'slice_signatures.json').read_text())}
    rows=[];catalogue=[];classification=[]
    for k,x in enumerate(fs,1):
        space=x['spaces'][0];c=space['certificate'];r=c['delays']
        assert space['scaling_and_shifts'] and c['positive'] and c['negative']
        hp,hm,period=c['positive']['h'],c['negative']['h'],c['labelled_tropical_seed_period']
        nv=slices[x['id']]['slice_vertices']
        rows.append(f'{k} & {x["id"]} & $({",".join(map(str,r))})$ & {hp} & {hm} & {period} & {nv}'+r'\\')
        if k==1 or (k-3)%3==0:
            if k>1:catalogue.append(r'\clearpage')
            end=2 if k==1 else min(k+2,len(fs))
            catalogue.append(r'\subsection*{'+('Class '+str(k) if k==end else 'Classes '+str(k)+'--'+str(end))+'}')
            catalogue.append(r'\noindent $B_\pm^{(k)}(z)=\operatorname{diag}(1+z^{r_1},\ldots,1+z^{r_4})-N_\pm^{(k)}(z)$.\par\medskip')
        catalogue.append(r'\noindent\begin{minipage}{\linewidth}')
        catalogue.append(r'\noindent\textbf{Class '+str(k)+r'}\quad\small (constant ID '+str(x['id'])+r')\hfill $\symbf{r}=('+','.join(map(str,r))+')$')
        catalogue.append(r'\par\smallskip\noindent\textcolor{rulegray}{\rule{\linewidth}{0.3pt}}')
        n0=sp.diag(*(1+z**v for v in r));pairs=[]
        for sign,field,cf in [('+','A_plus','N_plus_1'),('-','A_minus','N_minus_1')]:
            a=sp.Matrix([[sp.sympify(v,locals={'z':z}) for v in row] for row in c[field]])
            n=n0-a;a1=a.subs(z,1)
            assert a1==2*sp.eye(4)-sp.Matrix(constants[x['id']][cf])
            pairs.append((sign,n,a1))
        for which in [1,2]:
            if which==2:catalogue.append(r'\par\noindent')
            for pair in pairs:
                sign=pair[0]
                label=r'N_{'+sign+r'}^{('+str(k)+r')}(z)' if which==1 else r'A_{'+sign+r'}(1)'
                catalogue.append(r'\begin{minipage}[c]{0.495\linewidth}\centering\fontsize{9.7}{12}\selectfont')
                catalogue.append(r'\setlength{\abovedisplayskip}{4pt}\setlength{\belowdisplayskip}{4pt}')
                catalogue.append(r'\setlength{\abovedisplayshortskip}{4pt}\setlength{\belowdisplayshortskip}{4pt}')
                catalogue.append(r'\setlength{\arraycolsep}{3pt}\renewcommand{\arraystretch}{1.08}')
                catalogue.append(r'\['+label+'='+matrix(pair[which])+r'\]')
                catalogue.append(r'\end{minipage}%')
        catalogue.append(r'\par\smallskip\noindent\small $h_+='+str(hp)+r',\quad h_-='+str(hm)+r',\quad \Omega='+str(period)+r'.$\quad Vertices per slice: '+str(nv)+'.')
        catalogue.append(r'\end{minipage}\par\vfill')
        classification.append({'class':k,'constant_id':x['id'],'delays':r,'h_plus':hp,'h_minus':hm,'period':period,
                               'slice_vertices':nv,'A_plus':c['A_plus'],'A_minus':c['A_minus'],
                               'A_plus_1':[[int(v) for v in row] for row in pairs[0][2].tolist()],
                               'A_minus_1':[[int(v) for v in row] for row in pairs[1][2].tolist()]})
    table=[r'\begin{tabular}{@{}crcrrrr@{}}',r'\toprule',r'Class & ID & Delays & $h_+$ & $h_-$ & $\Omega$ & Slice size\\',r'\midrule']
    table+=rows
    table+=[r'\bottomrule',r'\end{tabular}']
    (HERE/'rank4-overview.tex').write_text('\n'.join(table)+'\n',encoding='utf-8',newline='\n')
    (HERE/'rank4-catalogue.tex').write_text('\n'.join(catalogue)+'\n',encoding='utf-8',newline='\n')
    (HERE/'classification.json').write_text(json.dumps(classification,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Generated 37 polynomial pairs and 74 checked constant matrices.')


if __name__=='__main__':main()
