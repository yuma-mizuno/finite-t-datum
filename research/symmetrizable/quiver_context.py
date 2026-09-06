"""Source-grounded geometric consequences of already certified mutation classes."""
def enrich(r,q):
    b=r['slice']['B'];n=len(b)
    if q['status']=='certified-finite-orbit' and n>6 and any(b[i][j]!=-b[j][i] for i in range(n) for j in range(n)):
        q['geometric_identification']={
            'label':'Orbifold type',
            'status':'deduced from the finite mutation classification',
            'explanation':f'The complete orbit proves mutation finiteness. This connected, non-skew-symmetric matrix has {n} vertices. The seven non-skew-symmetric exceptions in the finite mutation classification have at most six vertices, so its diagram is s-decomposable and comes from a triangulated weighted orbifold. The topology and marked-point data have not been identified; the stored orbit specifies the exact valued mutation class.',
            'references':[{'title':'Felikson–Shapiro–Tumarkin, Theorem 5.13 and the proof of Theorem 5.1 (orders of the seven exceptions)','url':'https://arxiv.org/abs/1006.4276'},
                          {'title':'Felikson–Shapiro–Tumarkin, Lemma 4.11 and Theorem 4.19 (weighted orbifolds)','url':'https://arxiv.org/abs/1111.3449'}]}
    return q
