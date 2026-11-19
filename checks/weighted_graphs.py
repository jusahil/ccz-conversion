from itertools import combinations,permutations
from random import Random
from math import inf

def xor(v):
    r=0
    for x in v:r^=x
    return r
rng=Random(7162)
active=0
for case in range(160):
    n=6;es=[e for e in combinations(range(n),2) if rng.random()<.55];m=len(es);eid={e:1<<i for i,e in enumerate(es)};label =[rng.randrange(4) for e in es];cost=[rng.randrange(6) for e in es]
    def ed(a,b):return eid.get(tuple(sorted((a,b))),0)
    def mask_path(p,cycle=False):
        ee=[ed(a,b) for a,b in zip(p,p[1:])]+([ed(p[-1],p[0])] if cycle else [])
        return sum(ee) if all(ee) else 0
    def lab(x):return xor(v for i,v in enumerate(label) if x>>i&1)
    def wt(x):return sum(v for i,v in enumerate(cost) if x>>i&1)
    cs=set()
    for k in range(3,n+1):
        for vs in combinations(range(n),k):
            for tail in permutations(vs[1:]):
                x=mask_path((vs[0],)+tail,True)
                if x:cs.add(x)
    A=min((wt(a|b) for a in cs for b in cs if lab(a&b)),default=inf)
    B=min((wt(c) for c in cs if lab(c)),default=inf)
    for u,v in combinations(range(n),2):
        r=[x for x in range(n) if x not in (u,v)];ps=[]
        for k in range(len(r)+1):
            for mid in permutations(r,k):
                x=mask_path((u,)+mid+(v,))
                if x and lab(x):ps.append((x,sum(1<<i for i in mid),lab(x)))
        for (a,x,i),(b,y,j),(c,z,k) in combinations(ps,3):
            if i==j==k and not ((a&b)|(a&c)|(b&c)|(x&y)|(x&z)|(y&z)):B=min(B,wt(a|b|c))
    assert A==B;active+=A<inf
print('160 seeded tensor-labelled graphs (two independent label coefficients, nonnegative costs); pair and minimal motif optimum agree; active',active)
