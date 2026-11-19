#!/usr/bin/env python3
"""Exact tensor and named-example checks for the incidence audit."""
from series_filtration import basis,series,lbl
from itertools import combinations,combinations_with_replacement,product
from fractions import Fraction
from random import Random

def fundamental(n,es):
    adj=[[] for _ in range(n)];B=[]
    for i,(u,v) in enumerate(es):
        stack=[(u,0)];seen=set();path=None
        while stack:
            a,p=stack.pop()
            if a==v:path=p;break
            if a in seen:continue
            seen.add(a);stack.extend((b,p|(1<<j)) for b,j in adj[a] if b not in seen)
        if path is not None:B.append(path|(1<<i))
        else:adj[u].append((v,i));adj[v].append((u,i))
    return B

def even_rank(es,n,L,B,dim=4):
    # Two independent alternating coefficients e123 and e124 (zero based).
    terms=[(0,1,2),(0,1,3)];h=len(B);m=len(es);coords=list(product(range(h),range(dim)))
    rows=[]
    for (b,v),(c,w) in combinations(coords,2):
        row=0
        for i,(a,u) in enumerate(coords):
            triple=tuple(sorted((u,v,w)))
            if triple in terms:
                val=(lbl(B[a]&B[b]&B[c],L)>>terms.index(triple))&1;row|=val<<i
        rows.append(row)
    direct=len(basis(rows));classrows=[]
    for mask in series(B,m).values():
        j=(mask&-mask).bit_length()-1;gamma=lbl(mask,L)
        for v,w in combinations(range(dim),2):
            row=0
            for i,(a,u) in enumerate(coords):
                triple=tuple(sorted((u,v,w)))
                if triple in terms and (B[a]>>j&1) and (gamma>>terms.index(triple)&1):row|=1<<i
            classrows.append(row)
    assert direct==len(basis(classrows));return direct

def cycle(v):return list(zip(v,v[1:]+v[:1]))
def bridgechain(lengths):
    es=[];n=0;last=None
    for l in lengths:
        v=list(range(n,n+l));es+=cycle(v)
        if last is not None:es.append((last,v[0]))
        last=v[1];n+=l
    return n,es

def phase_moment(es,B):
    # Contraction ranks for direct binary CCZ volumes on each cycle coordinate.
    h=len(B);D=3*h;co=[(a,u) for a in range(h) for u in range(3)]
    mats=[]
    for a,u in co:
        mat=[]
        for b,v in co:
            mat.append(sum(int(len({u,v,w})==3 and (B[a]&B[b]&B[c]).bit_count()%2)<<j for j,(c,w) in enumerate(co)))
        mats.append(mat)
    from collections import Counter
    hist=Counter()
    for x in range(1<<D):
        rows=[0]*D
        for i in range(D):
            if x>>i&1:rows=[a^b for a,b in zip(rows,mats[i])]
        hist[len(basis(rows))]+=1
    return D-(hist[0].bit_length()-1),sum(Fraction(a,1<<r) for r,a in hist.items())/(1<<D)

def main():
    rng=Random(9181)
    for i in range(120):
        n=6;es=[e for e in combinations(range(n),2) if rng.random()<.48];L=[rng.randrange(4) for _ in es];B=fundamental(n,es);even_rank(es,n,L,B)
    print('Heterogeneous tensor radical comparisons: 120 PASS')
    theta=[(0,1),(0,2),(2,3),(3,1),(0,4),(4,5),(5,1)]
    cactus=cycle([0,1,2])+cycle([0,3,4])+[(4,5)]
    for name,(n,es) in [('triangle-pentagon',bridgechain([3,5])),('two-squares',bridgechain([4,4])),('theta',(6,theta)),('two-triangle cactus',(6,cactus))]:
        B=fundamental(n,es);eta,m=phase_moment(es,B);print(name,'V',n,'E',len(es),'h',len(B),'eta',eta,'P2',m)
    assert phase_moment(theta,fundamental(6,theta))==(6,Fraction(79,1024))
    assert phase_moment(cactus,fundamental(6,cactus))==(6,Fraction(121,1024))
    for j in range(4):
        for G,mom in [(theta,Fraction(79,1024)),(cactus,Fraction(121,1024))]:
            es=list(G);n=6;last=2
            for k in range(j):
                es+=cycle([n,n+1,n+2])+[(last,n)];last=n+1;n+=3
            assert n==6+3*j and len(es)==7+4*j and len(fundamental(n,es))==j+2
        print('matched family j=',j,'V,E,h',6+3*j,7+4*j,j+2,'moments',Fraction(79,1024)*Fraction(11,32)**j,Fraction(121,1024)*Fraction(11,32)**j)
    # Threshold by explicit simple-cycle masks for the named models.
    for name,es,cs,expected in [('triangle',cycle([0,1,2]),[7],3),('square',cycle([0,1,2,3]),[15],None),('theta',theta,[15,113,126],4),('C9',cycle(list(range(9))),[511],9)]:
        val=min((max(a.bit_count(),b.bit_count(),c.bit_count()) for a,b,c in combinations_with_replacement(cs,3) if (a&b&c).bit_count()%2),default=None);assert val==expected;print(name,'support threshold',val)
if __name__=='__main__':main()
