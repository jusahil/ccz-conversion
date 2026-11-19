#!/usr/bin/env python3
"""Independent exact tensor and graph regression (standard library)."""
from itertools import combinations, product, permutations
from collections import Counter
from fractions import Fraction
from random import Random

def span(bs):
    v=[0]
    for b in bs:v += [a^b for a in v]
    return v

def basis(vs):
    p={}
    for x in vs:
        while x:
            k=x.bit_length()-1
            if k in p:x^=p[k]
            else:p[k]=x;break
    return list(p.values())

def subs(n):
    ss={frozenset([0])};pending=list(ss)
    while pending:
        w=pending.pop()
        for x in range(1<<n):
            if x not in w:
                t=w|frozenset(x^y for y in w)
                if t not in ss:ss.add(t);pending.append(t)
    return ss

def ev(terms,a,b,c):
    out=0
    for i,j,k in terms:
        for x,y,z in ((i,j,k),(i,k,j),(j,i,k),(j,k,i),(k,i,j),(k,j,i)):
            out^=((a>>x)&(b>>y)&(c>>z)&1)
    return out

def restrict_nonzero(fun,bs):return any(fun(*t) for t in combinations(bs,3))

def rank_contract(fun,a,bs):
    return len(basis(sum(fun(a,b,c)<<i for i,c in enumerate(bs)) for b in bs))

def mapv(cols,x):
    out=0
    for i,c in enumerate(cols):
        if (x>>i)&1:out^=c
    return out

def tensor_terms(n,mask):return [t for i,t in enumerate(combinations(range(n),3)) if mask>>i&1]

def pencil_tests():
    n=4;unit=[1<<i for i in range(n)];ss=subs(n);count=0
    for mask in range(16):
        terms=tensor_terms(n,mask)
        # cache original multilinear evaluations
        d=[[[ev(terms,a,b,c) for c in range(16)] for b in range(16)] for a in range(16)]
        for R in ss:
            rb=basis(R);bs=unit+[r<<n for r in rb]; u=span(bs)
            def T(a,b,c):
                a0,a1=a&15,a>>4;b0,b1=b&15,b>>4;c0,c1=c&15,c>>4
                return d[a0][b0][c0]^d[a0^a1][b0^b1][c0^c1]
            Z={a for a in range(16) if all(not d[a][r][b] for r in rb for b in unit)}
            radR={s for s in R if all(not d[s][b][c] for b,c in combinations(unit,2))}
            actual={a for a in u if not rank_contract(T,a,bs)}
            assert actual=={a|(s<<n) for a in Z for s in radR}
            for a in unit+[0]:
                rows=[sum(d[a][r][b]<<i for i,b in enumerate(unit)) for r in rb]
                assert rank_contract(T,a,bs)==2*len(basis(rows))
            count+=1
    print('Pencil/radical exhaustive comparisons:',count, 'PASS')

def hetero(cols1,cols2,terms1,terms2,rcols=()):
    n=3; rm=set(span(rcols)); col=cols1+cols2
    K=basis(x for x in range(64) if mapv(col,x) in rm)
    def T(a,b,c):return ev(terms1,a&7,b&7,c&7)^ev(terms2,a>>3,b>>3,c>>3)
    actual=restrict_nonzero(T,K)
    S=[set(span(cols1)),set(span(cols2))];V=[];W=[];P=[]
    for i,c in enumerate((cols1,cols2)):
        P.append(S[i]&{v^r for v in S[1-i] for r in rm});V.append(basis(x for x in range(8) if mapv(c,x) in P[-1]));W. append(basis(x for x in range(8) if not mapv(c,x)))
    mixed=any(ev(t,w,a,b) for i,t in enumerate((terms1,terms2)) for w in W[i] for a,b in combinations(V[i],2))
    if mixed:pred=True
    else:
        # Quotient tensor on reachable labelled charge pairs. Lift each pair of labels independently.
        h=basis(a|(b<<4) for a in P[0] for b in P[1] if a^b in rm)
        def lift(v):
            return next(x for x in range(8) if mapv(cols1,x)==(v&15))|(next(x for x in range(8) if mapv(cols2,x)==(v>>4)) <<3)
        pred=restrict_nonzero(T,[lift(x) for x in h])
    assert actual==pred

def hetero_tests():
    c=list(product(range(4),repeat=3)); t=[[],[(0,1,2)]]
    for a,b,i,j in product(c,c,t,t):hetero(a,b,i,j)
    rng=Random(20261117)
    for _ in range(400):hetero(tuple(rng.randrange(8) for _ in range(3)),tuple(rng.randrange(8) for _ in range(3)),rng. choice(t),rng.choice(t),basis(rng.randrange(8) for _ in range(rng.randrange(4))))
    print('Heterogeneous reachability:',len(c)**2*4,'exhaustive + 400 seeded with references PASS')

def fano_tests():
    for lab,triples in [('A',[(1,2,3),(4,5,6),(1,4,7),(2,5,7)]),('B',[(1,2,3),(1,4,5),(1,6,7),(2,4,6),(2,5,7),(3,4,7) ,(3,5,6)])]:
        tr=[tuple(i-1 for i in t) for t in triples];fun=lambda a,b,c:ev(tr,a,b,c);units=[1<<i for i in range(7)]
        hist=Counter(rank_contract(fun,a,units) for a in range(128));mom=lambda k:sum(Fraction(c,2**(h*(k-1))) for h,c in hist.items())/128
        H=[2,4,8,32,64] if lab=='A' else [1^2,1^8,4^16,4^32,1^64]
        assert len(basis(H))==5 and not restrict_nonzero(fun,H)
        assert hist[2]<63 and hist[0]==1 and mom(2)==Fraction(29,512)
        print(lab,'contraction histogram',dict(sorted(hist.items())),'P2',mom(2),'P3',mom(3),'kappa=2')

def graph_tests():
    allgraphs=0;activeall=0
    for n in range(1,7):
        edges=list(combinations(range(n),2));m=len(edges);ed={e:1<<i for i,e in enumerate(edges)}
        def pathmask(path):return sum(ed[tuple(sorted((a,b)))] for a,b in zip(path,path[1:]))
        cycles=set()
        for k in range(3,n+1):
            for v in combinations(range(n),k):
                for tail in permutations(v[1:]):
                    if tail[0]<tail[-1]:cycles.add(pathmask((v[0],)+tail+(v[0],)))
        # independent motif enumeration: all odd cycles and 3 internally disjoint odd paths.
        motifs={c for c in cycles if c.bit_count()%2}
        for a,b in edges:
            rem=[v for v in range(n) if v not in (a,b)];paths=[]
            for k in range(0,len(rem)+1,2): # k internal vertices -> odd length
                for p in permutations(rem,k):paths.append((sum(1<<v for v in p),pathmask((a,)+p+(b,))))
            for (i,x),(j,y),(k,z) in combinations(paths,3):
                if not ((i&j)|(i&k)|(j&k)|(x&y)|(x&z)|(y&z)):motifs.add(x|y|z)
        inf=m+1;best=[inf]*(1<<m);mot=[inf]*(1<<m)
        for c in cycles:
            for d in cycles:
                if (c&d).bit_count()%2:best[c|d]=(c|d).bit_count()
        for x in motifs:mot[x]=x.bit_count()
        for i in range(m):
            b=1<<i
            for x in range(1<<m):
                if x&b:best[x]=min(best[x],best[x^b]);mot[x]=min(mot[x],mot[x^b])
        assert mot==best
        # Graph coedge isolation + tensor activation and radical using independent cycle basis.
        if n<=5:
            for x in range(1<<m):
                cb=basis(c for c in cycles if c&x==c);C=span(cb)
                lam=[sum(((c>>i)&1)<<j for j,c in enumerate(cb)) for i in range(m) if x>>i&1]
                counts=Counter(a for a in lam if a);odd=[a for a,v in counts.items() if v%2]
                assert bool(odd)==(best[x]!=inf)
                # test trilinear overlap restricted to basis and repeated slots; and radical
                act=any((a&b&c).bit_count()%2 for a,b,c in product(cb,repeat=3));assert act==bool(odd)
                rad={a for a in C if not any((a&b&c).bit_count()%2 for b,c in product(cb,repeat=2))}
                assert len(rad)==1<<(len(cb)-len(basis(odd)))
        allgraphs+=1<<m;activeall+=sum(v!=inf for v in best)
        print('n',n,'graphs',1<<m,'active',sum(v!=inf for v in best),'motifs',len(motifs),'minimum-cost equality PASS')
    print('Graph motif minimum comparisons:',allgraphs,'active:',activeall,'PASS')
if __name__=='__main__':
    pencil_tests();hetero_tests();fano_tests();graph_tests()

def hetero_graph_tests():
    # Tensor labels represented by two independent coefficient bits, nonnegative costs.
    n=5; edges=list(combinations(range(n),2)); m=len(edges); ee={e:1<<i for i,e in enumerate(edges)}
    def pm(p):return sum(ee[tuple(sorted((a,b)))] for a,b in zip(p,p[1:]))
    cycles=set()
    for k in range(3,n+1):
        for vs in combinations(range(n),k):
            for tail in permutations(vs[1:]):
                if tail[0]<tail[-1]:cycles.add(pm((vs[0],)+tail+(vs[0],)))
    thetas=[]
    for a,b in edges:
        rem=[v for v in range(n) if v not in (a,b)];paths=[]
        for k in range(len(rem)+1):
            for p in permutations(rem,k):paths.append((sum(1<<v for v in p),pm((a,)+p+(b,))))
        for (i,x),(j,y),(k,z) in combinations(paths,3):
            if not ((i&j)|(i&k)|(j&k)|(x&y)|(x&z)|(y&z)):thetas.append((x,y,z))
    rng=Random(8871); tests=0
    for _ in range(160):
        tensor=[rng.randrange(4) for e in edges]; weights=[rng.randrange(6) for e in edges]; mask=rng.randrange(1<<m)
        def label(x):
            r=0
            for i in range(m):
                if x>>i&1:r^=tensor[i]
            return r
        def cost(x):return sum(weights[i] for i in range(m) if x>>i&1)
        cc=[c for c in cycles if c&mask==c];pairs=[c|d for c,d in product(cc,repeat=2) if label(c&d)]
        candidates=[c for c in cc if label(c)] + [x|y|z for x,y,z in thetas if (x|y|z)&mask==(x|y|z) and label(x)==label(y) ==label(z)!=0]
        assert min(map(cost,pairs),default=10**6)==min(map(cost,candidates),default=10**6)
        # Independent actual pair-overlap activation, via linear basis of cycles.
        cb=basis(cc); assert any(label(a&b) for a,b in product(cb,repeat=2))==bool(candidates)
        tests+=1
    print('Heterogeneous weighted graph motif tests:',tests,'seeded graphs with 2-bit tensor labels PASS')
if __name__=='__main__': hetero_graph_tests()
