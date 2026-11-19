#!/usr/bin/env python3
"""Independent exact tests of graph activation and filtration (standard library)."""
from itertools import combinations, combinations_with_replacement
from collections import Counter
from random import Random
from functools import reduce
from operator import xor

def basis(seq):
    piv={};out=[]
    for v in seq:
        x=v
        while x:
            i=x.bit_length()-1
            if i in piv:x^=piv[i]
            else:piv[i]=x;out.append(v);break
    return out

def lbl(x,L):return reduce(xor,(a for i,a in enumerate(L) if x>>i&1),0)
def wt(x,c):return sum(a for i,a in enumerate(c) if x>>i&1)
def cycles(n,es):
    # Enumerate connected even-incidence degree-two edge sets, not vertex permutations.
    m=len(es);out=[];d=[(1<<u)^(1<<v) for u,v in es];par=[0]*(1<<m)
    for x in range(1,1<<m):
        b=x&-x;par[x]=par[x^b]^d[b.bit_length()-1]
        if par[x] or x.bit_count()<3:continue
        adj={}
        for i,(u,v) in enumerate(es):
            if x>>i&1:adj.setdefault(u,[]).append(v);adj.setdefault(v,[]).append(u)
        if any(len(a)!=2 for a in adj.values()):continue
        seen=set();todo=[next(iter(adj))]
        while todo:
            u=todo.pop()
            if u not in seen:seen.add(u);todo.extend(adj[u])
        if len(seen)==len(adj):out.append(x)
    return sorted(out,key=lambda c:(c.bit_count(),c))
def theta(g,es):
    adj={}
    for i,(u,v) in enumerate(es):
        if g>>i&1:adj.setdefault(u,[]).append((v,i));adj.setdefault(v,[]).append((u,i))
    ends=[u for u in adj if len(adj[u])==3]
    if len(ends)!=2 or any(len(a) not in (2,3) for a in adj.values()):return None
    s,t=ends;paths=[]
    for u,i in adj[s]:
        mask=1<<i;prev=s
        while u!=t:
            if u==s:return None
            v,j=next((v,j) for v,j in adj[u] if v!=prev)
            if mask>>j&1:return None
            mask|=1<<j;prev,u=u,v
        paths.append(mask)
    if (paths[0]|paths[1]|paths[2])!=g or any(a&b for a,b in combinations(paths,2)):return None
    return paths

def motif(g,es,C,L):
    if g in C and lbl(g,L):return True
    ps=theta(g,es)
    return bool(ps and lbl(ps[0],L)!=0 and len({lbl(x,L) for x in ps})==1)
def series(B,m):
    S={}
    for i in range(m):
        col=sum(((b>>i)&1)<<j for j,b in enumerate(B))
        if col:S[col]=S.get(col,0)|(1<<i)
    return S

def main():
    total=active_total=classes=filtrations=0
    for n in range(1,7):
        es=list(combinations(range(n),2));m=len(es);C=cycles(n,es);Cset=set(C);N=1<<m;L=[1]*m
        act=[False]*N;ac=mc=0
        for g in range(N):
            cs=[c for c in C if c&g==c];B=basis(cs);S=series(B,m)
            for sm in S.values():
                containing=[c for c in cs if c&sm==sm]
                assert any((a&b)==sm for a,b in combinations_with_replacement(containing,2));classes+=1
            p=any((a&b).bit_count()&1 for a,b in combinations_with_replacement(B,2))
            t=any((a&b&c).bit_count()&1 for a,b,c in combinations_with_replacement(B,3))
            assert p==t==any(sm.bit_count()%2 for sm in S.values());act[g]=t;ac+=t
            isminimal=t and not any(act[g^(1<<i)] for i in range(m) if g>>i&1)
            assert isminimal==motif(g,es,Cset,L);mc+=isminimal
            for s in range(3,n+1):
                short=[b for b in B if b.bit_count()<=s];available=[c for c in cs if c.bit_count()<=s]
                assert len(short)==len(basis(available));filtrations+=1
        total+=N;active_total+=ac;print('n',n,'graphs',N,'active',ac,'minimal',mc,flush=True)
    print('census',total,active_total,'series-isolation checks',classes,'filtration ranks',filtrations,flush=True)
    rng=Random(47129);active=zero_period=weighted_filtrations=0;subsets=0
    n=5;es=list(combinations(range(n),2));m=len(es);C=cycles(n,es);Cset=set(C)
    for k in range(180):
        g=rng.randrange(1<<m);cost=[rng.randrange(6) for _ in es];phi=[rng.randrange(4) for _ in range(n)]
        L=[phi[a]^phi[b] for a,b in es] if k%2==0 else [rng.randrange(4) for _ in es]
        cs=[c for c in C if c&g==c];B=basis(cs);zero_period+=not any(lbl(c,L) for c in cs)
        best=motbest=None;x=g
        while True:
            Bs=basis(c for c in cs if c&x==c)
            t=any(lbl(a&b&c,L) for a,b,c in combinations_with_replacement(Bs,3));subsets+=1
            if t:best=min(best,wt(x,cost)) if best is not None else wt(x,cost)
            if motif(x,es,Cset,L):motbest=min(motbest,wt(x,cost)) if motbest is not None else wt(x,cost)
            if x==0:break
            x=(x-1)&g
        assert best==motbest;active+=best is not None
        for s in range(3,n+1):
            short=[c for c in cs if c.bit_count()<=s];Bs=[b for b in B if b.bit_count()<=s]
            assert any(lbl(a&b&c,L) for a,b,c in combinations_with_replacement(short,3))==any(lbl(a&b&c,L) for a,b,c in combinations_with_replacement(Bs,3));weighted_filtrations+=1
    print('weighted cases',180,'active',active,'zero-cycle-period',zero_period,'subsets',subsets,'weighted filtration comparisons',weighted_filtrations)
    for d in (3,4):
        rows=[sum(((v>>i)&1)<<(v-1) for v in range(1,1<<d)) for i in range(d)]
        byk=[[reduce(lambda a,b:a&b,t).bit_count()%2 for t in combinations(rows,k)] for k in (1,2,3)]
        assert not any(byk[0]+byk[1]);assert byk[2]==([1] if d==3 else [0]*4)
    print('non-graphic pair/triple counterexamples PASS')
if __name__=='__main__': main()
