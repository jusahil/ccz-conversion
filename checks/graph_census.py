"""Exhaustive labelled graphs through 6 vertices; independent pair and path-motif algorithms."""
from itertools import combinations, permutations
from math import inf

def transform(a,m):
    # Min over all edge subsets of current graph, in-place.
    for b in range(m):
        for x in range(1<<m):
            if x>>b&1:a[x]=min(a[x],a[x^(1<<b)])

def rank(vs):
    p={}
    for x in vs:
        while x:
            k=x.bit_length()-1
            if k in p:x^=p[k]
            else:p[k]=x;break
    return len(p)

def main(n):
    es=list(combinations(range(n),2));eid={e:1<<i for i,e in enumerate(es)};m=len(es);N=1<<m
    ed=lambda u,v:eid[tuple(sorted((u,v)))]
    cyc=set()
    for k in range(3,n+1):
        for vertices in combinations(range(n),k):
            for tail in permutations(vertices[1:]):
                c=(vertices[0],)+tail;cyc.add(sum(ed(c[i-1],c[i]) for i in range(k)))
    pair=[inf]*N;motif=[inf]*N
    C=sorted(cyc)
    for i,a in enumerate(C):
        for b in C[i:]:
            if (a&b).bit_count()%2:pair[a|b]=min(pair[a|b],(a|b).bit_count())
    for a in C:
        if a.bit_count()%2:motif[a]=a.bit_count()
    # independently enumerate all odd simple paths between endpoints; combine disjoint interiors and edges.
    for u,v in es:
        rest=[x for x in range(n) if x not in (u,v)];paths=[]
        for k in range(0,len(rest)+1,2):
            for mid in permutations(rest,k):
                seq=(u,)+mid+(v,);mask=sum(ed(a,b) for a,b in zip(seq,seq[1:]));inter=sum(1<<x for x in mid);paths.append(( mask,inter))
        for (a,x),(b,y),(c,z) in combinations(paths,3):
            if not ((x&y)|(x&z)|(y&z)|(a&b)|(a&c)|(b&c)):
                w=a|b|c;motif[w]=min(motif[w],w.bit_count())
    transform(pair,m);transform(motif,m);assert pair==motif
    # independently evaluate non-self-orthogonality of cycle space on a greedy spanning cycle basis.
    # Verify that pair/motif existence agrees in every graph, and sample exact minimum support using basis.
    active=0
    for g in range(N):
        B=[];p={}
        for c in C:
            if c&g!=c:continue
            x=c
            while x:
                k=x.bit_length()-1
                if k in p:x^=p[k]
                else:p[k]=x;B.append(c);break
        act=any((a&b).bit_count()%2 for a in B for b in B)
        assert act==(pair[g]<inf);active+=act
    print('labelled vertices',n,'graphs',N,'cycles in complete graph',len(C),'active graphs',active,flush=True)
    return N,active
if __name__=='__main__':
    rows=[main(n) for n in range(1,7)];print('total',sum(x[0] for x in rows),'active',sum(x[1] for x in rows))
