"""Exact finite cohomology and control regressions. Python 3, standard library only.
No random sampling is used except labelled-graph checks with fixed seed.
"""
from itertools import combinations, permutations, product
from collections import Counter
from fractions import Fraction
from random import Random

def rank(rows):
    piv={}
    for x in rows:
        while x:
            j=x.bit_length()-1
            if j in piv: x ^= piv[j]
            else: piv[j]=x; break
    return len(piv)

def mesh(tets):
    ts=sorted(set(tuple(sorted(t)) for t in tets))
    fs=Counter(f for t in ts for f in combinations(t,3))
    assert set(fs.values())=={2}, 'not closed mod-two chain'
    es=sorted(set(e for t in ts for e in combinations(t,2)))
    vs=sorted(set(v for t in ts for v in t))
    vi={v:i for i,v in enumerate(vs)}; ei={e:i for i,e in enumerate(es)}
    r1=rank((1<<vi[a])^(1<<vi[b]) for a,b in es)
    r2=rank(sum(1<<ei[e] for e in combinations(f,2)) for f in fs)
    return ts,vs,es,list(fs),len(es)-r1-r2

def cup(ts,a,b,c):
    return sum(a[t[:2]]*b[t[1:3]]*c[t[2:]] for t in ts)%2

def check(name,tets,cochains,expected,exact_shifts=False):
    ts,vs,es,fs,b1=mesh(tets)
    assert len(cochains)==b1
    assert all(sum(a[e] for e in combinations(f,2))%2==0 for a in cochains for f in fs)
    # Check independence modulo exact 1-cochains by row rank.
    rr=[sum((1<<i) for i,e in enumerate(es) if v in e) for v in vs]
    assert rank(rr+[sum(a[e]<<i for i,e in enumerate(es)) for a in cochains])==len(vs)-1+b1
    tests=0
    for i,j,k in product(range(b1),repeat=3):
        assert cup(ts,cochains[i],cochains[j],cochains[k])==expected(i,j,k)
        tests+=1
    shift_tests=0
    if exact_shifts:
        for v in vs:
            d={e:int(v in e) for e in es}
            for a,b in product(cochains,repeat=2):
                assert cup(ts,d,a,b)==cup(ts,a,d,b)==cup(ts,a,b,d)==0
                shift_tests+=3
    print(name, 'V,E,F,T,b1=',len(vs),len(es),len(fs),len(ts),b1,'cup checks',tests,'elementary exact shifts',shift_tests)
    return ts,vs,es

def torus(N=3):
    v=list(product(range(N),repeat=3)); ids={x:i for i,x in enumerate(v)}
    ts=[]
    for x in v:
        for p in permutations(range(3)):
            u=x;t=[ids[u]]
            for i in p:
                u=tuple((a+(j==i))%N for j,a in enumerate(u));t.append(ids[u])
            ts.append(tuple(t))
    ts,vs,es,fs,b1=mesh(ts)
    # wrap-crossing cocycles, evaluated with a globally consistent vertex order
    aa=[{(u,w):int({v[u][j],v[w][j]}=={0,N-1}) for u,w in es} for j in range(3)]
    return ts,aa

def s1xs2():
    # Simplicial circle (triangle boundary) x tetrahedron boundary, staircase product.
    ts=[]
    for (a,b), tri in product(combinations(range(3),2),combinations(range(4),3)):
        for p in [(0,1,1),(1,0,1),(1,1,0)]:
            i=j=0; t=[4*a+tri[0]]
            for k in p:
                if k==0:i+=1
                else:j+=1
                t.append(4*(a,b)[i]+tri[j])
            ts.append(tuple(t))
    ts,vs,es,fs,b1=mesh(ts)
    aa=[{(u,w):int({u//4,w//4}=={0,2}) for u,w in es}]
    return ts,aa

def connected_sum(A,B):
    ta,aa=A;tb,bb=B
    ta,va,ea,_,_=mesh(ta);tb,vb,eb,_,_=mesh(tb)
    x=ta[0]; y=tb[0]
    mapb=dict(zip(y,x)); base=max(va)+1
    mapb.update({v:base+i for i,v in enumerate(v for v in vb if v not in mapb)})
    ts=ta[1:]+[tuple(sorted(mapb[v] for v in t)) for t in tb[1:]]
    ts,vs,es,_,_=mesh(ts)
    def normalize(a, simplex):
        z={v:a[tuple(sorted((simplex[0],v)))] for v in simplex[1:]}
        return {e:b^z.get(e[0],0)^z.get(e[1],0) for e,b in a.items()}
    cc=[]
    for a in aa:
        a=normalize(a,x); cc.append({e:a.get(e,0) for e in es})
    for b in bb:
        b=normalize(b,y); b={tuple(sorted(mapb[v] for v in e)):t for e,t in b.items()};cc.append({e:b.get(e,0) for e in es })
    return ts,cc

def rp3():
    # Barycentric subdivision of boundary of 4D cross-polytope, quotient antipodal.
    # vertices of quotient are orbits of nonempty faces (40); chains lift uniquely up to sign.
    def canon(f):
        f=tuple(sorted(f,key=abs));return f if f[0]>0 else tuple(-v for v in f)
    faces=sorted({canon(tuple(s*i for i,s in zip(ids,signs))) for k in range(1,5) for ids in combinations(range(1,5),k) for signs in product((-1,1),repeat=k)},key=lambda f:(len(f),f))
    ids={f:i for i,f in enumerate(faces)};ts=[];a={}
    for signs in product((-1,1),repeat=4):
        f=tuple(s*(i+1) for i,s in enumerate(signs))
        for p in permutations(f):
            chain=[tuple(sorted(p[:k],key=abs)) for k in range(1,5)]
            vs=[ids[canon(c)] for c in chain];ts.append(tuple(vs))
            sgn=[int(c!=canon(c)) for c in chain]
            for i,j in combinations(range(4),2):
                e=(vs[i],vs[j]);val=sgn[i]^sgn[j]
                if e in a: assert a[e]==val
                a[e]=val
    return ts,[a]

A=torus(); check('T3',*A,lambda i,j,k:int(len({i,j,k})==3),True)
B=s1xs2();check('S1xS2',*B,lambda i,j,k:0,True)
C=connected_sum(A,B);check('T3#S1xS2',*C,lambda i,j,k:int({i,j,k}=={0,1,2}),True)
D=rp3();check('RP3',*D,lambda i,j,k:1,True)
# explicit topological reference costs from all triple vectors

def det(a,b,c):
    return sum(((a>>i)&1)*((b>>j)&1)*((c>>k)&1) for i,j,k in permutations(range(3)))%2

def costs(n,Ws):
    ts=[t for t in product(range(1<<n),repeat=3) if det(*t)]
    return [min(sum(v not in w for v,w in zip(t,ws)) for t in ts) for ws in Ws]
print('torus common restriction costs',costs(3,[(set(range(1<<(3-h))),)*3 for h in range(4)]))
Ws=[({0,1},{0,2},{0,4}),({0,1},{0,2},{0,1}),({0,1},)*3,({0,8},)*3]
assert costs(4,Ws)==[0,1,2,3];print('connected sum matched rank costs',costs(4,Ws))
# Matter chain and reference branch corrections with fully specified phase polynomial
f=lambda a,b,c:a*b*c
for s,t,u,a,b,c in product(range(2),repeat=6):
    F=f(a,b,s)^f(t,u,c)^f(a^t,b^u,s^c)
    Q=a*u*(s^c)^t*b*(s^c)^t*u*s
    assert F==f(a,b,c)^Q
for s,d2,d3,a,b,c in product(range(2),repeat=6):
    F=f(a^s,b,c)^f(s,b^d2,c^d3)
    Q=s*d3*b^s*d2*c^s*d2*d3
    assert F==f(a,b,c)^Q
    # reset directions: x=(a+s,b,c), y=(s,b+d2,c+d3), token tau -> tau+a
    for tau in range(2):
        assert (a^s^s^tau,b^(b^d2),c^(c^d3))==(0^a^tau,0^d2,0^d3)
print('chain and reference branch identities: 64 assignments each; charge conservation 128')
F=lambda x:int(x==7)
for m,n,x,y in product(range(8),repeat=4):
    before=F(x)^F(x^m^y)^F(y^n)
    correction=F(x^m^y)^F(x^y)^F(y^n)^F(y)
    assert before^correction==F(x)^F(x^y)^F(y)
print('matter half-link gluing: 4096 assignments')
# Exact fourth Pauli moment via real phase Walsh sums, no logarithms.
def moment(f,n):
    d=1<<n; num=0; null=0; N=[(-1)**f(x) for x in range(d)]
    for a in range(d):
        for b in range(d):
            w=sum(N[x]*N[x^a]*(-1)**((x&b).bit_count()%2) for x in range(d))
            num+=w**4;null+=abs(w)==d
    return Fraction(num,d**5), n-(null.bit_length()-1)
M3,nu3=moment(F,3);M6,nu6=moment(lambda z:F(z&7)^F(z>>3)^F((z&7)^(z>>3)),6)
assert (M3,nu3,M6,nu6)==(Fraction(11,32),3,Fraction(79,1024),6)
print('CCZ and three-matter chain fourth moments/nullities',M3,nu3,M6,nu6)
# Graph cycle-basis filtration vs full Eulerian-support triple enumeration.
def all_cycles(es):
    cs=[]
    for b in range(1,1<<len(es)):
        bd=0
        for i,(u,v) in enumerate(es):
            if b>>i&1:bd^=(1<<u)^(1<<v)
        if bd==0:cs.append(b)
    return sorted(cs,key=lambda c:(c.bit_count(),c))

def check_graph(es,labels):
    cs=all_cycles(es); bs=[];r=0
    for c in cs:
        t=rank(bs+[c])
        if t>r:bs.append(c);r=t
    def label(mask):
        z=0
        for j,v in enumerate(labels):
            if mask>>j&1:z^=v
        return z
    def threshold(ls):
        best=None
        for a,b,c in combinations_with_replacement(ls,3):
            if label(a&b&c):
                k=max(a.bit_count(),b.bit_count(),c.bit_count());best=k if best is None else min(best,k)
        return best
    x,y=threshold(cs),threshold(bs);assert x==y
    return x
from itertools import combinations_with_replacement
count=0
for n in range(1,6):
    E=list(combinations(range(n),2))
    for b in range(1<<len(E)):
        es=[e for i,e in enumerate(E) if b>>i&1];check_graph(es,[1]*len(es));count+=1
rng=Random(70506)
for _ in range(160):
    es=[e for e in combinations(range(6),2) if rng.random()<.55]
    check_graph(es,[rng.randrange(4) for e in es])
print('cycle filtration: exhaustive homogeneous graphs through 5 vertices',count,'and 160 seeded heterogeneous graphs through 6')
for name,es in [('triangle',[(0,1),(1,2),(0,2)]),('square',[(0,1),(1,2),(2,3),(0,3)]),('theta',[(0,1),(0,2),(2,3),(3,1) ,(0,4),(4,5),(5,1)]),('C9',[(i,(i+1)%9) for i in range(9)])]:
    print(name,'threshold',check_graph(es,[1]*len(es)))
# Edge-supply gluing: all syndrome/coset mappings for triangle, 3 bits per edge.
es=[(0,1),(1,2),(2,0)]
def syndrome(x):
    a=[0]*3
    for v,(i,j) in zip(x,es):a[i]^=v;a[j]^=v
    return tuple(a)
cosets={}
for x in product(range(8),repeat=3):cosets.setdefault(syndrome(x),[]).append(x)
assert len(cosets)==64 and {len(v) for v in cosets.values()}=={8}
for s,xs in cosets.items():
    z=xs[0]
    for x in xs:
        y=tuple(a^b for a,b in zip(x,z));assert syndrome(y)==(0,0,0) and y[0]==y[1]==y[2]
        before=F(x[0])^F(x[1])^F(x[2]);after=F(y[0])^F(y[1])^F(y[2])
        correction=0
        for a,b in zip(y,z):correction^=F(a^b)^F(a)
        assert before^correction==after==F(y[0])
print('edge supply gluing: 64 equiprobable syndromes, 512 assignments, exact CCZ')
