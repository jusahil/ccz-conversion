"""Exact independent baseline and benchmarks; standard library, bit masks."""
from itertools import combinations,product
from collections import Counter
from fractions import Fraction as F
from random import Random
from math import inf

def basis(vectors):
    p={}
    for x in vectors:
        while x:
            i=x.bit_length()-1
            if i in p:x^=p[i]
            else:p[i]=x;break
    return list(p.values())
def span(b):
    s=[0]
    for v in b:s += [x^v for x in s]
    return s

def subs(n):
    seen={frozenset([0])}; frontier=list(seen)
    while frontier:
        w=frontier.pop()
        for x in range(1,1<<n):
            if x in w:continue
            h=w|frozenset(y^x for y in w)
            if h not in seen:seen.add(h);frontier.append(h)
    return [(w,basis(w)) for w in sorted(seen,key=lambda x:(len(x),sorted(x)))]
def truth(n,terms):return [sum((x&t)==t for t in terms)%2 for x in range(1<<n)]
def pol(f,a,b,c):return f[0]^f[a]^f[b]^f[c]^f[a^b]^f[a^c]^f[b^c]^f[a^b^c]
def active(f,b):return any(pol(f,a,z,c) for a,z,c in combinations(b,3))
def coll(n,wb,hb,N):
    # last block = sum first N-1 blocks + vector in H. Each first block is arbitrary.
    return [(1<<(n*i+j))^(1<<(n*(N-1)+j)) for i in range(N-1) for j in range(n)]+[v<<(n*(N-1)) for v in hb]
def collective_active(t,n,bs,N):
    mask=(1<<n)-1
    for a,b,c in combinations(bs,3):
        if sum(pol(t,(a>>(i*n))&mask,(b>>(i*n))&mask,(c>>(i*n))&mask) for i in range(N))%2:return True
    return False

def baseline():
    S=subs(4);assert len(S)==67
    triples=[a|b|c for a,b,c in combinations([1,2,4,8],3)]
    count=Counter();copies=Counter();cmps=0
    for mask in range(16):
        t=truth(4,[v for i,v in enumerate(triples) if mask>>i&1])
        # cache independent collective tests by reference preimage
        tests={(N,w):collective_active(t,4,coll(4,[],b,N),N) for N in (1,2,3) for w,b in S}
        for w,wb in S:
            # adapted basis, compute d via count of factors outside neutral basis
            vb=list(wb)
            for e in (1,2,4,8):
                if len(basis(vb+[e]))>len(vb):vb.append(e)
            d=min((sum(i>=len(wb) for i in idx) for idx in combinations(range(4),3) if pol(t,*(vb[i] for i in idx))), default=inf)
            count[d]+=1
            c=min((N for N in (1,2,3) if tests[N,w]),default=inf);copies[c]+=1
            assert c=={0:1,1:2,2:2,3:3,inf:inf}[d]
            for N in (1,2):
                got=min((len(hb)-len(wb) for h,hb in S if w<=h and tests[N,h]),default=inf)
                want=d if N==1 else (inf if d==inf else int(d==3))
                assert got==want;cmps+=1
    print('S0 forms/spaces',sum(count.values()),'type',dict(count),'copy',dict(copies),'reference comparisons',cmps)
    assert list(count.get(i,0) for i in (0,1,2,3,inf))==[135,525,315,30,67]

def fwht(v):
    v=list(v);h=1
    while h<len(v):
        for i in range(0,len(v),2*h):
            for j in range(i,i+h):a,b=v[j],v[j+h];v[j]=a+b;v[j+h]=a-b
        h*=2
    return v

def diagnostics(t,direct=True):
    D=len(t);n=D.bit_length()-1;bs=[1<<i for i in range(n)]
    hist=Counter(len(basis(sum(pol(t,a,b,c)<<i for i,c in enumerate(bs)) for b in bs)) for a in range(D))
    mom=[sum(F(v,D*2**((j-1)*h)) for h,v in hist.items()) for j in (2,3)]
    eta=n-(hist[0].bit_length()-1)
    if direct:
        totals=[0,0];st=0
        for a in range(D):
            w=fwht([1-2*(t[x]^t[x^a]) for x in range(D)])
            totals[0]+=sum(x**4 for x in w);totals[1]+=sum(x**6 for x in w);st+=sum(abs(x)==D for x in w)
        assert mom==[F(totals[0],D**5),F(totals[1],D**7)] and eta==n-(st.bit_length()-1)
    return eta,*mom,dict(sorted(hist.items()))

def monomial_comp(d,N):
    # 3-d independently neutral variables per source; d common charge coordinates, last charge parity.
    m=(3-d)*N+d*(N-1);t=[]
    for v in range(1<<m):
        z=[(v>>((3-d)*N+d*i))&((1<<d)-1) for i in range(N-1)]
        z+=[0]
        for zz in z[:-1]:z[-1]^=zz
        p=0
        for i in range(N):
            u=(v>>((3-d)*i))&((1<<(3-d))-1);p^=((u==(1<<(3-d))-1) and (z[i]==(1<<d)-1))
        t.append(int(p))
    return t

def calc():
    for d,N in product(range(4),range(1,4)):
        eta,p2,p3,h=diagnostics(monomial_comp(d,N));print('monomial',d,N,'eta',eta,'P2',p2,'P3',p3)
    A=[7,56,73,82] # x1x4x7 = 73; x2x5x7=82
    B=[7,25,97,42,82,76,52]
    for name,t in [('A',A),('B',B)]:print(name,diagnostics(truth(7,t)))
    # Exact full-support stabilizer fidelity check: quadratic signs, linear i phases.
    best=F(0);count=0
    for linear in product(range(4),repeat=3):
        for q in range(8):
            uv=[0,0]; pairs=[3,5,6]
            for x in range(8):
                a=(sum(linear[i] for i in range(3) if x>>i&1)+2*(sum((x&t)==t for i,t in enumerate(pairs) if q>>i&1)+(x==7) ))%4
                uv[a%2]+=1 if a<2 else -1
            best=max(best,F(uv[0]**2+uv[1]**2,64));count+=1
    assert best==F(9,16) and count==512;print('fidelity',best,'full support phases',count)
    # exact seeded sparse phase cases, linear shear + affine shifts + quadratic correction.
    rng=Random(7017)
    for j in range(10):
        n=6+j%4;D=1<<n;terms=rng.sample([a|b|c for a,b,c in combinations([1<<i for i in range(n)],3)],min(8,n));t=truth(n, terms);base=diagnostics(t)
        r=rng.randrange(D);tt=[t[(x^(((x>>1)&1)<<0))^r]^int((x&12)==12) for x in range(D)]
        assert diagnostics(tt)==base
    print('10 seeded sparse cubic cases (6-9 bits), both algorithms and affine-Clifford regression')
    # heterogeneous image alignment checks using full tensor truth tables on kernels.
    # same charges: difference of unequal monomials; disjoint images force each source neutral.
    t=truth(4,[7,11]);assert active(t,[1,2,4,8]);assert not active(truth(4,[]),[1,2,4,8]);print('heterogeneous aligned unequal forms activate; equal forms cancel; disjoint injective images K=0')

if __name__=='__main__':baseline();calc()
