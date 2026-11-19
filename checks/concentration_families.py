#!/usr/bin/env python3
"""Exact finite checks for cubic restriction, independent of historical helpers.
No circuit-optimality claims are inferred from these finite domains.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from functools import reduce
from operator import xor

def rank(rows):
    piv={}
    for row in rows:
        while row:
            k=row.bit_length()-1
            if k in piv: row^=piv[k]
            else: piv[k]=row; break
    return len(piv)

def basis(rows):
    b=[]
    for x in rows:
        if rank(b+[x])>len(b):b.append(x)
    return b

def add(xs):return reduce(xor,xs,0)
def vol(a,b,c):
    return int(rank([a,b,c])==3)
def evalf(x,ts):return sum((x&t)==t for t in ts)&1

def form(x,y,z,ts):
    # Polarization by eight evaluations, including f(0)=0.
    return add(evalf(t,ts) for t in (0,x,y,z,x^y,x^z,y^z,x^y^z))

def hist(m,fn):
    # First extract coefficients in coordinate basis; construct linear pencil.
    E=[1<<i for i in range(m)]
    mats=[[sum(fn(a,b,c)<<k for k,c in enumerate(E)) for b in E] for a in E]
    out=Counter()
    for v in range(1<<m):
        inds=[i for i in range(m) if v>>i&1]
        out[rank(add(mats[i][j] for i in inds) for j in range(m))]+=1
    return dict(sorted(out.items()))

def moment(m,h,j): return sum((F(c,1<<(r*(j-1))) for r,c in h.items()),F(0))/(1<<m)

def truth_cubic(m,fn):
    """ANF via Mobius transform, exact cubic/higher terms."""
    tab=[fn(x) for x in range(1<<m)]
    for i in range(m):
        for x in range(1<<m):
            if x>>i&1:tab[x]^=tab[x^(1<<i)]
    return [x for x,c in enumerate(tab) if c and x.bit_count()>=3]

A=(7,56,73,82); B=(7,25,97,42,82,76,52)
for name,terms,h0,m3 in [('A',A,{0:1,2:3,4:76,6:48},F(383,32768)),('B',B,{0:1,2:7,4:56,6:64},F(107,8192))]:
    fn=lambda x,y,z:form(x,y,z,terms)
    h=hist(7,fn);assert h==h0 and moment(7,h,2)==F(29,512) and moment(7,h,3)==m3
    planes=[[2,4,8,32,64]] if name=='A' else [[3,9,20,36,65],[3,5,24,40,72]]
    for p in planes:
        assert rank(p)==5
        vals=[fn(*tr) for tr in combinations(p,3)]
        print(name,'five-plane',p,'values',vals);assert not any(vals)
    for l in range(1,128):
        H=basis(x for x in range(128) if (x&l).bit_count()%2==0)
        assert len(H)==6 and any(fn(*abc) for abc in combinations(H,3))
    print(name,'128 contractions, 127 hyperplanes excluded',h,'P2',moment(7,h,2),'P3',m3)
for o in (0,1):
    assert truth_cubic(6,lambda x:evalf(x|(o<<6),A))==[7,56]
print('A two x7 decoder outcomes verified')

# All-copy no-reference contraction spectra (coordinates are first N-1 slots).
contractions=0
for N in range(2,6):
    E=3*(N-1)
    def total(v):return add((v>>(3*i))&7 for i in range(N-1))
    def fn(x,y,z): return (sum(vol((x>>(3*i))&7,(y>>(3*i))&7,(z>>(3*i))&7) for i in range(N-1))+vol(total(x),total(y),total(z)))&1
    h=hist(E,fn); assert h=={0:8} if N==2 else h[0]==1 and min(r for r in h if r)>=(4)
    contractions+=1<<E
    print('unassisted N',N,'histogram',h)
print('total unassisted displacements',contractions)

# Decoder lower bounds: prescribe all measured bits and enumerate output cubic ANF.
checks=0
for N in range(2,5):
    d=3*(N-2)
    for a,b,c in product(range(2),repeat=3):
        def f(v):
            xs=[(v>>(3*i))&7 for i in range(N-2)]; s=add(xs)
            pen=(a|(b<<1)|((s&4)^(c<<2)))
            end=s^pen
            assert end>>2==c
            return sum(x==7 for x in xs+[pen,end])&1
        assert truth_cubic(d,f)==[7<<(3*i) for i in range(N-2)];checks+=1
    d=3*(N-1)
    for c in range(2):
        def f(v):
            xs=[(v>>(3*i))&7 for i in range(N-1)];end=(add(xs)&3)|(c<<2)
            return sum(x==7 for x in xs+[end])&1
        assert truth_cubic(d,f)==[7<<(3*i) for i in range(N-1)];checks+=1
print('decoder outcome ANFs',checks)
# Two-copy two-block orientation: x=(x1,...,x6), t in R.
for R,h0 in [((1,8),6),((1,2),5)]:
    def y(v):return (v&63)^add(r for i,r in enumerate(R) if v>>(6+i)&1)
    def fn(a,b,c):return form(a&63,b&63,c&63,(7,56))^form(y(a),y(b),y(c),(7,56))
    h=hist(8,fn);eta=8-(h[0].bit_length()-1);assert eta==h0;print('orientation',R,'nullity',eta,'histogram',h)

# Trace phase, s=2,3, over explicitly specified irreducible polynomials.
for s,poly in [(2,0b111),(3,0b1011)]:
    mask=(1<<s)-1
    def mul(a,b):
        acc=0
        while b:
            if b&1:acc^=a
            b>>=1;a<<=1
            if a>>s:a^=poly
        return acc
    def tr(a):
        ans=0
        for _ in range(s):ans^=a;a=mul(a,a)
        assert ans in (0,1);return ans
    def f(v):return tr(mul(mul(v&mask,(v>>s)&mask),(v>>(2*s))&mask))
    def fn(a,b,c):return add(f(v) for v in (0,a,b,c,a^b,a^c,b^c,a^b^c))
    h=hist(3*s,fn); assert h=={0:1,2*s:(1<<(3*s))-1}
    print('trace s',s,'histogram',h,'P2',moment(3*s,h,2),'P3',moment(3*s,h,3))
    if s==2:
        for x,y in product(range(8),repeat=2):
            a,b,c=[((x>>i)&1)|(((y>>i)&1)<<1) for i in range(3)]
            assert tr(mul(mul(a,b),c))==((x==7)^(y==7)^((x^y)==7))
        print('64 GF4 trace identities')
# Mixed monomial kernels type d=1,2,3; types counted on active cubic slots only.
for d in (1,2):
    for N in range(1,5):
        coords=[1<<(3*i+j) for i in range(N) for j in range(3-d)]
        coords += [(1<<(3*i+j))^(1<<(3*(N-1)+j)) for i in range(N-1) for j in range(3-d,3)]
        def emb(x):return add(v for i,v in enumerate(coords) if x>>i&1)
        terms=[7<<(3*i) for i in range(N)]
        def fn(x,y,z):return form(emb(x),emb(y),emb(z),terms)
        h=hist(len(coords),fn);eta=len(coords)-(h[0].bit_length()-1)
        expected=0 if N==1 else 3 if d==2 and N==2 else 3*N-d
        assert eta==expected
        print('mixed type',d,'N',N,'nullity',eta)
print('All exact finite checks passed.')
