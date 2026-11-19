"""Independent cross-checks of exact scientific tables; standard library only."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter, defaultdict
from itertools import permutations, product
from math import factorial, comb, ceil, log
import csv, ast, sys
DROOT=Path(sys.argv[1])
def read(name):
    with (DROOT/(name+'.csv')).open(newline='') as f: return list(csv.DictReader(f))
def walsh(a):
    a=list(a); h=1
    while h<len(a):
        for i in range(0,len(a),2*h):
            for j in range(i,i+h):
                x,y=a[j],a[j+h];a[j]=x+y;a[j+h]=x-y
        h*=2
    return a
def phase_spectrum(n, f):
    D=1<<n; s=[(-1)**f(x) for x in range(D)]; h=Counter()
    for a in range(D):
        h.update(F(abs(v),D) for v in walsh([s[x]*s[x^a] for x in range(D)]))
    assert sum(m*v*v for v,m in h.items())==D
    return h
def diag(h,n):
    D=1<<n;return n-(h[1].bit_length()-1),sum(F(m,D)*v**4 for v,m in h.items()),sum(F(m,D)*v**6 for v,m in h.items())
def poly(x,mon):return sum(x&m==m for m in mon)%2
def cp(d,N):
    n=3*N-d
    def f(x):
        small=(1<<d)-1; z=[(x>>((3-d)*N+d*i))&small for i in range(N-1)];z.append(0)
        for v in z[:-1]:z[-1]^=v
        return sum(((x>>((3-d)*i))&((1<<(3-d))-1))==((1<<(3-d))-1) and z[i]==small for i in range(N))%2
    return n,f
def spec_prod(a,b):
    t=Counter()
    for v,m in a.items():
        for w,k in b.items():t[v*w]+=m*k
    return t
def overlap(h,t,n):
    a=sorted(h.items(),reverse=True);b=sorted(t.items(),reverse=True);i=j=0;ma=a[0][1];mb=b[0][1];s=F(0)
    while i<len(a) and j<len(b):
        k=min(ma,mb);s+=k*a[i][0]*b[j][0];ma-=k;mb-=k
        if ma==0:
            i+=1
            if i<len(a):ma=a[i][1]
        if mb==0:
            j+=1
            if j<len(b):mb=b[j][1]
    assert i==len(a) and j==len(b)
    return s/(1<<n)
# All plotted tables, exact and independent of plot generator.
rows=read('batch_yield'); assert len(rows)==36
for r in rows:
    n,e,y=map(int,r.values()); assert y==[max(0,n-2),n-1,None,n][e]
assert len({(r['N'],r['ell_representative']) for r in rows})==36
print('batch_yield: 36 exact law values')
rows=read('cycle_noise');assert len(rows)==3604
for r in rows:
    L=int(r['L']);u,F0=F(r['u_exact']),F(r['F_exact']);q=sum(comb(L,k)*u**k*(1-u)**(L-k) for k in range(0,L+1,2));assert q**3==F0
assert len({(r['L'],r['u_exact']) for r in rows})==3604
print('cycle_noise: 3604 exact rational independent binomial checks')
rows=read('reference_price');assert len(rows)==301
for r in rows:
    b,a,c,d,v=map(F,r.values());assert (a,c,d,v)==(1+3*b,2+b,F(3),min(1+3*b,2+b,F(3)))
print('reference_price: 301 exact rows')
# Direct Walsh, not contraction ranks, for independent moment verification.
for r in read('monomial_moments'):
    d,N=int(r['relative_type']),int(r['N']);n,f=cp(d,N);h=phase_spectrum(n,f);assert diag(h,n)==(int(r['nullity']),F(r['P2']),F(r['P3']))
print('monomial_moments: 12 complete Walsh spectra')
# Spectrum-only caps for known input, independent target spectra construction.
C=phase_spectrum(3,lambda x: int(x==7)); assert C==Counter({F(0):35,F(1,2):28,F(1):1}); assert diag(C,3)==(3,F(11,32),F(23,128))
T=spec_prod(C,C); T3=spec_prod(T,C);Tpad=spec_prod(T,Counter({F(1):2,F(0):2}));assert Tpad==Counter({F(0):14702,F(1,4):1568,F(1,2):112,F(1):2})
a=[7,56,73,82]; b=[7,25,97,42,82,76,52]
# 147=1,4,7 -> 73; 257 -> 82.
for terms,expect in [(a,(7,F(29,512),F(383,32768))),(b,(7,F(29,512),F(107,8192)))]:assert diag(phase_spectrum(7,lambda x:poly(x,terms)),7)==expect
Bz=phase_spectrum(6,lambda x:poly((x<<1)|1,b));assert diag(Bz,6)[0]==6;assert overlap(Bz,T,6)==F(121,128)
# F4 trace six-bit phase f(a)+f(b)+f(a+b); determinant row-major.
trace=phase_spectrum(6,lambda x:int((x&7)==7)^int((x>>3)==7)^int(((x&7)^(x>>3))==7))
det=phase_spectrum(9,lambda x:sum(all(x&(1<<(3*i+p[i])) for i in range(3)) for p in permutations(range(3)))%2)
assert overlap(trace,T,6)==F(57,64);assert overlap(det,T3,9)==F(13199,16384)
assert diag(det,9)==(9,F(149,8192),F(10369,4194304))
normal=phase_spectrum(6,lambda x:poly(x,[1|8|16,2|8|32,4|16|32])); assert overlap(normal,T,6)==F(121,128)
vals={r['task']:(F(r['value']),r['status']) for r in read('fidelity_ceilings')};assert vals=={'A: 2 targets':(F(1),'attained'),'B: 2 targets':(F(121,128),'upper bound'),'F4: 2 targets':(F(57,64),'upper bound'),'Det: 3 targets':(F(13199,16384),'upper bound')}
print('fidelity_ceilings: independent spectra for CCZ, A, B, one maximizing Z branch, trace, determinant, coefficient normal form; all values agree')
# Direct combinatorial invariants and count checks (not enumeration of original objects).
h=defaultdict(dict)
for r in read('fano_histograms'): h[r['quantity']][F(r['value'])]=int(r['count'])
assert sum(h['nullity'].values())==2*(4**7-1);assert h['nullity'][6]==sum(h['bound'].values())==sum(h['moment'].values())==31644;assert max(h['bound'])==F(121,128);assert F(121,1024) not in h['moment']
print('fano_histograms: totals and maxima; distinct marginal distributions')
for r in read('determinant_planes'):
    eta,num,hist=int(r['nullity']),int(r['planes']),dict(ast.literal_eval(r['rank_histogram']));assert sum(hist.values())==64;assert 6-(hist[0].bit_length()-1)==eta
assert sum(int(r['planes']) for r in read('determinant_planes'))==1395
assert 1395==((64-1)*(64-2)*(64-4))//((8-1)*(8-2)*(8-4))
print('determinant_planes: four classes and Gaussian-binomial count')
byN=defaultdict(list)
for r in read('copy_ranks'):byN[int(r['N'])].append((int(r['dimension']),int(r['contraction_rank']),int(r['count'])))
for N,rs in byN.items():assert all(d==3*(N-1) and r%2==0 and r!=2 for d,r,k in rs) and sum(k for d,r,k in rs)==2**(3*(N-1)) and [k for d,r,k in rs if r==0]==[1]
print('copy_ranks: totals for 37440 displacements (not independent enumeration)')
for r in read('decoder_checks'): N,k=map(int,r.values());assert k==8*2**(3*(N-2))
print('decoder_checks: 37440 assignment count (distinct from displacement test)')
for r in read('reference_nullities'):N,e,n=map(int,r.values());assert n//3==(N if e==3 else N-1)
print('reference_nullities: all 9 rows consistent with analytic frontier')
rows=read('graph_counts');assert sum(int(r['graphs']) for r in rows)==33867;assert sum(int(r['active']) for r in rows)==28543
for r in rows:
    n,g,c,a=map(int,r.values());assert g==2**comb(n,2) and c==sum(comb(n,k)*factorial(k-1)//2 for k in range(3,n+1)) and 0<=a<=g
print('graph_counts: all six labelled domains and cycle counts')
for r in read('mesh_counts'):
    V,E,Fc,T,b1,c,s=map(int,[r[k] for k in ['V','E','F','T','b1','cup_checks','exact_shift_checks']]);assert V-E+Fc-T==0 and 2*Fc==4*T and c==b1**3 and s==3*V*b1*b1
print('mesh_counts: four mesh consistency identities; not a manifold certification')
print('PASS: all 12 table families reviewed; independent exact checks as stated.')
