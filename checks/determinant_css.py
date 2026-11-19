from itertools import combinations, permutations, product
from collections import Counter
from fractions import Fraction as F

def rank(xs):
    piv={}
    for x in xs:
        while x:
            i=x.bit_length()-1
            if i in piv:x^=piv[i]
            else:piv[i]=x;break
    return len(piv)
def span(bs):
    out=[0]
    for b in bs:out += [x^b for x in out]
    return out
def bit(x,i):return (x>>i)&1
def cubic(x,terms):return sum((x&m)==m for m in terms)%2
def polar(a,b,c,terms):return sum(cubic(t,terms) for t in (0,a,b,c,a^b,a^c,b^c,a^b^c))%2
def contractions(n,terms):
    # For one monomial xyz, contraction at a gives the alternating 2-form a_x yz + a_y xz + a_z xy.
    basis=[]
    for k in range(n):
        rows=[0]*n
        for m in terms:
            if bit(m,k):
                i,j=[i for i in range(n) if i!=k and bit(m,i)]
                rows[i]^=1<<j;rows[j]^=1<<i
        basis.append(rows)
    out=Counter()
    for a in range(1<<n):
        rows=[0]*n
        for k in range(n):
            if bit(a,k):rows=[r^s for r,s in zip(rows,basis[k])]
        out[rank(rows)]+=1
    return out
D=[sum(1<<(3*i+p[i]) for i in range(3)) for p in permutations(range(3))]
H=contractions(9,D); assert H=={0:1,4:49,6:294,8:168}
for a in range(512):
    # Independently calculate individual derivative matrices from evaluations, not from monomial structure.
    rows=[sum(polar(a,1<<i,1<<j,D)<<j for j in range(9)) for i in range(9)]
    assert rank(rows)==[0,4,6,8][rank([(a>>(3*i))&7 for i in range(3)])]
print('Determinant displacement/derivative checks: 512;',dict(sorted(H.items())))
assert sum(F(v,2**(9+r)) for r,v in H.items())==F(149,8192)
assert sum(F(v,2**(9+2*r)) for r,v in H.items())==F(10369,4194304)
def spectrum(n,H):
    return sorted([F(1,2**(r//2)) for r,c in H.items() for _ in range(c*2**r)],reverse=True)
def ceiling(n,A,B):return sum(x*y for x,y in zip(A,B))/2**n
ccz=[7]
T3=contractions(9,[7,7<<3,7<<6]);T2=contractions(6,[7,7<<3]);
assert ceiling(9,spectrum(9,H),spectrum(9,T3))==F(13199,16384)
print('Three-target determinant ceiling:',ceiling(9,spectrum(9,H),spectrum(9,T3)))
# all coefficient matrices, rows neutral coordinates, columns t2t3,t1t3,t1t2
coeffs={}; patterns={0:{0:64},1:{0:8,2:56},2:{0:2,2:30,4:32},3:{0:1,2:7,4:56}}
for L in range(512):
    rows=tuple((L>>(3*i))&7 for i in range(3));r=rank(rows)
    terms=[(1<<i)|((7^(1<<j))<<3) for i in range(3) for j in range(3) if bit(rows[i],j)]
    Hc=contractions(6,terms);assert Hc==patterns[r];coeffs[rows]=r
assert ceiling(6,spectrum(6,patterns[3]),spectrum(6,T2))==F(121,128)
print('All coefficient matrices: 512, 32768 contractions; rank classes',dict(sorted(Counter(coeffs.values()).items())))
# Deduplicate all independent bases via their sorted span; use a deterministic greedy basis per span.
planes=set()
for b in combinations(range(1,64),3):
    if rank(b)==3:planes.add(tuple(sorted(span(b))))
assert len(planes)==1395
hist=Counter()
for P in sorted(planes):
    b=[]
    for x in P:
        if rank(b+[x])>len(b):b.append(x)
    # neutral first column entries at 0,3,6; charged positions 1,2,4,5,7,8
    lift=lambda x:sum(bit(x,i)<<j for i,j in enumerate([1,2,4,5,7,8]))
    b=list(map(lift,b)); rows=[]
    for i in (0,3,6):
        rows.append(sum(polar(1<<i,b[(j+1)%3],b[(j+2)%3],D)<<j for j in range(3)))
    hist[rank(rows)]+=1
assert hist=={0:3,1:168,2:840,3:384}
print('Distinct 3-planes in F2^6: 1395; coefficient-rank counts',dict(sorted(hist.items())))
# Decoder: all outcomes of measuring omitted matching, direct phase residual degree <= 2.
keep=[0,4,8,1,5,6];omit=[2,3,7]
for s in range(8):
    vals=[]
    for y in range(64):
        x=sum(bit(y,i)<<j for i,j in enumerate(keep))|sum(bit(s,i)<<j for i,j in enumerate(omit))
        vals.append(cubic(x,D)^cubic(y,[7,56]))
    a=vals[:]
    for i in range(6):
        for x in range(64):
            if bit(x,i):a[x]^=a[x^(1<<i)]
    assert not any(a[x] for x in range(64) if x.bit_count()>2)
print('Determinant decoder: 8 branches, 512 labels')
# Relative cup examples: evaluate every ordered triple in the three fields of H1(T3#(S1xS2)).
def detcols(a,b,c):return int(rank([a&7,b&7,c&7])==3)
for ws,expected in [((1,2,4),0),((1,2,1),1),((1,1,1),2),((8,8,8),3)]:
    cost=min(sum(v not in (0,w) for v,w in zip(vs,ws)) for vs in product(range(16),repeat=3) if detcols(*vs))
    assert cost==expected
print('Relative kernel examples: 4 x 4096 ordered field triples, costs 0,1,2,3')
# Chain and token phases and full coherent reset coordinate maps.
for a,b,c,s,t,u in product(range(2),repeat=6):
    x=(a,b,s);y=(t,u,c); z=tuple(i^j for i,j in zip(x,y))
    f=lambda v:v[0]&v[1]&v[2]
    assert f(x)^f(z)^f(y)==(a*b*c ^ a*u*(s^c) ^ t*b*(s^c) ^ t*u*s)
for a,b,c,s,d2,d3,theta in product(range(2),repeat=7):
    x=(a^s,b,c);y=(s,b^d2,c^d3)
    assert f(x)^f(y)==(a*b*c ^ s*d3*b ^ s*d2*c ^ s*d2*d3)
    assert (x[0]^y[0]^theta)==((x[0]^a)^y[0]^(theta^a))
    assert x[1]^y[1]==(x[1]^b)^(y[1]^b) and x[2]^y[2]==(x[2]^c)^(y[2]^c)
print('Chain 64 phase identities; token 128 charge/phase identities (not full native circuit simulation)')
# Cubical CSS edge support: explicit periodic edge set, plaquettes and coordinate loops.
faces=loops=0
for L in range(3,9):
    V=list(product(range(L),repeat=3))
    def step(v,i):return tuple((x+(j==i))%L for j,x in enumerate(v))
    def edge(v,i):return (v,i)
    for k in range(3):
        membrane={edge(v,k) for v in V if v[k]==L-1};assert len(membrane)==L*L
        for v in V:
            for i,j in combinations(range(3),2):
                face={edge(v,i),edge(step(v,i),j),edge(step(v,j),i),edge(v,j)}
                assert len(face&membrane)%2==0;faces+=1
        for i in range(3):
            seen=set()
            for v in V:
                if v[i]:continue
                loop={edge(tuple((x+t*(j==i))%L for j,x in enumerate(v)),i) for t in range(L)}
                assert len(loop)==L and not (loop&seen)
                seen|=loop;assert len(loop&membrane)==int(i==k);loops+=1
print('Periodic CSS, L=3..8: plaquette parities',faces,'; loop-membrane pairings',loops)
print('PASS: all assertions. Analytic converses not replaced by finite enumeration.')
