"""Sparse exact amplitude tests. Amplitudes are integer signs; all nonzero branch
vectors have equal modulus and are normalized by their support size. Coherence is
preserved: these are basis-permutation and diagonal-phase simulations, not just
checks of probability distributions. All arithmetic exact."""
from itertools import product,permutations
from collections import defaultdict,Counter
from fractions import Fraction as F
from algebra import diagnostics,truth,pol,basis

def f(x):return int(x==7)
def phase(x):return 1-2*x

def put(d,k,v):
    assert k not in d;d[k]=v

def norm2(d):return sum(v*v for v in d.values())
def chain_check(z,l,es):
    q=[0]*len(z)
    for x,(u,v) in zip(l,es):q[u]^=x;q[v]^=x
    assert tuple(q)==z

def chain():
    es=[(0,1),(1,2)];branches=defaultdict(dict)
    for x,y in product(range(8),repeat=2):
        z=(x,x^y,y);chain_check(z,(x,y),es)
        put(branches[(x>>2,y&3)],(z,(x,y),0),phase(f(x)^f(x^y)^f(y)))
    assert len(branches)==8 and all(norm2(d)==8 for d in branches.values())
    for (s,tu),d in branches.items():
        t=tu&1;u=tu>>1;out={}
        for (z,(x,y),zero),v in d.items():
            a=x&1;b=(x>>1)&1;c=y>>2;k=a|b<<1|c<<2
            # explicit logical shifts on edge & both adjacent matter labels
            zz=list(z);ll=[x,y]
            for bit,i,r in [(a,0,1),(b,0,2),(c,1,4)]:
                if bit:
                    ll[i]^=r
                    for j in es[i]:zz[j]^=r
                chain_check(tuple(zz),tuple(ll),es)
            Q=(a*u*s)^(a*u*c)^(t*b*s)^(t*b*c)^(t*u*s) # cubic expansion: remaining x=(a,b,s), y=(t,u,c)
            put(out,(tuple(zz),tuple(ll),k),v*phase(Q))
        assert len({(z,l) for z,l,k in out})==1
        assert {k:v for (z,l,k),v in out.items()}=={k:phase(f(k)) for k in range(8)}
    print('matter-link chain: 8 branches, 64 coherent amplitudes, all intermediate dressed shifts Gauss-valid; exact target')

def gluing():
    # half links at v0,v1,v1,v2, local Gauss labels z0,z1,z2. No quotient coherence destroyed.
    branches=defaultdict(dict)
    for a,b,c,d in product(range(8),repeat=4):
        z=(a,b^c,d);h=(a,b,c,d);v=phase(f(z[0])^f(z[1])^f(z[2]));put(branches[(a^b,c^d)],(z,h),v)
    assert len(branches)==64 and all(norm2(v)==64 for v in branches.values())
    for (p,q),br in branches.items():
        out={}
        for (z,h),w in br.items():
            a,b,c,d=h;z0,z1,z2=z;b^=p;z1^=p;assert (a,b^c,d)==(z0,z1,z2)
            d^=q;z2^=q;assert (a,b^c,d)==(z0,z1,z2) and a==b and c==d
            Q=f(z1^p)^f(z1)^f(z2^q)^f(z2)
            put(out,((z0,z1,z2),(a,b,c,d)),w*phase(Q))
        assert out=={((x,x^y,y),(x,x,y,y)):phase(f(x)^f(x^y)^f(y)) for x,y in product(range(8),repeat=2)}
    print('half-link gluing: 64 full conditional amplitude maps, 4096 amplitudes; duplicates retained; all old Gauss laws satisfied')

def ref():
    d=defaultdict(dict)
    for x,y,t in product(range(8),range(8),range(2)):
        r=((x^y)>>1,y&1)
        put(d[r],(x,y,t),phase(f(x)^f(y)))
    assert len(d)==8 and all(norm2(v)==16 for v in d.values())
    # encoded physical map triple repetitions; each logical coordinate flip is X on all three qubits.
    enc=lambda x:sum(((x>>j&1)*7)<<(3*j) for j in range(3))
    def valid(v,n):return all(((v>>(3*j))&7) in (0,7) for j in range(n))
    for (delta,s),br in d.items():
        d2=delta&1;d3=delta>>1;out={};outminus={}
        for (x,y,t),v in br.items():
            aa=(x^y)&1;b=x>>1&1;c=x>>2&1;k=aa|b<<1|c<<2;p=x^y^t
            X,Y,T=enc(x),enc(y),7*t
            for (bit,xm,ym,tm) in [(aa,1,0,1),(b,2,2,0),(c,4,4,0)]:
                if bit:X^=enc(xm);Y^=enc(ym);T^=7*tm;x^=xm;y^=ym;t^=tm
                assert valid(X,3) and valid(Y,3) and valid(T,1) and x^y^t==p
            Q=s*d3*b^s*d2*c^s*d2*d3
            put(out,(X,Y,T,k),v*phase(Q));put(outminus,(X,Y,T,k),v*phase(Q)*phase(t^aa))
        assert len({(X,Y) for X,Y,T,k in out})==1
        assert {(T,k):v for (X,Y,T,k),v in out.items()}=={(7*t,k):phase(f(k)) for t,k in product(range(2),range(8))}
        assert {(T,k):v for (X,Y,T,k),v in outminus.items()}=={(7*t,k):phase(f(k)^t^(k&1)) for t,k in product(range(2), range(8))}
    print('dressed reference: 8 branches, 128 amplitudes each for +/- preparation; full token factorizes, minus branch adds neutral Z')

def graphcase(name,es,cycles,available,expected):
    # all coordinate directions on selected cycle supports; parity readout uses dual physical edges.
    h=len(cycles);m=len(es);n=3;phase_fn=lambda xs:sum(f(x) for x in xs)%2
    dual=[]
    for i in range(h):
        lam=next(b for b in range(1<<m) if all(((b&c).bit_count()%2)==(j==i) for j,c in enumerate(cycles)));dual.append(lam )
    for e in cycles:
        s=Counter(v for i,(a,b) in enumerate(es) if e>>i&1 for v in (a,b));assert all(c%2==0 for c in s.values())
    def enc(y):return tuple(__import__('functools').reduce(int.__xor__,(y[i] for i,c in enumerate(cycles) if c>>j&1),0) for j in range(m))
    def check(x):
        z=defaultdict(int)
        for a,(u,v) in zip(x,es):z[u]^=a;z[v]^=a
        assert not any(z.values())
    br=defaultdict(dict)
    for y in product(range(8),repeat=h):
        x=enc(y);check(x);q=y[available:];put(br[q],x,phase(phase_fn(x)))
    assert len(br)==8**(h-available)
    for q,d in br.items():
        anchor=enc((0,)*available+q);out={}
        for x,w in d.items():
            # coherent dual parity accumulation to neutral output (no per-edge measurements)
            y=tuple(__import__('functools').reduce(int.__xor__,(x[j] for j in range(m) if lam>>j&1),0) for lam in dual[: available]);xx=list(x)
            for a,c in zip(y,cycles[:available]):
                for bit in (1,2,4):
                    if a&bit:
                        for j in range(m):
                            if c>>j&1:xx[j]^=bit
                        check(tuple(xx))
            assert tuple(xx)==anchor
            p=phase_fn(enc(y+(0,)*(h-available)));diff=phase_fn(x)^p
            put(out,y,w*phase(diff))
        assert out=={y:phase(phase_fn(enc(y+(0,)*(h-available)))) for y in product(range(8),repeat=available)}
        # Check the correction has zero third derivative.
        vals=[phase_fn(enc(tuple((x>>(3*i))&7 for i in range(available))+q))^phase_fn(enc(tuple((x>>(3*i))&7 for i in range (available))+(0,)*(h-available))) for x in range(1<<(3*available))]
        assert not any(pol(vals,a,b,c) for a,b,c in __import__('itertools').combinations([1<<i for i in range(3*available) ],3))
    t=[phase_fn(enc(tuple((x>>(3*i))&7 for i in range(available))+(0,)*(h-available))) for x in range(1<<(3*available))]
    assert diagnostics(t)[:2]==expected;print(name,'quotient branches',len(br),'amplitudes',8**h,'accessible eta/P2', expected)

def det():
    ts=[sum(1<<(3*i+p[i]) for i in range(3)) for p in permutations(range(3))];t=truth(9,ts);r=diagnostics(t);assert r [:3]==(9,F(149,8192),F(10369,4194304));print('determinant',r)
    keep=[0,4,8,1,5,6];removed=[2,3,7]
    for s in range(8):
        y=[sum(((x>>i)&1)<<v for i,v in enumerate(keep))+sum(((s>>i)&1)<<v for i,v in enumerate(removed)) for x in range (64)]
        correction=[t[a]^int((i&7)==7)^int((i&56)==56) for i,a in enumerate(y)]
        from itertools import combinations
        assert not any(pol(correction,a,b,c) for a,b,c in combinations([1<<j for j in range(6)],3))
    print('determinant 8 affine branches extract 2CCZ; nullity saturation excludes 3')
    # 3-torus coordinate Wilson/membrane dual and face-flatness for L=3..8.
    for L in range(3,9):
        vs=list(product(range(L),repeat=3));sh=lambda v,i:tuple((a+(j==i))%L for j,a in enumerate(v))
        mem=[{(v,j) for v in vs if v[j]==L-1} for j in range(3)]
        for a in mem:
            assert len(a)==L*L
            for v in vs:
                for i,j in combinations(range(3),2):assert sum(e in a for e in [(v,i),(sh(v,i),j),(sh(v,j),i),(v,j)])%2==0
        for i in range(3):
            loop={(tuple(k if j==i else 0 for j in range(3)),i) for k in range(L)}
            assert [len(loop&a)%2 for a in mem]==[int(i==j) for j in range(3)]
    print('cubical T3: six sizes L=3..8, membrane flatness and coordinate loop duality')

def noise():
    # exact distribution of iid Z errors on an odd L cycle; event output Z parity zero.
    from math import comb
    for L in (3,5,9,21):
        p=F(1,100);ok=sum(comb(L,i)*p**i*(1-p)**(L-i) for i in range(0,L+1,2));assert ok==(1+(1-2*p)**L)/2
        print('cycle Z noise L',L,'F',float(ok**3))
if __name__=='__main__':
    chain();gluing();ref()
    graphcase('theta',[(0,1),(0,2),(2,3),(3,1),(0,4),(4,5),(5,1)],[15,113],2,(6,F(79,1024)))
    graphcase('triangle plus pentagon',[(0,1),(1,2),(2,0),(3,4),(4,5),(5,6),(6,7),(7,3)],[7,248],1,(3,F(11,32)))
    graphcase('odd cycle blocked',[(i,(i+1)%9) for i in range(9)],[511],0,(0,F(1)))
    det();noise()
