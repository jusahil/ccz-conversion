"""New exact controls; not a recovery of the unavailable C22/C22b files."""
from itertools import combinations, product
from collections import Counter
from fractions import Fraction as F
from math import comb, ceil, log
from algebra import basis, span, subs, truth, pol, diagnostics
from new_certificate import rankhist

def xorvec(bs,x):
    y=0
    for i,b in enumerate(bs):
        if x>>i&1: y^=b
    return y

def stabilizers():
    counts=Counter(); maxima={}; real=0
    for W,b in subs(3):
        d=len(b); pairs=[(1<<i)|(1<<j) for i,j in combinations(range(d),2)]
        for offset in sorted({min(x^w for w in W) for x in range(8)}):
            for linear in product(range(4),repeat=d):
                for q in range(1<<len(pairs)):
                    re=im=0
                    for y in range(1<<d):
                        x=offset^xorvec(b,y)
                        exponent=(sum(linear[i] for i in range(d) if y>>i&1)+2*(sum((y&t)==t for j,t in enumerate(pairs) if q>>j&1)+(x==7)))%4
                        if exponent%2: im+=1 if exponent==1 else -1
                        else: re+=1 if exponent==0 else -1
                    val=F(re*re+im*im,8*(1<<d));counts[d]+=1
                    maxima[d]=max(maxima.get(d,F(0)),val)
                    if d==3 and all(x%2==0 for x in linear): real+=1
    assert counts=={0:8,1:112,2:448,3:512} and real==64
    assert maxima=={0:F(1,8),1:F(1,4),2:F(1,2),3:F(9,16)}
    print('stabilizers',dict(counts),'maxima',maxima,'full-support real',real)

def vanishing():
    for name,terms,B in [('A',[7,56,73,82],[2,4,8,32,64]),('B',[7,25,97,42,82,76,52],[3,5,24,40,72])]:
        t=truth(7,terms); assert len(basis(B))==5 and not any(pol(t,*abc) for abc in combinations(B,3))
        count=0
        for l in range(1,128):
            H=basis(x for x in range(128) if (x&l).bit_count()%2==0)
            count+=not any(pol(t,*abc) for abc in combinations(H,3))
        assert count==0;print(name,'127 hyperplanes, 0 vanish; explicit five-plane verified',B)

def gf4():
    def mul(a,b):
        x=0
        for i in range(2):
            if b>>i&1:x^=a<<i
        return x^7 if x&4 else x
    def cc(x):return int(x==7)
    for x,y in product(range(8),repeat=2):
        a,b,c=[((x>>i)&1)|(((y>>i)&1)<<1) for i in range(3)]
        u=mul(mul(a,b),c); trace=u^mul(u,u)
        assert trace==cc(x)^cc(y)^cc(x^y)
    print('GF4 trace 64 assignments')

def normal_forms():
    h=Counter();nulls=Counter()
    for m in range(512):
        rows=[(m>>(3*i))&7 for i in range(3)];r=len(basis(rows));h[r]+=1
        # sum_i a_i (M_i0 b1b2 + M_i1 b1b3 + M_i2 b2b3)
        terms=[(1<<i)|(1<<(3+j))|(1<<(3+k)) for i in range(3) for v,(j,k) in enumerate(combinations(range(3),2)) if rows[i]>>v&1]
        hist=dict(rankhist(truth(6,terms)));eta=6-(hist[0].bit_length()-1);nulls[eta]+=1
        assert eta==[0,3,5,6][r]
    assert h=={0:1,1:49,2:294,3:168};assert nulls=={0:1,3:49,5:294,6:168}
    print('512 coefficient matrices',dict(sorted(h.items())),'nullities',dict(sorted(nulls.items())))

def spectrum(hist,n):
    out=Counter()
    for r,k in hist.items():
        assert r%2==0;out[F(1,2**(r//2))]+=k*2**r
    out[F(0)]=4**n-sum(out.values());return out

def overlap(lhs,rhs,n):
    a=sorted(lhs.items(),reverse=True);b=sorted(rhs.items(),reverse=True); i=j=0;tot=F(0)
    while i<len(a) and j<len(b):
        v,k=a[i];w,l=b[j];q=min(k,l);tot+=v*w*q
        a[i]=(v,k-q);b[j]=(w,l-q);i+=k==q;j+=l==q
    return tot/2**n

def ceilings():
    for name,h,n,k,expected in [('F4',{0:1,4:63},6,2,F(57,64)),('rank-three normal',{0:1,2:7,4:56},6,2,F(121,128)),('determinant',{0:1,4:49,6:294,8:168},9,3,F(13199,16384))]:
        target=Counter({0:1})
        for _ in range(k):
            nxt=Counter()
            for r,h0 in target.items():nxt[r]+=h0;nxt[r+2]+=7*h0
            target=nxt
        val=overlap(spectrum(h,n),spectrum(target,n),n);assert val==expected
        print(name,'sorted-spectrum bound',val)

def noise_and_prices():
    for L in [3,5,9,21]:
        for i in range(601):
            u=F(i,10000);q=sum(F(comb(L,j))*u**j*(1-u)**(L-j) for j in range(0,L+1,2));assert q==(1+(1-2*u)**L)/2
        print('noise',L,'u=.01',float(((1+(1-F(2,100))**L)/2)**3))
    # General unequal-rate identity, independently generated rational example.
    B=[1,2,3,1];u=[F(1,10),F(1,5),F(1,4),F(1,3)]
    prob=F(0)
    for e in product(range(2),repeat=4):
        z=0;p=F(1)
        for r,h,t in zip(B,e,u):z^=r*h;p*=t if h else 1-t
        if z==0:prob+=p
    char=F(0)
    for y in range(4):
        a=F(1)
        for r,t in zip(B,u):
            if (r&y).bit_count()%2:a*=1-2*t
        char+=a/4
    assert char==prob;print('unequal-rate independent example',prob)
    n=ceil((F(7,2)**2/(2*F(4,100)**2))*log(2/.05));assert n==14122;print('Hoeffding sufficient trials',n)
    for k in range(1,8):
        for i in range(12):
            w=F(i,4);cost=min(N+l*w for N in range(k,k+3) for l in range(4) if (max(0,N-2) if l==0 else N-1 if l<3 else N)>=k)
            assert cost==k+min(3*w,1+w,2)
    print('84 joint-block price checks; reference price ties 1/2 and 1')
    for lam in (F(0),F(1,2),F(1)):
        p=(1+lam)/2
        for T in range(1,5):
            shared=sum(w for e,w in [(0,p),(1,1-p)] if not any([e]*T))
            independent=sum((p**(T-sum(e))*(1-p)**sum(e)) for e in product(range(2),repeat=T) if not any(e))
            assert (shared,independent)==(p,p**T)
    print('common vs independent hidden sign: lambda=1/2, T=2: 3/4, 9/16')
def conditional_and_bridge():
    t=truth(6,[7]);E=[1<<i for i in range(6)]
    for W,d,N in [([1,2,4],0,1),([1,2,8],1,2),([8,16,32],3,3)]:
        got=min(j for j in range(4) if any(pol(t,*abc) for abc in product(*([W]*(3-j)+[E]*j))))
        assert got==d
        # Build collective kernel of common coordinate parity charge.
        n=6
        for copies in range(1,4):
            U=[v<<(n*i) for i in range(copies) for v in W]+[(v<<(n*i))^(v<<(n*(copies-1))) for i in range(copies-1) for v in E if v not in W]
            active=any(sum(pol(t,*[(v>>(n*i))&63 for v in abc]) for i in range(copies))%2 for abc in combinations(U,3))
            assert active==(copies>=N)
    print('conditional one-space cup types 0,1,3; copy thresholds 1,2,3')
    # Two triangles joined by bridge: every binary flow vanishes on bridge.
    es=[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(2,3)]
    def flows(es):
        a=[]
        for x in range(1<<len(es)):
            z=0
            for i,(u,v) in enumerate(es):
                if x>>i&1:z^=(1<<u)^(1<<v)
            if not z:a.append(x)
        return a
    assert len(flows(es))==4 and all(x>>6==0 for x in flows(es))
    split=es[:-1]+[(2,6),(6,3)]
    assert len(flows(split))==4 and all(x>>6==0 for x in flows(split))
    print('bridge and unattached subdivision invisible in all four flows')

def seeded_map_domain():
    # Reproduce the random stream in C12: forty radical cases precede map draws.
    from random import Random
    rng=Random(20261117)
    for _ in range(40):
        n=5
        [tri for tri in combinations(range(n),3) if rng.randrange(2)]
        [rng.randrange(1<<n) for j in range(rng.randrange(4))]
    counts=Counter();both=combined=0
    for _ in range(400):
        a=tuple(rng.randrange(8) for i in range(3));b=tuple(rng.randrange(8) for i in range(3))
        rng.choice([0,1]);rng.choice([0,1])
        ranks=[len(basis(a)),len(basis(b))];counts.update(ranks);both+=ranks==[3,3];combined+=len(basis(a+b))==3
    assert counts=={0:1,1:81,2:454,3:264} and (both,combined)==(45,347)
    print('C12 seeded map ranks',dict(sorted(counts.items())),'both full',both,'combined full',combined)

if __name__=='__main__':
    stabilizers();vanishing();gf4();normal_forms();ceilings();noise_and_prices();conditional_and_bridge();seeded_map_domain()
