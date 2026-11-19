"""S7 additional exact regression: two-copy pencil, heterogeneous charges and labels,
local-phase-averaged moments, and analytic scalable monomial ranks."""
from algebra import *
from itertools import permutations

def kernel(cols):
    # linear map over F2 with column vectors; kernel basis via elimination with witnesses
    p={};out=[]
    for i,v in enumerate(cols):
        t=1<<i
        while v:
            j=v.bit_length()-1
            if j in p:v^=p[j][0];t^=p[j][1]
            else:p[j]=(v,t);break
        if not v:out.append(t)
    return out

def pencil():
    rng=Random(7115);count=0
    for _ in range(36):
        n=rng.randrange(3,7);bs=[1<<i for i in range(n)];t=truth(n,rng.sample([sum(z) for z in combinations(bs,3)],rng. randrange(0,len(list(combinations(bs,3)))+1)))
        R=basis(rng.randrange(1<<n) for j in range(rng.randrange(4)));ell=len(R);D=1<<n
        rs=span(R);# span() ordering is binary index ordering
        # restricted difference phase direct truth.
        v=[t[x]^t[x^s] for s in rs for x in range(D)];U=[1<<i for i in range(n+ell)]
        def cv(a,b,c):return pol(t,a,b,c)
        rad=sum(not any(cv(a,b,c) for b,c in combinations(bs,2)) for a in range(D))
        Z=sum(not any(cv(a,b,c) for b in R for c in bs) for a in range(D))
        rr=sum(not any(cv(a,b,c) for b,c in combinations(bs,2)) for a in rs)
        directrad=0
        for a,s in product(range(D),range(len(rs))):
            actual=[sum(pol(v,a|s<<n,b,c)<<j for j,c in enumerate(U)) for b in U];directrad+=not any(actual)
            z=rs[s];slots=[(b,0) for b in bs]+[(0,b) for b in R]
            pred=[sum((cv(a,b,c)^cv(a^z,b^u,c^w))<<j for j,(c,w) in enumerate(slots)) for b,u in slots]
            assert actual==pred
            if s==0:assert len(basis(actual))==2*len(basis(sum(cv(a,b,c)<<j for j,c in enumerate(R)) for b in bs))
            count+=1
        assert directrad==Z*rr
    print('two-copy pencil: 36 seeded forms/refs,',count,'exact contraction matrices; radical checks')

def hetero():
    rng=Random(7219);num=0
    for _ in range(120):
        n1,n2=rng.randrange(3,5),rng.randrange(3,5);q1=[rng.randrange(16) for j in range(n1)];q2=[rng.randrange(16) for j in range(n2)]
        Q=lambda x,cols:__import__('functools').reduce(int.__xor__,[v for j,v in enumerate(cols) if x>>j&1],0)
        a,b=[truth(n,rng.sample([sum(z) for z in combinations([1<<i for i in range(n)],3)],rng.randrange(0,n-1))) for n in (n1,n2)]
        K=kernel(q1+q2);v=[a[x&((1<<n1)-1)]^b[x>>n1] for x in span(K)];direct=active(v,[1<<j for j in range(len(K))]);
        I=set(Q(x,q1) for x in range(1<<n1))&set(Q(y,q2) for y in range(1<<n2));L=[]
        for t,cols in ((a,q1),(b,q2)):
            src=[x for x in range(len(t)) if Q(x,cols) in I];w=[x for x in src if not Q(x,cols)];L.append(any(pol(t,x,y,z) for x in w for y,z in combinations(src,2)))
        if any(L):pred=True
        else:
            ib=basis(I);lift1=[next(x for x in range(len(a)) if Q(x,q1)==z) for z in ib];lift2=[next(x for x in range(len(b )) if Q(x,q2)==z) for z in ib];pred=any(pol(a,*(lift1[i] for i in ids))^pol(b,*(lift2[i] for i in ids)) for ids in combinations(range(len(ib)),3))
        assert pred==direct;num+=1
    print('heterogeneous charge images',num,'seeded exact overlap/kernel comparisons')

def phaseavg(t,axes):
    D=len(t);n=D.bit_length()-1;w={}
    for a,b in product(range(D),repeat=2):
        z=sum((1-2*(t[x]^t[x^a]^((x&b).bit_count()%2))) for x in range(D));k=(a&b).bit_count()%4;assert k%2==0 or z==0;w[a, b]=0 if k%2 else (-1 if k==2 else 1)*z
    av=F(0)
    for mask in range(1<<len(axes)):
        rr=sum(1<<i for j,i in enumerate(axes) if mask>>j&1);v=w.copy()
        for i in axes:
            if rr>>i&1:
                b0=1<<i
                for a,b in product(range(D),repeat=2):
                    if a&b0 and not b&b0:
                        x,y=v[a,b],v[a,b|b0];v[a,b]=x+y;v[a,b|b0]=y-x
        av+=sum(F(x**4,D**5*4**((a&rr).bit_count())) for (a,b),x in v.items())/(1<<len(axes))
    return av

def compare():
    c=truth(3,[7]);assert phaseavg([0,0],[0])==F(7,8)
    av=phaseavg(c,[0,1,2]);print('CCZ full independent local phase average P2',av,'6q with spectators',av*F(7,8)**3)
    for d in range(4):print('d',d,'charged-axes-only P2',phaseavg(c,list(range(d)))*F(7,8)**(3-d))
    # algebraic scalable copy checks; direct rank and phase truth up to N=4
    for d,N in product((1,2),range(2,5)):
        r=diagnostics(monomial_comp(d,N),direct=False)
        want=(3 if d==2 and N==2 else 3*N-d);assert r[0]==want
    print('monomial d=1,2 scalable nullity sanity N=2..4')
if __name__=='__main__':pencil();hetero();compare()
