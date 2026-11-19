"""Exact independent controls for robustness and resource budgets (Python standard library)."""
from fractions import Fraction as F
from itertools import product, combinations
from collections import Counter
from math import comb, ceil, log
from functools import reduce
from operator import xor
from random import Random

def X(xs): return reduce(xor,xs,0)
def parity(x): return x.bit_count()%2

def stabilizer_overlap():
    # Enumerate affine supports by subsets, then choose a coordinate basis.
    counts=Counter(); best=Counter(); real=0
    for d in range(4):
        for S in combinations(range(8),2**d):
            W={a^S[0] for a in S}
            if any(a^b not in W for a,b in product(W,repeat=2)): continue
            B=[]; span=[0]
            for a in sorted(W):
                if a not in span: B.append(a); span += [x^a for x in span]
            assert len(B)==d
            pairs=list(combinations(range(d),2))
            for lin in product(range(4),repeat=d):
                for quad in product(range(2),repeat=len(pairs)):
                    v=[0]*4
                    for y,x in enumerate(span):
                        bits=[(y>>j)&1 for j in range(d)]
                        phase=(sum(a*b for a,b in zip(lin,bits))+2*sum(q*bits[i]*bits[j] for q,(i,j) in zip(quad,pairs))+2*((x^S[0])==7))%4
                        v[phase]+=1
                    f=F((v[0]-v[2])**2+(v[1]-v[3])**2,8*2**d)
                    counts[d]+=1;best[d]=max(best[d],f)
                    real+=d==3 and all(a%2==0 for a in lin)
    assert list(counts.values())==[8,112,448,512] and list(best.values())==[F(1,8),F(1,4),F(1,2),F(9,16)] and real==64
    print('Stabilizer counts',dict(counts),'maxima',dict(best),'real full support',real)

def paulis(n,term):
    D=1<<n; v={}
    for a,b in product(range(D),repeat=2):
        w=sum((-1)**(int((x&term)==term)^int(((x^a)&term)==term)^parity(x&b)) for x in range(D))
        j=(a&b).bit_count()%4
        assert j%2==0 or w==0
        v[a,b]=w*(-1 if j==2 else 1) if j%2==0 else 0
    return v

def moment_quadrature(v,n,axes):
    # Traverse all angle settings. For transverse pairs an integer Hadamard
    # transform represents rotation by pi/4; track the omitted sqrt(2) factors.
    D=2**n; total=F(0)
    for sel in product((0,1),repeat=len(axes)):
        r=sum(1<<i for i,b in zip(axes,sel) if b); w=v.copy()
        for i in axes:
            if not r>>i&1: continue
            for a,b in v:
                if a>>i&1 and not b>>i&1:
                    pair=(a,b^(1<<i)); p,q=w[a,b],w[pair]; w[a,b]=p+q; w[pair]=q-p
        total+=sum(F(z**4,D**5*4**((a&r).bit_count())) for (a,b),z in w.items())
    return total/2**len(axes)

def diagnostics():
    v=paulis(3,7); hist=Counter(abs(x) for x in v.values());assert hist=={8:1,4:28,0:35}
    assert [sum(F(z,8)**(2*m) for z in v.values())/8 for m in (2,3)]==[F(11,32),F(23,128)]
    print('Pauli magnitude numerator histogram',dict(hist),'P2=11/32 P3=23/128')
    vals=[]
    for d in range(4):
        # d charged coordinates among the three CCZ sites, 3-d charged spectators
        vals.append(moment_quadrature(v,3,list(range(d)))*F(7,8)**(3-d))
    allavg=moment_quadrature(v,3,[0,1,2])*F(7,8)**3
    assert vals==[F(3773,16384),F(245,1024),F(133,512),F(151,512)] and allavg==F(51793,262144)
    print('Phase averages all axes',allavg,'charged axes',vals)
    # Direct estimator tested on all Z-shifted targets (orthogonal phase states).
    for e in range(8):
        obs=F(1,8)+F(7,4)*sum(F((1 if w>0 else -1)*w*(-1)**parity(a&e),8) for (a,b),w in v.items() if abs(w)==4)/28
        assert obs==int(e==0)
    print('8 exact Pauli estimator expectation checks; sufficient iid trials',ceil(F(7,2)**2*log(40)/(2*F(1,25)**2)))

def noise():
    count=0
    for L in [3,5,9,21]:
        for j in range(901):
            u=F(j,10000); q=sum(comb(L,i)*u**i*(1-u)**(L-i) for i in range(0,L+1,2));assert q==(1+(1-2*u)**L)/2;count+=1
        q=F(1,2)*(1+F(49,50)**L); cut=(1-(2*(9/16)**(1/3)-1)**(1/L))/2
        print('L',L,'F(.01)',float(q**3),'free overlap crossing',format(cut,'.9f'))
    print('exact cycle grid comparisons',count)
    rng=Random(1926);cases=[([1,2,3,1],[F(1,10),F(1,5),F(1,4),F(1,3)],2)]
    for k,m in product(range(1,4),range(3,8)):
        if m<k:continue
        for _ in range(8):
            rows=[1<<i for i in range(k)]+[rng.randrange(1<<k) for _ in range(m-k)]
            rng.shuffle(rows);rates=[F(rng.randrange(7),12) for _ in rows];cases.append((rows,rates,k))
    for i,(rows,rates,k) in enumerate(cases):
        # Dynamic convolution in the error syndrome group, vs Fourier characters.
        p=[F(1)]+[F(0)]*(2**k-1)
        for r,u in zip(rows,rates): p=[(1-u)*p[z]+u*p[z^r] for z in range(2**k)]
        char=[]
        for y in range(2**k):
            a=F(1)
            for r,u in zip(rows,rates):
                if parity(r&y):a*=1-2*u
            char.append(a)
        for z in range(2**k):assert p[z]==sum(a*(-1)**parity(y&z) for y,a in enumerate(char))/2**k
        if i==0:assert p[0]==F(239,600);print('unequal-rate example',p[0])
    print('exact full syndrome distributions, seeded finite encodings',len(cases))

def accounting():
    count=0
    for k in range(1,11):
        for a,b in product(range(5),range(9)):
            wS,wR=F(a,2),F(b,4)
            val=min(N*wS+l*wR for N in range(k,k+3) for l in range(4) if (max(0,N-2) if l==0 else N-1 if l<3 else N)>=k)
            assert val==k*wS+min(3*wR,wS+wR,2*wS);count+=1
    print('joint budget checks incl zero weights',count)
    cnt=0
    for T in range(1,9):
        for wR in [F(i,4) for i in range(33)]:
            # Enumerate one-shot menus across tasks, pay for maximal reference once.
            # counts of each option suffice; choices can use subspaces of one fixed token.
            val=min((a+2*b+3*(T-a-b))+wR*(3 if a else 1 if b else 0) for a in range(T+1) for b in range(T-a+1))
            assert val==min(T+3*wR,2*T+wR,3*T);cnt+=1
    print('separate-task menus incl mixed choices',cnt)
    for v,e in [(3,3),(6,7)]:print('zero-only fusion',v,e,'success',F(1,2**(3*(v-1))),'expected modules',e*2**(3*(v-1)))
    assert [9*N+3*l+3 for N,l in [(1,3),(2,1),(3,0)]]==[21,24,30]
    print('illustrative dressed widths 21, 24, 30')

def mixtures():
    c=F(9,16);checks=0
    for l,p in product([F(i,16) for i in range(17)],repeat=2):
        g=min(1-l,p);b=p-g;assert 0<=b<=l
        # linear objective at endpoints of feasible interval for good mass
        vals=[x+c*(p-x) for x in (max(F(0),p-l),min(p,1-l))]
        assert max(vals)==g+c*b;checks+=1
    print('readable flag LP grid',checks)
    for l in [F(i,10) for i in range(11)]:
        p=(1+l)/2
        assert max(c,p)<=c+(1-c)*l
        for T in range(1,9):
            # common sign has only two outcomes; independent has binomial weights
            assert sum(comb(T,j)*p**(T-j)*(1-p)**j for j in range(T+1))==1
            assert p>=p**T
    print('mixed reference checks',11*8)
    # Finite classical states saturate spectral-width trace bound: transfer mass
    # from accepted orthogonal failure to accepted target, holding p=1.
    for c in (F(9,16),F(121,128),F(57,64),F(13199,16384)):
        d=1-c;assert (1-c)*1==d
        print('ceiling',c,'symmetric deterministic allowance <',d/2)

if __name__=='__main__':
    stabilizer_overlap();diagnostics();noise();accounting();mixtures()
