"""Exact tensor and protocol checks. Requires Python 3 and NumPy. No random certificate of impossibility.
All random checks use a fixed seed; Fano one-drop enumeration is exhaustive.
"""
from itertools import combinations, product
from collections import Counter, defaultdict
from fractions import Fraction
from random import Random
from math import log2
import numpy as np
from time import perf_counter

def rank(rows):
    piv={}
    for a in rows:
        a=int(a)
        while a:
            i=a.bit_length()-1
            if i in piv:a^=piv[i]
            else:piv[i]=a;break
    return len(piv)
def span(B):
    S={0}
    for b in B:S|={x^b for x in S}
    return S
def basis(S):
    B=[];T={0}
    for x in sorted(S):
        if x not in T:B.append(x);T|={y^x for y in T}
    return B
def spaces(n):
    seen={frozenset({0})};todo=list(seen)
    while todo:
        S=todo.pop()
        for x in range(1<<n):
            if x not in S:
                T=S|{a^x for a in S}
                if T not in seen:seen.add(T);todo.append(T)
    return seen
def xorvec(B,x):
    z=0
    for i,b in enumerate(B):
        if x>>i&1:z^=b
    return z
def kernel(cols):
    # incremental elimination; records independent input relations
    piv={};K=[]
    for i,y in enumerate(cols):
        u=1<<i
        while y:
            j=y.bit_length()-1
            if j in piv:z,v=piv[j];y^=z;u^=v
            else:piv[j]=(y,u);break
        if not y:K.append(u)
    return K
def make(n,mons):return [sum((x&m)==m for m in mons)%2 for x in range(1<<n)]
def terms(strings):return [sum(1<<int(i) for i in s) for s in strings]
def pol(f,a,b,c):return f[0]^f[a]^f[b]^f[c]^f[a^b]^f[a^c]^f[b^c]^f[a^b^c]
def restrict(f,B):return [f[xorvec(B,x)] for x in range(1<<len(B))]
def nonzero(f,B):return any(pol(f,a,b,c) for a,b,c in combinations(B,3))
def hist(f):
    n=len(f).bit_length()-1
    return Counter(rank([sum(pol(f,a,1<<i,1<<j)<<j for j in range(n)) for i in range(n)]) for a in range(1<<n))
def moment(H,n):return sum(Fraction(v,2**(n+k)) for k,v in H.items())
def rad(f):
    n=len(f).bit_length()-1;E=[1<<i for i in range(n)]
    return {a for a in range(1<<n) if not any(pol(f,a,b,c) for b,c in combinations(E,2))}
def diffphase(f,R):return [f[x]^f[x^xorvec(R,t)] for t in range(1<<len(R)) for x in range(len(f))]
def radtest(f,R):
    n=len(f).bit_length()-1;E=[1<<i for i in range(n)];F=diffphase(f,R)
    K={a for a in range(1<<n) if not any(pol(f,t,a,b) for t in R for b in E)}
    rr=rad(f);pred={a|(t<<n) for a in K for t in range(1<<len(R)) if xorvec(R,t) in rr}
    assert pred==rad(F)
    return len(pred)
def phase_pair(f,g,n,x):return f[x&((1<<n)-1)]^g[x>>n]
def pairtest(q1,q2,f,g):
    n=3;B=kernel(q1+q2);I=span(q1)&span(q2)
    V1=[x for x in range(8) if xorvec(q1,x) in I];V2=[x for x in range(8) if xorvec(q2,x) in I]
    K1=[x for x in V1 if xorvec(q1,x)==0];K2=[x for x in V2 if xorvec(q2,x)==0]
    mixed=any(pol(f,a,b,c) for a in basis(K1) for b,c in combinations(basis(V1),2)) or any(pol(g,a,b,c) for a in basis(K2) for b,c in combinations(basis(V2),2))
    IB=basis(I)
    l1=[next(x for x in V1 if xorvec(q1,x)==i) for i in IB];l2=[next(x for x in V2 if xorvec(q2,x)==i) for i in IB]
    residual=any(pol(f,*[l1[i] for i in tri])^pol(g,*[l2[i] for i in tri]) for tri in combinations(range(len(IB)),3))
    direct=any(pol(f,a&7,b&7,c&7)^pol(g,a>>3,b>>3,c>>3) for a,b,c in combinations(B,3))
    assert direct==(mixed or residual)

def graphcheck(n,edges):
    B=kernel([(1<<u)|(1<<v) for u,v in edges]);k=len(B)
    patterns=Counter(sum(((b>>e)&1)<<j for j,b in enumerate(B)) for e in range(len(edges)))
    odd=[c for c,num in patterns.items() if c and num%2]
    # contraction rows of the cycle triple-overlap form (repeated indices included)
    rows=[sum(((a&b&c).bit_count()%2)<<i for i,c in enumerate(B)) for a in B for b in B]
    assert rank(rows)==rank(odd)
    direct=any((a&b&c).bit_count()%2 for a in B for b in B for c in B)
    assert direct==bool(odd)
    # Every coedge class is an exact intersection of two cycle vectors.
    cycles=list(span(B));
    for c in patterns:
        if c==0:continue
        S=sum(1<<e for e in range(len(edges)) if sum(((b>>e)&1)<<j for j,b in enumerate(B))==c)
        assert any(a&b==S for a in cycles if a&S==S for b in cycles if b&S==S)
    return k,rank(odd)

def pauli_integer(f):
    D=len(f);A=np.repeat(np.arange(D),D);B=np.tile(np.arange(D),D)
    weights=np.array([x.bit_count() for x in range(D)],dtype=np.int64)
    s=1-2*np.asarray(f,dtype=np.int64);w=[]
    for a in range(D):
        t=s*s[np.arange(D)^a];h=1
        while h<D:
            z=t.reshape(-1,2*h);u=z[:,:h].copy();v=z[:,h:].copy();z[:,:h]=u+v;z[:,h:]=u-v;h*=2
        w.extend(t)
    w=np.array(w,dtype=np.int64);wt=weights[A&B]
    assert np.all(w[wt%2==1]==0)
    w=w*np.array([1,0,-1,0],dtype=np.int64)[wt%4]
    return D,A,B,weights,w,wt

def direct_pauli(v,a,b):
    D=len(v);idx=np.arange(D);ph=(1j)**((a&b).bit_count()%4)
    signs=np.array([1-2*((int(x)&b).bit_count()%2) for x in idx]);out=np.empty(D,dtype=complex);out[idx^a]=ph*signs*v
    return out

def branches(f,seed=121):
    D,A,B,wts,w,wt=pauli_integer(f);n=D.bit_length()-1;C=Counter();H=Counter();matches=[];rng=Random(seed)
    assert np.count_nonzero(abs(w)==D)==1
    # direct check of 48 projected states, 30 Pauli entries each
    probes=set(rng.sample(range(1,D*D),48));v=np.array([1-2*x for x in f],dtype=np.int64)
    for p in range(1,D*D):
        a,b=int(A[p]),int(B[p]);idx=(A^a)*D+(B^b)
        comm=((wts[A&b]+wts[B&a])%2==0);phase=(wt+wt[p]-wt[idx]+2*wts[B&a])%4
        assert np.all(phase[comm]%2==0);sign=1-2*(phase//2)
        for s in (1,-1):
            den=int(D+s*w[p]);
            if den==0:continue
            nums=np.where(comm,w+s*sign*w[idx],0)
            if p in probes:
                u=v+s*direct_pauli(v,a,b);assert int(np.vdot(u,u).real)==2*den
                for j in rng.sample(range(D*D),30):
                    z=np.vdot(u,direct_pauli(u,int(A[j]),int(B[j])))
                    assert z.imag==0 and z.real==2*int(nums[j])
            ns=int(np.count_nonzero(abs(nums)==den));eta=n-(ns.bit_length()-1);assert ns==1<<(n-eta)
            C[eta]+=1
            if eta==6:
                z=int(np.sum(nums**4));val=Fraction(z,D*den**4);H[val]+=1
                if val==Fraction(121,1024):matches.append((a,b,s))
    return C,H,matches

def kerrow(n,p):
    j=p.bit_length()-1;return [(1<<i)^(((p>>i)&1)<<j) for i in range(n) if i!=j]

def verify_protocol(columns,f,R):
    n=len(columns);q=columns;R=basis(span(R));RB=span(R)
    # solve p(x) in R by brute force (only small test inputs)
    K=basis(x for x in range(1<<n) if xorvec(q,x) in RB)
    L=next((list(t) for t in combinations(K,3) if pol(f,*t)),None)
    if L is None:return False,0
    S=span(L);reps=sorted({min(x^y for y in S) for x in range(1<<n)})
    # check all outcomes and labels; use arbitrary extension of dual coordinates
    dual=[]
    for j in range(3):dual.append(next(l for l in range(1<<n) if all(((l&a).bit_count()%2)==(i==j) for i,a in enumerate(L)) ))
    count=0
    for x0 in reps:
        g=[f[x0^xorvec(L,y)]^int(y==7) for y in range(8)]
        assert pol(g,1,2,4)==0
        for y in range(8):
            x=x0^xorvec(L,y);out=sum((((l&(x^x0)).bit_count()%2)<<j) for j,l in enumerate(dual));assert out==y
            for ref in RB:
                shift=xorvec(L,out);xnew=x^shift;rnew=ref^xorvec(q,shift)
                assert xnew==x0 and rnew in RB and xorvec(q,x)^ref==xorvec(q,xnew)^rnew
                assert f[x]^g[out]==int(out==7);count+=1
    return True,count

if __name__=='__main__':
    rng=Random(20261117);start=perf_counter()
    # Two-copy radical identity: every 4D cubic and every charge subspace.
    trials=0
    S=spaces(4);mons=[sum(1<<i for i in tri) for tri in combinations(range(4),3)]
    for mask in range(16):
        f=make(4,[m for i,m in enumerate(mons) if mask>>i&1])
        for R in S:radtest(f,basis(R));trials+=1
    for _ in range(40):
        n=5;f=make(n,[sum(1<<i for i in tri) for tri in combinations(range(n),3) if rng.randrange(2)]);R=basis(rng. randrange(1<<n) for j in range(rng.randrange(4)));radtest(f,R)
    print('Two-copy radical identity:',trials,'exhaustive 4D cases + 40 seeded 5D cases',flush=True)
    # Heterogeneous: all maps F2^3 -> F2^2; additionally seeded F2^3 -> F2^3.
    f0=make(3,[]);f1=make(3,[7]);ct=0
    maps=list(product(range(4),repeat=3))
    for q1 in maps:
        for q2 in maps:
            for f,g in product([f0,f1],repeat=2):pairtest(q1,q2,f,g);ct+=1
    for _ in range(400):pairtest(tuple(rng.randrange(8) for i in range(3)),tuple(rng.randrange(8) for i in range(3)),rng. choice([f0,f1]),rng.choice([f0,f1]))
    print('Heterogeneous pair test:',ct,'exhaustive + 400 seeded rank-3 shared label cases',flush=True)
    # all small graph incidence constraints
    gc=0;active=0
    for n in range(1,6):
        E=list(combinations(range(n),2))
        for mask in range(1<<len(E)):
            edges=[e for i,e in enumerate(E) if mask>>i&1];k,o=graphcheck(n,edges);gc+=1;active+=o>0
    print('Graph parity/radical + exact cycle intersections:',gc,'graphs; active',active,flush=True)
    # Four types, N=1,2,3, selected reference supports; all parity branches and reference labels.
    checks=0
    for d in range(4):
        q=[0]*(3-d)+[1<<i for i in range(d)]
        for N in (1,2,3):
            n=3*N;F=make(n,[7<<(3*i) for i in range(N)]);qN=q*N
            opt={1:d,2:(1 if d==3 else 0),3:0}[N]
            R=[1<<i for i in range(opt)]
            yes,num=verify_protocol(qN,F,R);assert yes;checks+=num
            if opt:assert not verify_protocol(qN,F,R[:-1])[0]
    print('Witness circuit all branches/basis labels:',checks,'checks; 12 constructive cases + below-cost controls',flush= True)
    # Matched covectors example: block-specific directions, and all-independent selection.
    for d in range(1,5):
        if d>2:break # direct truth tables only for small support here
        f=make(3*d,[7<<(3*i) for i in range(d)])
        for r in range(d+1):
            R=[1<<(3*i) for i in range(r)];F=diffphase(f,R);H=hist(F)
            assert moment(H,3*d+r)==Fraction(11,32)**r
    print('Coordinate reference benchmark: d=1,2, all r<=d; exact moments agree',flush=True)
    # Vanishing codimension for two 6D target forms by all hyperplanes (codim1).
    F4=make(6,terms(['012','345'])) # overwrite with f(a)+f(b)+f(a+b)
    F4=[int((x&7)==7)^int((x>>3)==7)^int(((x&7)^(x>>3))==7) for x in range(64)]
    C2=make(6,terms(['012','345']))
    for name,f in [('CCZ^2',C2),('GF4',F4)]:
        H=hist(f);assert H[0]==1 and all(nonzero(f,kerrow(6,p)) for p in range(1,64));assert not nonzero(f,[1,2,8,16]); print(name,'rank histogram',dict(sorted(H.items())),'P2',moment(H,6),'vanishing codim 2',flush=True)
    # Matched scalar invariants with different exact yields, 7 qubits
    fA=make(7,terms(['012','345','036','146']));fB=make(7,terms(['012','034','056','135','146','236','245']))
    for name,f in [('A',fA),('Fano',fB)]:
        H=hist(f);assert H[0]==1 and moment(H,7)==Fraction(29,512);print(name,'rank histogram',dict(sorted(H.items())),'P2' ,moment(H,7),flush=True)
    assert restrict(fA,[1<<i for i in range(6)])==C2
    aff=sum(hist(restrict(fB,kerrow(7,p)))=={0:1,2:14,4:49} for p in range(1,128));assert aff==0
    print('Fano: 127 hyperplanes, zero CCZ^2 rank-spectrum matches',flush=True)
    C,H,M=branches(fB);assert not M and sum(C.values())==32766 and C[6]==31644
    print('Fano exhaustive one-drop search: nullity counts',dict(sorted(C.items())),flush=True)
    print('Fano eta=6 projected P2 histogram:',[(str(k),v) for k,v in sorted(H.items())],flush=True)
    print('Fano matches to CCZ^2 P2=121/1024:',len(M),'; all checks passed',flush=True)
