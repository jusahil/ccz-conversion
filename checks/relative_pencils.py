from itertools import combinations,product
from collections import Counter

def span(B):
    s={0}
    for b in B:s|={a^b for a in s}
    return frozenset(s)
def basis(S):
    S=list(S);B=[];s={0}
    for x in sorted(S):
        if x not in s:B.append(x);s|={a^x for a in s}
    assert s==set(S)
    return B
def spaces(n):
    S={frozenset([0])}
    for x in range(1,1<<n): S|={frozenset(set(s)|{a^x for a in s}) for s in S}
    return sorted(S,key=lambda s:(len(s),sorted(s)))
def rank(rows):
    d={}
    for a in rows:
        while a:
            i=a.bit_length()-1
            if i in d:a^=d[i]
            else:d[i]=a;break
    return len(d)
def makeform(n,mask):
    mons=[sum(1<<i for i in c) for k,c in enumerate(combinations(range(n),3)) if mask>>k&1]
    f=[sum((x&m)==m for m in mons)%2 for x in range(1<<n)]
    def t(a,b,c):return f[0]^f[a]^f[b]^f[c]^f[a^b]^f[a^c]^f[b^c]^f[a^b^c]
    return t

def nz(T,B):return any(T(*a) for a in combinations(B,3))
ss={n:spaces(n) for n in range(5)}
forms=0; acts=0; pencils=0; rads=0; lines=0
for m in range(16):
    T=makeform(4,m)
    for r in range(5):
        charge=lambda x:x & ((1<<r)-1)
        vals=[]
        for N in (1,2,3):
            Rset=ss[r] if N<3 else [frozenset([0])]
            possible=[]
            def total(x):
                y=0
                for j in range(N):y^=charge(x>>(4*j))
                return y
            def TT(a,b,c):return sum(T((a>>(4*j))&15,(b>>(4*j))&15,(c>>(4*j))&15) for j in range(N))%2
            for R in Rset:
                B=basis(x for x in range(1<<(4*N)) if total(x) in R)
                if nz(TT,B):possible.append(len(R).bit_length()-1)
                acts+=1
            vals.append(min(possible,default=99))
        d=vals[0];want={0:[0,0,0],1:[1,0,0],2:[2,0,0],3:[3,1,0],99:[99,99,99]}[d]
        assert vals==want,(m,r,vals,want)
        forms+=1
    # Full charge Q=F2^4. Direct difference phase contraction vs displayed block pencil.
    for R in ss[4]:
        J=basis(R); ell=len(J); k=4+ell
        emb=[__import__('functools').reduce(int.__xor__,(v for i,v in enumerate(J) if t>>i&1),0) for t in range(1<<ell)];embed=emb.__getitem__
        def TD(x,y,z):
            a,b,c=x&15,y&15,z&15
            return T(a,b,c)^T(a^embed(x>>4),b^embed(y>>4),c^embed(z>>4))
        E=[1<<i for i in range(k)]; Q=[1<<i for i in range(4)]
        rad=[]; moments=Counter()
        for v in range(1<<k):
            a=v&15;s=embed(v>>4)
            def predicted(i,j):
                if i<4 and j<4:return T(s,Q[i],Q[j])
                if i<4:return T(a^s,Q[i],J[j-4])
                if j<4:return T(a^s,J[i-4],Q[j])
                return T(a^s,J[i-4],J[j-4])
            rows=[sum(TD(v,x,y)<<j for j,y in enumerate(E)) for x in E]
            assert rows==[sum(predicted(i,j)<<j for j in range(k)) for i in range(k)]
            h=rank(rows);moments[h]+=1
            if not s: assert h==2*rank([sum(T(a,x,y)<<j for j,y in enumerate(J)) for x in Q])
            if not h:rad.append(v)
            pencils+=1
        pred_rad={a|(t<<4) for a in range(16) for t in range(1<<ell) if all(not T(a,x,y) for x in J for y in Q) and all(not T(embed(t),x,y) for x in Q for y in Q)}
        assert set(rad)==pred_rad;rads+=1
        if ell==1:
            from fractions import Fraction
            r=rank([sum(T(J[0],a,b)<<j for j,b in enumerate(Q)) for a in Q]);h=r//2
            actual=sum(Fraction(c,2**(k+j)) for j,c in moments.items())
            expected=Fraction(1,8)+Fraction(7,8*4**h) if h else 1
            assert actual==expected;lines+=1
print('activation cases (form, projection rank):',forms,'reference/copy zero tests:',acts)
print('pencil contractions:',pencils,'radical comparisons:',rads,'line moment comparisons:',lines)
