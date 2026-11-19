"""Exact Pauli spectra and rearrangement overlap for Boolean phase states."""
from fractions import Fraction as F
from collections import Counter
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
