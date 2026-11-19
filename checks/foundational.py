from itertools import combinations,product
from collections import Counter
from fractions import Fraction
from math import log2

def span(B):
 s={0}
 for b in B:s |= {x^b for x in s}
 return s

def basis(S):
 B=[];t={0}
 for x in sorted(S):
  if x not in t:B.append(x);t|={y^x for y in t}
 return B

def spaces(n):
 allS={frozenset({0})}; pending=list(allS)
 while pending:
  S=pending.pop()
  for x in range(1<<n):
   if x not in S:
    T=S|{x^y for y in S}
    if T not in allS:allS.add(T);pending.append(T)
 return allS

def val(mask,x):return (mask[x])
def pol(f,a,b,c):return f[0]^f[a]^f[b]^f[c]^f[a^b]^f[a^c]^f[b^c]^f[a^b^c]
def nonzero(f,B):return any(pol(f,a,b,c) for a,b,c in combinations(B,3))
def extend(W,n):
 B=basis(W);sp=span(B);C=[]
 for x in range(1<<n):
  if x not in sp:C.append(x);sp|={x^y for y in sp}
 return B,C

def kernbasis(W,n,N):
 B,C=extend(W,n)
 return [x<<(n*i) for i in range(N) for x in B]+[(x<<(n*i))^(x<<(n*(N-1))) for i in range(N-1) for x in C]

def combined(f,n,N):
 return [sum(f[(x>>(i*n))&((1<<n)-1)] for i in range(N))%2 for x in range(1<<(n*N))]

def restriction(f,B):
 return [f[sum(0 for _ in ()) ^ xorvec(B,x)] for x in range(1<<len(B))]

def xorvec(B,x):
 y=0
 for i,b in enumerate(B):
  if (x>>i)&1:y^=b
 return y

def fwht(a):
 a=a[:];h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):x,y=a[j],a[j+h];a[j]=x+y;a[j+h]=x-y
  h*=2
 return a

def moment(f):
 D=len(f);s=0
 for a in range(D):
  w=fwht([1-2*(f[x]^f[x^a]) for x in range(D)])
  s+=sum(t**4 for t in w)
 return Fraction(s,D**5)

def rank(rs):
 B={}
 for r in rs:
  while r:
   j=r.bit_length()-1
   if j in B:r^=B[j]
   else:B[j]=r;break
 return len(B)

def rhist(f):
 n=len(f).bit_length()-1
 return Counter(rank([sum(pol(f,a,1<<i,1<<j)<<j for j in range(n)) for i in range(n)]) for a in range(1<<n))

def mformula(h,n):return sum(Fraction(v,2**(n+k)) for k,v in h.items())

if __name__=='__main__':
 n=4;S=spaces(n);triples=list(combinations(range(n),3));cnt=Counter();ds=Counter();errors=[];total=0;refchecks=0
 for mask in range(16):
  mon=[sum(1<<i for i in t) for j,t in enumerate(triples) if mask>>j&1]
  f=[sum((x&t)==t for t in mon)%2 for x in range(16)]
  FF={N:combined(f,n,N) for N in (1,2,3)}
  for W in S:
   B=basis(W);eb=[1<<i for i in range(n)]
   d=next((j for j in range(4) if any(pol(f,*args) for args in product(*([B]*(3-j)+[eb]*j)))),99)
   pred={0:1,1:2,2:2,3:3,99:99}[d]
   actual=next((N for N in (1,2,3) if nonzero(FF[N],kernbasis(W,n,N))),99)
   if pred!=actual:errors.append((mask,W,pred,actual))
   cnt[actual]+=1;ds[d]+=1;total+=1
   # R ranges over subspaces of V/W, equivalently superspaces H containing W.
   # q_N^{-1}(R) spans ker q_N plus lifts of H/W in the final register.
   for N in (1,2):
    best=99
    for H in sorted((H for H in S if W<=H),key=len):
     BB=kernbasis(W,n,N)+[h<<(n*(N-1)) for h in basis(H)]
     BB=basis(span(BB))
     if nonzero(FF[N],BB):best=len(basis(H))-len(B);break
    target=d if N==1 else (99 if d==99 else (1 if d==3 else 0))
    if best!=target:errors.append(('ref',mask,W,N,best,target))
    refchecks+=1
 print('spaces',len(S),'pairs',total,'activation',dict(sorted(cnt.items())),'relative d',dict(sorted(ds.items())),'reference minima checks',refchecks,'errors',len(errors));assert not errors
 # Reduced active-coordinate polynomials; neutral spectator directions removed.
 A=lambda N:[sum(((x>>(3*i))&7)==7 for i in range(N))%2 for x in range(1<<(3*N))]
 def Bphase(N):
  if N==1:return [0]
  # first 2N coords u, last N-1 independent z, last z=sum
  D=3*N-1;out=[]
  for x in range(1<<D):
   zs=[(x>>(2*N+i))&1 for i in range(N-1)];zs.append(sum(zs)%2)
   out.append(sum(((x>>(2*i))&3)==3 and zs[i] for i in range(N))%2)
  return out
 def Cphase(N):
  if N<3:return [0]
  out=[]
  for x in range(1<<6):
   a=x&7;b=x>>3;out.append(int(a==7)^int(b==7)^int(a^b==7))
  return out
 for name,fn in [('A',A),('B',Bphase),('C',Cphase)]:
  for N in (1,2,3):
   f=fn(N);h=rhist(f);p=moment(f);assert p==mformula(h,len(f).bit_length()-1);print(name,N,'rank',dict(sorted(h.items())),'P2',p,'M2',f'{-log2(p):.9f}')
 # GF(4) trace in basis 1,t, t^2=t+1
 def mul(a,b):
  v=0
  for j in range(2):
   if b>>j&1:v^=a<<j
  if v&4:v^=7
  return v
 for a,b in product(range(8),repeat=2):
  z=[((a>>i)&1)|(((b>>i)&1)<<1) for i in range(3)]
  tr=mul(mul(*z[:2]),z[2])>>1
  assert tr==(int(a==7)^int(b==7)^int(a^b==7))
 print('GF4 trace: 64/64 basis labels matched')
 # Exact fidelity check over all full-support 3-qubit stabilizer phases i^lin (-1)^quad
 # Affine support <=4 gives bound <=1/2 by Cauchy-Schwarz, no enumeration needed.
 best=0
 for ell in product(range(4),repeat=3):
  for q in range(8):
   v=0j
   for x in range(8):
    b=[x>>i&1 for i in range(3)]
    quad=sum((q>>j&1)*b[a]*b[c] for j,(a,c) in enumerate(combinations(range(3),2)))%2
    v+=(-1)**(int(x==7)^quad)*(1j)**sum(ell[i]*b[i] for i in range(3))
   best=max(best,int(round(abs(v)**2)))
 print('512 full-support stabilizer phases: maximum squared overlap',Fraction(best,64));assert best==36

 # Reference star phases: exact Walsh and contraction-rank checks.
 for h in range(1,5):
  n=2*h+1
  f=[(x&1)*(sum(((x>>(2*i+1))&1)*((x>>(2*i+2))&1) for i in range(h))%2) for x in range(1<<n)]
  a=moment(f); b=mformula(rhist(f),n); expected=Fraction(1,8)+Fraction(7,8*4**h)
  assert a==b==expected
  print('Reference star',h,'P2',a)

# Additional exact checks. Full active + spectator coordinate spaces are used.
if __name__=='__main__':
 def f6(x): return int((x&7)==7)^int(((x>>3)&7)==7)
 for name,B,eta_expected,p_expected in [('A',[1,8],6,Fraction(121,1024)),('B',[1,2],5,Fraction(23,128))]:
  f=[f6(x&63)^f6((x&63)^xorvec(B,x>>6)) for x in range(256)]
  hist=rhist(f); eta=8-(hist[0].bit_length()-1); p=moment(f)
  assert hist[0]&(hist[0]-1)==0
  assert eta==eta_expected and p==p_expected==mformula(hist,8)
  for t in range(4):
   assert not nonzero(f[t*64:(t+1)*64],[1<<i for i in range(6)])
  print('Reference orientation',name,'rank',dict(sorted(hist.items())),'nullity',eta,'P2',p)
 for s,mod in [(1,3),(2,7),(3,11)]:
  D=1<<s
  def mul(a,b):
   z=0
   while b:
    if b&1:z^=a
    b>>=1;a<<=1
    if a&D:a^=mod
   return z
  def tr(a):
   v=0
   for _ in range(s):v^=a;a=mul(a,a)
   assert v in (0,1)
   return v
  f=[tr(mul(mul(x&(D-1),(x>>s)&(D-1)),x>>(2*s))) for x in range(D**3)]
  hist=rhist(f);p=moment(f);expected=Fraction(1,D**3)*(1+Fraction(D**3-1,D**2))
  assert hist==Counter({0:1,2*s:D**3-1}) and p==expected==mformula(hist,3*s)
  if s>=2: assert p!=Fraction(11,32)**s
  print('GF(2^%d) trace'%s,'rank',dict(sorted(hist.items())),'P2',p)
 assert Fraction(16,9)**2>3
 assert Fraction(1,2)/(1-Fraction(9,16)**2)==Fraction(128,175)
 print('Heralded rank-extent and flagged fidelity constants: exact checks passed')
