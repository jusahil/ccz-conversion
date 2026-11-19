"""Exact checks of charge alignment and coordinate-aligned Wilson restrictions."""
from itertools import product, permutations
def f(x): return int(x & 7 == 7)
for x in range(16):
    y=(x & 3) | ((x & 4)<<1) | ((x & 8)>>1)
    assert f(x)^f(y)==int(x&3==3)*(((x>>2)^(x>>3))&1)
def det(a,b,c):
    return sum(((a>>i)&1)*((b>>j)&1)*((c>>k)&1) for i,j,k in permutations(range(3)))%2
counts=[0,0]
for m in product(range(8),repeat=3):
    spaces=[[x for x in range(8) if x&w==x] for w in m]
    a=any(det(*v) for v in product(*spaces))
    b=any(all(m[r]>>p[r]&1 for r in range(3)) for p in permutations(range(3)))
    assert a==b;counts[int(a)]+=1
assert counts==[265,247]
print('16 charge-relabeling assignments; 512 coordinate restrictions;',counts)
