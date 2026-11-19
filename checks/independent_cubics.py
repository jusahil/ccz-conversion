from functools import cache
from itertools import combinations, product
from collections import Counter
from fractions import Fraction

def span_basis(vs):
    pivots = {}
    for v in vs:
        while v:
            pivot = v.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = v
                break
            v ^= pivots[pivot]
    return list(pivots.values())

def elements(basis):
    result = [0]
    for a in basis:
        result += [v ^ a for v in result]
    return result

@cache
def tensor(monomials, a, b, c):
    def f(x):
        return sum((x & m) == m for m in monomials) % 2
    return sum(f(x) for x in (0, a, b, c, a^b, a^c, b^c, a^b^c)) % 2

def cubic(monomials, coordinates):
    if len(coordinates) < 3:
        return ()
    n = len(coordinates[0])
    return tuple(sum(1 << i for i in triple)
                 for triple in combinations(range(len(coordinates)), 3)
                 if sum(tensor(monomials, coordinates[triple[0]][copy], coordinates[triple[1]][copy],
                               coordinates[triple[2]][copy]) for copy in range(n)) % 2)

def kernel_basis(n, w, k):
    extra = []
    curr = list(w)
    for i in range(n):
        candidate = 1 << i
        if candidate not in elements(curr):
            curr.append(candidate)
            extra.append(candidate)
    return [tuple(v if i == j else 0 for j in range(k)) for v in w for i in range(k)] + [
        tuple(v if (j == i or j == k - 1) else 0 for j in range(k)) for v in extra for i in range(k - 1)]

subs = {frozenset([0])}
for _ in range(4):
    subs |= {w | frozenset(x ^ a for x in w) for w in subs for a in range(16) if a not in w}
assert len(subs) == 67
cube4 = tuple(sum(1 << i for i in triple) for triple in combinations(range(4), 3))
counts = Counter()
for coeffs, sub in product(product((False, True), repeat=4), subs):
    M = tuple(v for v, use in zip(cube4, coeffs) if use)
    w = span_basis(sorted(sub))
    pred = (0 if not M else
            1 if any(tensor(M, a, b, c) for a, b, c in combinations(w, 3)) else
            2 if any(tensor(M, a, 1 << b, 1 << c) for a in w for b, c in combinations(range(4), 2)) else 3)
    exact = next((N for N in (1, 2, 3) if cubic(M, kernel_basis(4, w, N))), 0)
    assert pred == exact
    counts[exact] += 1
print('Copy classification:', dict(sorted(counts.items())))

def compress(M):
    used = 0
    for v in M:
        used |= v
    active = [i for i in range(used.bit_length()) if used & (1 << i)]
    return len(active), tuple(sum(1 << i for i, index in enumerate(active) if v & (1 << index)) for v in M)

def stats(m, M):
    rows = Counter()
    for a in range(2**m):
        B = [sum(tensor(M, a, 1 << i, 1 << j) << j for j in range(m)) for i in range(m)]
        rows[len(span_basis(B))] += 1
    rank_purity = sum(Fraction(count, 2**(m+r)) for r, count in rows.items())
    phase = [(-1)**sum((x & t) == t for t in M) for x in range(1 << m)]
    accum = 0
    for a in range(1 << m):
        walsh = [phase[x] * phase[x ^ a] for x in range(1 << m)]
        for k in range(m):
            stride = 1 << k
            for start in range(0, 1 << m, stride * 2):
                for j in range(start, start + stride):
                    left, right = walsh[j], walsh[j + stride]
                    walsh[j], walsh[j + stride] = left + right, left - right
        accum += sum(val**4 for val in walsh)
    assert rank_purity == Fraction(accum, 2**(5*m))
    return m, sorted(rows.items()), rank_purity

for state, M in zip('ABC', ([7], [11], [56])):
    for N in (1, 2, 3):
        print(state, N, *stats(*compress(cubic(tuple(M), kernel_basis(6, [1,2,4], N)))))

D=3
M=tuple(7 << (3*i) for i in range(D))
r=sum(1 << (3*i) for i in range(D))
B=cubic(M, [(1 << i, 1 << i) for i in range(3*D)] + [(0, r)])
m, B = compress(B)
print('two copies and one shared reference with disjoint cubic blocks:', D, m, B, stats(m, B)[1])

sub_bases = {w: span_basis(sorted(w)) for w in subs}
checked = 0
for coeffs, sub in product(product((False, True), repeat=4), subs):
    M = tuple(v for v, use in zip(cube4, coeffs) if use)
    if not M:
        continue
    w = sub_bases[sub]
    possible = sorted((s for s in subs if sub <= s), key=len)
    d = next(j for j in range(4) if any(
            tensor(M, *ws, *es)
            for ws in combinations(w, 3-j)
            for es in combinations((1,2,4,8), j)))
    lam1 = next(len(sub_bases[s]) - len(w) for s in possible if
                cubic(M, [(a,) for a in sub_bases[s]]))
    lam2 = next(len(sub_bases[s]) - len(w) for s in possible if
                cubic(M, [(a,0) for a in sub_bases[s]] + [(a,a) for a in (1,2,4,8)]))
    assert (lam1, lam2) == (d, int(d == 3))
    checked += 1
print('Reference frontiers checked:', checked)

X, Y, Z, I = 0, 1, 2, 3
# Each bit-pair gives a displacement and a character on three coordinates.
non_id_magnitudes = set()
for a,b in product(range(8), repeat=2):
    if not (a or b):
        continue
    non_id_magnitudes.add(abs(sum((-1)**((x==7)^((x^a)==7)^((b & x).bit_count()%2))
                                  for x in range(8)))/8)
assert non_id_magnitudes == {0,0.5}
print('CCZ nonidentity Pauli magnitudes:', sorted(non_id_magnitudes))
