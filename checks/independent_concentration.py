from collections import Counter
from itertools import combinations, permutations
from fractions import Fraction

def rank(rows):
    pivots = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)

def local_contraction(v):
    return (((v >> 2) & 1) << 1 | ((v >> 1) & 1) << 2,
            ((v >> 2) & 1) | (v & 1) << 2,
            ((v >> 1) & 1) | (v & 1) << 1)

expected = {
    3: {0: 1, 4: 63},
    4: {0: 1, 4: 49, 6: 294, 8: 168},
    5: {0: 1, 4: 70, 6: 420, 8: 1925, 10: 1680},
    6: {0: 1, 4: 105, 6: 840, 8: 4522, 10: 13860, 12: 13440},
}
for n in range(3, 7):
    blocks = n - 1
    hist = Counter()
    for label in range(1 << (3 * blocks)):
        values = [(label >> (3 * i)) & 7 for i in range(blocks)]
        total = 0
        for value in values:
            total ^= value
        shared = local_contraction(total)
        repeated = [sum(row << (3 * i) for i in range(blocks)) for row in shared]
        rows = []
        for i, value in enumerate(values):
            rows.extend(repeated[j] ^ (row << (3 * i))
                        for j, row in enumerate(local_contraction(value)))
        hist[rank(rows)] += 1
    assert hist == expected[n]
    print('N =', n, 'contraction histogram', dict(sorted(hist.items())))

planes = set()
for a,b,c in combinations(range(1,64),3):
    if rank([a,b,c]) == 3:
        planes.add(tuple(sorted((0,a,b,c,a^b,a^c,b^c,a^b^c))))
assert len(planes) == 1395
monomials = [sum(1 << (3 * column + row)
                 for column,row in enumerate(perm))
             for perm in permutations(range(3))]
def f(label):
    return sum((label & monomial) == monomial for monomial in monomials) & 1
def trilinear(a,b,c):
    return sum(f(t) for t in (0,a,b,c,a^b,a^c,b^c,a^b^c)) & 1
normal_ranks = Counter()
for plane in planes:
    basis = []
    for a in plane:
        if rank(basis + [a]) > len(basis):
            basis.append(a)
    pair_list = list(combinations(range(3),2))
    matrix = [sum(trilinear(1 << i, basis[j] << 3, basis[k] << 3) << position
                  for position,(j,k) in enumerate(pair_list)) for i in range(3)]
    normal_ranks[rank(matrix)] += 1
assert normal_ranks == {0:3,1:168,2:840,3:384}
print('Determinant three-plane normal-form ranks', dict(sorted(normal_ranks.items())))
