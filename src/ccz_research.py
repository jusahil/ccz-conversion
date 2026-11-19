"""Small public interface: exact cubic witnesses and analytic aligned-copy yield."""
from itertools import combinations
from .algebra import basis

def has_accessible_cubic(n, monomials, generators):
    """Whether the Boolean cubic polarization restricts nontrivially to the span.

    Masks encode binary vectors and monomials; repeated terms cancel in GF(2).
    Evaluates eight phase values per independent triple, without allocating 2**n.
    This is a witness test under the paper's primitive-operation assumptions.
    """
    if type(n) is not int or n < 0: raise ValueError('n must be a nonnegative integer')
    terms, vectors = list(monomials), list(generators)
    if any(type(x) is not int or x < 0 or x.bit_length() > n for x in terms + vectors):
        raise ValueError('mask outside dimension')
    if any(x.bit_count() > 3 for x in terms): raise ValueError('phase degree exceeds three')
    def f(x): return sum(x&t==t for t in terms)&1
    return any(sum(f(x) for x in (0,a,b,c,a^b,a^c,b^c,a^b^c))&1 for a,b,c in combinations(basis(vectors),3))

def aligned_yield(n, ell):
    """Evaluate the analytical frontier for N >= 1 and 0 <= ell <= 3."""
    if type(n) is not int or n<1 or type(ell) is not int or not 0<=ell<=3: raise ValueError('require N >= 1 and integer 0 <= ell <= 3')
    return max(0,n-2) if ell==0 else n-1 if ell<3 else n
