"""Minimal examples of accessible cubic polarization and the aligned frontier."""
from src.ccz_research import has_accessible_cubic, aligned_yield
print('Full cubic space:',has_accessible_cubic(3,[7],[1,2,4]))
print('Two-dimensional restriction:',has_accessible_cubic(3,[7],[1,2]))
print('Four aligned copies; reference dimensions 0,1,2,3:',[aligned_yield(4,e) for e in range(4)])
