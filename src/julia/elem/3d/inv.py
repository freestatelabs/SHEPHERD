""" Use sympy to print out the code for [A]^-1, where [A] is dim(3,3)
"""

from sympy import symbols
from sympy.matrices import Matrix

a11, a12, a13, a21, a22, a23, a31, a32, a33 = symbols('a11 a12 a13 a21 a22 a23 a31 a32 a33')

A = Matrix([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
print("Inverse:\n")
print(A.inv())

print("\nDeterminant:\n")
print(A.det())