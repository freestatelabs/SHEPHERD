using SparseArrays, LinearAlgebra 

# Create a random SPD array 
N = 100_000
p = 0.001

A = sprand(N, N, p)
A = A' * A
A = (A + A')/2

# Use SuiteSparse backend to calculate cholesky decomp
c = cholesky(A) 

