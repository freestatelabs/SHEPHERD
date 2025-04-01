
"""
Solution Process 

Receive K, f in: 
    K*x = f 

Solve for x: 

    1. compute LDL decomposition of K: 
        K = L*D*L' 
    2. Solve for y, use forward-substitution:
        L*y = f 
    3. Solve for x, use back-substitution: 
        D*L' * x = y 

"""

include("ldl.jl")
include("substitutions.jl")

function solve_direct(K::AbstractArray, f::AbstractArray) 

    L, D = LDL_decomposition(K) 
    DLt = D*transpose(L)

    y = forward_substitution(L, f) 
    x = back_substiution(DLt, y)

    return x
end