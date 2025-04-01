using LinearAlgebra, SparseArrays

"""
    generate_spd(N::Integer, p::AbstractFloat) 

    Generate a linear system to test pcg method on 
    Square symmetric positive definite matrix 

    Returns: (A, x, b)
        where A*x = b 
"""
function generate_spd_system(N::Integer; p=0.0) 

    # start with a random array A 
    if p > 0
        A = sprand(N,N,p) 
    else 
        A = rand(N,N) 
    end

    # A = randn(n,n); A = A'*A; A = (A + A')/2
    A = A' * A 
    A = (A + A')/2 

    return A
end 