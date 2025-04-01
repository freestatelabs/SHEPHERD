
using SparseArrays, LinearAlgebra 


"""
    LDL_decomposition(A::AbstractArray)

    Calculates the LDL decomposition of A: 
        A = L*D*L'
"""
function LDL_decomposition(A::AbstractArray) 

    N = size(A)[1]

    # allocate memory
    D = zeros(N)
    L = zeros(N,N)  # wasteful as only ~half entries are used

    # Compute the diagonal entries 
    # D(j) = A(j,j) - sum[k=1 to j-1](L(j,k)^2 * D(k))
    # L(i,j) = 1/D(j) * (A(i,j) - sum[k=1 to j-1](L(i,k) * L(j,k) * D(k))) (i>j)

    for j in 1:N 
        acc0 = 0.0
        for k in 1:j-1 
            acc0 += D[k] * L[j,k]^2
        end
        D[j] = A[j,j] - acc0 
        L[j,j] = 1

        for i = j+1:N
            acc0 = 0.0>
            for k in 1:j-1 
                acc0 += L[i,k] * L[j,k] * D[k] 
            end 
            L[i,j] = A[i,j]/D[j] - acc0 
        end
    end

    return D, L
end