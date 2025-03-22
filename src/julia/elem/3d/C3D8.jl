using LinearAlgebra

const gauss_pts = [-1/sqrt(3), 1/sqrt(3)]

const dShape111 = [ -0.311004   0.311004    0.0833334  -0.0833334  -0.0833334   0.0833334  0.0223291  -0.0223291;
                    -0.31108   -0.0832575   0.0832575   0.31108     0.0833537  -0.0223088  0.0223088   0.0833537;
                    -0.31108   -0.0832575  -0.0223088  -0.0833537   0.31108     0.0832575  0.0223088   0.0833537]

"""
    det3x3(A::AbstractArray)

Compute the determinant of a 3x3 matrix.
Roughly 30x faster than `det()` with no allocations 
"""
function det3x3(A::AbstractArray)
    return A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1]
end

"""
    inv3x3!(Ainv::AbstractArray, A::AbstractArray)

Compute the inverse of a 3x3 matrix; store in `Ainv`.
"""
function inv3x3!(Ainv::AbstractArray, A::AbstractArray) 
    # Roughly 45x faster than `inv()` with no allocations

    detA = det3x3(A)

    Ainv[1,1] = (A[2,2]*A[3,3] - A[2,3]*A[3,2])/detA
    Ainv[1,2] = (-A[1,2]*A[3,3] + A[1,3]*A[3,2])/detA
    Ainv[1,3] = (A[1,2]*A[2,3] - A[1,3]*A[2,2])/detA
    Ainv[2,1] = (-A[2,1]*A[3,3] + A[2,3]*A[3,1])/detA
    Ainv[2,2] = (A[1,1]*A[3,3] - A[1,3]*A[3,1])/detA
    Ainv[2,3] = (-A[1,1]*A[2,3] + A[1,3]*A[2,1])/detA
    Ainv[3,1] = (A[2,1]*A[3,2] - A[2,2]*A[3,1])/detA
    Ainv[3,2] = (-A[1,1]*A[3,2] + A[1,2]*A[3,1])/detA
    Ainv[3,3] = (A[1,1]*A[2,2] - A[1,2]*A[2,1])/detA

    return detA
end


function update_dShape!(dShape, xi1, xi2, xi3) 
# dShape .= (1/8) .*[
#     -(1-xi2)*(1-xi3) (1-xi2)*(1-xi3) (1+xi2)*(1-xi3) -(1+xi2)*(1-xi3) -(1-xi2)*(1+xi3) (1-xi2)*(1+xi3) (1+xi2)*(1+xi3) -(1+xi2)*(1+xi3);
#     -(1-xi1)*(1-xi3) -(1+xi1)*(1-xi3) (1+xi1)*(1-xi3) (1-xi1)*(1-xi3) -(1-xi1)*(1+xi3) -(1+xi1)*(1+xi3) (1+xi1)*(1+xi3) (1-xi1)*(1+xi3);
#     -(1-xi1)*(1-xi2) -(1+xi1)*(1-xi2) -(1+xi1)*(1+xi2) -(1-xi1)*(1+xi2) (1-xi1)*(1-xi2) (1+xi1)*(1-xi2) (1+xi1)*(1+xi2) (1-xi1)*(1+xi2)
# ]
    dShape[1,1] = -(1-xi2)*(1-xi3) 
    dShape[1,2] =  (1-xi2)*(1-xi3)
    dShape[1,3] =  (1+xi2)*(1-xi3)
    dShape[1,4] = -(1+xi2)*(1-xi3)
    dShape[1,5] = -(1-xi2)*(1+xi3)
    dShape[1,6] =  (1-xi2)*(1+xi3)
    dShape[1,7] =  (1+xi2)*(1+xi3)
    dShape[1,8] = -(1+xi2)*(1+xi3)
    dShape[2,1] = -(1-xi1)*(1-xi3)
    dShape[2,2] = -(1+xi1)*(1-xi3)
    dShape[2,3] =  (1+xi1)*(1-xi3)
    dShape[2,4] =  (1-xi1)*(1-xi3)
    dShape[2,5] =  (1-xi1)*(1+xi3)
    dShape[2,6] = -(1+xi1)*(1+xi3)
    dShape[2,7] =  (1+xi1)*(1+xi3)
    dShape[2,8] =  (1-xi1)*(1+xi3)
    dShape[3,1] = -(1-xi1)*(1-xi2)
    dShape[3,2] = -(1+xi1)*(1-xi2)
    dShape[3,3] = -(1+xi1)*(1+xi2)
    dShape[3,4] = -(1-xi1)*(1+xi2)
    dShape[3,5] =  (1-xi1)*(1-xi2)
    dShape[3,6] =  (1+xi1)*(1-xi2)
    dShape[3,7] =  (1+xi1)*(1+xi2)
    dShape[3,8] =  (1-xi1)*(1+xi2)
    dShape .*= 0.125
end

function update_dShape2!(dShape, xi1, xi2, xi3)
    # This is ~3x faster than calculating every loop

    if xi1 == 1 
        if xi2 == 1 
            if xi3 == 1 
                dShape .= dShape111
            else
                dShape .= dShape111 
            end

        else 
            if xi3 == 1
                dShape .= dShape111
            else
                dShape .= dShape111
            end
        end
    else
        if xi2 == 1 
            if xi3 == 1 
                dShape .= dShape111
            else 
                dShape .= dShape111
            end

        else 
            if xi3 == 1
                dShape .= dShape111
            else
                dShape .= dShape111
            end
        end
    end

end


"""
    K_C3D8!(K::AbstractArray, nodes::AbstractArray, C::AbstractArray,
            _K::AbstractArray, dShape::AbstractArray, J::AbstractArray, Jinv::AbstractArray,
            B::AbstractArray, aux::AbstractArray, BtC::AbstractArray, Bt::AbstractArray)

Calculate the 24x24 stiffness matrix for a 8-node isoparametric hexahedral 
finite element, store in `K`.

K = int(int(int(B' * C * B * |J| dr ds dt)))
"""
function K_C3D8!(K::AbstractArray, nodes::AbstractArray, C::AbstractArray,
                    _K::AbstractArray, dShape::AbstractArray, J::AbstractArray, Jinv::AbstractArray,
                    B::AbstractArray, aux::AbstractArray, BtC::AbstractArray, Bt::AbstractArray)

    K .= 0.0
    detJ = 0.0

    for xi1 in gauss_pts 
        for xi2 in gauss_pts
            for xi3 in gauss_pts

                # dShape .= (1/8) .*[
                #     -(1-xi2)*(1-xi3) (1-xi2)*(1-xi3) (1+xi2)*(1-xi3) -(1+xi2)*(1-xi3) -(1-xi2)*(1+xi3) (1-xi2)*(1+xi3) (1+xi2)*(1+xi3) -(1+xi2)*(1+xi3);
                #     -(1-xi1)*(1-xi3) -(1+xi1)*(1-xi3) (1+xi1)*(1-xi3) (1-xi1)*(1-xi3) -(1-xi1)*(1+xi3) -(1+xi1)*(1+xi3) (1+xi1)*(1+xi3) (1-xi1)*(1+xi3);
                #     -(1-xi1)*(1-xi2) -(1+xi1)*(1-xi2) -(1+xi1)*(1+xi2) -(1-xi1)*(1+xi2) (1-xi1)*(1-xi2) (1+xi1)*(1-xi2) (1+xi1)*(1+xi2) (1-xi1)*(1+xi2)
                # ]
                update_dShape!(dShape, xi1, xi2, xi3)

                mul!(J, dShape, nodes)
                detJ = inv3x3!(Jinv, J)
                mul!(aux, Jinv, dShape)

                # First 3 rows are normal strain 
                for i = 1:3 
                    for j = 0:7 
                        B[i,3*j+1+(i-1)] = aux[i,j+1]
                    end
                end

                # Next 3 rows are shear strains
                for j = 0:7 
                    B[4,3*j+1] = aux[2,j+1]
                    B[4,3*j+2] = aux[1,j+1]
                    B[5,3*j+3] = aux[2,j+1]
                    B[5,3*j+2] = aux[3,j+1]
                    B[6,3*j+1] = aux[3,j+1]
                    B[6,3*j+3] = aux[1,j+1]
                end

                # K += B' * C * B * det(J)
                Bt .= B'       
                mul!(BtC, Bt, C)
                mul!(_K, BtC, B)
                K .+= _K .* detJ
            
            end
        end
    end
end



# Test different shape function derivative methods
# dShape = zeros(3,8)
# xi1 = -0.57735
# xi2 = -0.57735 
# xi3 = -0.57735
# @btime update_dShape!($dShape, $xi1, $xi2, $xi3) 
# xi1 = 1
# xi2 = 1 
# xi3 = 1 
# @btime update_dShape2!($dShape, $xi1, $xi2, $xi3)
