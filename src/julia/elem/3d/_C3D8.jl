
"""
    K_C3D8!(Kelem::AbstractArray, nodes, dofs, C)

Calculate the 24x24 stiffness matrix for a C3D8 element

# Arguments
- Kelem: 24x24 array to modify 
- nodes: 8x3 array of node coordinates
- dofs: 
"""
function K_C3D8!(Kelem::AbstractArray, nodes::AbstractArray, C::AbstractArray)

    gauss_pts = [-1/sqrt(3), 1/sqrt(3)]

    for xi1 in gauss_pts 
        for xi2 in gauss_pts 
            for xi3 in gauss_pts

                # Shape function derivatives (3x8 matrix)
                dshape = (1/8)*[-(1-xi2)*(1-xi3) (1-xi2)*(1-xi3) (1+xi2)*(1-xi3) -(1+xi2)*(1-xi3) -(1-xi2)*(1+xi3) (1-xi2)*(1+xi3) (1+xi2)*(1+xi3) -(1+xi2)*(1+xi3);
                                -(1-xi1)*(1-xi3) -(1+xi1)*(1-xi3) (1+xi1)*(1-xi3) (1-xi1)*(1-xi3) -(1-xi1)*(1+xi3) -(1+xi1)*(1+xi3) (1+xi1)*(1+xi3) (1-xi1)*(1+xi3);
                                -(1-xi1)*(1-xi2) -(1+xi1)*(1-xi2) -(1+xi1)*(1+xi2) -(1-xi1)*(1+xi2) (1-xi1)*(1-xi2) (1+xi1)*(1-xi2) (1+xi1)*(1+xi2) (1-xi1)*(1+xi2)]
                
                jac = dshape * nodes 
                aux = inv(jac) * dshape 

                B = zeros(6, 24) 

                for i in 1:3 
                    for j in 0:7 
                        B[i, 3*j+1*(i-1)] = aux[i, j+1]
                    end 
                end 

                for j in 0:7 
                    B[4, 3*j+1] = aux[2, j+1]
                    B[4, 3*j+2] = aux[1, j+1]
                    B[5, 3*j+3] = aux[2, j+1]
                    B[5, 3*j+2] = aux[3, j+1]
                    B[6, 3*j+1] = aux[3, j+1]
                    B[6, 3*j+3] = aux[1, j+1]
                end 

                Kelem += B' * C * B * det(jac)
            end
        end
    end

end