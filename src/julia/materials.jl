""" SHEPHERD
    Material Library
"""

# Material supertype
abstract type Material end 


# Default material
struct NullMaterial <: Material 

end 


"""
    mutable struct LinearElastic <: Material 
"""
mutable struct LinearElastic <: Material
    E::Number 
    nu::Number 
end 


"""
    Cmatrix(mat::LinearElastic; dim=Dim3D)

Compute the constitutive matrix [C] for a linear-elastic material.
"""
function Cmatrix(E::Number, nu::Number; dim=Dim3D)

    if dim == Dim3D
        a = E / ((1 + nu) * (1 - 2nu))

        C = zeros(6,6)
        C[1,1] = C[2,2] = C[3,3] = 1-nu 
        C[4,4] = C[5,5] = C[6,6] = 0.5*(1-2nu)
        C[1,2] = C[1,3] = C[2,3] = C[2,1] = C[3,1] = C[3,2] = nu

        return a .* C
        
    elseif dim == Dim2D
        return (E/(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2]

    else 
        return E

    end

end


