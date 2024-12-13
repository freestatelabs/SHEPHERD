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


function Dmat(mat::LinearElastic, dim::Dim)


end

function Dmat(mat::LinearElastic, dim::Dim)

    if dim == Dim2D
        return (mat.E/(1-mat.nu^2)) * [1 mat.nu 0; mat.nu 1 0; 0 0 (1-mat.nu)/2]
    else 
        return mat.E
    end

end


