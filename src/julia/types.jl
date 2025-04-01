
mutable struct Nodes{INT_PRECISION, FP_PRECISION}
    numbers::Vector{INT_PRECISION}
    coordinates::Matrix{FP_PRECISION}
end

abstract type ElementSet end

mutable struct C3D8_Set{INT_PRECISION}<:ElementSet 
    material::String 
    element_numbers::Vector{INT_PRECISION}
    element_nodes::Matrix{INT_PRECISION}
end


mutable struct Model{INT_PRECISION, FP_PRECISION}
    nodes::Nodes{INT_PRECISION, FP_PRECISION}
    element_sets::Vector{ElementSet}
    # dofs::Matrix{Int32}
    # elements::Matrix{Int32}
    # Nsets::Dict 
    # materials::Dict 
    # boundaries::Dict 
    
    function Model()
        new()
    end
end