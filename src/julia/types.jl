

mutable struct Model 
    nodes::Matrix{Float32}
    dofs::Matrix{Int32}
    elements::Matrix{Int32}
    Nsets::Dict 
    materials::Dict 
    boundaries::Dict 
    
    function Model()
        new()
    end
end