

mutable struct Cloads 
    dofs::AbstractArray{Int64}
    forces::AbstractArray{Float64}

    function Cloads()
        new()
    end

    function Cloads(dofs::Vector{<:Integer}, forces::Vector{<:AbstractFloat})
        new(dofs, forces)
    end
end

mutable struct Model 
    nodes::AbstractArray{Float64}
    dofs::AbstractArray{Int64}
    elements::AbstractArray{Int64}
    cloads::Cloads
    
    function Model()
        new()
    end
end

