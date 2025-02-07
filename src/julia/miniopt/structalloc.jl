# Test why using structs sometimes allocs memory

mutable struct Model{S, T}

    nodes::AbstractArray{S}
    dofs::AbstractArray{T}

    function Model(nodes, dofs)

        new{Int32, Float32}(nodes, dofs)
    end 

end

using BenchmarkTools 
