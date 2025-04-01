

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


mutable struct Forces{T_int, T_float}
    nodes::Vector{T_int}
    directions::Vector{T_int}
    magnitudes::Vector{T_float}

    function Forces(Nforces)
        new{INT_PRECISION, FLOAT_PRECISION}(zeros(INT_PRECISION, Nforces), zeros(INT_PRECISION, Nforces), zeros(FLOAT_PRECISION, Nforces))
    end 
end



mutable struct Model{T_int, T_float}

    # Nnodes x 3: x,y,z coordinates of nodes
    nodes::Matrix{T_float}

    # Nelems x 8: node numbers in each element
    elements::Matrix{T_int}

    # Linear elastic material: [E, nu]
    material::Vector{T_float}

    # N_fixednodes: node numbers
    fixed_nodes::Vector{T_int}

    # N_forces: (node, direction, value)        
    forces::Forces{T_int, T_float}

    #
    # Compute after initialization
    #

    # Free DOFs: Nfree x 3: [dof #, node #, direction]
    free_dofs::Matrix{T_int}

    # Fixed DOFs: Nfixed x 3: [dof #, node #, direction]
    fixed_dofs::Matrix{T_int}

    function Model(nodes, elements, material, fixed_nodes, forces)

        new{INT_PRECISION, FLOAT_PRECISION}(
            convert.(FLOAT_PRECISION, nodes), 
            convert.(INT_PRECISION, elements), 
            convert.(FLOAT_PRECISION, material), 
            convert.(INT_PRECISION, fixed_nodes), 
            forces)
    end
end

