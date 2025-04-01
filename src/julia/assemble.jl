"""
"""

using SparseArrays


function assigndofs(model::Model)

    Ndofs = 3*size(model.nodes)[1]
    constrained_dofs = zeros(INT_PRECISION, 3*length(model.fixed_nodes))
    free_dofs = zeros(INT_PRECISION, Ndofs - length(constrained_dofs))
    i = 1

    for node in model.fixed_nodes 
        for dir in 1:3 
            constrained_dofs[i] = 3*(node - 1) + dir
            i += 1
        end
    end

    i = 1
    for j in 1:Ndofs
        if !(j in constrained_dofs) 
            free_dofs[i] = j
            i += 1
        end 
    end 

    return constrained_dofs, free_dofs
end


function assignforcedofs(model::Model)

    loaded_dofs = zeros(INT_PRECISION, length(model.forces.nodes))

    for i in eachindex(model.forces.nodes)
        loaded_dofs[i] = 3*(model.forces.nodes[i] - 1) + model.forces.directions[i]
    end

    return loaded_dofs
end


function reducesystem(model, K, free_dofs, force_dofs)
    println(free_dofs)

    Nfree = length(free_dofs) 
    Kr = zeros(FLOAT_PRECISION, Nfree, Nfree)
    F = zeros(FLOAT_PRECISION, Nfree) 

    for i in 1:Nfree 
        for j in 1:Nfree 
            Kr[i,j] = K[free_dofs[i],free_dofs[j]]
        end
    end

    for i in eachindex(force_dofs)
        F[i] = model.forces.magnitudes[i]
    end

    return Kr, F
end


"""
    assemble(nodes::AbstractArray{AbstractFloat}, dofs::AbstractArray{Integer}, 
                    elements::AbstractArray{Integer}, E::Number, nu::Number)

Assemble the global stiffness matrix for a linear-elastic finite-element mesh.
"""
function assemble(nodes::AbstractArray{<:AbstractFloat}, dofs::AbstractArray{<:Integer}, 
                    elements::AbstractArray{<:Integer}, E::Number, nu::Number)

    Nnodes = size(nodes)[1]
    Nelems = size(elements)[1]
    Ndof = 3*Nnodes
    Kglobal = zeros(Ndof, Ndof)
    C = Cmatrix(E, nu)

    # Initialize values for stiffness matrix 
    Kelement = zeros(24,24)
    dShape = zeros(3,8)
    J = zeros(3,3)
    B = zeros(6,24) 
    Bt = zeros(24,6)
    aux = zeros(3,8)
    BtC = zeros(24,6)
    _Ke = zeros(24,24)
    Jinv = zeros(3,3)

    enodes = zeros(8,3)
    edofs = zeros(Int64, 24)

    # Loop through all elements in model
    for e in 1:Nelems

        # Loop through all nodes in each element
        for i in 1:8
            enodes[i,:] .= nodes[elements[e,i],:]
            edofs[3*(i-1)+1:3*(i-1)+3] .= dofs[elements[e,i],:]
        end

        # Compute stiffness matrix for that element and store in `Kelement`
        K_C3D8!(Kelement, enodes, C, _Ke, dShape, J, Jinv, B, aux, BtC, Bt)

        # Add to the global stiffness matrix
        for i in 1:24
            for j in 1:24
                Kglobal[edofs[i],edofs[j]] += _Ke[i,j] 
            end
        end

    end

    return sparse(Kglobal)
end