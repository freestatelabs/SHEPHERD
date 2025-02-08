"""
"""

using SparseArrays


function reducedofs(dofs::AbstractArray{<:Integer}, constrained_dofs::AbstractArray{<:Integer}, 
                    cload_dofs, cload_forces)

    Ndofs = length(dofs)
    Nconstrained = length(constrained_dofs)
    Nreduced = Ndofs - Nconstrained
    reduced_dofs = zeros(Int64, Nreduced)
    f = zeros(Float64, Nreduced)
    i = 1 

        # Define loads 
        f = zeros(length(dofs))
        for i in eachindex(cload_dofs)
            f[cload_dofs[i]] = cload_forces[i]
        end

    for j in 1:Ndofs 
        if !(dofs[j] in constrained_dofs) 
            reduced_dofs[i] = dofs[j] 
            f[i] = cload_forces[cload_dofs[]]
            i += 1 
        end 
    end 

    return reduced_dofs
end


"""
    assemble(nodes::AbstractArray{AbstractFloat}, dofs::AbstractArray{Integer}, 
                    elements::AbstractArray{Integer}, E::Number, nu::Number)

Assemble the global stiffness matrix for a linear-elastic finite-element mesh.
"""
function assemble(nodes::AbstractArray{<:AbstractFloat}, dofs::AbstractArray{<:Integer}, 
                    elements::AbstractArray{<:Integer}, dof_index::AbstractArray{<:Integer}, 
                    E::Number, nu::Number)

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
                Kglobal[edofs[i],edofs[j]] += Ke[i,j] 
            end
        end

    end

    return sparse(Kglobal)
end



"""
    applyfixedbcs!(K, nodes)

    Apply fixed boundary conditions to a global stiffness matrix by
    zeroing out the rows and columns of each fixed DOF.
"""
function applyfixedbcs!(K, nodes)

    dofs = zeros(Int64, 3*length(nodes))

    for i in eachindex(nodes)
        dofs[(i-1)*3+1:(i-1)*3+3] .= [1,2,3] .+ 3*(nodes[i] - 1)
    end

    for dof in dofs 
        K[dof,:] .= 0.0 
        K[:,dof] .= 0.0 
    end

    return dofs
end