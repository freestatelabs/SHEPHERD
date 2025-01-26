"""
"""


using SparseArrays, BenchmarkTools
cd(@__DIR__)
include("types.jl")
include("elem/3d/C3D8.jl")
include("io.jl")

function assemble(nodes, dofs, elements)

    Nnodes = size(nodes)[1]
    Nelems = size(elements)[1]
    Ndof = 3*Nnodes
    K = zeros(Ndof, Ndof)
    C = Cmatrix(200e3, 0.3)

    # Initialize values for stiffness matrix 
    Ke = zeros(24,24)
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

        K_C3D8!(Ke, enodes, C, _Ke, dShape, J, Jinv, B, aux, BtC, Bt)

        for i in 1:24
            for j in 1:24
                K[edofs[i],edofs[j]] += Ke[i,j] 
            end
        end

    end

    return sparse(K)
end


function applybcs!(K, nodes)

    dofs = zeros(Int64, 3*length(nodes))

    for i in eachindex(nodes)
        dofs[(i-1)*3+1:(i-1)*3+3] .= [1,2,3] .+ 3*(nodes[i] - 1)
    end

    for dof in dofs 
        K[dof,:] .= 0.0 
        K[:,dof] .= 0.0 
    end

end

# Test 

# nodes = convert.(Float64, [
#     0 0 0;
#     1 0 0; 
#     1 1 0; 
#     0 1 0;
#     0 0 1; 
#     1 0 1;
#     1 1 1;
#     0 1 1
# ])

# nodes = hcat(nodes, [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21; 22 23 24])
# elements = [1 2 3 4 5 6 7 8]
# model = Model() 
# model.nodes = nodes 
# model.elements = elements 
# K = assemble(model)

# fn = "../../test/linear-2k.inp"
# fn = "../../test/single-elem.inp"
# nodes, dofs, elements = readinputfile(fn)
# # # K = assemble(model)
# K = assemble(nodes, dofs, elements)