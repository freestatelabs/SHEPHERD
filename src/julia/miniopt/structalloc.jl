# Test why using structs sometimes allocs memory

using Printf


precision = 32


struct Model{T_int, T_float}

    # Nnodes x 3: x,y,z coordinates of nodes
    nodes::Matrix{T_float}

    # Nelems x 8: node numbers in each element
    elements::Matrix{T_int}

    # Linear elastic material: [E, nu]
    material::Vector{T_float}

    # N_fixednodes: node numbers
    fixed_nodes::Vector{T_int}

    # N_forces: (node, direction, value)        
    forces::Vector{Tuple{T_int, T_int, T_float}}

    # # Free DOFs: Nfree x 3: [dof #, node #, direction]
    # free_dofs::Matrix{T_int}


    function Model(nodes, elements, material, fixed_nodes, forces, free_dofs)

        if precision == 32
            for i in 1:length(forces)
                forces[i] = (Int32(forces[i][1]), Int32(forces[i][2]), Float32(forces[i][3]))
            end
            new{Int32, Float32}(Float32.(nodes), Int32.(elements), Float32.(material), Int32.(fixed_nodes), forces)

        elseif precision == 64
            for i in 1:length(forces)
                forces[i] = (Int64(forces[i][1]), Int64(forces[i][2]), Float64(forces[i][3]))
            end
            new{Int64, Float64}(Float64.(nodes), Int64.(elements), Float64.(material), Int64.(fixed_nodes), forces)
        else
            @printf "Error. Precision `%i` not valid.\n" precision 
        end
    end 

end

using BenchmarkTools 

nodes = [1.0 2.0 3.0; 4.0 5.0 6.0]
elements = [1 2 3 4 5 6 7 8; 9 10 11 12 13 14 15 16]
material = [195e3, 0.3]
fixed_nodes = [1,2,3,4,5]
forces = [(1,1,3.0), (2,2,6)]
free_dofs = [1,2]

model = Model(nodes, elements, material, fixed_nodes, forces, free_dofs)


function f(model::Model{Int32, Float32})
    x = model.material[1]
    x = model.material[2]
    x = model.material[1]
    x = model.material[2]

    y = elements[1,5]
    y = elements[2,4]
    y = elements[1,5]
    y = elements[2,4]

end 


@benchmark f($model)


# struct mystruct
#     arr::AbstractArray
# end

# function f(s::mystruct)

#     x = s.arr[1]
#     x = s.arr[2]
#     x = s.arr[1] 
#     x = s.arr[2]

# end

# s = mystruct([1,2,3,4])

# @benchmark f($s)