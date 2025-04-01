
fn = "/Users/ryan/Github/SHEPHERD/test/linear-2k.inp"
include("types.jl")

using Printf

""" 
    scaninputfile(fn::AbstractString)

    Scans an input file to determine metadata
"""
function scaninputfile(fn::AbstractString)

    readingnodes = false 
    readingelems = false 
    Nnodes = 0 
    Nelems = 0
    lines =readlines(fn)

    for i in eachindex(lines)

        # Switch flags
        if (length(lines[i]) == 5) && (lines[i][1:5] == "*Node")
            readingnodes = true 
            continue
        end
        if (readingnodes) && (lines[i][1] == '*')
            readingnodes = false 
        end

        if (length(lines[i]) >= 8) && (lines[i][1:8] == "*Element")
            readingelems = true 
            continue
        end 
        if (readingelems) && (lines[i][1] == '*')
            readingelems = false 
        end 

        # Count nodes 
        if readingnodes 
            Nnodes += 1 
        end 

        # Count elements 
        if readingelems 
            Nelems += 1
        end

    end

    return Nnodes, Nelems
end

"""
    readinputfile(fn::AbstractString)

Reads a Calculix input file. 

Currently supports only the following keywords: 
    *Node 
    *Element
    *Material
# Returns 
(Dict, Dict, Dict) corresponding to nodes, elems, dofs

143ms
"""
function readinputfile(fn::AbstractString; verbose = false)

    if verbose
        @printf "Reading file: '%s'.\n" fn
    end

    Nnodes, Nelems = scaninputfile(fn)
    #Nnodes = 2025 
    #Nelems = 1280
    Ndofs = 3*Nnodes

    nodes = Matrix{Float32}(undef, Nnodes, 6)   # x, y, z, dof1, dof2, dof3     
    dofs= Matrix{Int32}(undef, Ndofs,2)  # node, direction
    elements = Matrix{Int32}(undef, Nelems, 8)       # Only supports C3D8
    dofs[1,:] = [1,1]

    model = Model()
    n = 1

    lines = readlines(fn) 
    imax = length(lines)
    i = 1   # line number 
    d = 1
    while i <= imax

        # Comment line (saves time if we just skip to the next loop)
        if lines[i][1:2] == "**"
            i += 1
            continue 
        end
             
        # 
        # Read nodes
        # 
        if (length(lines[i]) == 5) && (lines[i][1:5] == "*Node")
            for j in 1:Nnodes
                nodes[j,1:3] = parse.(Float32, split(lines[i+j], ", "))[2:end]
                dofs[d,:] = [j,1]
                dofs[d+1,:] = [j,2]
                dofs[d+2,:] = [j,3] 
                d+= 3
            end
            i += Nnodes
        end

        #
        # Read elem 
        #
        if (length(lines[i]) > 8) && lines[i][1:8] == "*Element"
            for j in 1:Nelems
                # Read an element from the file 
                elements[j,:] = parse.(Int32, [x for x in split(lines[i+j], ", ")])[2:end] 
            end
            i += Nelems
        end 

        i += 1
    end

    model.nodes = nodes 
    model.elements = elements 
    model.dofs = dofs
    return model
end

model = readinputfile(fn);
println()

# It's pretty slow, lots of allocs
using BenchmarkTools
@btime readinputfile($fn)


"""
    loadinputfile(fn::String)

    Load a Calculix-formatted input file into memory.
"""
function loadinputfile(fn::String)

    model = Model()

    lines = readlines(fn)
    Nlines = length(lines)

    line_number = 1 

    while line_number <= Nlines 

        if length(lines[line_number]) >= 5 && lines[line_number] == "*Node" 
            
        end

        line_number += 1 
    end

end

