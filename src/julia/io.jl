

using Printf

""" 
    scaninputfile(fn::AbstractString)

    Scans an input file to determine metadata
"""
function scaninputfile(fn::AbstractString)

    readingnodes = false 
    readingelems = false 
    readingcloads = false
    Nnodes = 0 
    Nelems = 0
    Ncloads = 0
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

        if (length(lines[i]) == 6) && (lines[i][1:6] == "*Cload")
            readingcloads = true 
            continue
        end 
        if (readingcloads) && (lines[i][1] == '*')
            readingcloads = false 
        end 

        # Count nodes 
        if readingnodes 
            Nnodes += 1 
        end 

        # Count elements 
        if readingelems 
            Nelems += 1
        end

        # Count Cloads 
        if readingcloads 
            Ncloads += 1
        end
    end

    return Nnodes, Nelems, Ncloads
end

"""
    readinputfile(fn::AbstractString)

Reads a Calculix input file. 

Currently supports only the following keywords: 
    *Node 
    *Element
    *Elastic
    *Cload
# Returns 


143ms
"""
function readinputfile(fn::AbstractString; verbose = false)

    if verbose
        @printf "Reading file: '%s'.\n" fn
    end

    Nnodes, Nelems, Ncloads = scaninputfile(fn)
    @printf "\tNnodes: %i, Nelems: %i, Ncloads: %i\n" Nnodes Nelems Ncloads
    #Nnodes = 2025 
    #Nelems = 1280
    #Ncloads = 25
    Ndofs = 3*Nnodes

    nodes = zeros(Float64, Nnodes, 3)           # x, y, z, dof1, dof2, dof3     
    dofs = zeros(Int64, Nnodes, 3)              # dof1, dof2, dof3
    elements = zeros(Int64, Nelems, 8)          # Only supports C3D8
    cload_dofs = zeros(Int64, Ncloads)
    cload_forces = zeros(Float64, Ncloads)

    n = 1

    lines = readlines(fn) 
    imax = length(lines)
    i = 1   # line number 
    d = 1
    E = 0.0
    nu = 0.0 

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
                nodes[j,1:3] = parse.(Float64, split(lines[i+j], ", "))[2:end]
                dofs[j,:] = [d, d+1, d+2]
                d += 3
            end
            i += Nnodes
        end

        #
        # Read elem 
        #
        if (length(lines[i]) > 8) && lines[i][1:8] == "*Element"
            for j in 1:Nelems
                # Read an element from the file 
                elements[j,:] = parse.(Int64, [x for x in split(lines[i+j], ", ")])[2:end] 
            end
            i += Nelems
        end 

        # Read Cload
        if (length(lines[i]) == 6) && lines[i][1:6] == "*Cload"
            for j in 1:Ncloads
                # Read a load from the file 
                line = parse.(Float64, [x for x in split(lines[i+j], ", ")])
                node = convert(Int64, line[1]) 
                k = convert(Int64, line[2]) 
                # f = line[3]
                cload_dofs[j] = (node - 1)*3 + k 
                cload_forces[j] = line[3]
            end
            i += Ncloads
        end 

        if (length(lines[i]) == 8) && (lines[i] == "*Elastic")
            i += 1
            line = parse.(Float64, [x for x in split(lines[i], ", ")])
            E = line[1]
            nu = line[2] 
        end

        i += 1
    end

    return nodes, dofs, elements, cload_dofs, cload_forces, E, nu
end


# fn = "/Users/ryan/Github/SHEPHERD/test/linear-2k.inp"
# model = readinputfile(fn);
# println()

# # It's pretty slow, lots of allocs
# using BenchmarkTools
# @btime readinputfile($fn)
