

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

    nodes = zeros(FLOAT_PRECISION, Nnodes, 3)           # x, y, z, dof1, dof2, dof3     
    elements = zeros(INT_PRECISION, Nelems, 8)          # Only supports C3D8
    material = [0.0, 0.0]
    fixed_nodes = [0]
    forces = Forces(Ncloads)

    E = 0.0
    nu = 0.0 

    lines = readlines(fn) 
    linenumber = 1  

    while linenumber <= length(lines)

        # Comment line (saves time if we just skip to the next loop)
        if lines[linenumber][1:2] == "**"
            linenumber += 1
            continue 
        end
             
        # 
        # Read nodes
        # 
        if (length(lines[linenumber]) == 5) && (lines[linenumber][1:5] == "*Node")
            for j in 1:Nnodes
                nodes[j,1:3] = parse.(FLOAT_PRECISION, split(lines[linenumber+j], ", "))[2:end]
            end
            linenumber += Nnodes
        end

        #
        # Read elem 
        #
        if (length(lines[linenumber]) > 8) && lines[linenumber][1:8] == "*Element"
            for j in 1:Nelems
                # Read an element from the file 
                elements[j,:] = parse.(INT_PRECISION, [x for x in split(lines[linenumber+j], ", ")])[2:end] 
            end
            linenumber += Nelems
        end 

        # Read Cload
        if (length(lines[linenumber]) == 6) && lines[linenumber][1:6] == "*Cload"
            for j in 1:Ncloads
                # Read a load from the file 
                line = parse.(FLOAT_PRECISION, [x for x in split(lines[linenumber+j], ", ")])

                forces.nodes[j] = convert(INT_PRECISION, line[1]) 
                forces.directions[j] = convert(INT_PRECISION, line[2]) 
                forces.magnitudes[j] = convert(FLOAT_PRECISION, line[3])
            end
            linenumber += Ncloads
        end 

        if (length(lines[linenumber]) == 8) && (lines[linenumber] == "*Elastic")
            linenumber += 1
            line = parse.(FLOAT_PRECISION, [x for x in split(lines[linenumber], ", ")])
            E = line[1]
            nu = line[2] 
        end

        linenumber += 1
    end



    return Model(nodes, elements, material, fixed_nodes, forces)
end


# fn = "/Users/ryan/Github/SHEPHERD/test/linear-2k.inp"
# model = readinputfile(fn);
# println()

# # It's pretty slow, lots of allocs
# using BenchmarkTools
# @btime readinputfile($fn)
