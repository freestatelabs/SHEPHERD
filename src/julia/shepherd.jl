""" SHEPHERD

"""

# module Shepherd 

PRECISION = 32


function SETPRECISION(PRECISION::Integer) 
    global INT_PRECISION, FLOAT_PRECISION
    if PRECISION == 32 
        INT_PRECISION = Int32 
        FLOAT_PRECISION = Float32
    else 
        INT_PRECISION = Int64 
        FLOAT_PRECISION = Float64 
    end
end

export SETPRECISION

SETPRECISION(PRECISION)


include("types.jl")
include("utils.jl")
include("materials.jl")
include("io.jl")
include("elem/elem.jl")
include("assemble.jl")
include("solve.jl")

# end # module