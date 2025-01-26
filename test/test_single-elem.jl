
include("../src/julia/shepherd.jl")

fn = "two_elems.inp"
constraints = [1, 2, 3, 10, 12, 13, 14, 22] 

K, f, q = solve(fn; constraints=constraints)