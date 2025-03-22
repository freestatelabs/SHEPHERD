cd(@__DIR__)
include("../src/julia/shepherd.jl")

fn = "two-elem.inp"
constraints = [1, 2, 3, 10, 12, 13, 14, 22] 
fixed_nodes = [1,4,5,8]

K, Kr, f, q = solve(fn; fixed_nodes=fixed_nodes)

qr(Matrix(Kr), Val(true))\f