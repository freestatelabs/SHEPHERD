
cd(@__DIR__)
include("../src/julia/shepherd.jl")

fn = "linear-2k.inp"
fixed_nodes = [3, 4, 6, 8, 91, 92, 93, 255, 256, 257, 340, 341, 342, 346, 347, 348, 
                1297, 1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305]

K, Kr, f, q = solve(fn; solver=" ", fixed_nodes=fixed_nodes);