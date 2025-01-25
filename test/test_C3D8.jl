""" Test Functionality of C3D8 element
"""

cd(@__DIR__)
include("../src/julia/elem/3d/C3D8.jl")
using Printf

# Simple cube element with node 1 at (0,0,0), all other nodes in + xyz space
nodes = convert.(Float64, [
    0 0 0;
    1 0 0; 
    1 1 0; 
    0 1 0;
    0 0 1; 
    1 0 1;
    1 1 1;
    0 1 1
])

# Material properties
E = 200e3; nu = 0.3
C = Cmatrix(E, nu)

# Boundary conditions
# Minimally constrained at -x face
constraints = [1, 2, 3, 10, 12, 13, 14, 22]  

# Nodal force loads at +x face 
dx = 1
F = 0.25*(E*dx) .* [0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0]

# Output stiffness matrix 
K = zeros(24,24)

# Intermediate calculation matrices
dShape = zeros(3,8)
J = zeros(3,3)
B = zeros(6,24) 
Bt = zeros(24,6)
aux = zeros(3,8)
BtC = zeros(24,6)
_K = zeros(24,24)
Jinv = zeros(3,3)

# Calculate the stiffness matrix 
K_C3D8!(K, nodes, C, _K, dShape, J, Jinv, B, aux, BtC, Bt)

# Apply constraints
for constraint in constraints
    K[constraint,:] .= 0.0 
    K[:,constraint] .= 0.0 
end 

# Solve the system of equations
q = qr(K, Val(true)) \ F 

@printf "Results: \n"
for i in 1:8
  @printf "Node %i: " i
  for j in 1:3 
      @printf "%.3f " q[(i-1)*3+j]
  end
  @printf "\n"
end

report = @benchmark K_C3D8!($K, $nodes, $C, $_K, $dShape, $J, $Jinv, $B, $aux, $BtC, $Bt)
display(report)
# Should be zero allocs, median ~8.82 us on M1 Macbook Pro