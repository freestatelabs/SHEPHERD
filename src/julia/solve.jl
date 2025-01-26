""" SHEPHERD

    Solve a finite element problem
"""


using Printf, Pardiso


"""
    solve(fn::AbstractString)

Solve a finite element problem by loading a Calculix input file.
"""
function solve(fn::AbstractString; fixed_nodes=[], constraints=[], verbose=false)
    
    # Read the input file, load data structures
    # Needs to support *NODE, *ELEMENT, *MATERIAL, *Cload

    if verbose @printf "\nReading input file from '%s' \n" fn end

    tstart = time()
    nodes, dofs, elements, cload_dofs, cload_forces, E, nu = readinputfile(fn)
    @printf "File read elapsed time:       %.3f s\n" time() - tstart

    # Assemble the stiffness matrix 
    t1 = time()
    K = assemble(nodes, dofs, elements, E ,nu)
    @printf "Matrix assembly elapsed time: %.3f s\n" time() - t1

    # Apply boundary conditions 
    # fixed_nodes = [3, 4, 6, 8, 91, 92, 93, 255, 256, 257, 340, 341, 342, 346, 347, 348, 
    #                 1297, 1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305]
    t1 = time()
    if length(fixed_nodes) > 0
        applyfixedbcs!(K, fixed_nodes)
    end 

    if length(constraints) > 0 
        for constraint in constraints
            K[constraint,:] .= 0.0 
            K[:,constraint] .= 0.0 
        end 
    end

    @printf "BC elapsed time:              %.3f s\n" time() - t1


    # Define loads 
    f = zeros(length(dofs))
    for i in eachindex(cload_dofs)
        f[cload_dofs[i]] = cload_forces[i]
    end

    # Send to solver 

    t1 = time()
    q = qr(Matrix(K), Val(true)) \ f
    @printf "Solver elapsed time:          %.3f s\n" time() - t1
    
    ps = MKLPardisoSolver()
    solve!(ps, q, K, f)

    @printf "Total SHEPHERD time:          %.3f s\n" time() - tstart
    return K, f, q
end



"""
    conjugate_gradient!(A, b, x)

Return the solution to `A * x = b` using the conjugate gradient method.
"""
function conjugate_gradient!(
    A::AbstractMatrix, b::AbstractVector, x::AbstractVector; tol=eps(eltype(b)), itmax=10000
)
    # Initialize residual vector
    residual = b - A * x
    # Initialize search direction vector
    search_direction = copy(residual)
    # Compute initial squared residual norm
    norm(x) = sqrt(sum(x.^2))
    old_resid_norm = norm(residual)

    it = 0

    # Iterate until convergence
    while (old_resid_norm > tol) && (it < itmax)
        A_search_direction = A * search_direction
        step_size = old_resid_norm^2 / (search_direction' * A_search_direction)
        # Update solution
        @. x = x + step_size * search_direction
        # Update residual
        @. residual = residual - step_size * A_search_direction
        new_resid_norm = norm(residual)
        
        # Update search direction vector
        @. search_direction = residual + 
            (new_resid_norm / old_resid_norm)^2 * search_direction
        # Update squared residual norm for next iteration
        old_resid_norm = new_resid_norm
        it += 1
    end

    @printf "Iterations: %i. residual: %.3e\n\n" it old_resid_norm
    return x
end
