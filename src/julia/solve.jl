""" SHEPHERD

    Solve a finite element problem
"""


using Printf, Pardiso


"""
    solve(fn::AbstractString)

Solve a finite element problem by loading a Calculix input file.
"""
function solve(fn::AbstractString; solver="CG", fixed_nodes=[], constraints=[], verbose=false)
    
    # Read the input file, load data structures
    # Needs to support *NODE, *ELEMENT, *MATERIAL, *Cload
    if verbose @printf "\nReading input file from '%s' \n" fn end
    tstart = time()
    model = readinputfile(fn)
    @printf "File read elapsed time:       %.3f s\n" time() - tstart

    # Find constrained, free, and loaded dofs
    model.fixed_nodes = fixed_nodes
    constrained_dofs, free_dofs = assigndofs(model)
    force_dofs = assignforcedofs(model)

    # Assemble the global stiffness matrix 
    t1 = time()
    dofs = [x for x in 1:3*length(model.nodes)]
    K = assemble(model.nodes, dofs, model.elements, model.material[1], model.material[2])
    @printf "Matrix assembly elapsed time: %.3f s\n" time() - t1

    # Apply boundary conditions 
    # fixed_nodes = [3, 4, 6, 8, 91, 92, 93, 255, 256, 257, 340, 341, 342, 346, 347, 348, 
    #                 1297, 1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305]
    t1 = time()
    # Kr = copy(K)
    # if length(fixed_nodes) > 0
    #     applyfixedbcs!(Kr, fixed_nodes)
    # end 

    # if length(constraints) > 0 
    #     for constraint in constraints
    #         Kr[constraint,:] .= 0.0 
    #         Kr[:,constraint] .= 0.0 
    #     end 
    # end

    Kr, F = reducesystem(model, K, free_dofs, force_dofs)

    @printf "BC elapsed time:              %.3f s\n" time() - t1

    # Send to solver 

    t1 = time()
    q = zeros(size(F))
    if solver == "CG"
        conjugate_gradient!(Matrix(Kr), F, q)
    elseif solver == "pardiso"
        solve!(MKLPardisoSolver(), q, sparse(Kr), F)
    else 
        q = qr(Matrix(Kr), Val(true)) \ F
    end
    @printf "Solver elapsed time:          %.3f s\n" time() - t1

    @printf "Total SHEPHERD time:          %.3f s\n" time() - tstart
    return K, Kr, F, q
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

        @printf "Iteration: %i, max residual: %.3e\n" it old_resid_norm
    end

    @printf "Iterations: %i. residual: %.3e\n\n" it old_resid_norm
    return x
end
