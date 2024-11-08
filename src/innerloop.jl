"""Inner-loop problem"""


function solve_innerloop!(
    problem::TwoStageShootingProblem,
    maxiter::Int,
    eps_inner::Float64,
    verbosity::Int = 1;
    set_to_problem::Bool = false,
)
    return solve_innerloop!(
        problem,
        problem.nodes,
        maxiter,
        eps_inner,
        verbosity;
        set_to_problem = set_to_problem,
    )
end



"""Inner-loop solve algorithm with Newton's method"""
function solve_innerloop!(
    problem::TwoStageShootingProblem,
    nodes_copy,
    maxiter::Int,
    eps_inner::Float64,
    verbosity::Int = 1;
    set_to_problem::Bool = false,
)
    if verbosity > 0
        @printf("  iter  |  res-norm  |\n")
    end

    success_flag = false
    residuals = zeros(problem.nx, problem.N)
    sols = []

    for it in 1:maxiter
        # propagate nodes
        sols, residuals, J_inner = propagate_nodes(problem, nodes_copy)
        res_norms = zeros(problem.Nseg)

        # update nodes sequentially by updating the velocity
        for iseg in 1:problem.Nseg
            _J = J_inner[iseg,:,:]
            nodes_copy[4:6,iseg] -= inv(_J) * residuals[1:3,iseg+1]
            res_norms[iseg] = norm(residuals[1:3,iseg+1])
        end
        
        # print information
        if verbosity > 0
            @printf("  %4.0f  | %1.4e |\n", it, maximum(res_norms))
        end

        if maximum(res_norms) < eps_inner
            success_flag = true
            break
        end
    end
    if set_to_problem == true
        problem.nodes[:,:] = nodes_copy
    end
    return sols, success_flag, residuals
end