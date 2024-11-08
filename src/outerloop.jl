"""Outerloop problem"""


"""Outer-loop solve algorithm with Newton's method"""
function solve_outerloop!(
    problem::TwoStageShootingProblem,
    maxiter::Int,
    eps_outer::Float64;
    verbosity::Int = 1,
    maxiter_inner::Int = 20,
    eps_inner::Float64 = 1e-12,
    verbosity_inner::Int = 0,
)
    sols = []
    residuals = zeros(problem.nx, problem.N)
    status = :NotStarted

    if verbosity >= 1
        @printf(" iter  |    res-norm   |\n")
    end
    
    sols = []
    Fnorm_last = 1e16

    for it in 1:maxiter
        # solve inner-loop
        sols, success_inner, residuals = solve_innerloop!(
            problem, maxiter_inner, eps_inner, verbosity_inner;
            set_to_problem = true
        )
        if success_inner == false
            status = :InnerLoopFailed
        end

        # function to compute Jacobian
        if verbosity >= 2
            @printf("    Computing Jacobian...\n")
        end
        _solve_innerloop_for_jacobian = function (r_flattened_intermediate)
            _nodes = vcat(
                hcat(problem.nodes[1:3,1], reshape(r_flattened_intermediate, (3, problem.N-2)), problem.nodes[1:3,end]),
                problem.nodes[4:6,:]
            )
            _, _, _residuals = solve_innerloop!(
                problem,
                _nodes,
                maxiter_inner,
                eps_inner,
                verbosity_inner;
                set_to_problem = false
            )
            return reshape(_residuals[4:6,:], 3*problem.N)
        end

        # compute outer-loop jacobian
        r_flattened_intermediate = reshape(problem.nodes[1:3,2:end-1], 3*(problem.N-2))
        DF = FiniteDiff.finite_difference_jacobian(
            _solve_innerloop_for_jacobian,
            r_flattened_intermediate,
        )

        # iterate least-squares update
        F = reshape(residuals[4:6,:],3*problem.N)
        new_r_flat = r_flattened_intermediate - inv(transpose(DF)*DF)*transpose(DF) * F

        if abs(norm(F) - Fnorm_last) < eps_outer
            status = :Converged
            break
        else
            Fnorm_last = norm(F)
        end

        # store into nodes
        problem.nodes[1:3,2:end-1] = reshape(new_r_flat, 3, problem.N-2)

        if verbosity >= 1
            @printf("  %3.0f  | %1.7e |\n",
                it,
                Fnorm_last)
        end

        if it == maxiter
            status = :MaxIterReached
        end
    end

    if verbosity >= 1
        println("")
        @printf("    Status : %s\n", string(status))
        println("")
    end
    return status, sols, residuals
end