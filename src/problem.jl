"""Two-stage shooting problem class"""


mutable struct TwoStageShootingProblem

    N::Int
    Nseg::Int
    nx::Int
    rv0
    rvf
    times::Vector
    nodes::Union{Matrix{Float64}, Adjoint{Float64, Matrix{Float64}}}
    rhs_stm!::Function
    params

    method
    reltol::Union{Nothing,Float64}
    abstol::Union{Nothing,Float64}
    ode_dt::Union{Nothing,Real}

    function TwoStageShootingProblem(
        rv0,
        rvf,
        times::Vector,
        nodes::Union{Matrix{Float64}, Adjoint{Float64, Matrix{Float64}}},
        rhs_stm!::Function,
        params;
        method = Tsit5(),
        reltol = 1e-12,
        abstol = 1e-12,
        ode_dt = nothing,
    )
        Nseg = length(times) - 1
        nx = size(nodes, 1)
        @assert size(nodes, 2) == Nseg + 1 "nodes should have size (nx, Nseg+1)!"

        new(
            Nseg + 1,
            Nseg,
            nx,
            rv0,
            rvf,
            times,
            nodes,
            rhs_stm!,
            params,
            method,
            reltol,
            abstol,
            ode_dt,
        )
    end
end


function update_velocity_initial_guess!(problem::TwoStageShootingProblem)
    for i in 1:problem.Nseg
        vi_guess = (problem.nodes[1:3,i+1] - problem.nodes[1:3,i])/(problem.times[i+1] - problem.times[i])
        problem.nodes[4:6,i] = vi_guess
    end
end


"""
Overload method for showing OptimalControlSCPProblem
"""
function Base.show(io::IO, problem::TwoStageShootingProblem)
    println("Two-stage shooting problem")
    @printf("  number of nodes    : %d\n", problem.N)
    @printf("  number of segments : %d\n", problem.Nseg)
end

# function propagate_nodes(problem::TwoStageShootingProblem)
#     return propagate_nodes(problem, problem.nodes)
# end


function propagate_nodes(problem::TwoStageShootingProblem, nodes)
    sols = [] #Vector{ODESolution}[]
    residuals = zeros(problem.nx, problem.N)
    J_inner = zeros((problem.Nseg, 3, 3))

    # initial node residuals
    residuals[1:6,1] = nodes[1:6,1] - problem.rv0
    nodes[:,end] = problem.rvf

    for iseg in 1:problem.Nseg
        # set ODE problem and solve
        x0_aug = [nodes[:,iseg]; reshape(I(problem.nx), problem.nx^2)]
        tspan  = [problem.times[iseg], problem.times[iseg+1]]
        ode_problem = ODEProblem(
            problem.rhs_stm!,
            x0_aug,
            tspan,
            problem.params,
        )
        if isnothing(problem.ode_dt)           # use adaptive integration
            sol = solve(
                ode_problem,
                problem.method;
                reltol = problem.reltol,
                abstol = problem.abstol
            )
        else                                    # use fixed-step integration
            sol = solve(
                ode_problem,
                problem.method;
                dt = problem.ode_dt,
            )
        end
        stmf = reshape(sol.u[end][problem.nx+1:end], (problem.nx, problem.nx))'
        J_inner[iseg,:,:] = -stmf[1:3,4:6]    # sensitivity of final position w.r.t. initial velocity

        # store
        push!(sols, sol)
        residuals[:,iseg+1] = nodes[:,iseg+1] - sol.u[end][1:problem.nx]
    end
    return sols, residuals, J_inner
end