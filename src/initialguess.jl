"""Methods to construct initial guess"""


"""
Construct linear guess by connecting initial and final position with straight line
"""
function initialguess_linear(rv0::Vector, rvf::Vector, times)
    @assert length(rv0) == length(rvf) == 6 "rv0 and rvf should have length 6"

    Nseg = length(times) - 1
    nodes = hcat([[el for el in LinRange(rv0[i], rvf[i], Nseg+1)] for i in 1:6]...)'

    # overwrite each velocity by connecting with straight line
    for i in 1:Nseg
        vi_guess = (nodes[1:3,i+1] - nodes[1:3,i])/(times[i+1] - times[i])
        nodes[4:6,i] = vi_guess
    end
    return nodes
end


function initialguess_gradual_transit(rv0, rvf, times, rhs!::Function, params_ode;
    method = Tsit5(), reltol = 1e-12, abstol = 1e-12)
    # propagate initial state
    sol_fwd = solve(
        ODEProblem(rhs!, rv0, [times[1] times[end]], params_ode),
        method; reltol = reltol, abstol = abstol
    )
    sol_bck = solve(
        ODEProblem(rhs!, rvf, [times[end] times[1]], params_ode),
        method; reltol = reltol, abstol = abstol
    )

    # evaluate solution at times
    rv_fwd = [sol_fwd(t) for t in times]
    rv_bck = [sol_bck(t) for t in times]

    # take gradual average
    alphas = LinRange(1, 0, length(times))
    nodes = hcat([alpha*rv_fwd[i] + (1-alpha)*rv_bck[i] for (i, alpha) in enumerate(alphas)]...)
    return nodes
end