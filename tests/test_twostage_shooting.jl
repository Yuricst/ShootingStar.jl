"""
Test two-stage shooting method
"""

using LinearAlgebra
using DifferentialEquations
using BenchmarkTools
using Plots

pyplot()

include("../../julia-R3BP/R3BP/src/R3BP.jl")
include("../src/ShootingStar.jl")


## Initialize ODE settings
reltol = 1.e-12
abstol = 1.e-12
method = Tsit5()

## define parameters
params = R3BP.get_cr3bp_param(399, 301)
mu = params.mu
println("mu: $mu")
lp = R3BP.lagrangePoints(mu)

# initial condition of halo
X0 = [1.176924090973164, 0.0, -0.060210863312217, 0.0, -0.173836346247689, 0.0]

T = 3.385326412831325

X2 = [1.176924090973164, 0.0, -0.060210863312217, 0.0, -0.173836346247689, 0.0] + 1.e-12*rand(6)

tspan = (0.0, 0.5T)
prob = ODEProblem(R3BP.rhs_cr3bp_sv!, X0, tspan, (mu))
sol = solve(prob, method, reltol=reltol, abstol=abstol)
X1 = sol.u[end] + 1.e-12*rand(6)  #+ [0.0, 0.0,-1.e-5, 0.0, -1.e-6, 0.0]


## prepare for two-stage shooting
n_sv = 6
tofs = [T/2, T/2]
svs = [X0, X1, X2]
n = length(svs)

prob_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, vcat(X0, reshape(I(6), (36,)))[:], T, (mu))
nodes_conv = ShootingStar.twostage_shooting(n_sv, svs, tofs, prob_stm)

# propagate result
prob = ODEProblem(R3BP.rhs_cr3bp_sv!, X0, T, (mu))
function prob_func(prob, i, repeat)
    remake(prob, u0=nodes_conv[i], tspan=(0.0, tofs[i]))
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob, method, EnsembleThreads(), trajectories=n-1, reltol=reltol, abstol=abstol);

display(plot(sim, linealpha=0.4, vars=(1,2), flip=false, aspect_ratio=:equal, size=(800,650), c=:orangered,
     frame_style=:box, gridalpha=0.4, xlabel="x", ylabel="y"))
println("Done!")

