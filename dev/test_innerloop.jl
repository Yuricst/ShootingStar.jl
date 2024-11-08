"""Test for inner loop"""


using GLMakie
using LinearAlgebra
using Printf
using OrdinaryDiffEq

include(joinpath(@__DIR__, "..", "src", "ShootingStar.jl"))

# ODE parameters
μ = 1.215058560962404e-02
params_ode = [μ,]

# create initial guess
rv0 = [1.0809931218390707E+00,
    0.0000000000000000E+00,
    -2.0235953267405354E-01,
    1.0157158264396639E-14,
    -1.9895001215078018E-01,
    7.2218178975912707E-15]
period_0 = 2.3538670417546639E+00

rvf = [1.1648780946517576,
    0.0,
    -1.1145303634437023E-1,
    0.0,
    -2.0191923237095796E-1,
    0.0]
period_f = 3.3031221822879884

# initial & final LPO
sol_lpo0 = solve(
    ODEProblem(ShootingStar.rhs_cr3bp_sv!, rv0, [0.0, period_0], params_ode),
    Tsit5(); reltol = 1e-12, abstol = 1e-12
)
sol_lpof = solve(
    ODEProblem(ShootingStar.rhs_cr3bp_sv!, rvf, [0.0, period_f], params_ode),
    Tsit5(); reltol = 1e-12, abstol = 1e-12
)

# create problem
tf = 2.5
Nseg = 20
times = [el for el in LinRange(0.0, tf, Nseg+1)]
nodes = ShootingStar.initialguess_gradual_transit(rv0, rvf, times, 
    ShootingStar.rhs_cr3bp_sv!, params_ode)

# create problem
prob = ShootingStar.TwoStageShootingProblem(
    rv0,
    rvf,
    times,
    nodes,
    ShootingStar.rhs_cr3bp_svstm!,
    params_ode,
)

# solve innerloop
# println("Solving inner-loop!")
# sols, success_flag, residuals = ShootingStar.solve_innerloop!(prob, 10, 1e-12; set_to_problem = true)
# _, residuals, J_inner = ShootingStar.propagate_nodes(prob, prob.nodes)

# solve outerloop
maxiter = 10
status, sols, residuals = ShootingStar.solve_outerloop!(
    prob,
    maxiter,
    1e-6;
    verbosity = 1,
)

# create plot
fig = Figure()
ax3d = Axis3(fig[1,1],
    aspect = :data,
    xlabel = "x",
    ylabel = "y",
    zlabel = "z",
    title = "Solve test",
)
lines!(ax3d, Array(sol_lpo0)[1,:], Array(sol_lpo0)[2,:], Array(sol_lpo0)[3,:], color = :blue)
lines!(ax3d, Array(sol_lpof)[1,:], Array(sol_lpof)[2,:], Array(sol_lpof)[3,:], color = :green)
for sol in sols
    lines!(ax3d, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:], color = :black)
end
# scalar_DV = 1.0
# arrows!(ax3d, prob.nodes[1,:], prob.nodes[2,:], prob.nodes[3,:],
#     scalar_DV*residuals[4,:], scalar_DV*residuals[5,:], scalar_DV*residuals[6,:], 
#     linewidth = 0.1, arrowsize = Vec3f(0.01, 0.01, 0.01), color = :red)
display(fig)