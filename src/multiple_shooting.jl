"""
Propagation module
"""

# struct for handling information of the problem
struct ShootingProblem
	n::Int
	n_sv::Int
	prob_stm::ODEProblem
	tofs::Array
	r0::Array    # fixed values
	rf::Array    # fixed values
	method
	reltol::Float64
	abstol::Float64
end


"""
Multiple-shooting 
"""
function two_stage_shooting(n_sv::Int, x0s::Array, tofs::Array, prob_stm::ODEProblem; kwargs...)

	reltol = 1.e-12
	abstol = 1.e-12
	method = Tsit5()

	n  = length(x0s)  # number of nodes
	nr = n_sv ÷ 2   # length of positions (== length of velocities)

	# initialize
	rs, vs = [], []
	for (i,x0) in enumerate(x0s)
		push!(rs, x0[1:nr])
		push!(vs, x0[nr+1:end])
	end

	# fixed initial and final position
	r0 = rs[1]
	rf = rs[end]

	# construct ShootingProblem
	Shoot = ShootingProblem(n, n_sv, prob_stm, tofs, r0, rf, method, reltol, abstol)

end


function outer_loop_shooting!(Shoot::ShootingProblem)
end


function inner_loop_shooting!(Shoot::ShootingProblem)

	# propagate n nodes
	for i = 1:Shoot.n-1
		x0i =  vcat(r_i, v_i, reshape(I(Shoot.n_sv), (Shoot.n_sv*Shoot.n_sv)))[:]
		_prob = remake(Shoot.prob_stm; tspan=(0.0, Shoot.tofs[i]), u0=x0_stm)
        sol = DifferentialEquations.solve(_prob, Shoot.method, reltol=Shoot.reltol, abstol=Shoot.abstol)
	end
end