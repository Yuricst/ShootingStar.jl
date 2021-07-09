"""
Propagation module
"""

# struct for handling information of the problem
struct ShootingSettings
	n::Int
	nr::Int
	n_sv::Int
	prob_stm::ODEProblem
	tofs::Array
	r0::Array    # fixed vector
	rf::Array    # fixed vector
	v0::Array
	vf::Array
	method
	reltol::Float64
	abstol::Float64
	tolDC::Float64
	tolConv::Float64
	maxiter::Int
	verbosity::Int
end


"""
Two-stage shooting algorithm for n-impulsive trajectory design
"""
function twostage_shooting(n_sv::Int, svs::Array, tofs::Array, prob_stm::ODEProblem; kwargs...)

	# FIXME -- take from kwargs
	reltol = 1.e-12
	abstol = 1.e-12
	method = Tsit5()
	tolDC = 1.e-11
	tolConv = 1.e-6
	maxiter = 15
	verbosity = 1

	n  = length(svs)  # number of nodes
	nr = n_sv ÷ 2     # number of elements of positions (== length of velocities)

	println("Using $n nodes")
	for (i,sv) in enumerate(svs)
		print("Node $i: \n")
		println(sv)
	end

	# initial and final position and velocity vectors
	r0 = svs[1][1:nr]
	v0 = svs[1][nr+1:end]
	rf = svs[end][1:nr]
	vf = svs[end][nr+1:end]

	# inner-loop decision vector: velocities at first n-1 nodes
	x_inner = zeros(nr*(n-1))
	for i = 1:n-1
		x_inner[1+(i-1)*nr:i*nr] = svs[i][nr+1:end]
	end

	# outer-loop decision vector
	x_outer = zeros((n-2)*nr)
	for j = 1:n-2
		x_outer[1+(j-1)*nr:j*nr] = svs[j+1][1:nr]
	end

	println("\nx_outer: $x_outer")
	println("x_inner: $x_inner\n\n")

	# construct ShootingSettings
	Settings = ShootingSettings(n, nr, n_sv, prob_stm, tofs, r0, rf, v0, vf, 
		method, reltol, abstol, tolDC, tolConv, maxiter, verbosity)

	# initialize storage
	J_inner = zeros(nr*(n-1), nr*(n-1))

	# outer-loop shooting to minimize velocity discontinuity
	outer_loop_shooting!(x_outer, x_inner, J_inner, Settings)

	println("\nx_outer: $x_outer")
	println("x_inner: $x_inner\n\n")

	# construct final solution
	return
end



"""
Outer-loop finds least-square solution to minimize velocity discontinuity.
Function mutates `x0_outer`
"""
function outer_loop_shooting!(x0_outer, x0_inner, J_inner, Settings::ShootingSettings)

	# initialize old velocities
	old_vs = vcat(Settings.v0, x0_inner, Settings.vf)

	for i_outer = 1:Settings.maxiter

		# inner-loop shooting to ensure position continuity
		inner_loop_shooting!(x0_outer, x0_inner, J_inner, Settings, true)

		# compute new velocities
		new_vs = vcat(Settings.v0, x0_inner, Settings.vf)

		# break if norm of cost stops improving
		if (norm(old_vs) - norm(new_vs) < Settings.tolConv)
			println("Outer-loop achieved tolerance!")
			break
		else  # update old velocities
			old_vs[:] = vcat(Settings.v0, x0_inner, Settings.vf)
		end

		# compute J_outer via finite difference
		function g(x)
			inner_loop_shooting!(x, x0_inner, J_inner, Settings, false)
			return x
		end
		J_outer = FiniteDiff.finite_difference_jacobian(g, x0_outer)
		
		# least-square update
		print("x0_outer: "); println(length(x0_outer))
		print("J_outer: "); println(size(J_outer))
		print("new_vs: "); println(length(new_vs))
		x0_outer[:] = x0_outer -inv(transpose(J_outer)*J_outer) * transpose(J_outer) * new_vs
	end

	return
end



"""
Inner-loop corrects velocity vectors to ensure position continuity.
Function mutates `x0_inner`
"""
function inner_loop_shooting!(x0_outer, x0_inner, J_inner, Settings::ShootingSettings, verbose::Bool=false)

	r_i2_prop = zeros(Settings.nr * (Settings.n-1))

	for i_inner = 1:Settings.maxiter

		# propagate n nodes
		for i = 1:Settings.n-1
			# re-construct array of positions from node 1 ~ n-1
			rs_vec = vcat(Settings.r0, x0_outer)
			#
			r_i = rs_vec[1+(i-1)*Settings.nr : i*Settings.nr]
			v_i = x0_inner[1+(i-1)*Settings.nr : i*Settings.nr]
			x0i =  vcat(r_i, v_i, reshape(I(Settings.n_sv), (Settings.n_sv*Settings.n_sv)))[:]
			_prob = remake(Settings.prob_stm; tspan=(0.0, Settings.tofs[i]), u0=x0i)
	        sol = DifferentialEquations.solve(_prob, Settings.method, reltol=Settings.reltol, abstol=Settings.abstol)

	        # store results
	        r_i2_prop[1+(i-1)*Settings.nr : i*Settings.nr] = sol.u[end][1:Settings.nr]

	    	# fill-in upper-left submatrix of STM into Jacobian
	    	STM = transpose( reshape(sol.u[end][Settings.n_sv+1:end], (Settings.n_sv, Settings.n_sv)) )
	    	J_inner[1+(i-1)*Settings.nr:i*Settings.nr, 1+(i-1)*Settings.nr:i*Settings.nr] = -STM[1:Settings.nr, Settings.nr+1:end]  # FIXME
		end

		# compute final position offset: δr_i2 = r_i2_guess - r_i2_prop
		δr_i2 = vcat(x0_outer, Settings.rf) - r_i2_prop

		# check breaking condition
		err = norm(δr_i2)
		if verbose==true
			println("Iteration $i_inner ... err: $err")
			println("x0_inner: $x0_inner")
		end
		if err < Settings.tolDC
			if verbose==true
				println("Inner-loop achieved tolerance!")
			end
			break
		end

		# correct velocity: v_i = v_i_guess - J_inner^-1 * δr_i2
		x0_inner[:] = x0_inner - inv(J_inner) * δr_i2

	end

	return
end
