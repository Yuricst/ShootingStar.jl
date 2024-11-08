module ShootingStar

    using FiniteDiff
    using LinearAlgebra
    using OrdinaryDiffEq
    using Printf

    include("eom.jl")
    include("problem.jl")
    include("innerloop.jl")
    include("outerloop.jl")

end # module
