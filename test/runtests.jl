"""Run tests for ShootingStar.jl"""

using Test

@testset "Solve" begin
    include("test_solve.jl")
end