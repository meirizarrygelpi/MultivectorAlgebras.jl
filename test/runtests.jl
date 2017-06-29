using MultivectorAlgebras
using Base.Test: @testset

@testset "1-dimensional" begin include("one-dimensional_tests.jl") end
@testset "2-dimensional" begin include("two-dimensional_tests.jl") end
