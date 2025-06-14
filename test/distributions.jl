using RMA, Test
using Distances
using Distributions
using DifferentialEquations

##
##      Lorenz system
function lorenz!(du, u, p, dt)
    x, y, z = u
    σ, ρ, β = p
        
    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
end

prob = ODEProblem(lorenz!, rand(Uniform(0, 1), 3), (0, 5000), [10.0, 28.0, 8.0/3.0])
sol = solve(prob, dt = 0.00001)

##
##      Data
uniform = rand(Uniform(0, 1), 5000)

##
##      Motif size
motif_n = [2, 3, 4]

## ----------------------------------------------------------------------------------------
##          Square motifs
@testset "square motif" begin
    for n in motif_n
        th, s = find_parameters(uniform, n)

        ##  Sampling mode: random
        dist_1 = distribution(uniform, th, n; num_samples = 0.5, sampling_mode = :random)
        dist_2 = distribution(uniform, th, n; num_samples = 0.5, sampling_mode = :random)
        @test length(dist_1) == length(dist_2) == 2^(n * n)
        @test sqrt(js_divergence(dist_1, dist_2)) < 0.1

        ##  Sampling mode: full
        dist_1 = distribution(uniform, th, n; sampling_mode = :full)
        dist_2 = distribution(uniform, th, n; sampling_mode = :full)
        @test length(dist_1) == length(dist_2) == 2^(n * n)
        @test dist_1 == dist_2

        ##  Sampling mode: columnwise
        dist_1 = distribution(uniform, th, n; num_samples = 0.95, sampling_mode = :columnwise)
        dist_2 = distribution(uniform, th, n; num_samples = 0.95, sampling_mode = :columnwise)
        sz = size(dist_1)
        @test sz == size(dist_2)
        for i in 1:sz[2]
            @test sqrt(js_divergence(dist_1[:, i], dist_2[:, i])) < 0.25
        end

        ##  Sampling mode: columnwise_full
        dist_1 = distribution(uniform, th, n; num_samples = 0.8, sampling_mode = :columnwise_full)
        dist_2 = distribution(uniform, th, n; num_samples = 0.8, sampling_mode = :columnwise_full)
        sz = size(dist_1)
        @test sz == size(dist_2)
        for i in 1:sz[2]
            @test dist_1[:, i] == dist_2[:, i]
        end

        ##  Sampling mode: triangle up
        dist_1 = distribution(uniform, th, n; num_samples = 0.5, sampling_mode = :triangleup)
        dist_2 = distribution(uniform, th, n; num_samples = 0.5, sampling_mode = :triangleup)
        @test length(dist_1) == length(dist_2) == 2^(n * n)
        @test sqrt(js_divergence(dist_1, dist_2)) < 0.1

        ##  Test using "prepare"
        dist_1 = distribution(sol, 2.7, n)
        dist_2 = distribution(sol, 2.7, n)
        @test length(dist_1) == length(dist_2) == 2^(n * n)
        @test sqrt(js_divergence(dist_1, dist_2)) < 0.1
    end
end
## ----------------------------------------------------------------------------------------
##          Triangle motifs
@testset "triangle motif" begin
    for n in motif_n
        th, s = find_parameters(uniform, n; shape = :triangle)

        ##  Sampling mode: random
        dist_1 = distribution(uniform, th, n; num_samples = 0.5, sampling_mode = :random, shape = :triangle)
        dist_2 = distribution(uniform, th, n; num_samples = 0.5, sampling_mode = :random, shape = :triangle)
        @test sqrt(js_divergence(dist_1, dist_2)) < 0.1

        ##  Sampling mode: full
        dist_1 = distribution(uniform, th, n; sampling_mode = :full, shape = :triangle)
        dist_2 = distribution(uniform, th, n; sampling_mode = :full, shape = :triangle)
        @test sqrt(js_divergence(dist_1, dist_2)) < 0.1
    end
end
## ----------------------------------------------------------------------------------------
##          