using RMA
using CairoMakie
using Statistics
using Distributions
using ProgressMeter

x = rand(Uniform(0, 1), 1, 10000)
thres = range(0, 1, 100) * sqrt(3)

results = []

@showprogress for t in 1:100
    dist = distribution(x, thres[t], 4; sampling_mode = :full)
    a = findall(x -> x != 0, dist)
    push!(results, length(a))
end

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Threshold", ylabel = "Number of motifs not null")

lines!(ax, thres, results, color = :orange, linewidth=2)

println(string("Max N = ", findmax(results)))

fig