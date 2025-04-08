using RMA
using CairoMakie
using Statistics
using RecurrenceAnalysis
using DifferentialEquations

include("../systems.jl")

##
##      Function sin(x * pi)
x = range(0, 12, 500)
y = sin.(x * pi)

th_1, _ = RMA.find_parameters(y, 3)
rp_1 = RecurrenceMatrix(y, th_1)

##
##      Lorenz
problem = ODEProblem(lorenz!, rand(Float64, 3), (0, 5000), [10.0, 28.0, 8.0/3.0])
solution = solve(problem)

lorenz = prepare(solution, 0.2; transient = 10000, max_length = 500)

norm_data = lorenz
norm_data[1, :] .= (norm_data[1, :] .- mean(norm_data[1, :])) ./ std(norm_data[1, :])
norm_data[2, :] .= (norm_data[2, :] .- mean(norm_data[2, :])) ./ std(norm_data[2, :])
norm_data[3, :] .= (norm_data[3, :] .- mean(norm_data[3, :])) ./ std(norm_data[3, :])
lorenz = norm_data

data =  lorenz[1, :]
th_2, _ = RMA.find_parameters(data, 3; threshold_max = 10.0)
rp_2 = RecurrenceMatrix(data, th_2)

##
##      Make graphs
fig = Figure(size = (720, 720))

ax_1 = Axis(fig[1, 1], yautolimitmargin = (0, 0), ylabel = "Time", xgridvisible = false, ygridvisible = false)
lines!(ax_1, y, x, color = :Black)
hidexdecorations!(ax_1)

ax_2 = Axis(fig[1, 2])
heatmap!(ax_2, rp_1, colormap = :binary)
hidedecorations!(ax_2)

ax_3 = Axis(fig[2, 2], xautolimitmargin = (0, 0), xlabel = "Time", xgridvisible = false, ygridvisible = false)
lines!(ax_3, x, y, color = :Black)
hideydecorations!(ax_3)

colsize!(fig.layout, 1, Relative(0.2))
rowsize!(fig.layout, 2, Relative(0.2))

save("test/rp/sin.png", fig)

fig = Figure(size = (720, 720))

ax_1 = Axis(fig[1, 1], yautolimitmargin = (0, 0), ylabel = "Time", xgridvisible = false, ygridvisible = false)
lines!(ax_1, data, 1:500, color = :Black)
hidexdecorations!(ax_1)

ax_2 = Axis(fig[1, 2])
heatmap!(ax_2, rp_2, colormap = :binary)
hidedecorations!(ax_2)

ax_3 = Axis(fig[2, 2], xautolimitmargin = (0, 0), xlabel = "Time", xgridvisible = false, ygridvisible = false)
lines!(ax_3, 1:500, data, color = :Black)
hideydecorations!(ax_3)

colsize!(fig.layout, 1, Relative(0.2))
rowsize!(fig.layout, 2, Relative(0.2))

save("test/rp/lorenz.png", fig)