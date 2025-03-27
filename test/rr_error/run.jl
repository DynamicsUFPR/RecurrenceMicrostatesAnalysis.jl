include("../../src/RMA.jl")
using .RMA
using CairoMakie
using Statistics
using Distributions
using RecurrenceAnalysis
using DifferentialEquations

using ProgressMeter

##
##
include("../systems.jl")

##
##      Settings
data_len = 3000
samples_percent = 0.05

num_tests = 100

##
##      Compute the errors
function compute()
    result = zeros(Float64, 5, 5, num_tests)

    ##
    ##      Uniform Distribution.
    println("Running: Uniform Distribution")
    @showprogress for i in 1:num_tests
        ##
        ##      Generate the data.
        uniform = rand(Uniform(0, 1), data_len)
        thres, _ = RMA.find_parameters(uniform, 3)
        
        rr = rqa(RecurrenceMatrix(uniform, thres))[:RR]
        
        rr_square = rrate(uniform, thres, 3; num_samples = samples_percent, shape = :square)
        rr_triangle = rrate(uniform, thres, 3; num_samples = samples_percent, shape = :triangle)
        rr_diagonal = rrate(uniform, thres, 3; num_samples = samples_percent, shape = :diagonal)
        rr_column = rrate(uniform, thres, 3; num_samples = samples_percent, shape = :column)
        rr_pair = rrate(uniform, thres, 3; num_samples = samples_percent, shape = :pair)

        result[1, 1, i] = (abs(rr_square - rr) / rr)
        result[2, 1, i] = (abs(rr_triangle - rr) / rr)
        result[3, 1, i] = (abs(rr_diagonal - rr) / rr)
        result[4, 1, i] = (abs(rr_column - rr) / rr)
        result[5, 1, i] = (abs(rr_pair - rr) / rr)

        # println(string(rr, ", s: ", result[1, 1, i], ", t: ", result[2, 1, i], ", d: ", result[3, 1, i], ", c: ", result[4, 1, i], ", p: ", result[5, 1, i]))
    end

    ##
    ##      Lorenz System
    println("Running: Lorenz System")
    @showprogress for i in 1:num_tests
        ##
        ##      Generate the data.
        problem = ODEProblem(lorenz!, rand(Float64, 3), (0, 5000), [10.0, 28.0, 8.0/3.0])
        solution = solve(problem)

        ##      Normalize the data.
        means = mean.(eachrow(hcat(solution.u...)))
        stds = std.(eachrow(hcat(solution.u...)))

        normalized_sol_u = [(u .- means) ./ stds for u in solution.u]
        sol_normalized = deepcopy(solution)
        sol_normalized.u .= normalized_sol_u

        ##      Prepare the solution.
        lorenz = prepare(sol_normalized, 0.2; transient = 10000, max_length = data_len)

        ##      Get threshold.
        thres, _ = RMA.find_parameters(lorenz, 3; threshold_max = 10.0)
        
        rr = rqa(RecurrenceMatrix(StateSpaceSet((hcat(lorenz...))'), thres))[:RR]
        
        rr_square = rrate(lorenz, thres, 3; num_samples = samples_percent, shape = :square)
        rr_triangle = rrate(lorenz, thres, 3; num_samples = samples_percent, shape = :triangle)
        rr_diagonal = rrate(lorenz, thres, 3; num_samples = samples_percent, shape = :diagonal)
        rr_column = rrate(lorenz, thres, 3; num_samples = samples_percent, shape = :column)
        rr_pair = rrate(lorenz, thres, 3; num_samples = samples_percent, shape = :pair)

        result[1, 2, i] = (abs(rr_square - rr) / rr)
        result[2, 2, i] = (abs(rr_triangle - rr) / rr)
        result[3, 2, i] = (abs(rr_diagonal - rr) / rr)
        result[4, 2, i] = (abs(rr_column - rr) / rr)
        result[5, 2, i] = (abs(rr_pair - rr) / rr)

        # println(string("thres = ", thres, " | ", rr, ", s: ", result[1, 2, i], ", t: ", result[2, 2, i], ", d: ", result[3, 2, i], ", c: ", result[4, 2, i], ", p: ", result[5, 2, i]))
    end

    ##
    ##      Logistic Map
    println("Running: Logistic Map")
    @showprogress for i in 1:num_tests
        ##
        ##      Generate the data.
        logistic = logistic_map(4.0, data_len)
        thres, _ = RMA.find_parameters(logistic, 3)
        
        rr = rqa(RecurrenceMatrix(logistic, thres))[:RR]
        
        rr_square = rrate(logistic, thres, 3; num_samples = samples_percent, shape = :square)
        rr_triangle = rrate(logistic, thres, 3; num_samples = samples_percent, shape = :triangle)
        rr_diagonal = rrate(logistic, thres, 3; num_samples = samples_percent, shape = :diagonal)
        rr_column = rrate(logistic, thres, 3; num_samples = samples_percent, shape = :column)
        rr_pair = rrate(logistic, thres, 3; num_samples = samples_percent, shape = :pair)

        result[1, 3, i] = (abs(rr_square - rr) / rr)
        result[2, 3, i] = (abs(rr_triangle - rr) / rr)
        result[3, 3, i] = (abs(rr_diagonal - rr) / rr)
        result[4, 3, i] = (abs(rr_column - rr) / rr)
        result[5, 3, i] = (abs(rr_pair - rr) / rr)

        # println(string("thres = ", thres, " | ", rr, ", s: ", result[1, 3, i], ", t: ", result[2, 3, i], ", d: ", result[3, 3, i], ", c: ", result[4, 3, i], ", p: ", result[5, 3, i]))
    end

    ##
    ##      Lorenz System
    println("Running: Rössler System")
    @showprogress for i in 1:num_tests
        ##
        ##      Generate the data.
        problem = ODEProblem(rossler!, rand(Float64, 3), (0, 10000), [0.2, 0.2, 18.0])
        solution = solve(problem, Tsit5())
        
        ##      Normalize the data.
        means = mean.(eachrow(hcat(solution.u...)))
        stds = std.(eachrow(hcat(solution.u...)))

        normalized_sol_u = [(u .- means) ./ stds for u in solution.u]
        sol_normalized = deepcopy(solution)
        sol_normalized.u .= normalized_sol_u

        rossler = prepare(sol_normalized, 0.18; transient = 20000, max_length = data_len)
        thres, _ = RMA.find_parameters(rossler, 3; threshold_max = 10.0)
        
        rr = rqa(RecurrenceMatrix(StateSpaceSet((hcat(rossler...))'), thres))[:RR]
        
        rr_square = rrate(rossler, thres, 3; num_samples = samples_percent, shape = :square)
        rr_triangle = rrate(rossler, thres, 3; num_samples = samples_percent, shape = :triangle)
        rr_diagonal = rrate(rossler, thres, 3; num_samples = samples_percent, shape = :diagonal)
        rr_column = rrate(rossler, thres, 3; num_samples = samples_percent, shape = :column)
        rr_pair = rrate(rossler, thres, 3; num_samples = samples_percent, shape = :pair)

        result[1, 4, i] = (abs(rr_square - rr) / rr)
        result[2, 4, i] = (abs(rr_triangle - rr) / rr)
        result[3, 4, i] = (abs(rr_diagonal - rr) / rr)
        result[4, 4, i] = (abs(rr_column - rr) / rr)
        result[5, 4, i] = (abs(rr_pair - rr) / rr)

        # println(string("thres = ", thres, " | ", rr, ", s: ", result[1, 4, i], ", t: ", result[2, 4, i], ", d: ", result[3, 4, i], ", c: ", result[4, 4, i], ", p: ", result[5, 4, i]))
    end

    ##
    ##      Bernoulli Shifted Generalized (BSG / Beta-X)
    println("Running: BSG/Beta-X")
    @showprogress for i in 1:num_tests
        ##
        ##      Generate the data.
        bsg = beta_x(4.99, data_len)
        thres, _ = RMA.find_parameters(bsg, 3)
        
        rr = rqa(RecurrenceMatrix(bsg, thres))[:RR]
        
        rr_square = rrate(bsg, thres, 3; num_samples = samples_percent, shape = :square)
        rr_triangle = rrate(bsg, thres, 3; num_samples = samples_percent, shape = :triangle)
        rr_diagonal = rrate(bsg, thres, 3; num_samples = samples_percent, shape = :diagonal)
        rr_column = rrate(bsg, thres, 3; num_samples = samples_percent, shape = :column)
        rr_pair = rrate(bsg, thres, 3; num_samples = samples_percent, shape = :pair)

        result[1, 5, i] = (abs(rr_square - rr) / rr)
        result[2, 5, i] = (abs(rr_triangle - rr) / rr)
        result[3, 5, i] = (abs(rr_diagonal - rr) / rr)
        result[4, 5, i] = (abs(rr_column - rr) / rr)
        result[5, 5, i] = (abs(rr_pair - rr) / rr)

        # println(string("thres = ", thres, " | ", rr, ", s: ", result[1, 3, i], ", t: ", result[2, 3, i], ", d: ", result[3, 3, i], ", c: ", result[4, 3, i], ", p: ", result[5, 3, i]))
    end

    return result
end

function graph()
    dataset = compute()

    ##
    ##      Create the figure.
    fig = Figure(size = (1024, 1080))

    ##
    ##      Categories for boxplot.
    cat = ones(Int, size(dataset, 3))
    for i = 2:size(dataset, 1)
        cat = vcat(cat, ones(Int, size(dataset, 3)) * i)
    end

    ##
    ##      Uniform distribution
    values = dataset[1, 1, :]
    for i = 2:size(dataset, 1)
        values =vcat(values, dataset[i, 1, :])
    end

    ax_u = Axis(fig[1, 1],
        xticks = (1:5, [":square", ":triangle", ":diagonal", ":column", ":pair"]),
        xgridvisible = false,
        ygridvisible = false,
        title = "a) Uniform Distribution",
        ylabel = "ΔError (%)")

    ax_u.ytickformat = y -> string.(round.(y, digits=2), "%")
    boxplot!(ax_u, cat, values .* 100, color = :orange)

    ##
    ##      Lorenz System
    values = dataset[1, 2, :]
    for i = 2:size(dataset, 1)
        values =vcat(values, dataset[i, 2, :])
    end

    ax_lz = Axis(fig[1, 2],
        xticks = (1:5, [":square", ":triangle", ":diagonal", ":column", ":pair"]),
        xgridvisible = false,
        ygridvisible = false,
        title = "b) Lorenz System",
        ylabel = "ΔError (%)")

    ax_lz.ytickformat = y -> string.(round.(y, digits=2), "%")
    boxplot!(ax_lz, cat, values .* 100, color = :orange)

    ##
    ##      Logistic Map
    values = dataset[1, 3, :]
    for i = 2:size(dataset, 1)
        values =vcat(values, dataset[i, 3, :])
    end

    ax_lm = Axis(fig[2, 1],
        xticks = (1:5, [":square", ":triangle", ":diagonal", ":column", ":pair"]),
        xgridvisible = false,
        ygridvisible = false,
        title = "c) Logistic Map",
        ylabel = "ΔError (%)")

    ax_lm.ytickformat = y -> string.(round.(y, digits=2), "%")
    boxplot!(ax_lm, cat, values .* 100, color = :orange)

    ##
    ##      Rössler System
    values = dataset[1, 4, :]
    for i = 2:size(dataset, 1)
        values =vcat(values, dataset[i, 4, :])
    end

    ax_r = Axis(fig[2, 2],
        xticks = (1:5, [":square", ":triangle", ":diagonal", ":column", ":pair"]),
        xgridvisible = false,
        ygridvisible = false,
        title = "d) Rössler System",
        ylabel = "ΔError (%)")

    ax_r.ytickformat = y -> string.(round.(y, digits=2), "%")
    boxplot!(ax_r, cat, values .* 100, color = :orange)

    ##
    ##      BSG
    values = dataset[1, 5, :]
    for i = 2:size(dataset, 1)
        values =vcat(values, dataset[i, 5, :])
    end

    ax_b = Axis(fig[3, 1],
        xticks = (1:5, [":square", ":triangle", ":diagonal", ":column", ":pair"]),
        xgridvisible = false,
        ygridvisible = false,
        title = "e) Bernoulli Shifted Generalized",
        ylabel = "ΔError (%)")

    ax_b.ytickformat = y -> string.(round.(y, digits=2), "%")
    boxplot!(ax_b, cat, values .* 100, color = :orange)

    ##
    ##      Overview
    values = reshape(dataset, size(dataset, 1), size(dataset, 2) * size(dataset, 3))
    cat = ones(Int, size(values, 2))
    for i = 2:size(dataset, 1)
        cat = vcat(cat, ones(Int, size(values, 2)) * i)
    end

    values = reshape(values, size(values, 1) * size(values, 2))
    cat = reshape(cat, size(cat, 1) * size(cat, 2))

    ax_o = Axis(fig[3, 2],
        xticks = (1:5, [":square", ":triangle", ":diagonal", ":column", ":pair"]),
        xgridvisible = false,
        ygridvisible = false,
        title = "f) Overview",
        ylabel = "ΔError (%)")

    ax_o.ytickformat = y -> string.(round.(y, digits=2), "%")
    boxplot!(ax_o, cat, values .* 100)

    ##
    save("test/rr_error/fig.png", fig)
end

graph()