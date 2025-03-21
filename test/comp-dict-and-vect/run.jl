#
#
#       This code tests the performance between the run modes :dict and :vect.
#       How it generate some graphics, you will need the libraries CairoMakie and JLD2 installed
#   to run this code.
#
#       It code define two functions: compute and graph, compute generate the data and save as a
#   JLD2 file that will be load by the graph function to make a visual representation of that data.
#
using RMA
using JLD2
using CairoMakie
using ProgressMeter
using BenchmarkTools
#
const max_samples = 10000
const intervals = range(0.01, 1.0, 20)
const n_values = [2, 3, 5]
#
#       Generate a random dataset.
const data = rand(Float64, 1, max_samples)
#
#       Alloc memory to store the results...
result = zeros(Float64, (length(intervals), length(n_values), 2))
#
#       I need to declare some global variables to benchmark works.
n_applied = 0
i_applied = 0
#
#       Calculate the mean time of process.
function calc_time(n, i)
    global n_applied = n_values[n]
    global i_applied = intervals[i]

    dict_test = @benchmark distribution(data, 0.1, n_applied; run_mode = :dict, num_samples = i_applied)
    result[i, n, 1] = dict_test.times |> mean

    vect_test = @benchmark distribution(data, 0.1, n_applied; run_mode = :vect, num_samples = i_applied)
    result[i, n, 2] = vect_test.times |> mean
end
#
#
function compute()
    @showprogress for i in eachindex(intervals)
        for n in eachindex(n_values)
            calc_time(n, i)
        end
    end

    result ./= 1e6

    save_object("test/comp-dict-and-vect/data.jld2", result)
end