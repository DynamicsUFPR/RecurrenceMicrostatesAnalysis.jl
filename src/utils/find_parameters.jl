#
#       RMA Utils
#
#       TODO - Improve it. The current version is only for tests.
function find_parameters(x::AbstractArray, n::Int; threshold_min::Float64 = 0.0, threshold_max::Float64 = maximum(pairwise(Euclidean(), x, x)), large_prec = 20, small_prec = 50, num_samples = 0.05, shape = :square)
    ##
    ##      Define a range of thresholds and alloc memory to store the entropy.
    threshold_range = range(threshold_min, threshold_max, large_prec)
    s_max = 0.0
    t_indeces = zeros(Int, 3)

    ##
    ##      Compute the max entropy for a "large precision"
    for i in eachindex(threshold_range)
        dist = distribution(x, threshold_range[i], n; num_samples = num_samples, shape = shape)
        s = rentropy(dist)

        if (s > s_max)
            t_indeces[1] = t_indeces[2]
            t_indeces[2] = i
            s_max = s
        elseif (s < s_max)
            t_indeces[3] = i
            break
        end
    end

    if (t_indeces[3] == 0)
        throw(ArgumentError("The max entropy is not between the values of 'threshold_min' and 'threshold_max'."))
    end

    ##
    ##      Increase the precision.
    threshold_range = range(threshold_range[t_indeces[1]], threshold_range[t_indeces[3]], small_prec)
    s_max = 0.0
    t_indeces = zeros(Int, 3)

    for i in eachindex(threshold_range)
        dist = distribution(x, threshold_range[i], n)
        s = rentropy(dist)

        if (s > s_max)
            t_indeces[1] = t_indeces[2]
            t_indeces[2] = i
            s_max = s
        elseif (s <= s_max)
            t_indeces[3] = i
            break
        end
    end

    return threshold_range[t_indeces[2]], s_max
end


"""
***THIS IS MERELY A TEMPLATE THAT MAY OR MAY NOT BE USED***

This function calculates the maximum microstate entropy for the RP of the input
time series given the microstate size `k` and the input ratio `r` of the
total number of microstates.

Input:
* `x`: time series with `N` data points.
* `k` **(kwarg)**: number of columns of the microsatate (default to `k = 2`).
* `r` **(kwarg)**: ratio of the total number of microstates to be sampled
for the histogram (default to `r = 0.05`)

Output:
* `Smax`: maximum recurrence microstates entropy for the input time series.
* `εopt`: the value of the vicinity parameter that maximizes the recurrence
microstates entropy.
"""
# function find_parameters_test(
#         x::VecOrMat{Float64};
#         k::Int=2,
#         r::Float64=0.05,
#     )

#     fraction_1::UInt8 = 5
#     fraction_2::UInt8 = 1 * fraction_1
    
#     ε::Float64 = 1e-6
#     εopt::Float64 = maximum(pairwise(Euclidean(), eachrow(x))) * (0.5 - ε)
#     Δε::Float64 = (εopt - ε) / fraction_2

#     Smax::Float64 = 0.0

#     for _ ∈ 1:fraction_1
#         for _ ∈ 1:fraction_2
#             S = # CALCULATE THE ENTROPY

#             if S > Smax
#                 Smax = S
#                 εopt = ε
#             end

#             ε += Δε
#         end

#         ε = εopt - Δε
#         Δε *= 2 / fraction_2
#     end

#     return Smax, εopt
# end
