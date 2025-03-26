#
#           RMA Analysis
#
"""
    rrate([probs])

Compute the approximated recurrence rate of a Recurrence Plot from a probability distribution of recurrence microstates.
"""
function rrate(probs::Vector{Float64})
    result = 0.0
    hv = Int(log2(length(probs)))

    ##
    ##       Compute the recurrence rate.
    for i in eachindex(probs)
        rr = sum(digits(i - 1, base = 2)) / hv
        result += rr * probs[i]
    end

    return result
end

"""
    rrate([x], parameters, n::Int; {shape}, {sampling_mode}})

Compute the approximated recurrence rate of a Recurrence Plot from a dataset [x] using a probability distribution of recurrence microstates computed
from it.
"""
function rrate(x::AbstractArray, parameters, n::Int; shape::Symbol = :square, sampling_mode::Symbol = :random, num_samples::Union{Float64, Int} = 0.05)
    ##
    ##      Compute the probabilities.
    probs = distribution(x, parameters, n; shape = shape, sampling_mode = sampling_mode, num_samples = num_samples)

    return rrate(probs)
end