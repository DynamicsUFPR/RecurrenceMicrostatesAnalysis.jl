#
#
"""
"""
function recurrence_rate(probs::Vector{Float64})
    result = 0.0
    hv = Int(log2(length(probs)))

    #       Compute the recurrence rate.
    for i in eachindex(probs)
        rr = sum(digits(i, base = 2)) / hv
        result += rr * probs[i]
    end

    return result
end