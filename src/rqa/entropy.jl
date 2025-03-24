#
#           RMA Analysis
#
"""
    rentropy([probs]; {ignore_motifs})

Compute the recurrence entropy, proposed by [Corso2018](@cite). `ignore_motifs::Vector{Int}` defines a set of
motif's indeces to be ignored while the function compute the entropy.
"""
function rentropy(probs::Vector{Float64}; ignore_motifs::Vector{Int} = [-1])
    s::Float64 = 0.0

    for i in eachindex(probs)
        if (any(x -> x == i, ignore_motifs))
            continue
        end

        if (probs[i] > 0)
            s += (-1) * probs[i] * log(probs[i])
        end
    end

    return s
end

function rentropy(probs::Dict{Int, Float64}; ignore_motifs::Vector{Int} = [-1])
    s::Float64 = 0.0

    for i in collect(keys(probs))
        if (any(x -> x == i, ignore_motifs))
            continue
        end

        if (probs[i] > 0)
            s += (-1) * probs[i] * log(probs[i])
        end
    end

    return s
end