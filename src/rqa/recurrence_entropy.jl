#
#
#
function recurrence_entropy(probs::Vector{Float64}; ignore_motifs::AbstractVector = [])
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
#
#
#
function recurrence_entropy(probs::Dict{Int, Float64}; ignore_motifs::AbstractVector = [])
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