
"""
    rentropy([probs]; [ignore_motifs])

Compute the recurrence entropy, as proposed by [Corso2018](@cite).

Input:
* `[probs]`: a `Vector{Float64}` returned by the function `distribution(...)`.
* `[ignore_motifs]` **(kwarg)**: list of motifs to ignore.

Output: return the recurrence entropy as a `Float64`.
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