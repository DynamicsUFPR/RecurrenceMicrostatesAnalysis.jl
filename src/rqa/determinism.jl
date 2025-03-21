#
#
"""
"""
function determinism(probs::Vector{Float64}, rr::Float64)
    if (length(probs) != 512)
        throw(ArgumentError("Determinism must be computed using motifs with n = 3 and sampling mode :triangleup."))
    end

    #       Compute the indexes...
    values = zeros(Int, 64)
    v_idx = 1
    
    for a1 in 0:1, a2 in 0:1, a3 in 0:1, a4 in 0:1, a5 in 0:1, a6 in 0:1
        I = 2^4 + 2 * a1 + 4 * a2 + 8 * a3 + 32 * a4 + 64 * a5 + 128 * a6
        values[v_idx] = I
        v_idx += 1
    end

    #       Sum the probabilities that we need...
    pl = 0.0
    for i in values
        pl += probs[i]
    end

    #       Return the determinism.
    return 1 - ((1/rr) * pl)
end

function determinism(data::AbstractMatrix, parameters; num_samples::Union{Int, Float64} = 0.05, metric::Metric = euclidean_metric)
    #
    #       Compute the distribution.
    dist = distribution(data, parameters, 3; sampling_mode = :triangleup, num_samples = num_samples, metric = metric)
    #
    #       Get the recurrence rate.
    rr = recurrence_rate(dist)
    #
    #       Compute the determinism.
    return determinism(dist, rr)
end