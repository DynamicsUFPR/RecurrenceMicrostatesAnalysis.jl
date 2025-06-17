
"""
    find_parameters([x], n::Int; r::Float64 = 0.05, ε_max_range = 0.5)

This function calculates the maximum microstate entropy for the RP of the input
time series given the microstate size `n` and the input ratio `r` of the
total number of microstates.

Input:
* `[x]`: input data.
* `n`: microstate size.
* `r` **(kwarg)**: ratio of the total number of microstates to be sampled for the histogram (default `r = 0.05`)
* `ε_max_range` **(kwarg)**: percentage of the maximum distance to be used as the range in the process. (default `ε_max_range = 0.5`)
* `fraction` **(kwarg)**: iteration fraction. (default `fraction = 5`)
* `shape` **(kwarg)**: motif shape. (default `shape = :square`)

Output: (is a `Tuple{Float64, Float64}`)
* `εopt`: the value of the vicinity parameter that maximizes the recurrence microstates entropy.
* `Smax`: maximum recurrence microstates entropy for the input time series.
"""
function find_parameters(
        x::AbstractArray,
        n::Int;
        r::Float64 = 0.05,
        ε_max_range::Float64 = 0.5,
        fraction::Int = 5,
        shape::Symbol = :square
    )
    
    ε::Float64 = 1e-6
    εopt::Float64 = maximum(pairwise(Euclidean(), eachrow(x))) * (ε_max_range - ε)
    Δε::Float64 = (εopt - ε) / fraction

    Smax::Float64 = 0.0

    ##      Just to allocate memory.
    dist::Vector{Float64} = distribution(x, 0.01, n; num_samples = 2, shape = shape)

    for _ ∈ 1:fraction
        for _ ∈ 1:fraction
            ##      Compute the recurrence entropy.
            dist .= distribution(x, ε, n; num_samples = r, shape = shape)
            S = rentropy(dist)

            if S > Smax
                Smax = S
                εopt = ε
            end

            ε += Δε
        end

        ε = εopt - Δε
        Δε *= 2 / fraction
    end

    return εopt, Smax
end