
"""
    rrate([probs])

Compute the approximated recurrence rate of a RP from a probability distribution of recurrence microstates. Here, we use a relation between
the mean recurrence rate of each motif and the desired value. It can be written as

```math
RR \\approx \\sum_{I = 0}^N \\mathbf{p}_I^{(k)}\\left(\\frac{1}{k^2}\\sum_{i=1}^k\\sum_{j=1}^k \\mathbf{M}_{ij}^{(I)}\\right),
```

where \$\\mathbf{M}_{ij}^{(I)}\$ is the motif structure.

Input:
* `[probs]`: the vector of probabilities \$\\mathbf{p}^{(k)}\$ computed using `distribution(...)`.

Output: returns the recurrence rate as a `Float64`.
"""
function rrate(
        probs::Vector{Float64}
    )

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
    rrate([x], [parameters], n::Int; shape::Symbol, sampling_mode::Symbol, r::Float64})

Compute the approximated recurrence rate of a Recurrence Plot from a data `[x]` using a probability
distribution of recurrence microstates computed from it.

Input:
* `[x]`: input data.
* `[parameter]`: set of parameters used to compute the recurrence microstate distribution.
* `n`: microstate size.
* `shape` **(kwarg)**: shape of the used motifs. `:square` by default, it can be: `:square, :triangle, :pair, :diagonal, :line`.
* `sampling_mode` **(kwarg)**: sampling mode used. `:random` by default, it can be: `:square, :triangle, :pair, :diagonal, :line`.
* `r` **(kwarg)**: ratio of the total number of microstates to be sampled for the histogram. (default `r = 0.05`)

Output: returns the recurrence rate as a `Float64`.
"""
function rrate(
        x::AbstractArray, 
        parameters, 
        n::Int; 
        shape::Symbol = :square, 
        sampling_mode::Symbol = :random, 
        r::Float64= 0.05
    )
    
    if (sampling_mode == :columnwise || sampling_mode == :columnwise_full)
        throw("The sampling $(sampling_mode) is invalid.")
    end

    probs = distribution(x, parameters, n; shape = shape, sampling_mode = sampling_mode, num_samples = r)
    return rrate(probs)
end