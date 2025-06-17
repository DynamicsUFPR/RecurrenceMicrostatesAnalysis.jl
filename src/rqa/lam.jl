
"""
    laminarity(rr::Float64, [probs])

Estimate the laminarity from a distribution. If the distribution has 512 elements, the function will consider square motifs, computing laminarity using

```math
I^{(\\beta)} = \\frac{1}{RR}(2 + 8\\mathbf{M}_{2,1}^{(\\beta)} + 16\\mathbf{M}_{2,2}^{(\\beta)} + 32\\mathbf{M}_{2,3}^{(\\beta)} + 64\\mathbf{M}_{3,1}^{(\\beta)} + 128\\mathbf{M}_{3,2}^{(\\beta)} + 256\\mathbf{M}_{3,3}^{(\\beta)}),
```

where \$\\mathbf{M}\$ is the motif structure. It defines a class of motifs \$(C_L) \\ni I^{(\\beta)}\$ that we use to estimate LAM:

```math
LAM \\approx 1 - \\frac{1}{RR}\\sum_{I\\in (C_L)} \\mathbf{p}^{(3)}_I.
```

```julia
dist = distribution(data, th, 3)
rr = rrate(dist)
lam = laminarity(rr, dist)
```

If the distribution has 8 elements, this function will consider line motifs, which makes the process simpler. In this case, we just need the motif with \$I = 2\$:

```math
LAM \\approx 1 - \\frac{\\mathbf{p}^{(3)}_2}{RR}.
```

```julia
dist = distribution(data, th, 3; shape = :line)
rr = rrate(dist)
lam = laminarity(rr, dist)
```

Input:
* `rr`: recurrence rate.
* `[probs]`: a `Vector{Float64}` returned by the function `distribution(...)`.

Output: returns the laminarity as a `Float64`.
"""
function laminarity(
        rr::Float64, 
        probs::Vector{Float64}
    )

    if (length(probs) != 512 && length(probs) != 8)
        throw(ArgumentError(string("Determinism must be computed using square motifs with n = 3. Actual value results in ", length(probs))))
    end

    function __default_mode()
        ##
        ##      Compute the indeces...
        values = zeros(Int, 64)
        v_idx = 1

        for a1 in 0:1, a2 in 0:1, a3 in 0:1, a4 in 0:1, a5 in 0:1, a6 in 0:1
            I_1 = 2 + 8 * a1 + 16 * a2 + 32 * a3 + 64 * a4 + 128 * a5 + 256 * a6
            values[v_idx] = I_1 + 1
            v_idx += 1
        end

        ##
        ##      Sum the probabilities that we need.
        pl = 0.0
        for i in values
            pl += probs[i]
        end

        ##
        ##      Return the determinism.
        return 1 - ((1/rr) * pl)
    end

    function  __line_mode()
        ##
        ##      With line mode we need only a specific probability, the motif (0 1 0) with index 2 (3 in Julia).
        return 1 - ((1/rr) * probs[3])
    end

    ##
    ##      Compute the determinism.
    return length(probs) != 512 ? __line_mode() : __default_mode()
end
"""
    laminarity([x], threshold::Float64; r::Float64)

Estimate the laminarity from a data `[x]`` using a probability distribution and a RR computed from it. 

Input:
* `[x]`: input data.
* `[parameter]`: set of parameters used to compute the recurrence microstate distribution.
* `r` **(kwarg)**: ratio of the total number of microstates to be sampled for the histogram. (default `r = 0.05`)

Output: returns a `Tuple{Float64, Float64}`.
* `lam`: laminarity as `Float64`.
* `rr`: recurrence rate as `Float64`.
"""
function laminarity(
        x::AbstractArray, 
        threshold::Float64; 
        r::Float64 = 0.05
    )
    
    ##
    ##      Compute the probabilities.
    probs = distribution(x, (0.0, threshold), 3; shape = :line, num_samples = r)
    ##
    ##      Compute the recurrence rate.
    rr = rrate(probs)
    ##
    ##      Return the determinism.
    return laminarity(rr, probs), rr
end