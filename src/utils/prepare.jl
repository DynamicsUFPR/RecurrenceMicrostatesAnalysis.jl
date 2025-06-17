
"""
    prepare([solution], σ::Union{Float64, Int}; transient::Int, K::Int)

Prepare a problem solved by the library `DifferentialEquations.jl` to be used in `RecurrenceMicrostatesAnalysis.jl`. 
This function applies the sampling parameter (σ) to discretize the continuous time series, as proposed by [Thiago2024](@cite).

Input:
*   `[solution]`: solution returned by the library `DifferentialEquations.jl`.
*   `σ`: sampling parameter; it defines the time resolution of discretized data.
*   `transient` **(kwarg)**: number of points, without application of sampling, that will be ignored.
*   `K` **(kwarg)**: maximum length of the returned data series.

Output:
*   `data`: returns the prepared data in the format of a `Matrix{Float64}`. Each row represents a system component, and each column represents a time step.
"""
function prepare(
        solution, σ::Union{Float64, Int}; 
        transient::Int = 0, 
        K::Int = 0
    )

    time = solution.t

    if (transient >= length(time))
        throw(ArgumentError("Transient value is greater than solution size."))
    end
    
    pos = (transient+1):length(time)
    new_pos = [1]

    time = time[pos]
    for i in eachindex(time)
        if (i == 1)
            continue
        end

        if (time[new_pos[length(new_pos)]] + σ <= time[i])
            push!(new_pos, i)
        end
    end

    pos =  transient .+ new_pos
    data = (solution[:, :])[:, pos]

    if (size(data, 2) < K)
        println("Warning: the result data series is shorter than the maximum length.")
        return data
    elseif (K == 0)
        return data
    end

    return (solution[:, :])[:, pos[1:K]]
end