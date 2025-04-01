#
#       RMA Utils
#
"""
    prepare([solution], vicinity::Union{Float64, Int}; transient::Int, max_length::Int)

Prepare a problem solved by the library DifferencialEquations.jl to be used in RMA.jl.
It can apply the vicinity parameter to discretize the "continuous" time serie, like proposed
by [Thiago2023](@cite).
"""
function prepare(solution, vicinity::Union{Float64, Int}; transient::Int = 0, max_length::Int = 0)
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

        if (time[new_pos[length(new_pos)]] + vicinity <= time[i])
            push!(new_pos, i)
        end
    end

    pos =  transient .+ new_pos
    data = (solution[:, :])[:, pos]

    if (max_length == 0 || size(data, 2) < max_length)
        return data
    end

    return (solution[:, :])[:, pos[1:max_length]]
end