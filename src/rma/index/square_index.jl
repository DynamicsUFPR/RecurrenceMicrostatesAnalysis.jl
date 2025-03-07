#
#
#
"""
    
"""
function compute_square_index(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, fixed::Vector{Int}, recursive::Vector{Int}, power_vector::Vector{Int}, metric::Metric) where {F}
    index = 0

    #       Reset the recursive register.
    copy!(recursive, fixed)

    for m in 1:length(power_vector)
        if  @inline func(data_x, data_y, parameters, recursive, dim, metric)
            index += power_vector[m]
        end

        recursive[1] += 1
        for k in 1:length(structure) - 1
            if (recursive[k] >= fixed[k] + structure[k])
                recursive[k] = fixed[k]
                recursive[k + 1] += 1
            else
                break
            end
        end
    end

    return 1 + index
end