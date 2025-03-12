#
#
#
"""

"""
function compute_triangle_index(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, fixed::Vector{Int}, recursive::Vector{Int}, power_vector::Vector{Int}, metric::Metric) where {F}
    index = 0
    pvec_index = 1

    #       Reset the recursive register.
    copy!(recursive, fixed)

    for j in 0:structure[2] - 1
        recursive[2] = fixed[2] + j
        for i in j:structure[1] - 1
            recursive[1] = fixed[1] + i
            if  @inline func(data_x, data_y, parameters, recursive, dim, metric)
                index += power_vector[pvec_index]
            end
            pvec_index += 1
        end
    end

    return 1 + index
end