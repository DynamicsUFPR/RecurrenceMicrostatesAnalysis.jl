#
#           
#

function square_columnwise(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int}, space_size::AbstractVector{Int}, num_samples::Int, func::F, dim::AbstractVector{Int}, hypervolume::Int, metric::Metric) where {F}
    #
    #       Alloc memory for histogram and the index list
    hg = zeros(Int, 2^hypervolume, space_size[1])
    idx = zeros(Int, length(space_size))
    recursive_index = zeros(Int, length(structure))
    #
    #       Compute the Power Vector
    p_vect = zeros(Int, hypervolume)
    for i in 1:hypervolume
        p_vect[i] = 2^(i-1)
    end

    #
    #       Do the process...
    @inbounds for i in 1:space_size[1]
        idx[1] = i
        for _ in 1:num_samples
            idx[2] = rand(1:space_size[2])
            @fastmath hg[compute_square_index(data_x, data_y, parameters, structure, func, dim, idx, recursive_index, p_vect, metric), i] += 1
        end
    end
    #
    #
    return hg
end