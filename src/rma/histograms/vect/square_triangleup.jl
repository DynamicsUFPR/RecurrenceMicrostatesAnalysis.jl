#
#           
#

function square_triangleup(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int}, space_size::AbstractVector{Int}, num_samples::Int, func::F, dim::AbstractVector{Int}, hypervolume::Int, metric::Metric) where {F}
    #
    #       Verify the data length.
    if (space_size[1] != space_size[2])
        throw(ArgumentError("To use the sampling mode :triangleup the data_x and data_y must have same length."))
    end
    #
    #       Alloc memory for histogram and the index list
    hg = zeros(Int, 2^hypervolume)
    idx = zeros(Int, length(space_size))
    recursive_index = zeros(Int, length(structure))
    #
    #       Compute the Power Vector
    p_vect = zeros(Int, hypervolume)
    for i in 1:hypervolume
        p_vect[i] = 2^(i-1)
    end

    segment = (structure[2] + 1):(space_size[2] - (structure[2] - 1))

    #
    #       Do the process...
    @inbounds for _ in 1:num_samples
        idx[2] = rand(segment)
        idx[1] = rand(1:(idx[2] - structure[1]))

        @fastmath hg[compute_square_index(data_x, data_y, parameters, structure, func, dim, idx, recursive_index, p_vect, metric)] += 1
    end
    #
    #
    return hg   
end