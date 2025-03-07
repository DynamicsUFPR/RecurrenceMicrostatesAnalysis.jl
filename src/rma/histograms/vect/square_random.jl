#
#           
#

function square_random(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int}, space_size::AbstractVector{Int}, num_samples::Int, func::F, dim::AbstractVector{Int}, hypervolume::Int, metric::Metric) where {F}
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
    
    #
    #       Do the process...
    @inbounds for _ in 1:num_samples
        for s in eachindex(space_size)
            idx[s] = rand(1:space_size[s])
        end

        @fastmath hg[compute_square_index(data_x, data_y, parameters, structure, func, dim, idx, recursive_index, p_vect, metric)] += 1
    end
    #
    #
    return hg
end