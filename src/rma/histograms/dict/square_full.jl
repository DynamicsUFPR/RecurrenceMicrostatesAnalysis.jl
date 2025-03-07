#
#           
#

function dict_square_full(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int}, space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hypervolume::Int, total_microstates::Int, metric::Metric) where {F}
    #
    #       Alloc memory for histogram and the index list
    hg = Dict{Int, Int}()
    idx = ones(Int, length(space_size))
    recursive_index = zeros(Int, length(structure))
    #
    #       Compute the Power Vector
    p_vect = zeros(Int, hypervolume)
    for i in 1:hypervolume
        p_vect[i] = 2^(i-1)
    end
    #
    #       Do the process...
    @inbounds for _ in 1:total_microstates
        p = @fastmath compute_square_index(data_x, data_y, parameters, structure, func, dim, idx, recursive_index, p_vect, metric)
        hg[p] = get(hg, p, 0) + 1

        idx[1] += 1
        for k in 1:length(idx) - 1
            if (idx[k] > space_size[k])
                idx[k] = 1
                idx[k+1] += 1
            else
                break
            end
        end
    end
    #
    #
    return hg
end