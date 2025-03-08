#
#           
#

function square_columnwise_async(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int}, space_size::AbstractVector{Int}, num_samples::Int, func::F, dim::AbstractVector{Int}, hypervolume::Int, metric::Metric) where {F}
    #
    #       Compute the Power Vector
    p_vect = zeros(Int, hypervolume)
    for i in 1:hypervolume
        p_vect[i] = 2^(i-1)
    end

    #
    #       Creates a function to work as our async task.
    function square_columnwise_task(segment::UnitRange{Int})
        hg = zeros(Int, 2^hypervolume, length(segment))
        idx = zeros(Int, length(space_size))
        recursive_index = zeros(Int, length(structure))

        @inbounds for i in 1:length(segment)
            idx[1] = i
            for _ in 1:num_samples
                idx[2] = rand(1:space_size[2])
                @fastmath hg[compute_square_index(data_x, data_y, parameters, structure, func, dim, idx, recursive_index, p_vect, metric), i] += 1
            end
        end
        
        return hg, segment
    end

    #       Split the samples between the number of available threads
    int_columns_value = floor(Int, space_size[1] / Threads.nthreads())
    rest_columns_value = space_size[1] % Threads.nthreads()

    #
    #       Initialize our tasks...
    tasks = []
    start_value = 1
    for _ in 1:Threads.nthreads()
        incrementor = int_columns_value + (rest_columns_value > 0 ? 1 : 0)
        segment = start_value:start_value+incrementor - 1

        push!(tasks, Threads.@spawn square_columnwise_task(segment))

        start_value += incrementor
        if (rest_columns_value > 0)
            rest_columns_value -= 1
        end
    end

    results = fetch.(tasks)
    res = zeros(Int, 2^hypervolume, space_size[1])
    for r in results
        hg, segment = r
        res[:, segment] .= hg
    end

    return res
end