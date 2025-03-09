#
#           
#

function square_triangleup_async(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int}, space_size::AbstractVector{Int}, num_samples::Int, func::F, dim::AbstractVector{Int}, hypervolume::Int, metric::Metric) where {F}
    #
    #       Verify the data length.
    if (space_size[1] != space_size[2])
        throw(ArgumentError("To use the sampling mode :triangleup the data_x and data_y must have same length."))
    end
    #
    #       Compute the Power Vector
    p_vect = zeros(Int, hypervolume)
    for i in 1:hypervolume
        p_vect[i] = 2^(i-1)
    end

    rand_segment = (structure[2] + 1):(space_size[2] - (structure[2] - 1))

    #
    #       Creates a function to work as our async task.
    function square_triangleup_task(segment::UnitRange{Int})
        hg = zeros(Int, 2^hypervolume)
        idx = zeros(Int, length(space_size))
        recursive_index = zeros(Int, length(structure))

        @inbounds for _ in segment
            idx[2] = rand(rand_segment)
            idx[1] = rand(1:(idx[2] - structure[1]))

            @fastmath hg[compute_square_index(data_x, data_y, parameters, structure, func, dim, idx, recursive_index, p_vect, metric)] += 1
        end

        return hg
    end

    #       Split the samples between the number of available threads
    int_sampling_value = floor(Int, num_samples / Threads.nthreads())
    rest_sampling_value = num_samples % Threads.nthreads()

    #
    #       Initialize our tasks...
    tasks = []
    start_value = 1
    for _ in 1:Threads.nthreads()
        incrementor = int_sampling_value + (rest_sampling_value > 0 ? 1 : 0)
        segment = start_value:start_value+incrementor - 1

        push!(tasks, Threads.@spawn square_triangleup_task(segment))

        start_value += incrementor
        if (rest_sampling_value > 0)
            rest_sampling_value -= 1
        end
    end

    results = fetch.(tasks)

    res = zeros(Int, 2^hypervolume)
    for r in results
        res .+= r
    end

    return res
end