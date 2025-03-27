#
#       RMA Core
#
"""
    vect_pair_columnwise([x], [y], parameters, [structure], space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, samples::Int, metric::Metric)

Get a histogram of a random set of microstates available on each column of a recurrence space. The result is a vector with a probability distribution.
It is only available for 2D recurrence plot !!
"""
function vect_pair_columnwise(x::AbstractVector, y::AbstractVector, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, samples::Int, metric::Metric) where {F}

    ##
    ##      Alloc memory for the histogram and the indeces list.
    hg = zeros(Int, 4, space_size[1])
    idx = ones(Int, length(space_size))
    itr = zeros(Int, length(space_size))

    ##
    ##      Compute the power vector.
    p_vect = zeros(Int, 2)
    for i in 1:2
        p_vect[i] = 2^(i - 1)
    end

    ##
    ##      Get the samples and compute the histogram.
    @inbounds for i in 1:space_size[1]
        idx[1] = i

        for _ in 1:samples
            idx[2] = rand(1:space_size[2])
            p = @fastmath compute_index_pair(x, y, parameters, structure, func, dim, idx, itr, metric)
            hg[p, i] += 1
        end
    end

    ##
    ##      Return the histogram.
    return hg
end

"""
    vect_pair_columnwise_async([x], [y], parameters, [structure], space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, samples::Int, metric::Metric)

Get a histogram of a random set of microstates available on each column of a recurrence space. The result is a vector with a probability distribution.
It is only available for 2D recurrence plot !!
"""
function vect_pair_columnwise_async(x::AbstractVector, y::AbstractVector, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, samples::Int, metric::Metric) where {F}

    ##
    ##      Compute the power vector.
    p_vect = zeros(Int, 2)
    for i in 1:2
        p_vect[i] = 2^(i - 1)
    end

    ##
    ##      Define a task to compute the histograms.
    function func_task(segment)
        ##
        ##      Alloc memory for the histogram and the indeces list.
        hg = zeros(Int, 4, length(segment))
        idx = ones(Int, length(space_size))
        itr = zeros(Int, length(space_size))

        @inbounds for i in 1:length(segment)
            idx[1] = i
            
            for _ in 1:samples
                idx[2] = rand(1:space_size[2])

                ##
                ##      Compute the index and register the motif.
                p = @fastmath compute_index_pair(x, y, parameters, structure, func, dim, idx, itr, metric)
                hg[p, i] += 1
            end
        end

        ##
        ##      Return the partial histogram.
        return hg, segment
    end

    ##
    ##      Split the samples between the number of available threads
    int_columns_value = floor(Int, space_size[1] / Threads.nthreads())
    rest_columns_value = space_size[1] % Threads.nthreads()

    ##
    ##      Initialize our tasks...
    tasks = []
    start_value = 1
    for _ in 1:Threads.nthreads()
        incrementor = int_columns_value + (rest_columns_value > 0 ? 1 : 0)
        segment = start_value:start_value+incrementor - 1

        push!(tasks, Threads.@spawn func_task(segment))

        start_value += incrementor
        if (rest_columns_value > 0)
            rest_columns_value -= 1
        end
    end

    ##
    ##      Wait the result.
    result = fetch.(tasks)

    ##
    ##      Get the results
    res = zeros(Int, 4, space_size[1])
    for r in result
        hg, segment = r
        res[:, segment] .= hg
    end

    ##
    ##      Return the histogram.
    return res
end