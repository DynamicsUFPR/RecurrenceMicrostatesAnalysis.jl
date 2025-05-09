#
#       RMA Core - Get a random set of microstates available on each column of a recurrence space, using the square shape.
#
"""
    vect_square_columnwise([x], [y], parameters, [structure], space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, samples::Int, metric::Metric)

Get a histogram of a random set of microstates available on each column of a recurrence space. The result is a vector with a probability distribution.
It is only available for 2D recurrence plot !!
"""
function vect_square_columnwise_full(x::Matrix{Float64}, y::Matrix{Float64}, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, metric) where {F}
    
    ##
    ##      Alloc memory for the histogram and the indeces list.
    hg = zeros(Int, 2^hv, space_size[1])
    idx = ones(Int, length(space_size))
    itr = zeros(Int, length(space_size))

    ##
    ##      Compute the power vector.
    p_vect = zeros(Int, hv)
    for i in 1:hv
        p_vect[i] = 2^(i - 1)
    end

    ##
    ##      Get the samples and compute the histogram.
    @inbounds for i in 1:space_size[1]
        idx[1] = i
        for j in 1:space_size[2]
            idx[2] = j
            p = @fastmath compute_index_square(x, y, parameters, structure, func, dim, idx, itr, p_vect, metric)
            hg[p, i] += 1
        end
    end

    ##
    ##      Return the histogram.
    return hg
end

"""
    vect_square_columnwise([x], [y], parameters, [structure], space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, samples::Int, metric::Metric)

Get a histogram of a random set of microstates available on each column of a recurrence space, using an async structure. The result is a vector with a probability distribution.
It is only available for 2D recurrence plot !!
"""
function vect_square_columnwise_full_async(x::Matrix{Float64}, y::Matrix{Float64}, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, metric) where {F}

    ##
    ##      Compute the power vector.
    p_vect = zeros(Int, hv)
    for i in 1:hv
        p_vect[i] = 2^(i - 1)
    end

    ##
    ##      Define a task to compute the histograms.
    function func_task(segment)
        ##
        ##      Alloc memory to the partial histogram, and the indeces.
        hg = zeros(Int, 2^hv, length(segment))
        idx = zeros(Int, length(space_size))
        itr = zeros(Int, length(space_size))

        ##
        ##      Get the samples and compute the histogram.
        @inbounds for i in eachindex(segment)
            idx[1] = segment[i]

            for j in 1:space_size[2]
                idx[2] = j
                p = @fastmath compute_index_square(x, y, parameters, structure, func, dim, idx, itr, p_vect, metric)
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
    res = zeros(Int, 2^hv, space_size[1])
    for r in result
        hg, segment = r
        res[:, segment] .= hg
    end

    ##
    ##      Return the histogram.
    return res
end
