#
#       RMA Core - Get all microstates available, using the square shape.
#
"""
    dict_square_full([x], [y], parameters, [structure], space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, total_microstates::Int, metric::Metric)

Get a histogram of all microstates available on a recurrence space. The result is a dict with a probability distribution.
"""
function dict_square_full(x::AbstractArray, y::AbstractArray, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, total_microstates::Int, metric::Metric) where {F}
    
    ##
    ##      Alloc memory for the histogram and the indeces list.
    hg = Dict{Int, Int}()
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
    @inbounds for _ in 1:total_microstates
        ##
        ##      Compute the index and register the motif.
        p = @fastmath compute_index_square(x, y, parameters, structure, func, dim, idx, itr, p_vect, metric)
        hg[p] = get(hg, p, 0) + 1

        ##
        ##      Move to the next index.
        idx[1] += 1
        for k in 1:length(idx) - 1
            if (idx[k] > space_size[k])
                idx[k] = 1
                idx[k + 1] += 1
            else
                break
            end
        end
    end

    ##
    ##      Return the histogram.
    return hg
end

##      TODO : adapt to the spatial generalization.
"""
    dict_square_full_async([x], [y], parameters, [structure], space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, metric::Metric)

Get a histogram of all microstates available on a recurrence space, using an async structure. The result is a dict with a probability distribution.
It is only available for 2D recurrence plot !!
"""
function dict_square_full_async(x::Matrix{Float64}, y::Matrix{Float64}, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, metric::Metric) where {F}

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
        hg = Dict{Int, Int}()
        idx = zeros(Int, length(space_size))
        itr = zeros(Int, length(space_size))

        @inbounds for i in segment
            ##
            ##      Fix the first index =D
            idx[1] = i
            for j in 1:space_size[2]
                idx[2] = j

                ##
                ##      Compute the index and register the motif.
                p = @fastmath compute_index_square(x, y, parameters, structure, func, dim, idx, itr, p_vect, metric)
                hg[p] = get(hg, p, 0) + 1
            end
        end

        ##
        ##      Return the partial histogram.
        return hg
    end

    ##
    ##      Split the rows between the number of available threads.
    int_rows_value = floor(Int, space_size[1] / Threads.nthreads())
    rest_rows_value = space_size[1] % Threads.nthreads()

    ##
    ##      Initialize our tasks...
    tasks = []
    start_value = 1
    for _ in 1:Threads.nthreads()
        incrementor = int_rows_value + (rest_rows_value > 0 ? 1 : 0)
        segment = start_value:start_value + incrementor - 1

        push!(tasks, Threads.@spawn func_task(segment))

        start_value += incrementor
        rest_rows_value -= 1
    end

    ##
    ##      Wait the result.
    result = fetch.(tasks)

    ##
    ##      Get the results
    res = Dict{Int, Int}()
    for r in result
        for (k, v) in r
            res[k] = get(res, k, 0) + v
        end
    end

    ##
    ##      Return the histogram.
    return res
end

"""
    vect_square_full([x], [y], parameters, [structure], space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, total_microstates::Int, metric::Metric)

Get a histogram of all microstates available on a recurrence space. The result is a vect with a probability distribution.
"""
function vect_square_full(x::AbstractArray, y::AbstractArray, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, total_microstates::Int, metric::Metric) where {F}

    ##
    ##      Alloc memory for the histogram and the indeces list.
    hg = zeros(Int, 2^hv)
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
    @inbounds for _ in 1:total_microstates
        ##
        ##      Compute the index and register the motif.
        p = @fastmath compute_index_square(x, y, parameters, structure, func, dim, idx, itr, p_vect, metric)
        hg[p] += 1

        ##
        ##      Move to the next index.
        idx[1] += 1
        for k in 1:length(idx) - 1
            if (idx[k] > space_size[k])
                idx[k] = 1
                idx[k + 1] += 1
            else
                break
            end
        end
    end

    ##
    ##      Return the histogram.
    return hg
end

##      TODO : adapt to the spatial generalization.
"""
    vect_square_full_async([x], [y], parameters, [structure], space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, metric::Metric)

Get a histogram of all microstates available on a recurrence space, using an async structure. The result is a vector with a probability distribution.
It is only available for 2D recurrence plot !!
"""
function vect_square_full_async(x::Matrix{Float64}, y::Matrix{Float64}, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, metric::Metric) where {F}

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
        hg = zeros(Int, 2^hv)
        idx = zeros(Int, length(space_size))
        itr = zeros(Int, length(space_size))

        @inbounds for i in segment
            ##
            ##      Fix the first index =D
            idx[1] = i
            for j in 1:space_size[2]
                idx[2] = j

                ##
                ##      Compute the index and register the motif.
                p = @fastmath compute_index_square(x, y, parameters, structure, func, dim, idx, itr, p_vect, metric)
                hg[p] += 1
            end
        end

        ##
        ##      Return the partial histogram.
        return hg
    end

    ##
    ##      Split the rows between the number of available threads.
    int_rows_value = floor(Int, space_size[1] / Threads.nthreads())
    rest_rows_value = space_size[1] % Threads.nthreads()

    ##
    ##      Initialize our tasks...
    tasks = []
    start_value = 1
    for _ in 1:Threads.nthreads()
        incrementor = int_rows_value + (rest_rows_value > 0 ? 1 : 0)
        segment = start_value:start_value + incrementor - 1

        push!(tasks, Threads.@spawn func_task(segment))

        start_value += incrementor
        rest_rows_value -= 1
    end

    ##
    ##      Wait the result.
    result = fetch.(tasks)

    ##
    ##      Get the results
    res = zeros(Int, 2^hv)
    for r in result
        res .+= r
    end

    ##
    ##      Return the histogram.
    return res
end