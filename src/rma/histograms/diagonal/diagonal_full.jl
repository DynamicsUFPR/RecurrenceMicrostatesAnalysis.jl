#
#       RMA Core - Get all microstates available, using the square shape.
#
function vect_diagonal_full_async(x::Matrix{Float64}, y::Matrix{Float64}, parameters, structure::AbstractVector{Int},
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
                p = @fastmath compute_index_diagonal(x, y, parameters, structure, func, dim, idx, itr, p_vect, metric)
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