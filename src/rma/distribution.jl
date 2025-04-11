#
#           RMA Core - Public function to compute the probability distribution of some datasets.
#
# """
# ## Keywords arguments
# - `shape::Symbol`: can be `:square`, `:triangle`, `:pair`, `:diagonal` and `:line`. The value `:square` refers to 
# default square format of motifs, based on the work of [Corso2018](@cite) with a generalization for spatial data based 
# on the work of [Marwan2006](@cite). On  the other hand, the value `:triangle` is only available for 2-dimensional 
# recurrence spaces (i.e., time series) and is based on the work of [Hirata2021](@cite). `:pair`, `:diagonal` and 
# `:line` are experimental shapes and there are not an work about them yet.

# - `run_mode::Symbol`: can be `:default`, `:dict`, or `:vect`. `:dict` and `:vect` set the return format to vector and
# dictionary respectively. The mode `:default` uses a vector for motifs with an hypervolume lesser than 28 and 
# a dictionary otherwise.

# - `sampling_mode::Symbol`: can be `:full`, `:random`, `:columnwise` or `:triangleup`. The sampling mode `:full` retrieves all 
# available microstates in the recurrence space sequentially; it is not yet available for multi-threading. The `:random`
# mode retrieves a number of samples based on the value of `num_samples` and creates a distribution for the entire 
# recurrence space. Finally, `:columnwise` creates a distribution for each column of the recurrence space. It is only 
# available for `:vect` run mode because it returns an array where each column is a probability distribution of motifs.

# - `num_samples::Union{Int, Float64}`: the number of samples used to create the probability distribution. CIt can be a 
# percentage of the total possible motifs or the desired number of samples.

# - `metric::Metric`: a structure used to compute distances from the [Distances.jl](https://github.com/JuliaStats/Distances.jl) 
# library.

# - `func = (x, y, p, ix, dim, metric) -> recurrence(x, y, p, ix, dim, met)`: this is the recurrence function that returns true when 
# there is a recurrence and false otherwise. This function takes 5 parameters:
#     - x::AbstractArray : it is a reference to the parameter `data_x`.
#     - y::AbstractArray : it is a reference to the parameter `data_y`.
#     - p : it is a reference to the parameter `parameters`.
#     - [ix] : it is a vector with `length(structure)` elements that indicates the elements used to compute the 
# recurrence.
#     - [dim] : it is a vector with two elements: the number of dimenions of `data_x` and `data_y`. It is used only for
# high-dimensional problems.
#     - metric::Metric : a metric from the [Distances.jl](https://github.com/JuliaStats/Distances.jl) library. It used
#     the same value as the `metric` parameter of distribution function.

# For reference, you can see the code [`src/rma/recurrence.jl`](src/rma/recurrence.jl).

# - `threads::Bool`: defines whether the library will use threads in the computational process. If set to true, the 
# library will use the number of available threads defined by `JULIA_NUM_THREADS`.
"""
### Based on Recurrence Plot

        distribution([x], parameters, n::Int; kwords...)

Compute the distribution of recurrence microstates probabilities from the dataset `x`. The input `parameters`
consists of the constant values used to calculate the recurrence between two points. `n` is an integer that
represents the length of motifs side.

This function can return a vector, an array or a dictionary based on the number of possible microstates and
the setting of `run_mode` or `sampling_mode`.
"""
function distribution(x::AbstractArray, parameters, n::Int;
    shape::Symbol = :square, run_mode::Symbol = :default, sampling_mode::Symbol = :random,
    num_samples::Union{Int, Float64} = 0.05, threads::Bool = Threads.nthreads() > 1, 
    metric = euclidean_metric, func = (x, y, p, ix, metric, dim) -> recurrence(x, y, p, ix, metric, dim))
    
    if (ndims(x) == 1)
        x = Matrix(x')
    end
    
    dim = ndims(x) - 1
    structure = ones(Int, 2 * dim) .* n

    return distribution(x, x, parameters, structure; shape = shape, run_mode = run_mode, sampling_mode = sampling_mode, num_samples = num_samples, threads = threads, metric = metric, func = func)
end
"""
---
### Based on Cross-Recurrence Plot

    distribution([x], [y], parameters, n::Int; kwords...)

Compute the distribution of recurrence microstates probabilities from the datasets `x` and `y`. The input `parameters`
consists of the constant values used to calculate the recurrence between two points. `n` is an integer that
represents the length of motifs side.

This function can return a vector, an array or a dictionary based on the number of possible microstates and
the setting of `run_mode` or `sampling_mode`.
"""
function distribution(x::AbstractArray, y::AbstractArray, parameters, n::Int;
    shape::Symbol = :square, run_mode::Symbol = :default, sampling_mode::Symbol = :random,
    num_samples::Union{Int, Float64} = 0.05, threads::Bool = Threads.nthreads() > 1, 
    metric = euclidean_metric, func = (x, y, p, ix, metric, dim) -> recurrence(x, y, p, ix, metric, dim))

    if (ndims(x) == 1)
        x = Matrix(x')
    end
    if (ndims(y) == 1)
        y = Matrix(y')
    end

    dim_x = ndims(x) - 1
    dim_y = ndims(y) - 1
    structure = ones(Int, dim_x + dim_y) .* n

    return distribution(x, y, parameters, structure; shape = shape, run_mode = run_mode, sampling_mode = sampling_mode, num_samples = num_samples, threads = threads, metric = metric, func = func)
end
"""
---
### Using DifferencialEquations.jl

    distribution([solution], parameters, n::Int, vicinity::Union{Int, Float64}; kwords...)

Compute the distribution of recurrence microstates probabilities from the `solution` of a differencial equation solved by 
the library DifferencialEquations.jl. The input `parameters` consists of the constant values used to calculate the recurrence 
between two points. `n` is an integer that represents the length of motifs side. `vicinity` is the time separation used to
discretize a continuous problem.

It presents some new kwords:
- `transient::Int`: defines an interval of time that will be ignored, and taked as a transient.
- `max_length::Int`: defines the maximum size of the result time series.

This function can return a vector, an array or a dictionary based on the number of possible microstates and
the setting of `run_mode` or `sampling_mode`.
"""
function distribution(solution, parameters, n::Int, vicinity::Union{Int, Float64} = 0.0;
    shape::Symbol = :square, run_mode::Symbol = :default, sampling_mode::Symbol = :random,
    num_samples::Union{Int, Float64} = 0.05, threads::Bool = Threads.nthreads() > 1, 
    metric = euclidean_metric, func = (x, y, p, ix, metric, dim) -> recurrence(x, y, p, ix, metric, dim),
    transient::Int = 0, max_length::Int = 0)
    
    ##
    ##      Prepare the dataset.
    data = prepare(solution, vicinity; transient = transient, max_length = max_length)


    dim = ndims(data) - 1
    structure = ones(Int, 2 * dim) .* n

    return distribution(data, data, parameters, structure; shape = shape, run_mode = run_mode, sampling_mode = sampling_mode, num_samples = num_samples, threads = threads, metric = metric, func = func)
end
"""
---
### Main Version

    distribution([x], [y], parameters, [structure]; kwords...)

Compute the distribution of recurrence microstates probabilities from the datasets `x` and `y`. The input `parameters`
consists of the constant values used to calculate the recurrence between two points. Meanwhile, the input `structure`
is a vector where each element represents a side of the motif.

This function can return a vector, an array or a dictionary based on the number of possible microstates and
the setting of `run_mode` or `sampling_mode`.
"""
function distribution(x::AbstractArray, y::AbstractArray, parameters, structure::AbstractVector{Int};
    shape::Symbol = :square, run_mode::Symbol = :default, sampling_mode::Symbol = :random,
    num_samples::Union{Int, Float64} = 0.05, threads::Bool = Threads.nthreads() > 1, 
    metric = euclidean_metric, func = (x, y, p, ix, metric, dim) -> recurrence(x, y, p, ix, metric, dim))

    if (ndims(x) == 1)
        x = Matrix(x')
    end
    if (ndims(y) == 1)
        y = Matrix(y')
    end

    ##
    ##      Get the number of dimenions.
    d_x = ndims(x) - 1
    d_y = ndims(y) - 1

    ##
    ##      ** It is need that d_x + d_y = length(structure) !!
    if ((length(structure) != d_x + d_y))
        throw(ArgumentError("The structure and the given data are not compatible."))
    end 

    ##
    ##      Compute the size of all recurrence space.
    total_microstates = 1
    space_size::Vector{Int} = []

    ##  From [x]
    for d in 1:d_x
        len = size(x, d + 1) - (sampling_mode == :triangleup ? 0 : (structure[d] - 1))
        total_microstates *= len
        push!(space_size, len)
    end
    ##  From [y]
    for d in 1:d_y
        len = size(y, d + 1) - (sampling_mode == :triangleup ? 0 : (structure[d_x + d] - 1))
        total_microstates *= len
        push!(space_size, len)
    end

    ##
    ##      Prepare the number of samples that we will be using.
    if (num_samples isa Float64 || num_samples == 1)
        if (num_samples <= 0 || num_samples > 1)
            throw(ArgumentError("num_samples must be in the range (0, 1]."))
        end
        if (sampling_mode == :columnwise)
            num_samples = Int(round(num_samples * space_size[2] * structure[1]))
        else
            num_samples = Int(round(num_samples * total_microstates))
        end
    else
        if (num_samples <= 0 || num_samples > total_microstates)
            throw(ArgumentError(string("num_samples must be in the range (1, ", total_microstates,"] for the given data.")))
        end
    end

    ##
    ##      Compute the hypervolume (it is not applied to triangle shape)
    ##      It represents the number of elements inside a motif.
    hv = (shape == :line || shape == :diagonal) ? structure[1] : reduce(*, structure)
    ##      Verify if need to use dictionary or not.
    use_dict = run_mode == :dict
    use_dict = hv > 28 ? true : use_dict
    use_dict = run_mode == :vect ? false : use_dict

    ##
    ##      Verify if hypervolume is compatible with Julia.
    if (hv >= 64 && shape != :pair)
        throw(ArgumentError("Due to memory limitations imposed by Julia, the hyper-volume of a microstate cannot exceed 64."))
    end
    
    ##
    ##      Call the method to compute the histogram and return it.
    ##  i. Verify about the shape.
    ## =====================================================================================================================
    ##          * Shape: square
    if (shape == :square)
        ##
        ##  ii. Sampling mode
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: full
        if (sampling_mode == :full)
            histogram = use_dict ? (
                threads ? (     #   --- Run Mode: dictionary
                    dict_square_full_async(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, metric)) : (
                    dict_square_full(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, total_microstates, metric))
                ) : (
                threads ? (     #   --- Run Mode: vector
                    vect_square_full_async(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, metric)) : (
                    vect_square_full(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, total_microstates, metric)))

            ##      iii. Return the distribution
            return histogram isa Dict{Int, Int} ? (
                total = sum(values(histogram));
                dist = Dict(k => v / total for (k, v) in histogram);
                return dist
            ) : histogram ./ sum(histogram)
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: random
        elseif (sampling_mode == :random)
            histogram = use_dict ? (
                threads ? (     #   --- Run Mode: dictionary
                    dict_square_random_async(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, num_samples, metric)) : (
                    dict_square_random(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, num_samples, metric))
                ) : (
                threads ? (     #   --- Run Mode: vector
                    vect_square_random_async(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, num_samples, metric)) : (
                    vect_square_random(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, num_samples, metric)))

            ##      iii. Return the distribution
            return histogram isa Dict{Int, Int} ? (
                total = sum(values(histogram));
                dist = Dict(k => v / total for (k, v) in histogram);
                return dist
            ) : histogram ./ sum(histogram)
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: columnwise
        elseif (sampling_mode == :columnwise)
            histogram = threads ? (     #   --- Run Mode: vector
                    vect_square_columnwise_async(x, y, parameters, structure, space_size, num_samples, func, [d_x, d_y], hv, metric)) : (
                    vect_square_columnwise(x, y, parameters, structure, space_size, num_samples, func, [d_x, d_y], hv, metric))

            ##      iii. Return the distribution
            return histogram ./ num_samples
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: triangleup
        elseif (sampling_mode == :triangleup)
            histogram = use_dict ? (
                threads ? (     #   --- Run Mode: dictionary
                    dict_square_triangleup_async(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, num_samples, metric)) : (
                    dict_square_triangleup(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, num_samples, metric))
                ) : (
                threads ? (     #   --- Run Mode: vector
                    vect_square_triangleup_async(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, num_samples, metric)) : (
                    vect_square_triangleup(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, num_samples, metric)))

            ##      iii. Return the distribution
            return histogram isa Dict{Int, Int} ? (
                total = sum(values(histogram));
                dist = Dict(k => v / total for (k, v) in histogram);
                return dist
            ) : histogram ./ sum(histogram)
        end
        ## =================================================================================================================
        ##      * Shape: triangle
    elseif (shape == :triangle)
        ##
        ##  ii. Sampling mode
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: full
        if (sampling_mode == :full)
            histogram = use_dict ? (
                threads ? (     #   --- Run Mode: dictionary
                    dict_triangle_full_async(x, y, parameters, structure[1], space_size, func, [d_x, d_y], metric)) : (
                    dict_triangle_full(x, y, parameters, structure[1], space_size, func, [d_x, d_y], total_microstates, metric))
                ) : (
                threads ? (     #   --- Run Mode: vector
                    vect_triangle_full_async(x, y, parameters, structure[1], space_size, func, [d_x, d_y], metric)) : (
                    vect_triangle_full(x, y, parameters, structure[1], space_size, func, [d_x, d_y], total_microstates, metric)))

            ##      iii. Return the distribution
            return histogram isa Dict{Int, Int} ? (
                total = sum(values(histogram));
                dist = Dict(k => v / total for (k, v) in histogram);
                return dist
            ) : histogram ./ sum(histogram)
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: random
        elseif (sampling_mode == :random)
            histogram = use_dict ? (
                threads ? (     #   --- Run Mode: dictionary
                    dict_triangle_random_async(x, y, parameters, structure[1], space_size, func, [d_x, d_y], num_samples, metric)) : (
                    dict_triangle_random(x, y, parameters, structure[1], space_size, func, [d_x, d_y], num_samples, metric))
                ) : (
                threads ? (     #   --- Run Mode: vector
                    vect_triangle_random_async(x, y, parameters, structure[1], space_size, func, [d_x, d_y], num_samples, metric)) : (
                    vect_triangle_random(x, y, parameters, structure[1], space_size, func, [d_x, d_y], num_samples, metric)))

            ##      iii. Return the distribution
            return histogram isa Dict{Int, Int} ? (
                total = sum(values(histogram));
                dist = Dict(k => v / total for (k, v) in histogram);
                return dist
            ) : histogram ./ sum(histogram)
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: columnwise
        elseif (sampling_mode == :columnwise)
            throw(ArgumentError("The sampling mode ':columnwise' is not implemented to shape ':triangle'.'"))
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: triangleup
        elseif (sampling_mode == :triangleup)
            throw(ArgumentError("The sampling mode ':triangleup' is not implemented to shape ':triangle'.'"))
        end
        ## =================================================================================================================
        ##      * Shape: pair
    elseif (shape == :pair)
        ##
        ##  ii. Sampling mode
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: full
        if (sampling_mode == :full)
            throw(ArgumentError("The sampling mode ':full' is not implemented to shape ':pair'.'"))
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: random
        elseif (sampling_mode == :random)
            histogram = threads ? (     #   --- Run Mode: vector
                vect_pair_random_async(x, y, parameters, structure, space_size, func, [d_x, d_y], num_samples, metric)) : (
                vect_pair_random(x, y, parameters, structure, space_size, func, [d_x, d_y], num_samples, metric))

            ##      iii. Return the distribution
            return histogram ./ sum(histogram)
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: columnwise
        elseif (sampling_mode == :columnwise)
            histogram = threads ? (     #   --- Run Mode: vector
                    vect_pair_columnwise_async(x, y, parameters, structure, space_size, func, [d_x, d_y], num_samples, metric)) : (
                    vect_pair_columnwise(x, y, parameters, structure, space_size, func, [d_x, d_y], num_samples, metric))

            ##      iii. Return the distribution
            return histogram ./ num_samples
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: triangleup
        elseif (sampling_mode == :triangleup)
            throw(ArgumentError("The sampling mode ':triangleup' is not implemented to shape ':pair'.'"))
        end
        ## =================================================================================================================
        ##      * Shape: diagonal
    elseif (shape == :diagonal)
        ##
        ##  ii. Sampling mode
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: full
        if (sampling_mode == :full)
            throw(ArgumentError("The sampling mode ':full' is not implemented to shape ':diagonal'.'"))
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: random
        elseif (sampling_mode == :random)
            histogram = threads ? (     #   --- Run Mode: vector
                vect_diagonal_random_async(x, y, parameters, space_size, func, [d_x, d_y], hv, num_samples, metric)) : (
                vect_diagonal_random(x, y, parameters, space_size, func, [d_x, d_y], hv, num_samples, metric))

            ##      iii. Return the distribution
            return histogram ./ sum(histogram)
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: columnwise
        elseif (sampling_mode == :columnwise)
            throw(ArgumentError("The sampling mode ':columnwise' is not implemented to shape ':diagonal'.'"))
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: triangleup
        elseif (sampling_mode == :triangleup)
            throw(ArgumentError("The sampling mode ':triangleup' is not implemented to shape ':diagonal'.'"))
        end
        ## =================================================================================================================
        ##      * Shape: line
    elseif (shape == :line)
        ##
        ##  ii. Sampling mode
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: full
        if (sampling_mode == :full)
            throw(ArgumentError("The sampling mode ':full' is not implemented to shape ':line'.'"))
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: random
        elseif (sampling_mode == :random)
            histogram = threads ? (     #   --- Run Mode: vector
                vect_line_random_async(x, y, parameters, space_size, func, [d_x, d_y], hv, num_samples, metric)) : (
                vect_line_random(x, y, parameters, space_size, func, [d_x, d_y], hv, num_samples, metric))

            ##      iii. Return the distribution
            return histogram ./ sum(histogram)
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: columnwise
        elseif (sampling_mode == :columnwise)
            throw(ArgumentError("The sampling mode ':columnwise' is not implemented to shape ':line'.'"))
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: triangleup
        elseif (sampling_mode == :triangleup)
            throw(ArgumentError("The sampling mode ':triangleup' is not implemented to shape ':line'.'"))
        end
        ## =================================================================================================================
    else
        throw(ArgumentError(string("The shape mode '", string(shape), "' is not valid. Please, use ':square', ':triangle', ':pair', ':diagonal', or ':line'.")))
    end
end
