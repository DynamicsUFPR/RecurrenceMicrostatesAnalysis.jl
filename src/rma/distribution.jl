
"""
### Based on Recurrence Plot

        distribution([x], [parameters], n::Int; kwargs...)

Compute the distribution of recurrence microstates probabilities from the dataset `x`. The input `parameters`
consists of the constant values used to calculate the recurrence between two points. `n` is an integer that
represents the length of motifs side.

Input:
* `[x]`: input dataset.
* `[parameter]`: set of parameters used to compute the recurrence microstate distribution.
* `n`: microstate size.

Output: this function can return a vector, an array or a dictionary based on the number of possible microstates and
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
    return distribution(x, x, parameters, ones(Int, 2 * dim) .* n; shape = shape, run_mode = run_mode, sampling_mode = sampling_mode, num_samples = num_samples, threads = threads, metric = metric, func = func)
end
"""
### Based on Cross-Recurrence Plot

    distribution([x], [y], parameters, n::Int; kwords...)

Compute the distribution of recurrence microstates probabilities from the datasets `x` and `y`. The input `parameters`
consists of the constant values used to calculate the recurrence between two points. `n` is an integer that
represents the length of motifs side.

Input:
* `[x]`: input dataset.
* `[y]`: input dataset.
* `[parameter]`: set of parameters used to compute the recurrence microstate distribution.
* `n`: microstate size.

Output: this function can return a vector, an array or a dictionary based on the number of possible microstates and
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
    return distribution(x, y, parameters, ones(Int, dim_x + dim_y) .* n; shape = shape, run_mode = run_mode, sampling_mode = sampling_mode, num_samples = num_samples, threads = threads, metric = metric, func = func)
end

"""
### Using DifferencialEquations.jl

    distribution([solution], parameters, n::Int, vicinity::Union{Int, Float64}; kwords...)

Compute the distribution of recurrence microstates probabilities from the `solution` of a differencial equation solved by 
the library DifferencialEquations.jl. The input `parameters` consists of the constant values used to calculate the recurrence 
between two points. `n` is an integer that represents the length of motifs side. `vicinity` is the time separation used to
discretize a continuous problem.

Input:
* `[solution]`: solution returned by the library `DifferentialEquations.jl`.
* `[parameter]`: set of parameters used to compute the recurrence microstate distribution.
* `n`: microstate size.
* `σ`: sampling parameter; it defines the time resolution of discretized data.

Specific **kwargs**:
* `transient`: defines an interval of time that will be ignored, and taked as a transient.
* `K`: defines the maximum size of the result time series.

Output: this function can return a vector, an array or a dictionary based on the number of possible microstates and
the setting of `run_mode` or `sampling_mode`.
"""
function distribution(solution, parameters, n::Int, σ::Union{Int, Float64} = 0.0;
    shape::Symbol = :square, run_mode::Symbol = :default, sampling_mode::Symbol = :random,
    num_samples::Union{Int, Float64} = 0.05, threads::Bool = Threads.nthreads() > 1, 
    metric = euclidean_metric, func = (x, y, p, ix, metric, dim) -> recurrence(x, y, p, ix, metric, dim),
    transient::Int = 0, K::Int = 0)
    
    ##
    ##      Prepare the dataset.
    data = prepare(solution, σ; transient = transient, K = K)

    dim = ndims(data) - 1
    structure = ones(Int, 2 * dim) .* n

    return distribution(data, data, parameters, structure; shape = shape, run_mode = run_mode, sampling_mode = sampling_mode, num_samples = num_samples, threads = threads, metric = metric, func = func)
end

"""
### Main

    distribution([x], [y], parameters, [structure]; kwords...)

Compute the distribution of recurrence microstates probabilities from the datasets `x` and `y`. The input `parameters`
consists of the constant values used to calculate the recurrence between two points. Meanwhile, the input `structure`
is a vector where each element represents a side of the motif.

Input:
* `[x]`: input dataset.
* `[y]`: input dataset.
* `[parameter]`: set of parameters used to compute the recurrence microstate distribution.
* `[structure]`: microstate structure.

**kwargs**:
* `shape`: microstate shape. Can be `:square`, `:triangle`, `:pair`, `:diagonal` or `:line`. (default `:square`)
* `run_mode`: define the output format. It can be `:vect` for a `Vector{Float64}`, or `:dict` for a `Dict{Int, Float64}`. If you are you sampling_mode `:columnwise``
* `sampling_mode`: define how the library will take motifs in a RP. Can be `:full`, `:random`, `:triangleup`, `:columnwise` or `:columnwise_full`. (default `:random`)
* `num_samples`: number of samples used to compute the distribution. Can be an `Int` value or a `Float64`, which will be interpretad as a proportion of the total population of microstates in a RP. (This is not required for `:full` and `:columnwise_full` sampling modes)
* `threads`: set if library will use asyncronous jobs or not.
* `metric`: metric defined using the library `Distances.jl`.
* `func`: recurrence function.

Output: this function can return a vector, an array or a dictionary based on the number of possible microstates and
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
        ##          * Mode: columnwise full
        elseif (sampling_mode == :columnwise_full)
            histogram = threads ? (     #   --- Run Mode: vector
                    vect_square_columnwise_full_async(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, metric)) : (
                    vect_square_columnwise_full(x, y, parameters, structure, space_size, func, [d_x, d_y], hv, metric))

            ##      iii. Return the distribution

            res = Float64.(histogram)
            for i in 1:space_size[1]
                res[:, i] ./= space_size[2] 
            end

            return res
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
        ## -----------------------------------------------------------------------------------------------------------------
        else
            throw(string("The sampling mode '", sampling_mode, "' is invalid."))
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
        ## -----------------------------------------------------------------------------------------------------------------
        else
            throw(string("The sampling mode '", sampling_mode, "' is invalid."))
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
        ## -----------------------------------------------------------------------------------------------------------------
        else
            throw(string("The sampling mode '", sampling_mode, "' is invalid."))
        end
        ## =================================================================================================================
        ##      * Shape: diagonal
    elseif (shape == :diagonal)
        ##
        ##  ii. Sampling mode
        ## -----------------------------------------------------------------------------------------------------------------
        ##          * Mode: full
        if (sampling_mode == :full)
            histogram = vect_diagonal_full_async(x, y, parameters, space_size, func, [d_x, d_y], hv, metric)
            return histogram ./ sum(histogram)
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
        ## -----------------------------------------------------------------------------------------------------------------
        else
            throw(string("The sampling mode '", sampling_mode, "' is invalid."))
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
        ## -----------------------------------------------------------------------------------------------------------------
        else
            throw(string("The sampling mode '", sampling_mode, "' is invalid."))
        end
        ## =================================================================================================================
    else
        throw(ArgumentError(string("The shape mode '", string(shape), "' is not valid. Please, use ':square', ':triangle', ':pair', ':diagonal', or ':line'.")))
    end
end
