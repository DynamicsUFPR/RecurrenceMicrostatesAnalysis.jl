"""
    microstate(data_x::AbstractArray, data_y::AbstractArray, parameters, [structure]; kwargs...)

Compute the distribution of recurrence microstates probabilities from the datasets `data_x` and `data_y`.
The input `parameters` consists of the constant values used to calculate the recurrence between two points, 
with a default value of Float64 for default. Meanwhile, the input `structure` is a vector where each element 
represents a side of the motif.

This function can return a vector, an array or a dictionary based on the number of possible microstates and
the setting of `run_mode` or `sampling_mode`.

## Keywords arguments
- `shape::Symbol`: can be `:square` or `:triangle`. The value `:square` refers to the default square format of motifs,
    based on the work of [Corso2018](@cite) with a generalization for spatial data based on the work of [Marwan2006](@cite).
    On  the other hand, the value `:triangle` is only available for 2-dimensional recurrence spaces (i.e., time series) 
    and is based on the work of [Hirata2021](@cite).

- `run_mode::Symbol`: can be `:default`, `:dict`, or `:vect`. `:dict` and `:vect` set the return format to vector and
    dictionary respectively. The mode `:default` uses a vector for motifs with an hypervolume lesser than 28 and 
    a dictionary otherwise.

- `sampling_mode::Symbol`: can be `:full`, `:random`, or `:columnwise`. The sampling mode `:full` retrieves all 
available microstates in the recurrence space sequentially; it is not yet available for multi-threading. The `:random`
mode retrieves a number of samples based on the value of `num_samples` and creates a distribution for the entire 
recurrence space. Finally, `:columnwise` creates a distribution for each column of the recurrence space and is only 
available for `:vect` run mode.

- `num_samples::Union{Int, Float64}`: the number of samples used to create the probability distribution. CIt can be a 
percentage of the total possible motifs or the desired number of samples.

- `metric::Metric`: a structure used to compute distances from the [Distances.jl](https://github.com/JuliaStats/Distances.jl) 
library.

- `func = (x, y, p, ix, dim, metric) -> recurrence(x, y, p, ix, dim, met)`: this is the recurrence function that returns true when 
there is a recurrence and false otherwise. This function takes 5 parameters:
    - x::AbstractArray : it is a reference to the parameter `data_x`.
    - y::AbstractArray : it is a reference to the parameter `data_y`.
    - p : it is a reference to the parameter `parameters`.
    - [ix] : it is a vector with `length(structure)` elements that indicates the elements used to compute the 
recurrence.
    - [dim] : it is a vector with two elements: the number of dimenions of `data_x` and `data_y`. It is used only for
high-dimensional problems.
    - metric::Metric : a metric from the [Distances.jl](https://github.com/JuliaStats/Distances.jl) library. It used
    the same value as the `metric` parameter of distribution function.

For reference, you can see the code [`src/rma/recurrence.jl`](src/rma/recurrence.jl).

- `use_threads::Bool`: defines whether the library will use threads in the computational process. If set to true, the 
library will use the number of available threads defined by `JULIA_NUM_THREADS`.
"""
function distribution(data_x::AbstractArray, data_y::AbstractArray, parameters, structure::AbstractVector{Int};
        shape::Symbol = :square, run_mode::Symbol = :default, sampling_mode::Symbol = :random, 
        num_samples::Union{Int, Float64} = 1.0, use_threads::Bool = Threads.nthreads() > 1, 
        metric::Metric = euclidean_metric, func = (x, y, p, ix, dim, metric) -> recurrence(x, y, p, ix, dim, metric))

    #       Verify the arguments
    #   1. Structure must have at least two elements.
    if (length(structure) < 2)
        throw(ArgumentError("The microstate structure required at least two values."))
    end
    #   2. Structure must have same number of dimenions that the recurrence space.
    d_x = ndims(data_x) - 1
    d_y = ndims(data_y) - 1
    if (length(structure) != d_x + d_y)
        throw(ArgumentError("The structure and the given data are not compatible."))
    end
    #
    #       Number of samples
    #   1. Take the recurrence space size that we can use.
    total_microstates = 1
    space_size::Vector{Int} = []
    for d in 1:d_x
        len = size(data_x, d + 1) - structure[d]

        total_microstates *= len
        push!(space_size, len)
    end
    for d in 1:d_y
        len = size(data_y, d + 1) - structure[d_x + d]

        total_microstates *= len
        push!(space_size, len)
    end
    #   2. Prepare the number of samples that we will be using.
    if (num_samples isa Float64 || num_samples == 1)
        if (num_samples <= 0 || num_samples > 1)
            throw(ArgumentError("num_samples must be in the range (0, 1]."))
        end
        num_samples = Int(round(num_samples * total_microstates))
    else
        if (num_samples <= 0 || num_samples > total_microstates)
            throw(ArgumentError(string("num_samples must be in the range (1, ", total_microstates,"] for the given data.")))
        end
    end
    #
    #       Compute the hypervolume, which represents the number of elements inside a motif.
    #   TODO - Adapt it to triangle (just for 2D RP)
    hv = reduce(*, structure)
    #       Verify if need to use dictionary or not.
    use_dict = run_mode == :dict
    use_dict = hv > 28 ? true : use_dict
    use_dict = run_mode == :vect ? false : use_dict
    #
    #       Verify again the settings...
    if (!use_dict && hv >= 64)
        throw(ArgumentError("Due to memory limitations imposed by Julia, the hyper-volume of a microstate cannot exceed 64."))
    end

    #   TODO - Move to line 60 and organize
    if (shape == :triangle && length(structure) > 2)
        throw(ArgumentError("The shape mode `:triangle` is only available for a recurrence space with two dimensions."))
    end
    if (sampling_mode == :columnwise && length(structure) > 2)
        throw(ArgumentError("The sampling mode `:columnwise` is only available for a recurrence space with two dimensions."))
    end
    if (sampling_mode != :columnwise && sampling_mode != :full && sampling_mode != :random)
        throw(ArgumentError("Invalid sampling mode. Use :full, :random or :columnwise"))
    end

    #       Call the process...
    if (shape == :square)
        if (sampling_mode == :full)
            histogram = use_dict ? (
                throw("Dictionaries not implemented yet")) : (
                use_threads ? (
                    throw("Threads is not yet implemented to :full sampling mode.")) : (
                    square_full(data_x, data_y, parameters, structure, space_size, func, [d_x, d_y], hv, total_microstates, metric)))
            #
            #   Compute the distribution from the histogram.
            return histogram ./ sum(histogram)
        elseif (sampling_mode == :random)
            histogram = use_dict ? (
                throw("Dictionaries not implemented yet")) : (
                use_threads ? (
                    square_random_async(data_x, data_y, parameters, structure, space_size, num_samples, func, [d_x, d_y], hv, metric)) : (
                    square_random(data_x, data_y, parameters, structure, space_size, num_samples, func, [d_x, d_y], hv, metric)))
            #
            #   Compute the distribution from the histogram.
            return histogram ./ sum(histogram)
        elseif (sampling_mode == :columnwise)
            throw("Sampling mode not implemented yet.")
        end
    elseif (shape == :triangle)
        throw("Shape mode not implemented yet.")
    else
        throw(ArgumentError("Invalid shape. Use :square or :triangle"))
    end
end