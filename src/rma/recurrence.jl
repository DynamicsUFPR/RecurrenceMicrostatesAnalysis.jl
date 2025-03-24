#
#           RMA Core - Recurrence functions
#
#       Here we define some recurrence functions of the bibligraphy. It is a
#   low code to compute the recurrence between two objects and define the
#   recurrence of a specific position of a recurrence space.
#
#       This code have two different ways to compute a same recurrence value: one for
#   the standard recurrence and other to the spatial recurrence.
#
# ----------------------------------    2D RECURRENCE SPACE      ----------------------------------------------------
"""
    recurrence([x], [y], parameters, idx::AbstractVector{Int}, metric::Metric, dim::AbstractVector{Int})

Compute the recurrence between two position defined by `idx` of the datasets `x` and `y`. `parameters` defines
which type of recurrence it will use, like the standard recurrence and the recurrence with corridor threshold.
`metric` defines the norm applied to the datasets, it uses the library Distances.jl to compute the metric. `dim` is
just used for high-dimensional problems, such as images.

    #       Examples
    #   For a 2D recurrence space.
    @inline function recurrence(...)
        return @inbounds evaluate(metric, x[idx[1]], y[idx[2]]) <= threshold
    end

    #   For a recurrence tensor space (generalization to spatial data)
    @inline function recurrence(...)
        return @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2]))) <= threshold
    end

These functions return `true` when we have a recurrence, and `false` otherwise.
"""
@inline function recurrence(x::Matrix{Float64}, y::Matrix{Float64}, threshold::Float64, idx::AbstractVector{Int}, metric::Metric, _::AbstractVector{Int})
    return @inbounds evaluate(metric, x[idx[1]], y[idx[2]]) <= threshold
end

@inline function recurrence(x::Matrix{Float64}, y::Matrix{Float64}, thresholds::Tuple{Float64, Float64}, idx::AbstractVector{Int}, metric::Metric, _::AbstractVector{Int})
    distance = @inbounds evaluate(metric, x[idx[1]], y[idx[2]])
    return (distance > thresholds[1] && distance <= thresholds[2])
end

# ----------------------------------      RECURRENCE TENSOR SPACE      ----------------------------------------------
@inline function recurrence(x::AbstractArray, y::AbstractArray, threshold::Float64, idx::AbstractVector{Int}, metric::Metric, dim::AbstractVector{Int})
    return @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2]))) <= threshold
end

@inline function recurrence(x::AbstractArray, y::AbstractArray, thresholds::Float64, idx::AbstractVector{Int}, metric::Metric, dim::AbstractVector{Int})
    distance = @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2])))
    return (distance > thresholds[1] && distance <= thresholds[2])
end

#######################
#   jensen-shannon
#       sqrt dela.
