
"""
    recurrence([x], [y], parameters, idx::AbstractVector{Int}, metric::Metric, dim::AbstractVector{Int})

Compute the recurrence between two position defined by `idx` of the datasets `x` and `y`. `parameters` defines
which type of recurrence it will use, like the standard recurrence and the recurrence with corridor threshold.
`metric` defines the norm applied to the datasets, it uses the library `Distances.jl` to compute the metric. `dim` is
just used for high-dimensional problems, such as images. A recurrence function is a function of the form

```math
R_{ij} = \\Theta(\\varepsilon - \\|\\mathbf{x}_i - \\mathbf{y}_j\\|),
```

where \$\\Theta\$ is the Heaviside function, \$\\|\\cdot\\|\$ denotes an appropriate norm, and \$\\varepsilon\$ is 
a threshold parameter that defines the maximum distance between two points for them to be considered \$\\varepsilon\$-recurrent
to each other.

    #       Examples
    #   For a 2D recurrence space.
    @inline function recurrence(...)
        return @inbounds evaluate(metric, x[idx[1]], y[idx[2]]) <= threshold
    end

    #   For a recurrence tensor space (generalization to spatial data)
    @inline function recurrence(...)
        return @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2]))) <= threshold
    end

Input:
* `[x]`: a dataset.
* `[y]`: a dataset.
* `[parameter]`: set of parameters used to compute the recurrence microstate distribution, i.e., the value of \$\\varepsilon\$.
* `[idx]`: vector of indeces from `[x]` and `[y]` to calculate the recurrence between them.
* `metric`: metric from `Distances.jl` used to compute the recurrence. It defines the norm \$\\|\\cdot\\|\$.
* `[dim]`: number of dimensions derived from `[x]` and `[y]`. If you are using a time series it is usually `[1,1]`.

Output: Recurrence functions return `true` when we have a recurrence, and `false` otherwise.
"""
@inline function recurrence(x::Matrix{Float64}, y::Matrix{Float64}, threshold::Float64, idx::AbstractVector{Int}, metric, _::AbstractVector{Int})
    return @inbounds evaluate(metric, view(x, :, idx[1]), view(y, :, idx[2])) <= threshold
end

@inline function recurrence(x::Matrix{Float64}, y::Matrix{Float64}, thresholds::Tuple{Float64, Float64}, idx::AbstractVector{Int}, metric, _::AbstractVector{Int})
    distance = @inbounds evaluate(metric, view(x, :, idx[1]), view(y, :, idx[2]))
    return (distance > thresholds[1] && distance <= thresholds[2])
end
    
# ----------------------------------      RECURRENCE TENSOR SPACE      ----------------------------------------------
@inline function recurrence(x::AbstractArray, y::AbstractArray, threshold::Float64, idx::AbstractVector{Int}, metric, dim::AbstractVector{Int})
    return @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2]))) <= threshold
end

@inline function recurrence(x::AbstractArray, y::AbstractArray, thresholds::Tuple{Float64, Float64}, idx::AbstractVector{Int}, metric, dim::AbstractVector{Int})
    distance = @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2])))
    return (distance > thresholds[1] && distance <= thresholds[2])
end

# ----------------------------------      JRP      ----------------------------------------------
@inline function jrp(x::Matrix{Float64}, y::Matrix{Float64}, threshold::Float64, idx::AbstractVector{Int}, metric, _::AbstractVector{Int})
    return @inbounds (evaluate(metric, view(x, :, idx[1]), view(x, :, idx[2])) <= threshold) && (evaluate(metric, view(y, :, idx[1]), view(y, :, idx[2])) <= threshold)
end


#######################
#   jensen-shannon
#       sqrt dela.
