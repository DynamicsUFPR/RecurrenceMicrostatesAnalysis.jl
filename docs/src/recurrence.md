#   Recurrence functions
Given a trajectory $\mathbf{x}_t \in \mathbb{R}^m$, with $t \in {1, 2, ..., K}$ and $K$ being the length of the analyzed time series, we say that two points $i$ and $j$ from that trajectory are or not recurrent based on a recurrence function. Originally, the recurrence function used for RPs are:

```math
\mathbf{R}_{i,j} = \Theta(\varepsilon - \|\mathbf{x}_i - \mathbf{x}_j\|),
```

where $\Theta$ is the Heaviside function, $\|\cdot\|$ denotes an appropriate norm, and $\varepsilon$ is 
a threshold parameter that defines the maximum distance between two points for them to be considered $\varepsilon$-recurrent
to each other. However, there are several variations that can be used as recurrence functions, as compiled in [Marwan2007](@cite), and `RecurrenceMicrostatesAnalysis.jl` aims to provide support for their implementation.

The library include two recurrence functions: the standard, presented above, and a recurrence function with corridor threshold, proposed in [Iwanski1998](@cite). The standard version code is

```julia
@inline function recurrence(x::Matrix{Float64}, y::Matrix{Float64}, threshold::Float64, idx::AbstractVector{Int}, metric, _)
    return @inbounds evaluate(metric, view(x, :, idx[1]), view(y, :, idx[2])) <= threshold
end
```

The metric is computed using the `Distances.jl` library. Here, you can see that the function receives six arguments: the datasets `x` and `y`, the threshold, the iterator `idx` (which is computed by the shape and sampling functions), a metric from `Distances.jl`, and an ignored parameter that we will discuss later.

The structure of a recurrence function is consistent, so the recurrence function with a corridor threshold — given by the equation below

```math
\mathbf{R}_{ij}=\Theta(\|\mathbf{x}_i-\mathbf{x}_j\|-\varepsilon_{min})\cdot\Theta(\varepsilon_{max} -\|\mathbf{x}_i-\mathbf{x}_j\|), 
```

was implemented as follows in Julia:
```julia
@inline function recurrence(x::Matrix{Float64}, y::Matrix{Float64}, thresholds::Tuple{Float64, Float64}, idx::AbstractVector{Int}, metric, _)
    distance = @inbounds evaluate(metric, view(x, :, idx[1]), view(y, :, idx[2]))
    return (distance > thresholds[1] && distance <= thresholds[2])
end
```

Note that the only difference between them is the `thresholds` parameter, which is a `Tuple{Float64, Float64}` instead of a single `Float64`. This is because it is a free parameter that can be adapted to your specific use case. So, if you implement a different recurrence function, this approach allows you to easily pass constant parameters.

It is also possible to adapt this for the spatial generalization of recurrence plots, as proposed in [Marwan2006](@cite). For this purpose, we use the previously ignored parameter, which represents the number of dimensions of `x` and `y` — more precisely, `[dims(x) - 1, dims(y) - 1]`. This allows us to access the iterator correctly and to unpack each component vector from a spatial problem using a `view`:

```julia
@inline function recurrence(x::AbstractArray, y::AbstractArray, threshold::Float64, idx::AbstractVector{Int}, metric, dim::AbstractVector{Int})
    return @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2]))) <= threshold
end

@inline function recurrence(x::AbstractArray, y::AbstractArray, thresholds::Tuple{Float64, Float64}, idx::AbstractVector{Int}, metric, dim::AbstractVector{Int})
    distance = @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2])))
    return (distance > thresholds[1] && distance <= thresholds[2])
end
```