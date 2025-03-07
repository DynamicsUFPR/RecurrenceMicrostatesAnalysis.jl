#   TODO - Documentation
@inline function recurrence(x::AbstractArray, y::AbstractArray, threshold::Float64, idx::AbstractVector{Int}, dim::AbstractVector{Int}, metric::Metric)
    return @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2]))) <= threshold
end

#   TODO - Documentation
@inline function recurrence(x::AbstractArray, y::AbstractArray, threshold::Tuple{Float64, Float64}, idx::AbstractVector{Int}, dim::AbstractVector{Int}, metric::Metric)
    distance = @inbounds evaluate(metric, view(x, :, view(idx, 1:dim[1])), view(y, :, view(idx, dim[1]+1:dim[1] + dim[2])))
    return ( distance >= threshold[1] && distance <= threshold[2])
end