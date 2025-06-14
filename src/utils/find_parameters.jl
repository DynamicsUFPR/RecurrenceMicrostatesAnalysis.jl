#
#       RMA Utils
#
#       TODO - Improve it. The current version is only for tests.
function find_parameters(x::AbstractArray, n::Int; threshold_min::Float64 = 0.0, threshold_max::Float64 = maximum(pairwise(Euclidean(), x, x)), large_prec = 20, small_prec = 50, num_samples = 0.05, shape = :random)
    ##
    ##      Define a range of thresholds and alloc memory to store the entropy.
    threshold_range = range(threshold_min, threshold_max, large_prec)
    s_max = 0.0
    t_indeces = zeros(Int, 3)

    ##
    ##      Compute the max entropy for a "large precision"
    for i in eachindex(threshold_range)
        dist = distribution(x, threshold_range[i], n; num_samples = num_samples, shape = shape)
        s = rentropy(dist)

        if (s > s_max)
            t_indeces[1] = t_indeces[2]
            t_indeces[2] = i
            s_max = s
        elseif (s < s_max)
            t_indeces[3] = i
            break
        end
    end

    if (t_indeces[3] == 0)
        throw(ArgumentError("The max entropy is not between the values of 'threshold_min' and 'threshold_max'."))
    end

    ##
    ##      Increase the precision.
    threshold_range = range(threshold_range[t_indeces[1]], threshold_range[t_indeces[3]], small_prec)
    s_max = 0.0
    t_indeces = zeros(Int, 3)

    for i in eachindex(threshold_range)
        dist = distribution(x, threshold_range[i], n)
        s = rentropy(dist)

        if (s > s_max)
            t_indeces[1] = t_indeces[2]
            t_indeces[2] = i
            s_max = s
        elseif (s <= s_max)
            t_indeces[3] = i
            break
        end
    end

    return threshold_range[t_indeces[2]], s_max
end