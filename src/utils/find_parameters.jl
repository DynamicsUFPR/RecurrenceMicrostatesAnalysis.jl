#
#       RMA Utils
#
#       TODO - Improve it. The current version is only for tests.
function find_parameters(x::AbstractArray, n::Int; threshold_min::Float64 = 0.0, threshold_max::Float64 = 2 * std(x), large_prec = 20, small_prec = 50)
    ##
    ##      Define a range of thresholds and alloc memory to store the entropy.
    threshold_range = range(threshold_min, threshold_max, large_prec)
    s_max = 0.0
    t_indeces = zeros(Int, 3)

    ##
    ##      Compute the max entropy for a "large precision"
    for i in eachindex(threshold_range)
        dist = distribution(x, threshold_range[i], n)
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

function find_parameters(solution, n::Int, transient::Int, len::Int; threshold_min::Float64 = 0.0, threshold_max::Float64 = -1.0, vicinity_min::Float64 = 0.0, vicinity_max::Float64 = 0.5, large_prec = 20, small_prec = 50)
    ##
    ##      Define a range of vicinity and alloc memory to store the entropy.
    vicinity_range = range(vicinity_min, vicinity_max, large_prec)
    s_max = 0.0
    threshold = 0.0
    v_indeces = zeros(Int, 3)

    ##
    ##      Compute the max entropy for a "large precision"
    for i in eachindex(vicinity_range)
        data = prepare(solution, vicinity_range[i]; transient = transient, max_length = len)
        _, s = find_parameters(data, n; threshold_min = threshold_min, threshold_max = threshold_max < 0 ? 2 * std(data) : threshold_max, large_prec = large_prec, small_prec = small_prec)

        if (s > s_max)
            v_indeces[1] = v_indeces[2]
            v_indeces[2] = i
            s_max = s
        elseif (s <= s_max)
            v_indeces[3] = i
        end
    end

    if (v_indeces[3] == 0)
        throw(ArgumentError("The max entropy is not between the values of 'vicinity_min' and 'vicinity_max'."))
    end

    ##
    ##      Increase the precision.
    vicinity_range = range(vicinity_range[v_indeces[1]], vicinity_range[v_indeces[3]], small_prec)
    s_max = 0.0
    v_indeces = zeros(Int, 3)

    for i in eachindex(vicinity_range)
        data = prepare(solution, vicinity_range[i]; transient = transient, max_length = len)
        m_threshold, s = find_parameters(data, n; threshold_min = threshold_min, threshold_max = threshold_max < 0 ? 2 * std(data) : threshold_max, large_prec = large_prec, small_prec = small_prec)

        if (s > s_max)
            v_indeces[1] = v_indeces[2]
            v_indeces[2] = i
            s_max = s
            threshold = m_threshold
        elseif (s <= s_max)
            v_indeces[3] = i
        end
    end

    return threshold, vicinity_range[v_indeces[2]], s_max
end