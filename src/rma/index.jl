#
#           RMA Core - Low level functions to convert a motif to an index.
#
#       We have a different function to each motif shape.
#
"""
These functions compute and convert a motif to a decimal value that will be used as an index
for store the motif in the system memory.

It have three types:
* Type-A: Using a generic structure, applied to square (or rectangle) motifs and time.
* Type-B: For a 2D space, such as the triangle motif, proposed by [Hirata2021](@cite). It must not be applied to
spatial data using the recurrence plot generalization, proposed by [Marwan2006](@cite).
* Type-C: Using fixed structures, like diagonal, and line shapes.

Each type have a specific method to get the fixed indeces.
"""
# -----------------------------------------------------------  TYPE A:  GENERIC STRUCTURE ------------------------
function compute_index_square(x::AbstractArray, y::AbstractArray, parameters, structure::AbstractVector{Int}, func::F, 
    dim::AbstractVector{Int}, fixed::Vector{Int}, itr::Vector{Int}, power_vector::Vector{Int}, metric) where {F}

    ##      Let a variable to store the index.
    index = 0

    ##      Copy the values of the fixed indeces to the vector of iterative indeces.
    copy!(itr, fixed)

    ##      Iterate to compute the index.
    for m in 1:length(power_vector)
        ##      Calculate the recurrence between two positions.
        if @inline func(x, y, parameters, itr, metric, dim)
            index += power_vector[m]
        end

        ##      Move the iterative indeces to next position.
        itr[1] += 1
        for k in length(structure) - 1
            if (itr[k] >= fixed[k] + structure[k])
                itr[k] = fixed[k]
                itr[k + 1] += 1
            else
                break
            end
        end
    end
    
    ##      Return the computed index adapted to Julia's indexing.
    return index + 1
end

function compute_index_pair(x::AbstractArray, y::AbstractArray, parameters, structure::AbstractVector{Int}, func::F,
    dim::AbstractVector{Int}, fixed::Vector{Int}, itr::Vector{Int}, metric) where {F}
    
    ##  Let a variable to store the index.
    index = 0

    ##  Copy the values of the fixed indeces to the vector of iterative indeces.
    copy!(itr, fixed)
    ##  Add the structure values to itr index.
    itr .+= structure

    ##  Compute the index.
    if @inline func(x, y, parameters, fixed, metric, dim)
        index += 1
    end

    if @inline func(x, y, parameters, itr, metric, dim)
        index += 2
    end

    ##      Return the computed index adapted to Julia's indexing.
    return index + 1
end

# -----------------------------------------------------------  TYPE B:   FOR 2D SPACE     -----------------------
function compute_index_triangle(x::Matrix{Float64}, y::Matrix{Float64}, parameters, len::Int, func::F, dim::AbstractVector{Int},
    fixed::Vector{Int}, itr::Vector{Int}, power_vector::Vector{Int}, metric) where {F}

    ##      Let a variable to store the result, and another to store the index to access the power vector.
    index = 0
    m = 1

    ##      Copy the values of the fixed indeces to the vector of iterative indeces.
    copy!(itr, fixed)

    ##      Iterate to compute the index.
    for j in 0:len - 1
        itr[2] = fixed[2] + j
        for i in j:len - 1
            itr[1] = fixed[1] + i

            ##      Calculate the recurrence between two positions.
            if @inline func(x, y, parameters, itr, metric, dim)
                index += power_vector[m]
            end

            m += 1
        end
    end
    
    ##      Return the computed index adapted to Julia's indexing.
    return index + 1
end

# -----------------------------------------------------------  TYPE C:   FOR FIXED STRUCTURE   --------------------
function compute_index_diagonal(x::AbstractArray, y::AbstractArray, parameters, func::F, dim::AbstractVector{Int},
    fixed::Vector{Int}, itr::Vector{Int}, power_vector::Vector{Int}, metric) where {F}
    
    ##      Let a variable to store the index.
    index = 0

    ##      Copy the values of the fixed indeces to the vector of iterative indeces.
    copy!(itr, fixed)

    ##      Iterate to compute the index.
    for m in eachindex(power_vector)
        ##      Calculate the recurrence between two positions.
        if @inline func(x, y, parameters, itr, metric, dim)
            index += power_vector[m]
        end

        ##      Move to the next position.
        itr .+= 1
    end

    ##      Return the computed index adapted to Julia's indexing.
    return index + 1
end

function compute_index_line(x::AbstractArray, y::AbstractArray, parameters, func::F, dim::AbstractVector{Int},
    fixed::Vector{Int}, itr::Vector{Int}, power_vector::Vector{Int}, metric) where {F}
    
    ##      Let a variable to store the index.
    index = 0

    ##      Copy the values of the fixed indeces to the vector of iterative indeces.
    copy!(itr, fixed)

    ##      Iterate to compute the index.
    for m in eachindex(power_vector)
        ##      Calculate the recurrence between two positions.
        if @inline func(x, y, parameters, itr, metric, dim)
            index += power_vector[m]
        end

        itr[2] += 1
    end

    ##      Return the computed index adapted to Julia's indexing.
    return index + 1
end