#   Motifs: shapes and sampling
The library supports different motif shapes and sampling modes. Here, we present a brief explanation of how these mechanisms work and how you can create your own shape or sampling mode.

##  Shapes
By default, the library includes five predefined motif shapes: `:square`, `:triangle`, `:pair`, `:line`, and `:diagonal`. These shapes are defined in the file `src/rma/index.jl`, and they determine how a motif is drawn by the library and how it is converted to a decimal value used as an index. Therefore, when constructing a motif shape, it is important to consider how it will be converted into a decimal value.

For example, a square motif can be converted into a decimal value using the following equation (without spatial generalization).

```math
I = \sum_{r = 0}^{n - 1}\sum_{c = 0}^{n - 1} 2^{rn+c}~\mathbf{R}_{i+r, j+c},
```

where $2^{rn+c}$ is responsible for associating each position of the motif with a power of 2, converting the binary structure into a decimal value. A Julia's function to compute it can be written as:

```julia
function compute_index_square(x::AbstractArray, y::AbstractArray, parameters, structure::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, fixed::Vector{Int}, itr::Vector{Int}, metric) where {F}

    ##      Let a variable to store the index.
    I = 0

    ##      Copy the values of the fixed indeces to the vector of iterative indeces.
    copy!(itr, fixed)       ##  We do it to avoid memory allocations =D

    ##      Iterate to compute the index.
    for r in 0:(structure[1] - 1)
        for c in 0:(structure[2] - 1)
            ##      Change the iterator.
            itr[1] = fixed[1] + r
            itr[2] = fixed[2] + c

            ##      Calculate the recurrence between two positions.
            if @inline func(x, y, parameters, itr, metric, dim)
                index += 2^((r * structure[1]) + c)
            end
        end
    end

    return I + 1 ##     It is necessary for Julia indexing! i = I + 1
end
```

Knowing an algebraic expression to convert a motif into a decimal value is not strictly necessary, but it is recommended — especially considering the importance of understanding how this process will work for any value of $n$ (if applicable). For example, consider a motif with an X-shape (for $n = 3$):

```math
\begin{pmatrix}
\xi_1   &       &  \xi_2    \\
        & \xi_3 &           \\
\xi_4   &       &  \xi_5
\end{pmatrix}
```

It is easy to convert this shape into an index using something like:
```julia
function compute_index_x(x::AbstractArray, y::AbstractArray, parameters, func::F, dim::AbstractVector{Int}, fixed::Vector{Int}, itr::Vector{Int}, metric) where {F}

    ##      Let a variable to store the index.
    I = 0

    ##      Copy the values of the fixed indeces to the vector of iterative indeces.
    copy!(itr, fixed)       ##  We do it to avoid memory allocations =D

    ##  1. \xi_1
    if @inline func(x, y, parameters, itr, metric, dim)
        I += 1
    end

    ##  2. \xi_2
    itr[2] = fixed[2] + 2
    if @inline func(x, y, parameters, itr, metric, dim)
        I += 2
    end

    ##  3. \xi_3
    itr[1] = fixed[1] + 1
    itr[2] = fixed[2] + 1
    if @inline func(x, y, parameters, itr, metric, dim)
        I += 4
    end

    ##  4. \xi_4
    itr[1] = fixed[1] + 2
    itr[2] = fixed[2]
    if @inline func(x, y, parameters, itr, metric, dim)
        I += 8
    end

    ##  5. \xi_5
    itr[2] = fixed[2] + 2
    if @inline func(x, y, parameters, itr, metric, dim)
        I += 16
    end

    return I + 1 ##     It is necessary for Julia indexing! i = I + 1
end
```

But it is also clear that this code is not scalable, as it relies on a fixed structure. Transforming this into a scalable version that can handle any motif structure is one of the key computational challenges in the process of defining a new motif shape.

The application of a shape is associated with a sampling mode, so it is necessary to create a sampling function for each defined shape. This is because the sampling mode determines how the `fixed` parameter is set.

## Sampling
The sampling defines how each motif will be extracted from an RP, without computing the full RP explicitly. The library includes five predefined sampling modes: `:full`, `:random`, `:triangleup`, `:columnwise` and `:columnwise_full`. These modes are defined in the folder `src/rma/histograms/`, being called by the function `distribution(...)`.

The difficulty of creating a sampling function is proportional to the complexity of your problem. For example, if you are sampling motifs randomly from the RP, the sampling function is simple. However, if the motifs need to be extracted from a specific structure, the difficulty increases.

The sampling code for a `:random` sampling using square motifs is

```julia
function vect_square_random(x::AbstractArray, y::AbstractArray, parameters, structure::AbstractVector{Int},
    space_size::AbstractVector{Int}, func::F, dim::AbstractVector{Int}, hv::Int, samples::Int, metric) where {F}

    ##
    ##      Alloc memory for the histogram and the index list.
    hg = zeros(Int, 2^hv)   ##  Our histogram!
    fixed = ones(Int, length(space_size))
    itr = zeros(Int, length(space_size))

    ##
    ##      Compute the power vector.
    p_vect = zeros(Int, hv)
    for i in 1:hv
        p_vect[i] = 2^(i - 1)
    end

    ##
    ##      Get the samples and compute the histogram.
    @inbounds for _ in 1:samples
        ##
        ##      Take a random index.
        for s in eachindex(space_size)
            fixed[s] = rand(1:space_size[s])
        end

        ##
        ##      Compute the index and register the motif.
        p = @fastmath compute_index_square(x, y, parameters, structure, func, dim, fixed, itr, p_vect, metric)
        hg[p] += 1
    end

    ##
    ##      Return the histogram.
    return hg
end
```

In this code, we compute the power vector (`p_vect`) before calling the `compute_index_square` function, so it is not necessary to calculate powers of 2 inside the loop. Moreover, we allocate memory for the iterator `itr` and the fixed indices `fixed` outside the loop to avoid unnecessary calls to the garbage collector. The `for`...

```julia
for s in eachindex(space_size)
    fixed[s] = rand(1:space_size[s])
end
```

... retrieves an index set $(i,j)$ to define the first recurrence $R_{ij}$, , from which a motif is constructed using the shape function `compute_index_square`.