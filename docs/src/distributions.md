#   Recurrence Motifs Probability Distributions
`RecurrenceMicrostatesAnalysis.jl` aims to be a user-friendly library with powerful capabilities. It can be used through simple function calls or more advanced configurations that offer greater control. We will begin with the simpler usage, explaining its arguments and settings, and gradually move toward more complex configurations throughout this discussion.

##  One-dimensional data
This section presents a run similar to the one shown on the [quick start](quickstart.md) page, but with a more detailed explanation. For one-dimensional problems, such as the logistic map or the generalized Bernoulli shift map (Beta-X), you can use a vector of positions along the trajectory as input. To illustrate this, let's consider a uniform distribution:
```@repl dist_one
using Distributions
data = rand(Uniform(0, 1), 3000)
```

Computing the recurrence motif distribution is straightforward once the `threshold` and `n` (motif size) parameters are defined. A good value for `threshold` can be estimated using the `find_parameters` function, which we recommend using in most cases.
```@repl dist_one
using RecurrenceMicrostatesAnalysis
th, s = find_parameters(data, 3)
dist = distribution(data, th, 3)
```

!!! warning
    We do not recommend the use of `find_parameters` inside a loop, as it needs to compute several distributions to find the `threshold` value that maximizes recurrence entropy, which can significantly reduce the library's performance. For this reason, we have not created an overload of the `distribution` function that automatically calculates the `threshold`. Instead, we suggest using an average `threshold` value computed from a few representative snippets of your dataset using the `find_parameters` function.

The `distribution` function includes several keyword arguments for configuration. Before moving on to the next section, we will discuss these arguments, as they apply to every call of the distribution function.

### Motif constrained shape
There are variations in motif constraint shapes proposed in the literature, such as the triangular motif. Supporting these shape generalizations is one of the goals of `RecurrenceMicrostatesAnalysis.jl`, and it is also a computational challenge. Adapting the conversion of motifs with a generic shape from a binary structure to a decimal value can be a very complex problem, and to support this in the library, we need to adapt the pipeline that converts a motif for each specific shape.

Currently, `RecurrenceMicrostatesAnalysis.jl` supports five shapes: square, triangle, diagonal, line, and pair. The way the library converts these motifs constrained shapes to decimal values is detailed on the [motifs](motifs.md) page. You can change the shape using the **kword** `shape`, which can be set to `:square`, `:triangle`, `:diagonal`, `:line`, or `:pair`. By default, the library uses `:square` as the default shape.
```@repl dist_one
dist = distribution(data, th, 3; shape = :triangle)
```

!!! note
    The shape `:pair` doesn't require a value of `n`, since it always uses `n=2`. However, it is still necessary to informe a value to this parameter, that will be interpreted as the separation between two points in a diagonal.
    ```@repl dist_one
    dist = distribution(data, th, 6; shape = :pair)
    ```

    When workign with shape `:pair`, we recommend you to use the full structure of `distribution` function.
    ```@repl dist_one
    structure = [3, 9]
    dist = distribution(data, data, th, structure; shape = :pair)
    ```
    Here, `structure` defines the position of the second element based on the random position of the first element.

### Motifs sampling
The sampling mode defines how `RecurrenceMicrostatesAnalysis.jl` selects motifs from a recurrence space. Currently, the library supports four sampling modes: full, random, columnwise and triangle up. You can learn more about them on the [motifs](motifs.md) page, where we discuss how each mode works. The sampling mode can be configured using the keyword argument `sampling_mode`, which can be set to `:full`, `:random`, `:columnwise`, `:columnwise_full`, or `:triangleup`. By default, the library uses `:random` as the default sampling mode.

```@repl dist_one
dist = distribution(data, th, 3; sampling_mode = :full)
```

!!! compat
    Not all sampling modes are compatible with certain motif constrained shapes, and the following table illustrates the compatibility between them.

    |             | `:full`    | `:random`  | `:columnwise` | `:columnwise_full` | `:triangleup` |
    |-------------|:----------:|:----------:|:-------------:|:------------------:|:-------------:|
    | `:square`   |$\checkmark$|$\checkmark$|$\checkmark$   |$\checkmark$        |$\checkmark$   |
    | `:triangle` |$\checkmark$|$\checkmark$|               |                    |               |
    | `:diagonal` |            |$\checkmark$|               |                    |               |
    | `:time`     |            |$\checkmark$|               |                    |               |
    | `:pair`     |            |$\checkmark$|$\checkmark$   |                    |               |

### Run mode
`RecurrenceMicrostatesAnalysis.jl` has two run modes that results in a different output type. The run mode `:vect` allocates all required memory in beginning of the process, and return the distribution as a vector. This is the default configuration of the library for $n < 6$.
```@repl dist_one
dist = distribution(data, th, 4; run_mode = :vect)
```

The run mode `:dict` uses dictionaries to allocate memory just when needed. The total allocation of dictionary mode can be greater than when using vectors, but the real memory allocation is smaller.
```@repl dist_one
dist = distribution(data, th, 4; run_mode = :dict)
```

!!! compat
    It is important to note that the shapes `:diagonal`, `:line`, and `:pair` are not compatible with run mode `:dict`. Additionally, sampling modes `:columnwise` and `:columnwise_full` return a matrix in which each column represents a probability distribution for a specif time value.
    ```@repl dist_one
    dist = distribution(data, th, 2; sampling_mode = :columnwise)
    dist = distribution(data, th, 2; sampling_mode = :columnwise_full)
    ```

### Number of samples
With exception of sampling modes `:full` and `:columnwise_full`, all sampling modes take motifs randomly in a recurrence space. The **kword** `num_samples` defines the number of samples that will be used by the library, it can be either an integer value that specifies the exact number, or a decimal value interpreted as the percentage of samples taken from the entire available population. By default, `RecurrenceMicrostatesAnalysis.jl` uses $5\%$.
```@repl dist_one
dist = distribution(data, th, 3; num_samples = 0.1)
```
```@repl dist_one
dist = distribution(data, th, 3; num_samples = 50000)
```

### Threads
`RecurrenceMicrostatesAnalysis.jl` is highly compatible with CPU asynchronous jobs, that can increase significantly the computational performance of the library. The **kword** `threads` defines if the library will use threads or not, being `true` by default. The number of threads used is equal to the number of threads available to Julia, being it configured by the environment variable `JULIA_NUM_THREADS`, or by the running argument `--threads T` in Julia initiation: For example, using `julia --threads 8`.
```@repl dist_one
using BenchmarkTools
@benchmark distribution(data, th, 4; sampling_mode = :full, threads = false)
@benchmark distribution(data, th, 4; sampling_mode = :full, threads = true)
```

!!! warning
    `RecurrenceMicrostatesAnalysis.jl` allocates memory for each thread, so how many threads you use, more memory the library will allocate. It is done to increase the performance, and avoid the memory concurrency.

### Metrics
`RecurrenceMicrostatesAnalysis.jl` uses the library [Distances.jl](https://github.com/JuliaStats/Distances.jl) to simplify the configuration of metrics, and increase the computation performance. With it, modify the metric is a easy process that can be done with the **kword** `metric`.
```@repl dist_one
using Distances
my_metric = KLDivergence()
dist = distribution(data, th, 2; metric = my_metric)
```

!!! warning
    The default recurrence functions were configured to metrics with two arguments, like `euclidean(x, y)`, so if you need to use another type of metric, it is needed to define a new recurrence function, see [Recurrence functions](recurrence.md) page to know more about it.

### Recurrence functions
A recurrence function defines if two points of a trajectory recurr or not. Actually the library have two recurrence functions available
1. Standard recurrence: $R(\mathbf{x}, \mathbf{y})=\Theta(\varepsilon - \|\mathbf{x}-\mathbf{y}\|)$
2. Recurrence with corridor threshold: $R(\mathbf{x}, \mathbf{y})=\Theta(\|\mathbf{x}-\mathbf{y}\| - \varepsilon_{min}) \cdot \Theta(\varepsilon_{max} - \|\mathbf{x}-\mathbf{y}\|)$

`RecurrenceMicrostatesAnalysis.jl` automatically change between them with the type of parameters, so if you use as parameter a `Float64`, the library will apply the standard recurrence, or, if you use a `Tuple{Float64, Float64}`, the library will apply the recurrence with corridor threshold.

```@repl dist_one
dist = (distribution(data, th, 2))'
dist = (distribution(data, (0.0, th), 2))'
```

It is possible to write your own recurrece function, we talk more about it in the [Recurrence functions](recurrence.md) page.

##  High-dimensionality data
If you are working with a dynamical system or a data time serie with two or more dimensions, it is important to note that `RecurrenceMicrostatesAnalysis.jl` effectively not works with vectors, but matrices. In this situation, each row of the matrix will represent a coordinate, and each column a set of coordinates along a trajectory. For example, if we want a uniform distribution with three dimension and 3,000 points, we will have something like:
```@repl dist_high
using Distributions
data = rand(Uniform(0, 1), 3, 3000)
```

This format of data is effectvely what the library uses. In the case of previous section, when we are working with vectors, `RecurrenceMicrostatesAnalysis.jl` converts it to a matrix $1\times 3,000$ but when we are working which data with a dimensionality different than one, it is necessary to use the proper format.
```@repl dist_high
using RecurrenceMicrostatesAnalysis
th, s = find_parameters(data, 3)
dist = distribution(data, th, 3)
```

##  Continuous problems
Continuous problems means numerically integrate a differential equation problem and take the values as input to `RecurrenceMicrostatesAnalysis.jl`. Thinking in it, we make the library compatible with a powerful tool to solve these problems in Julia: the library [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/). The way to apply this kind of data in the library is similar with the other two cases discussed before, as we will demonstrate in this section. 

!!! info
    The code of Lorenz system used in these examples was get from Example 2 of [DifferentialEquations.jl documentation](https://docs.sciml.ai/DiffEqDocs/stable/getting_started/)

```@repl continuous
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

using DifferentialEquations
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 1000.0)
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob)
```


With the data computed, it is easy to apply to `RecurrenceMicrostatesAnalysis.jl`, with a simply memory access given by `DifferentialEquations.jl`.
```@repl continuous
using RecurrenceMicrostatesAnalysis
data = sol[:, :]
th, s = find_parameters(data, 3)
dist = distribution(data, th, 3)
```

!!! warning
    Although it is possible to compute the distribution as demonstrated above, we strongly advise against doing so in this way.

We recommend you to apply a transient into your data and take a correct time resolution while doing the process of discretization, it is needed to maximize the information available. `RecurrenceMicrostatesAnalysis.jl` has a utilitary function to help with this process.
```@repl continuous
prepared_data = prepare(sol, 0.2; transient = 10000, max_length = 1000)
th, s = find_parameters(prepared_data, 3)
dist = distribution(prepared_data, th, 3)
```

If you have the threshold parameter, it is also possible to simplify the call using:
```@repl continuous
dist = distribution(sol, th, 3, 0.2; transient = 10000, max_length = 1000)
```

##  Spatial data
`RecurrenceMicrostatesAnalysis.jl` is compatible with generalised recurrence plot analysis for spatial data proposed by Marwan, Kurths and Saparin at 2006 [Marwan2006](@cite). It allow the library to calculate a probability distribution of motifs in a tensorial recurrence space, for example, to images the recurrence space have four dimensions.

!!! todo
    Since this is an open research field, the library is designed to support exploration and estimation for research purposes. We don’t recommend applying it in production environments 😉

The application of `RecurrenceMicrostatesAnalysis.jl` to spatial data is very similar to the others presented before, but the input format is more complex. Instead to matrices we need to use abstract arrays with dimension $D$, where the first dimension will be interpreted as a coordinate dimension (such as for high-dimensionaly data), and rest of the dimensions will be the spatial data dimensionality. To illustrate it, let an image with RGB. It can be represented as an abstract array with 3 dimensions, where the first dimension will have a length 3, being each element a color value (red, blue and green), and the others two dimensions are relative to each pixel that compose the image. We will demonstrate it using a uniform distribution, where each position can be interpreted as a RGB pixel for an image 100x100.
```@repl spatial
using Distributions
data = rand(Uniform(0, 1), 3, 100, 100)
```

When we work with spatial data is necessarity to use the complete structure of `distribution` function, defining a vector structure where each value represents the length of a motif constrained side. For example, to a square tensorial motif constrained with side 2, we can use:
```@repl spatial
using RecurrenceMicrostatesAnalysis
dist = distribution(data, data, 0.5, [2, 2, 2, 2])
```

Since the recurrence space has four dimensions, in this examples, it is necessary for `structure` has the same number of elements, where each element will represent the motif' side lenght for each dimension.

!!! warning
    The `find_parameters` function is not compatible with spatial data.

!!! compat
    It is important to note that this functionality is only available to motif shapes `:square`, `:diagonal`, `:line` and `:pair`, for `:random` sampling mode.