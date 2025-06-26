#   A Julia library for analyzing dynamical systems with recurrence microstates
Recurrence Microstates Analysis (RMA) is an advanced approach that generalizes the analysis of recurrence structures by capturing the statistical properties of recurrence motifs. `RecurrenceMicrostatesAnalysis.jl` is an efficient Julia package for performing RMA, offering support for a wide range of motif shapes, flexible sampling strategies, and comprehensive distribution computation capabilities. Furthermore, the library features an optimized pipeline for estimating standard RQA quantifiers, with significantly reduced memory and computational requirements, making it particularly well-suited for large-scale datasets.

##  Installation
Download [Julia 1.8](https://julialang.org/) or later, preferably the current stable release. 

You can add `RecurrenceMicrostatesAnalysis.jl` using Julia's package manager. In the Julia prompt, you can use the following code snippets:
```julia
using Pkg
Pkg.add("RecurrenceMicrostatesAnalysis")
```
or, in `Pkg REPL` mode write:
```julia
] add RecurrenceMicrostatesAnalysis
```

!!! todo "GitHub"
    `RecurrenceMicrostatesAnalysis.jl` is an open-source library available at GitHub repository [DynamicsUFPR/RMA.jl](https://github.com/DynamicsUFPR/RMA.jl). If you have found this library useful, please consider starring it on [GitHub](https://github.com/DynamicsUFPR/RMA.jl) 😉.
```

##  Learning RMA
If you have worked with recurrence microstates analysis before, the [Quick Start](quickstart.md) page offers a brief guide on how to apply the `RecurrenceMicrostatesAnalysis.jl` to time series data and dynamical systems.

If you haven't, then you might prefer the [Theoretical Overview](theory.md) page, which provides a quick and simple introduction about the recurrence microstates field. The rest of the guide explains how to use the library to compute recurrence motifs probability distributions and calculate common recurrence quantifiers, along with descriptions of all available configuration options. We also include the [Utils](utils.md) page, which covers utility functions to simplify the use of `RecurrenceMicrostatesAnalysis.jl`, and the [Performance Tips](performance.md) page, where we discuss how to improve the library’s usage performance.