# RecurrenceMicrostatesAnalysis.jl

![GitHub repo size](https://img.shields.io/github/repo-size/DynamicsUFPR/RMA.jl)
![GitHub license](https://img.shields.io/github/license/DynamicsUFPR/RMA.jl)
[![Julia](https://img.shields.io/badge/Julia-1.8%2B-blue?logo=julia)](https://julialang.org/)

![Library logo](doc/logo.png)

A simple and faster library for recurrence microstates analysis.

### 📦 Dependencies

The library uses the package Distances.jl, that is installed along with the library when you use Julia's Pkg.

[![Julia](https://img.shields.io/badge/Julia-Package-red?logo=julia)](https://juliahub.com/ui/Packages/Distances)


### ⚙️ Installation

1. Using `Pkg.add`:
  With the Julia terminal open, type:

```julia
using Pkg
Pkg.add("RecurrenceMicrostatesAnalysis")
```

2. Using the Pkg REPL mode:

```julia
] add RecurrenceMicrostatesAnalysis
```

## Library usage guide

Access the library documentation [here](https://dynamicsufpr.github.io/RecurrenceMicrostatesAnalysis.jl/).

##  Citation
```bibtex
@article{10.1063/5.0293708,
    author = {Vinicius Ferreira, Gabriel and Lopes da Cruz, Felipe Eduardo and Marghoti, Gabriel and de Lima Prado, Thiago and Roberto Lopes, Sergio and Marwan, Norbert and Kurths, Jürgen},
    title = {RecurrenceMicrostatesAnalysis.jl: A Julia library for analyzing dynamical systems with recurrence microstates},
    journal = {Chaos: An Interdisciplinary Journal of Nonlinear Science},
    volume = {35},
    number = {11},
    pages = {113123},
    year = {2025},
    month = {11},
    issn = {1054-1500},
    doi = {10.1063/5.0293708},
    url = {https://doi.org/10.1063/5.0293708}
}
```
