#   Recurrence Motifs Distributions in One Minute
The following code is an example of how to use `RecurrenceMicrostatesAnalysis.jl` to compute the motif distribution of a uniform distribution. Try pasting it into the REPL prompt 😉.

```julia
##  Install everything that we need
using Pkg; Pkg.add("Distributions"); Pkg.add(url="https://github.com/DynamicsUFPR/RMA.jl")
using Distributions                             #   For generate our uniform distribution
import RecurrenceMicrostatesAnalysis as RMA     #   !! Import RecurrenceMicrostatesAnalysis.jl

##  Generate our data
data = rand(Uniform(0, 1), 1000)

##  Square motif side
n = 3

##  Compute the threshold that maximize the recurrence entropy.
th, s = find_parameters(data, n)

##  Compute the recurrence motif probabilities distribution.
dist = distribution(data, th, n)
```

The output will be a set of 512 probabilities. The `distribution` function samples $5\%$ of all motifs available in the RP, regardless of overlap or repetition. Each value represents the probability of encountering a motif with a decimal representation $I$ within the RP.

It is not necessary to compute the Recurrence Plot (RP), as the library calculates the recurrences internally without needing to construct one explicitly.

!!! note
    It is important to remember that Julia's indexing starts at $1$ instead of $0$. Therefore, in the library, we define $I = i - 1$, where $i$ is the Julia's index.

###   Easy RQA Estimation
One of the main goals of `RecurrenceMicrostatesAnalysis.jl` is to provide high performance when estimating typical RQA quantifiers, such as determinism (DET) and laminarity (LAM). This process is very simple, as you can see in the following code, give it a try 😁.
```julia
entropy = rentropy(dist)        #   Recurrence entropy
rr = rrate(dist)                #   Recurrence rate
det = determinism(rr, dist)     #   Determinism
lam = laminarity(rr, dist)      #   Laminarity
```

It is also possible to skip computing the recurrence distribution by using an alternative overload of the `determinism` and `laminarity` functions.
```julia
det = determinism(data, th)
lam = laminarity(data, th)
```