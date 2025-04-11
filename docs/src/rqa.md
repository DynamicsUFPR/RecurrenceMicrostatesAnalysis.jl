#   Recurrence Quantification Analysis
The recurrence microstate analysis allows us to estimate values of tipical RQA measures, such as determinism and laminarity, with a good precision, and define some novels quantifier. We will demonstate in this page how compute these quantifiers using a uniform distribution as input.
```@repl rqa
using RMA, Distributions
data = rand(Uniform(0, 1), 3000)
```
##  Recurrence Entropy

##  Recurrence Rate

##  Determinism

##  Laminarity