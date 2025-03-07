#
#           RecurrenceMicrostates.jl
#       DynamicsUFPR - https://github.com/DynamicsUFPR
#       
#
#           Abstract: The idea of this project is to create a standard library that
#   can be used in the application and development of Recurrence Microstates Theory
#   and its associated quantifiers. Here we try to generalize the computational process
#   so that it can be applied to any data situation, such as time series, images, or 
#   high-dimensional data.
#
#       GitHub - Julia Version: https://github.com/DynamicsUFPR/RMA.jl
#
#       ----- BEGIN CODE
#module RMA
    #
    #       Libraries needed for the code to work.
    using Distances

    #       Metrics
    const euclidean_metric = Euclidean()

    #
    #       Import the code source.
    #   - Probabilities (RMA)
    include("rma/distribution.jl")
    include("rma/index/square_index.jl")
    include("rma/histograms/vect/square_full.jl")
    include("rma/histograms/vect/square_random.jl")
    include("rma/histograms/vect/square_random_async.jl")
    #   - Utils
    include("rma/recurrence.jl")

    #
    #       Export some functions.
    #   - Probabilities (RMA)
    export distribution
#end
#       ----- END CODE