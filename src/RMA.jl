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
module RMA
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
    include("rma/index/triangle_index.jl")
    #       - Vectors
    include("rma/histograms/vect/square_full.jl")
    include("rma/histograms/vect/square_random.jl")
    include("rma/histograms/vect/square_triangleup.jl")
    include("rma/histograms/vect/square_columnwise.jl")
    include("rma/histograms/vect/triangle_full.jl")
    include("rma/histograms/vect/triangle_random.jl")
    #       - Vectors async.
    include("rma/histograms/vect/square_random_async.jl")
    include("rma/histograms/vect/square_triangleup_async.jl")
    include("rma/histograms/vect/square_columnwise_async.jl")
    include("rma/histograms/vect/triangle_random_async.jl")
    #       - Dict
    include("rma/histograms/dict/square_full.jl")
    include("rma/histograms/dict/square_random.jl")
    include("rma/histograms/dict/square_triangleup.jl")
    include("rma/histograms/dict/triangle_full.jl")
    include("rma/histograms/dict/triangle_random.jl")
    #       - Dict async.
    include("rma/histograms/dict/square_random_async.jl")
    include("rma/histograms/dict/square_triangleup_async.jl")
    include("rma/histograms/dict/triangle_random_async.jl")
    #   - Utils
    include("rma/recurrence.jl")
    include("utils/prepare.jl")

    #   - RQA
    include("rqa/recurrence_entropy.jl")
    include("rqa/recurrence_rate.jl")
    include("rqa/determinism.jl")

    #
    #       Export some functions.
    #   - Probabilities (RMA)
    export distribution
    #   - Quantifiers (RQA)
    export recurrence_entropy
end
#       ----- END CODE