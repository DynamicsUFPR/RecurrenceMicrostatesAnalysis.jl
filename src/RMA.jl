#
#           RMA.jl
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
    # ======================================================================================================= #
    #           * RMA Core
    #   i. Import the recurrence functions.
    include("rma/recurrence.jl")
    #   ii. Import the functions that convert a motif to an index.
    include("rma/index.jl")
    #   iii. Import the distribution function.
    include("rma/distribution.jl")
    #   iv. Import histogram functions
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: square
    include("rma/histograms/square/square_full.jl")
    include("rma/histograms/square/square_random.jl")
    include("rma/histograms/square/square_triangleup.jl")
    include("rma/histograms/square/square_columnwise.jl")
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: triangle
    include("rma/histograms/triangle/triangle_full.jl")
    include("rma/histograms/triangle/triangle_random.jl")
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: time pair
    include("rma/histograms/timepair/timepair_random.jl")
    include("rma/histograms/timepair/timepair_columnwise.jl")
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: diagonal
    include("rma/histograms/diagonal/diagonal_random.jl")
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: column
    include("rma/histograms/column/column_random.jl")
    # ======================================================================================================= #
    #           * RMA Analysis
    #       - Recurrence Rate (RR)
    include("rqa/rr.jl")
    #       - Determinism (DET)
    include("rqa/det.jl")
    #       - Laminarity (LAM)
    include("rqa/lam.jl")
    #       - Recurrence Entropy (RETR)
    include("rqa/entropy.jl")
    # ======================================================================================================= #
    #           * RMA Utils
    include("utils/prepare.jl")
    # ======================================================================================================= #
    #
    #       Export some functions to the main scope.
    # ======================================================================================================= #
    export rrate
    export prepare
    export rentropy
    export laminarity
    export determinism
    export distribution
    # ======================================================================================================= #
end
#       ----- END CODE