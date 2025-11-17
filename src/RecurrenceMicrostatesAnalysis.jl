#
#           RecurrenceMicrostatesAnalysis.jl
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
module RecurrenceMicrostatesAnalysis
    #
    #       Libraries needed for work.
    using Distances
    using Random
    using Statistics

    #       Basic metric
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
    include("rma/histograms/square/square_columnwise_full.jl")
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: triangle
    include("rma/histograms/triangle/triangle_full.jl")
    include("rma/histograms/triangle/triangle_random.jl")
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: time pair
    include("rma/histograms/pair/pair_random.jl")
    include("rma/histograms/pair/pair_columnwise.jl")
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: diagonal
    include("rma/histograms/diagonal/diagonal_full.jl")
    include("rma/histograms/diagonal/diagonal_random.jl")
    # ------------------------------------------------------------------------------------------------------- #
    #       - Shape: line
    include("rma/histograms/line/line_random.jl")
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
    #       - Disorder (Ξ)
    include("rqa/disorder.jl")
    # ======================================================================================================= #
    #           * RMA Utils
    include("utils/prepare.jl")
    include("utils/find_parameters.jl")
    # ======================================================================================================= #
    using .Disorder
    # ======================================================================================================= #
    #
    #       Export some functions to the main scope (PublicAPI)
    # ======================================================================================================= #
    export rrate
    export prepare
    export rentropy
    export disorder
    export laminarity
    export recurrence
    export determinism
    export distribution
    export find_parameters

    export jrp
    # ======================================================================================================= #
end
#       ----- END CODE