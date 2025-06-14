using RMA, Test
using Distances
using Distributions
using RecurrenceAnalysis

##
##      Data
uniform = rand(Uniform(0, 1), 5000)

##
##      Base
rp = RecurrenceMatrix(uniform, 0.27)
based = rqa(rp)
det_base = based[:DET]
lam_base = based[:LAM]

##
##      RMA
det_rma = RMA.determinism(uniform, 0.27)
lam_rma = RMA.laminarity(uniform, 0.27)

error_det = abs(det_rma - det_base) / det_base
error_lam = abs(lam_rma - lam_base) / lam_base

@testset "rqa" begin
    @test error_det < 0.01
    @test error_lam < 0.01
end