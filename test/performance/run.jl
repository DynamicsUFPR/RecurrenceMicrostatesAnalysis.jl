using RMA
using JLD2
using Colors
using CairoMakie
using Distributions
using ProgressMeter
using BenchmarkTools
using RecurrenceAnalysis

##
##      Settings
data_len = round.(Int, range(500, 10000, 30))
samp_range = range(0.05, 0.5, 5)

data = rand(Uniform(0, 1), data_len[1])
thres = 0.0
samp_percent = 0.0

##
##      Compute using recurrence analysis.
function classic_compute()
    rp = RecurrenceMatrix(data, thres)
    det = RecurrenceAnalysis.determinism(rp; l_min = 2)
    lam = RecurrenceAnalysis.laminarity(rp; l_min = 2)

    return det, lam
end

##
##      Compute using RMA
function rma_compute()
    det = RMA.determinism(data, thres; num_samples = samp_percent)
    lam = RMA.laminarity(data, thres; num_samples = samp_percent)

    return det, lam
end

##
##      Compute the performance.
function compute()
    results = zeros(Float64, 1 + length(samp_range), 2, length(data_len))

    @showprogress for i in eachindex(data_len)
        global data = rand(Uniform(0, 1), data_len[i])

        th, _ = RMA.find_parameters(data, 3)
        global thres = th

        ##
        ##      First, using Recurrence RecurrenceAnalysis.
        b_c = @benchmark classic_compute()
        results[1, 1, i] = (b_c.times |> mean) * (10^(-9))             # In s
        results[1, 2, i] = (b_c.memory) / (1024 * 1024 * 1024)         # In GB

        ##
        ##      Now, for RMA.
        for s in eachindex(samp_range)
            global samp_percent = samp_range[s]

            b_r = @benchmark rma_compute()
            results[s + 1, 1, i] = (b_r.times |> mean) * (10^(-9))
            results[s + 1, 2, i] = (b_r.memory) / (1024 * 1024 * 1024)
        end
    end

    save_object("data.jld2", results)
end

function graph()
    dat = load_object("test/performance/data.jld2")

    ##
    ##      Performance - CPU
    fig_cpu = Figure()
    rowsize!(fig_cpu.layout, 1, Relative(0.1))

    ax_cpu = Axis(fig_cpu[2, 1], xlabel = "Data length", ylabel = "Mean calculation time [sec]")
    ##  - RMA
    cpu_rma_sn_1 = scatter!(ax_cpu, data_len, dat[2, 1, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#7fc97f", strokewidth = 1)
    cpu_rma_sn_2 = scatter!(ax_cpu, data_len, dat[3, 1, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#beaed4", strokewidth = 1)
    cpu_rma_sn_3 = scatter!(ax_cpu, data_len, dat[4, 1, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#fdc086", strokewidth = 1)
    cpu_rma_sn_4 = scatter!(ax_cpu, data_len, dat[5, 1, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#eaea72", strokewidth = 1)
    cpu_rma_sn_5 = scatter!(ax_cpu, data_len, dat[6, 1, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#386cb0", strokewidth = 1)

    ##  - RecurrenceAnalysis
    cpu_rma_ca = scatter!(ax_cpu, data_len, dat[1, 1, :], marker = :xcross, color = colorant"#f0027f",)

    Legend(fig_cpu[1, 1], 
        [cpu_rma_ca, cpu_rma_sn_1, cpu_rma_sn_2, cpu_rma_sn_3, cpu_rma_sn_4, cpu_rma_sn_5], 
        ["Standard approach", 
        string("RMA, sn = ", round(samp_range[1] * 100; digits = 1), "%"),
        string("RMA, sn = ", round(samp_range[2] * 100; digits = 1), "%"),
        string("RMA, sn = ", round(samp_range[3] * 100; digits = 1), "%"),
        string("RMA, sn = ", round(samp_range[4] * 100; digits = 1), "%"),
        string("RMA, sn = ", round(samp_range[5] * 100; digits = 1), "%")]; 
        tellheight = false, tellwidth = false, orientation = :horizontal, nbanks = 2, labelsize = 12)

    save("test/performance/fig_cpu.png", fig_cpu)

    ##
    ##      Performance - CPU
    fig_ram = Figure()
    rowsize!(fig_ram.layout, 1, Relative(0.1))

    ax_ram = Axis(fig_ram[2, 1], xlabel = "Data length", ylabel = "Total Allocated Memory [GiB]")
    ##  - RMA
    ram_rma_sn_1 = scatter!(ax_ram, data_len, dat[2, 2, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#7fc97f", strokewidth = 1)
    ram_rma_sn_2 = scatter!(ax_ram, data_len, dat[3, 2, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#beaed4", strokewidth = 1)
    ram_rma_sn_3 = scatter!(ax_ram, data_len, dat[4, 2, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#fdc086", strokewidth = 1)
    ram_rma_sn_4 = scatter!(ax_ram, data_len, dat[5, 2, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#eaea72", strokewidth = 1)
    ram_rma_sn_5 = scatter!(ax_ram, data_len, dat[6, 2, :], marker = :rect, markersize = 10, color = :white, strokecolor = colorant"#386cb0", strokewidth = 1)

    ##  - RecurrenceAnalysis
    ram_rma_ca = scatter!(ax_ram, data_len, dat[1, 2, :], marker = :xcross, color = colorant"#f0027f",)

    Legend(fig_ram[1, 1], 
        [ram_rma_ca, ram_rma_sn_1, ram_rma_sn_2, ram_rma_sn_3, ram_rma_sn_4, ram_rma_sn_5], 
        ["Standard approach", 
        string("RMA, sn = ", round(samp_range[1] * 100; digits = 1), "%"),
        string("RMA, sn = ", round(samp_range[2] * 100; digits = 1), "%"),
        string("RMA, sn = ", round(samp_range[3] * 100; digits = 1), "%"),
        string("RMA, sn = ", round(samp_range[4] * 100; digits = 1), "%"),
        string("RMA, sn = ", round(samp_range[5] * 100; digits = 1), "%")]; 
        tellheight = false, tellwidth = false, orientation = :horizontal, nbanks = 2, labelsize = 12)

    save("test/performance/fig_ram.png", fig_ram)
end

graph()