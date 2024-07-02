using Plots
using DataFrames
using CSV
using Statistics
using StatsPlots
using StatsBase
using ColorSchemes
using Measures

# Bar graph with error
function bar_with_err(ptitle, subtitle, xlabs, ylab, vals, errs)
    pallete = cgrad(:blues, 3)
    p= groupedbar(
        xlabs,
        vals,
        size=(1000,800),
        yerr=(zeros(12),errs),  # Bottom error = 0
        group = repeat(["1-Low", "2-Med", "3-High"], outer=4),
        color = [pallete[i] for i ∈ repeat([1,2,3],outer=4)],
        legendtitle = "Supply Strategy",
        margin=10mm,
        xtickfontsize = 12,
        ytickfontsize = 12,
        legendfontsize = 12,
        legend = (occursin("Total Time", ptitle) || occursin("Drug Overage", ptitle)  ? :bottomright : :best) # Fix legend position for graphs where all bars are near the top
        
    )
    title!("\n"*ptitle*"\n"*subtitle)
    ylabel!(ylab)
    xlabel!("IRT Specification")

    return p
end

# 3x4 image of all time histograms
function time_hists(times::Array, fname::String)
    plts = []
    pallete = cgrad(:blues, 3)

    titles = ["F0 Low Supply", "F0 Medium Supply", "F0 High Supply","F1 Low Supply", "F1 Medium Supply", "F1 High Supply",
    "F2a Low Supply", "F2a Medium Supply", "F2a High Supply","F2b Low Supply", "F2b Medium Supply", "F2b High Supply"]

    for i in 1:length(eachrow(times))
        p = histogram(times[i,:], legend=false, bins=20, color=pallete[(i-1)%3+1])
        ylabel!("Count")
        xlabel!("Time to complete recruitment (days)")
        title!(titles[i])
        annotate!([
            (xlims()[2]-0.05*xlims()[2], ylims()[2], text(string("Min: ",round(minimum(times[i,:]),digits=2)), :black, :right, 10)),
            (xlims()[2]-0.05*xlims()[2], ylims()[2]-0.05*ylims()[2], text(string("Q25: ",round(quantile(times[i,:], 0.25),digits=2)), :black, :right, 10)),
            (xlims()[2]-0.05*xlims()[2], ylims()[2]-0.1*ylims()[2], text(string("Mean: ",round(mean(times[i,:]),digits=2)), :black, :right, 10)),
            (xlims()[2]-0.05*xlims()[2], ylims()[2]-0.15*ylims()[2], text(string("Median: ",round(median(times[i,:]),digits=2)), :black, :right, 10)),
            (xlims()[2]-0.05*xlims()[2], ylims()[2]-0.2*ylims()[2], text(string("Q75: ",round(quantile(times[i,:], 0.75),digits=2)), :black, :right, 10)),
            (xlims()[2]-0.05*xlims()[2], ylims()[2]-0.25*ylims()[2], text(string("Max: ",round(maximum(times[i,:]),digits=2)), :black, :right, 10)),
        ])

        push!(plts, p)
    end

    t_plot = plot(plts..., layout = (4,3), size=(2000,2000))
    savefig(t_plot, fname)
end

# Export all statistics for time distribution to csv
function export_time_stats_csv(times::Array, fname::String)
    irt_approaches = String["F0","F0","F0","F1","F1","F1","F2a","F2a","F2a","F2b","F2b","F2b"]
    resupp_strategies = String["Low","Medium","High","Low","Medium","High","Low","Medium","High","Low","Medium","High"]
    mat = [irt_approaches, resupp_strategies, minimum.(eachrow(times)), quantile.(eachrow(times), 0.25), median.(eachrow(times)) ,
            quantile.(eachrow(times), 0.75), maximum.(eachrow(times)), mean.(eachrow(times)), std.(eachrow(times)), quantile.(eachrow(times), 0.9)]
    df = DataFrame(mat, ["IRT Approach","Re-supply Strategy","Min", "Q25", "Median", "Q75", "Max", "Mean", "SD", "Q90"])
    CSV.write(fname, df)
end

# Export "stat" of all characterstics to csv
function export_csv_stat(chars::Array, labels::Array, fname::String, stat::Function)
    irt_approaches = ["F0" "F0" "F0" "F1" "F1" "F1" "F2a" "F2a" "F2a" "F2b" "F2b" "F2b"]
    resupp_strategies = ["Low" "Medium" "High" "Low" "Medium" "High" "Low" "Medium" "High" "Low" "Medium" "High"]
    mat = vcat(irt_approaches, resupp_strategies, reduce(hcat, stat(chars[:,1,:],dims=2)), reduce(hcat, stat(chars[:,2,:],dims=2)), reduce(hcat, stat(chars[:,3,:],dims=2)), reduce(hcat, stat(chars[:,4,:],dims=2)), 
    reduce(hcat, stat(chars[:,5,:],dims=2)), reduce(hcat, stat(chars[:,6,:],dims=2)), reduce(hcat, stat(chars[:,7,:],dims=2)), reduce(hcat, stat(chars[:,8,:],dims=2)))
    df = DataFrame(permutedims(mat,(2,1)), labels)
    CSV.write(fname, df)
end
