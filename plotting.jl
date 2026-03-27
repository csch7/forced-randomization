using Plots
using DataFrames
using CSV
using Statistics
using StatsPlots
using StatsBase
using ColorSchemes
using Measures
using Distributions
using LinearAlgebra

# Grouped bar chart with one-sided upper error bars.
function grouped_bar_with_error(ptitle, subtitle, xlabs, ylab, vals, errs)
    palette = cgrad(:blues, 3)
    p = groupedbar(
        xlabs,
        vals,
        size          = (1000, 800),
        yerr          = (zeros(12), errs),
        group         = repeat(["1-Low", "2-Med", "3-High"], outer=4),
        color         = [palette[i] for i in repeat([1, 2, 3], outer=4)],
        legendtitle   = "Supply Strategy",
        margin        = 10mm,
        xtickfontsize = 12,
        ytickfontsize = 12,
        legendfontsize = 12,
        legend        = (occursin("Total Time", ptitle) || occursin("Drug Overage", ptitle) ? :bottomright : :best),
    )
    title!("\n" * ptitle * "\n" * subtitle)
    ylabel!(ylab)
    xlabel!("IRT Specification")
    return p
end

# Returns a single line-plot panel for one statistic over patient index m.
# dmzs: (1, simulations, patients). func is applied over the simulations dimension.
# cutoff: optional x-axis limit (trims the series at that patient index).
function imbalance_line_panel(
    dmzs::Array,
    ylabel_str::String,
    func::Function;
    max_value = :auto,
    cutoff::Union{Int,Nothing} = nothing,
    start::Int = 1,
    logy::Bool = false
)
    color     = cgrad(:blues)[0.8]
    stat_vals = vec(func(dmzs, dims=2))
    n         = cutoff === nothing ? length(stat_vals) : min(cutoff, length(stat_vals))
    xs        = start:n
    vals      = stat_vals[xs]
    p = logy ?
        plot(xs, vals, legend=false, color=color, yaxis=(:log10, [0.01, max_value])) :
        plot(xs, vals, legend=false, color=color, ylims=(:auto, max_value))
    ylabel!(p, ylabel_str)
    xlabel!(p, "m")
    title!(p, ylabel_str)
    return p
end

# Saves a grid of line panels with an overall title and generous subplot spacing.
function save_line_summary(panels::Vector, nrows::Int, ncols::Int, filepath::String, title::String)
    t_plot = plot(panels...,
        layout      = (nrows, ncols),
        size        = (600 * ncols, 450 * nrows),
        plot_title  = title,
        margin      = 8mm,
        left_margin = 12mm,
    )
    savefig(t_plot, filepath)
end

# Grid of histograms of treatment imbalance at several patient indices m.
# dmz1s/dmz2s: (1, simulations, patients) arrays for each stratum.
# Shows 2 rows (one per stratum) x 4 columns (patient positions near end of trial).
function plot_imbalance_histograms(dmz1s, dmz2s, min_m1, min_m2, filepath, title)
    plts = []
    for j in 1:4
        push!(plts, imbalance_histogram(
            dmz1s[1, :, j + min_m1],
            title * " (m=" * string(j + min_m1) * ", z=1)"
        ))
    end
    for j in 1:4
        push!(plts, imbalance_histogram(
            dmz2s[1, :, j + min_m2],
            title * " (m=" * string(j + min_m2) * ", z=2)"
        ))
    end
    t_plot = plot(plts..., layout=(2, 4), margin = 8mm, left_margin = 12mm, size=(2000, 1300))
    savefig(t_plot, filepath)
end

function imbalance_histogram(dmzs::Array, title)
    lo = floor(Int, minimum(dmzs))
    hi = ceil(Int, maximum(dmzs))
    bins = (lo - 0.5):1:(hi + 0.5)
    p = histogram(dmzs, legend=false, bins=bins)
    ylabel!("Count")
    xlabel!("dm(z)")
    title!(title)
    return p
end

# 4×3 grid of recruitment-time histograms, one per scenario.
function plot_recruitment_time_histograms(times::Array, filepath::String)
    palette = cgrad(:blues, 3)
    titles  = [
        "F0a Low Supply",  "F0a Medium Supply",  "F0a High Supply",
        "F0b Low Supply",  "F0b Medium Supply",  "F0b High Supply",
        "F1a Low Supply",  "F1a Medium Supply",  "F1a High Supply",
        "F1b Low Supply",  "F1b Medium Supply",  "F1b High Supply",
    ]
    plts = []
    for i in 1:size(times, 1)
        row = times[i, :]
        p   = histogram(row, legend=false, bins=20, color=palette[(i-1)%3+1])
        ylabel!("Count")
        xlabel!("Time to complete recruitment (days)")
        title!(titles[i])
        annotate!([
            (xlims()[2]*0.95, ylims()[2]*1.00, text("Min:    " * string(round(minimum(row),    digits=2)), :black, :right, 10)),
            (xlims()[2]*0.95, ylims()[2]*0.95, text("Q25:    " * string(round(quantile(row, 0.25), digits=2)), :black, :right, 10)),
            (xlims()[2]*0.95, ylims()[2]*0.90, text("Mean:   " * string(round(mean(row),     digits=2)), :black, :right, 10)),
            (xlims()[2]*0.95, ylims()[2]*0.85, text("Median: " * string(round(median(row),   digits=2)), :black, :right, 10)),
            (xlims()[2]*0.95, ylims()[2]*0.80, text("Q75:    " * string(round(quantile(row, 0.75), digits=2)), :black, :right, 10)),
            (xlims()[2]*0.95, ylims()[2]*0.75, text("Max:    " * string(round(maximum(row),    digits=2)), :black, :right, 10)),
        ])
        push!(plts, p)
    end
    t_plot = plot(plts..., layout=(4, 3), size=(2000, 2000))
    savefig(t_plot, filepath)
end

# Export summary statistics for the recruitment-time distribution to CSV.
function export_time_statistics(times::Array, filepath::String)
    irt_approaches      = repeat(["F0a", "F0b", "F1a", "F1b"], inner=3)
    resupply_strategies = repeat(["Low", "Medium", "High"], outer=4)
    df = DataFrame(
        "IRT Approach"      => irt_approaches,
        "Resupply Strategy" => resupply_strategies,
        "Min"               => minimum.(eachrow(times)),
        "Q25"               => quantile.(eachrow(times), 0.25),
        "Median"            => median.(eachrow(times)),
        "Q75"               => quantile.(eachrow(times), 0.75),
        "Max"               => maximum.(eachrow(times)),
        "Mean"              => mean.(eachrow(times)),
        "SD"                => std.(eachrow(times)),
        "Q90"               => quantile.(eachrow(times), 0.9),
    )
    CSV.write(filepath, df)
end

# Export a chosen summary statistic (e.g. mean, std, quantile_90) for all
# simulation characteristics to CSV.
function export_scenario_statistics(chars::Array, labels::Vector{String}, filepath::String, stat::Function)
    irt_approaches      = ["F0a" "F0a" "F0a" "F0b" "F0b" "F0b" "F1a" "F1a" "F1a" "F1b" "F1b" "F1b"]
    resupply_strategies = ["Low" "Medium" "High" "Low" "Medium" "High" "Low" "Medium" "High" "Low" "Medium" "High"]
    cols = [reduce(hcat, stat(chars[:, k, :], dims=2)) for k in 1:8]
    mat  = vcat(irt_approaches, resupply_strategies, cols...)
    df   = DataFrame(permutedims(mat, (2, 1)), labels)
    CSV.write(filepath, df)
end

# Chi-squared Q-Q plot of squared Mahalanobis distances for joint bivariate normality.
# d500z1s/d500z2s: end-of-trial imbalance vectors (one value per simulation) for strata 1 and 2.
# Under bivariate normality, D² = (x - μ)ᵀ Σ⁻¹ (x - μ) ~ χ²(2).
# The histogram is density-normalised and overlaid with the χ²(2) PDF.
function plot_joint_normality_mahalanobis(
    d500z1s::AbstractVector,
    d500z2s::AbstractVector,
    filepath::String,
    title::String,
)
    n = length(d500z1s)
    X = hcat(Float64.(d500z1s), Float64.(d500z2s))  # n × 2

    μ     = vec(mean(X, dims=1))
    Σ     = cov(X)
    Σ_inv = inv(Σ)

    D2 = map(1:n) do i
        x = X[i, :] .- μ
        dot(x, Σ_inv * x)
    end

    chi2  = Chisq(2)
    x_max = max(quantile(chi2, 0.999), maximum(D2) * 1.02)
    xs    = range(0, x_max, length=400)

    p = histogram(D2;
        normalize  = :pdf,
        bins       = 40,
        label      = "Empirical D²",
        alpha      = 0.55,
        color      = cgrad(:blues)[0.65],
        xlabel     = "Mahalanobis Distance² (D²)",
        ylabel     = "Density",
        title      = title,
        size       = (800, 600),
        margin     = 10mm,
        legend     = :topright,
        titlefont  = font(13),
    )
    plot!(p, xs, pdf.(chi2, xs);
        label     = "χ²(2) density",
        linewidth = 2.5,
        color     = :red,
    )

    savefig(p, filepath)
end
