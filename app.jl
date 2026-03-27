### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ cell-imports
begin
    using PlutoUI
    using Plots
    using Statistics
    using LinearAlgebra
    using Printf
    using StatsBase
    using DataFrames
    using CSV
    using Distributions
    cd(@__DIR__)
    include("funcs.jl")
    include("Simulation.jl")
    include("plotting.jl")
    include("stratifiedMain.jl")
    mkpath("plots")
end

# ╔═╡ cell-defaults
_defaults = default_params()

# ╔═╡ cell-trial-design
md"""
## Trial Design
Sample size: $(@bind sample_size_val Slider(100:100:5000, default=_defaults.sample_size, show_value=true))

Centers: $(@bind centers_val Slider(1:200, default=_defaults.centers, show_value=true))

Block size: $(@bind block_size_val Select([2, 4, 6, 8], default=_defaults.block_size))

*Treatment arms and allocation ratio are fixed at 2 arms, 1:1.*
"""

# ╔═╡ cell-resupply
md"""
## Resupply Logistics
Resupply period (days): $(@bind resupply_period_val Slider(1:30, default=_defaults.resupply_period, show_value=true))

Resupply time (days): $(@bind resupply_time_val Slider(1:21, default=_defaults.resupply_time, show_value=true))

Kit cost: $(@bind kit_cost_val NumberField(0:1000000, default=_defaults.kit_cost))

Shipping cost: $(@bind ship_cost_val NumberField(0:100000, default=_defaults.ship_cost))
"""

# ╔═╡ cell-recruitment
md"""
## Recruitment Rate (Gamma(α, 1/β))
α: $(@bind alpha_val NumberField(0.1:0.1:10.0, default=_defaults.alpha))

β values (comma-separated): $(@bind beta_options_str TextField(default=join(Int.(_defaults.beta_options), ", ")))
"""

# ╔═╡ cell-strata-fr
md"""
## Stratification & Forced Randomization
Stratum 1 probability: $(@bind strata_prob_val Slider(0.01:0.01:0.99, default=_defaults.strata_assignment_probability, show_value=true))

Initial cap (fraction of sample size): $(@bind initial_cap_frac_val Slider(0.0:0.05:1.0, default=_defaults.initial_cap / _defaults.sample_size, show_value=true))
"""

# ╔═╡ cell-supply-strategies
md"""
## Supply Strategies

**Low supply**
Resupply amount: $(@bind low_resupply_val NumberField(1:50, default=_defaults.low_resupply))
Init arm 1: $(@bind low_init_1_val NumberField(1:50, default=_defaults.low_init[1]))
Init arm 2: $(@bind low_init_2_val NumberField(1:50, default=_defaults.low_init[2]))
Critical point: $(@bind low_critical_val NumberField(0:20, default=_defaults.low_critical))

**Medium supply**
Resupply amount: $(@bind med_resupply_val NumberField(1:50, default=_defaults.med_resupply))
Init arm 1: $(@bind med_init_1_val NumberField(1:50, default=_defaults.med_init[1]))
Init arm 2: $(@bind med_init_2_val NumberField(1:50, default=_defaults.med_init[2]))
Critical point: $(@bind med_critical_val NumberField(0:20, default=_defaults.med_critical))

**High supply**
Resupply amount: $(@bind high_resupply_val NumberField(1:50, default=_defaults.high_resupply))
Init arm 1: $(@bind high_init_1_val NumberField(1:50, default=_defaults.high_init[1]))
Init arm 2: $(@bind high_init_2_val NumberField(1:50, default=_defaults.high_init[2]))
Critical point: $(@bind high_critical_val NumberField(0:20, default=_defaults.high_critical))
"""

# ╔═╡ cell-nsims
md"""
## Simulation
Number of simulations: $(@bind number_simulations_val Slider(100:100:5000, default=_defaults.number_simulations, show_value=true))

$(number_simulations_val > 2000 ? md"⚠️ **Warning:** values above 2000 may take several minutes." : md"")
"""

# ╔═╡ cell-build-params
current_params = SimParams(
    sample_size                   = sample_size_val,
    treatment_arms                = _defaults.treatment_arms,
    allocation_ratio              = _defaults.allocation_ratio,
    block_size                    = block_size_val,
    centers                       = centers_val,
    resupply_period               = resupply_period_val,
    resupply_time                 = resupply_time_val,
    kit_cost                      = Float64(kit_cost_val),
    ship_cost                     = Float64(ship_cost_val),
    initial_cap                   = round(Int, initial_cap_frac_val * sample_size_val),
    alpha                         = alpha_val,
    beta_options                  = parse.(Float64, strip.(split(beta_options_str, ","))),
    strata_assignment_probability = strata_prob_val,
    number_simulations            = number_simulations_val,
    low_resupply                  = low_resupply_val,
    low_init                      = [low_init_1_val, low_init_2_val],
    low_critical                  = low_critical_val,
    med_resupply                  = med_resupply_val,
    med_init                      = [med_init_1_val, med_init_2_val],
    med_critical                  = med_critical_val,
    high_resupply                 = high_resupply_val,
    high_init                     = [high_init_1_val, high_init_2_val],
    high_critical                 = high_critical_val,
)

# ╔═╡ cell-run-button
@bind run_button Button("▶  Run Simulation")

# ╔═╡ cell-run-simulation
sim_results = let
    run_button
    run_simulation(current_params)
end

# ╔═╡ cell-imbalance-hists
let
    r = sim_results[1]
    plot_imbalance_histograms(
        r.dm1s[1:1, :, :], r.dm2s[1:1, :, :],
        r.max_z1 - 4, r.max_z2 - 4,
        joinpath("plots", "dm_hists_f1a.png"), "F1a Low"; num_panels=1
    )
    plot_imbalance_histograms(
        r.dm1s[2:2, :, :], r.dm2s[2:2, :, :],
        r.max_z1 - 4, r.max_z2 - 4,
        joinpath("plots", "dm_hists_f1b.png"), "F1b Low"
    )
    md"""
    ### Imbalance Histograms
    **F1a** $(PlutoUI.LocalResource(joinpath("plots", "dm_hists_f1a.png")))
    **F1b** $(PlutoUI.LocalResource(joinpath("plots", "dm_hists_f1b.png")))
    """
end

# ╔═╡ cell-joint-normality
let
    r = sim_results[1]
    plot_joint_normality_mahalanobis(
        r.d500z1s[1, :], r.d500z2s[1, :],
        joinpath("plots", "joint_normality_f1a.png"),
        "F1a Low Supply — Joint Normality of End-of-Trial Imbalances",
    )
    md"""
    ### Joint Normality (F1a)
    $(PlutoUI.LocalResource(joinpath("plots", "joint_normality_f1a.png")))
    """
end

# ╔═╡ cell-varcov
let
    r = sim_results[1]
    lines = String[]
    for (label, idx) in [("F1a Low", 1), ("F1b Low", 2)]
        A = cov(hcat(r.d500z1s[idx, :], r.d500z2s[idx, :]))
        evs = eigvals(A)
        push!(lines, "**$label** — Cov: $(round.(A, digits=4)), Eigenvalues: $(round.(evs, digits=4))")
    end
    Markdown.parse("### Variance-Covariance\n\n" * join(lines, "\n\n"))
end

# ╔═╡ cell-fa-rates
let
    r = sim_results[1]
    rows = [
        ("F1a low, z=1",  mean(r.characteristics[1, 1, 1, :])),
        ("F1a low, z=2",  mean(r.characteristics[1, 1, 2, :])),
        ("F1a low, total",mean(r.characteristics[1, 1, 3, :])),
        ("F1b low, z=1",  mean(r.characteristics[2, 1, 1, :])),
        ("F1b low, z=2",  mean(r.characteristics[2, 1, 2, :])),
        ("F1b low, total",mean(r.characteristics[2, 1, 3, :])),
    ]
    table = join(["| $(row[1]) | $(round(row[2], digits=4)) |" for row in rows], "\n")
    Markdown.parse("### FA Rate Summary\n\n| Scenario | Mean FA Rate |\n|---|---|\n$table")
end

# ╔═╡ cell-recruitment-times
let
    r = sim_results[1]
    p1 = histogram(r.recruitment_times[1, :]; bins=30, title="F1a Low — Recruitment Time",
                   xlabel="Days", ylabel="Count", legend=false)
    p2 = histogram(r.recruitment_times[2, :]; bins=30, title="F1b Low — Recruitment Time",
                   xlabel="Days", ylabel="Count", legend=false)
    combined = plot(p1, p2, layout=(1, 2), size=(1200, 500))
    savefig(combined, joinpath("plots", "recruitment_times.png"))
    md"""
    ### Recruitment Times
    $(PlutoUI.LocalResource(joinpath("plots", "recruitment_times.png")))
    """
end

# ╔═╡ cell-csv-exports
let
    r = sim_results[1]

    time_path = joinpath("plots", "time_stats.csv")
    time_df = DataFrame(
        Scenario = ["F1a Low", "F1b Low"],
        Min      = minimum.(eachrow(r.recruitment_times)),
        Q25      = quantile.(eachrow(r.recruitment_times), 0.25),
        Median   = median.(eachrow(r.recruitment_times)),
        Q75      = quantile.(eachrow(r.recruitment_times), 0.75),
        Max      = maximum.(eachrow(r.recruitment_times)),
        Mean     = mean.(eachrow(r.recruitment_times)),
        SD       = std.(eachrow(r.recruitment_times)),
    )
    CSV.write(time_path, time_df)

    fa_path = joinpath("plots", "fa_rates.csv")
    fa_df = DataFrame(
        Scenario = ["F1a Low", "F1a Low", "F1a Low", "F1b Low", "F1b Low", "F1b Low"],
        Stratum  = ["z=1", "z=2", "Total", "z=1", "z=2", "Total"],
        Mean_FA  = [
            mean(r.characteristics[1, 1, 1, :]),
            mean(r.characteristics[1, 1, 2, :]),
            mean(r.characteristics[1, 1, 3, :]),
            mean(r.characteristics[2, 1, 1, :]),
            mean(r.characteristics[2, 1, 2, :]),
            mean(r.characteristics[2, 1, 3, :]),
        ],
    )
    CSV.write(fa_path, fa_df)

    md"""
    ### CSV Exports
    $(DownloadButton(read(time_path), "time_stats.csv"))
    $(DownloadButton(read(fa_path), "fa_rates.csv"))
    """
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
"""
