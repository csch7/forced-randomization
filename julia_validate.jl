# Run from the project root in a Julia REPL:
#   PlutoRunner = true   # skip the auto-run block in stratifiedMain.jl
#   include("julia_validate.jl")

using Random, Statistics

PlutoRunner = true
include("stratifiedMain.jl")

p = default_params()
p2 = SimParams(
    sample_size          = 500,
    treatment_arms       = p.treatment_arms,
    allocation_ratio     = p.allocation_ratio,
    block_size           = p.block_size,
    centers              = p.centers,
    resupply_period      = p.resupply_period,
    resupply_time        = p.resupply_time,
    kit_cost             = p.kit_cost,
    ship_cost            = p.ship_cost,
    initial_cap          = 250,
    alpha                = p.alpha,
    beta_options         = Float64[3.0],
    num_strata           = p.num_strata,
    strata_probabilities = p.strata_probabilities,
    number_simulations   = 200,
    run_label            = "",
    low_resupply         = p.low_resupply,
    low_init             = copy(p.low_init),
    low_critical         = p.low_critical,
    med_resupply         = p.med_resupply,
    med_init             = copy(p.med_init),
    med_critical         = p.med_critical,
    high_resupply        = p.high_resupply,
    high_init            = copy(p.high_init),
    high_critical        = p.high_critical,
)

Random.seed!(42)
results = run_simulation(p2)
r  = results[1]
ns = p2.num_strata

println("=== F1a Low (scenario 7) ===")
for k in 1:ns
    println("  Mean FA rate  (z=$k): ", round(mean(r.characteristics[1, 1, k, :]), digits=4))
end
println("  Mean FA rate  (total): ", round(mean(r.characteristics[1, 1, ns+1, :]), digits=4))
for k in 1:ns
    dm = r.d500s[1, :, k]
    println("  Mean end imbalance (z=$k): ", round(mean(dm), digits=4),
            "  |  Var: ", round(var(dm), digits=4))
end
for k in 1:ns
    println("  Avg slots skipped/FA (z=$k): ",
            round(mean(r.avg_slots_skipped_f1a[k, :]), digits=4))
end

println("\n=== F1b Low (scenario 10) ===")
for k in 1:ns
    println("  Mean FA rate  (z=$k): ", round(mean(r.characteristics[2, 1, k, :]), digits=4))
end
println("  Mean FA rate  (total): ", round(mean(r.characteristics[2, 1, ns+1, :]), digits=4))
for k in 1:ns
    dm = r.d500s[2, :, k]
    println("  Mean end imbalance (z=$k): ", round(mean(dm), digits=4),
            "  |  Var: ", round(var(dm), digits=4))
end

println("\n=== Recruitment times ===")
println("  F1a: mean=", round(mean(r.recruitment_times[1, :]), digits=1),
        "  F1b: mean=", round(mean(r.recruitment_times[2, :]), digits=1))
