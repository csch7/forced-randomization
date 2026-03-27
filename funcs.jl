using Random
using Statistics
using StatsBase

include("constants.jl")

# Generate inter-arrival times for each center using an exponential process,
# then return all patients sorted by arrival time.
# Returns an N×3 matrix: [patient_id, arrival_time, center_id]
function generate_patient_arrivals(
    rates::Vector{Float64},
    center_activations::AbstractVector{<:Integer},
    num_centers::Int,
    num_patients::Int
)
    all_arrivals = Matrix{Float64}(undef, 0, 3)
    for i in 1:num_centers
        # Inverse-transform sampling: u ~ Uniform(0,1) → Exponential inter-arrivals
        u = rand(num_patients)
        inter_arrivals = -log.(1 .- u) ./ rates[i]
        arrival_times  = cumsum(inter_arrivals) .+ center_activations[i] * 30
        block = hcat(
            Float64.((i-1)*num_patients .+ (1:num_patients)),
            arrival_times,
            fill(Float64(i), num_patients)
        )
        all_arrivals = vcat(all_arrivals, block)
    end
    return all_arrivals[sortperm(all_arrivals[:, 2]), :]
end

# Assign each patient to a stratum: stratum 1 with probability p, stratum 2 otherwise.
function generate_strata_assignments(num_patients::Int, p::Float64)
    return [u ≤ p ? Int8(1) : Int8(2) for u in rand(num_patients)]
end

# Generate a permuted block design treatment sequence.
# Returns a 1×N matrix. For stratified use, the caller reshapes and transposes
# to get a (num_strata × positions_per_stratum) matrix.
function generate_treatment_blocks(
    ratio::Tuple,
    sample_size::Int,
    treat_arms::Int,
    block_size::Int
)
    block_template = [t for t in 1:treat_arms for _ in 1:(ratio[t] * block_size ÷ sum(ratio))]
    n_blocks = ceil(Int, 2 * sample_size / block_size)
    return reduce(hcat, [reshape(shuffle(block_template), 1, :) for _ in 1:n_blocks])
end

# Return supply and FR parameters for a given scenario (1–12).
#
# Scenario layout:
#   Rows (IRT spec):     F0a (1–3), F0b (4–6), F1a (7–9), F1b (10–12)
#   Columns (supply):    Low, Medium, High
#
# Returns: (resupply_amount, init_supply, critical_point, fr_allowed, backfill_enabled, cap)
function scenario_params(scenario::Int, cap::Int)
    supply_configs = [
        (LOW_RESUPPLY,  copy(LOW_INIT),  LOW_CRITICAL),
        (MED_RESUPPLY,  copy(MED_INIT),  MED_CRITICAL),
        (HIGH_RESUPPLY, copy(HIGH_INIT), HIGH_CRITICAL),
    ]
    fr_configs = [
        (false, false, 0),    # F0a: no forced randomization; send home if supply missing
        (true,  false, 0),    # F0b: FR enabled but cap = 0 (never triggers)
        (true,  false, cap),  # F1a: FR enabled, skip forward without backfilling
        (true,  true,  cap),  # F1b: FR enabled, skip forward with backfilling
    ]
    supply_idx = ((scenario - 1) % 3) + 1
    fr_idx     = ((scenario - 1) ÷ 3) + 1
    resupply, init_supply, crit = supply_configs[supply_idx]
    fr_allowed, backfill, eff_cap = fr_configs[fr_idx]
    return resupply, init_supply, crit, fr_allowed, backfill, eff_cap
end

# 90th quantile approximation: mean + 1.645 * std, computed along a given dimension.
function quantile_90(values; dims::Int)
    return mean(values, dims=dims) .+ 1.645 .* std(values, dims=dims)
end
