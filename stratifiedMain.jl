using DataFrames
using Distributions
using Random
using Printf
using Statistics
using StatsBase
using ProgressBars
using LinearAlgebra
include("funcs.jl")
include("plotting.jl")
include("Simulation.jl")

mkpath("plots")

const HEADINGS = [
    "IRT Approach", "Resupply Strategy",
    "Treatment_Imbalance", "Pct_FAs", "Pct_Patients_Sent_Home",
    "Drug_Overage", "Cost", "Time", "Patients_NA", "Patients_Waitlisted",
]

# Run the stratified simulation for all beta values in params.beta_options.
# Returns a Vector of NamedTuples, one per beta, containing all result arrays.
function run_simulation(params::SimParams)
    all_results = []

    for BETA in params.beta_options

        max_z1 = Int(params.strata_assignment_probability * params.sample_size)
        max_z2 = params.sample_size - max_z1

        # characteristics[scenario_offset, characteristic, stratum_or_total, simulation]
        characteristics   = zeros(Float64, 2, 1, 3, params.number_simulations)

        dm1s              = fill(NaN, 2, params.number_simulations, max_z1)
        dm2s              = fill(NaN, 2, params.number_simulations, max_z2)

        dlgz1s            = fill(NaN, 2, params.number_simulations, max_z1)
        dlgz2s            = fill(NaN, 2, params.number_simulations, max_z2)

        d500z1s           = zeros(Float64, 2, params.number_simulations)
        d500z2s           = zeros(Float64, 2, params.number_simulations)

        unfilled_slots_z1 = zeros(Float64, 1, params.number_simulations, max_z1)
        unfilled_slots_z2 = zeros(Float64, 1, params.number_simulations, max_z2)

        # recruitment_times[scenario_offset, simulation] — arrival time of last patient
        recruitment_times = zeros(Float64, 2, params.number_simulations)

        for sim in tqdm(1:params.number_simulations)

            center_rates = rand(Gamma(params.alpha, 1/BETA), params.centers)
            center_acts  = rand(0:4, params.centers)
            patients     = generate_patient_arrivals(center_rates, center_acts, params.centers, params.sample_size)

            raw_blocks       = generate_treatment_blocks(params.allocation_ratio, params.sample_size, params.treatment_arms, params.block_size)
            treatment_blocks = Matrix(transpose(reshape(raw_blocks, (params.sample_size, 2))))

            strata = generate_strata_assignments(params.sample_size * 2, params.strata_assignment_probability)

            for (si, scenario) in enumerate([7, 10])

                resupply_amt, init_supply, critical_pt, fr_allowed, backfill_enabled, cap =
                    scenario_params(scenario, params)

                center_supplies = Dict{Int, Vector{Int}}()
                for i in 1:params.centers
                    center_supplies[i] = copy(init_supply)
                end

                total_drugs = params.centers * sum(init_supply)
                total_cost  = params.centers * sum(init_supply) * params.kit_cost

                S = Simulation(
                    [Int16[] for _ in 1:2],
                    zeros(Int16, 2),
                    Set{Int16}(),
                    [Int16[] for _ in 1:2],
                    zeros(Int16, 2),
                    zeros(Int16, 2),
                    [Int16[] for _ in 1:2],
                    params.sample_size,
                    center_supplies,
                    critical_pt,
                    cap,
                    treatment_blocks,
                )

                next_supply_check  = params.resupply_period
                next_resupply      = next_supply_check + params.resupply_time
                sent_supply        = Dict{Int, Vector{Int}}()
                tot_delayed        = zeros(Int16, 2)
                num_waitlisted     = zeros(Int16, 2)
                patients_per_stratum = zeros(Int16, 2)
                break_loop         = false

                i = 1
                while i <= S.num_patients
                    center = Int16(patients[i, 3])

                    if patients[i, 2] >= next_supply_check
                        for j in S.need_supply
                            new_supply = zeros(Int, params.treatment_arms)
                            for k in eachindex(S.supplies[j])
                                if S.supplies[j][k] <= critical_pt
                                    top_up         = resupply_amt - S.supplies[j][k]
                                    total_cost    += top_up * params.kit_cost
                                    total_drugs   += top_up
                                    new_supply[k]  = resupply_amt
                                end
                            end
                            sent_supply[j] = new_supply
                            total_cost    += params.ship_cost
                        end
                        next_supply_check += params.resupply_period
                        empty!(S.need_supply)
                    end

                    if patients[i, 2] >= next_resupply
                        for (cent, new_supply) in sent_supply
                            for k in eachindex(new_supply)
                                if new_supply[k] != 0
                                    S.supplies[cent][k] = new_supply[k]
                                end
                            end
                        end

                        num_delayed   = [length(S.delayed_patients[1]), length(S.delayed_patients[2])]
                        tot_delayed  .+= num_delayed

                        for z in 1:2
                            cts = countmap(S.delayed_patients[z])
                            for (_, num) in cts
                                num_waitlisted[z] += min(num, resupply_amt * params.treatment_arms)
                            end
                        end

                        flattened_delayed = vcat(S.delayed_patients[1], S.delayed_patients[2])

                        for j in eachindex(flattened_delayed)
                            new_delayed   = length(S.delayed_patients[1]) + length(S.delayed_patients[2]) - sum(num_delayed)
                            patient_index = i - sum(S.patients_sent_home) - sum(tot_delayed) - new_delayed
                            stratum       = strata[patient_index]
                            center_delayed = Int16(flattened_delayed[j])
                            stratum_index = trunc(Int16,
                                (i - patients_per_stratum[abs(stratum - 3)]) -
                                S.patients_sent_home[stratum] -
                                tot_delayed[stratum] -
                                (length(S.delayed_patients[stratum]) - num_delayed[stratum])
                            )

                            if !fr_allowed
                                allocate_f0a!(S, center_delayed, stratum, stratum_index)
                            elseif cap == 0
                                allocate_f0b!(S, center_delayed, stratum, stratum_index)
                            elseif !backfill_enabled
                                fr_allowed = allocate_f1a!(S, center_delayed, stratum, stratum_index)
                                if stratum == 1
                                    unfilled_slots_z1[si, sim, min(max_z1, stratum_index)] = S.treatments_skipped[1]
                                else
                                    unfilled_slots_z2[si, sim, min(max_z2, stratum_index)] = S.treatments_skipped[2]
                                end
                            else
                                fr_allowed = allocate_f1b!(S, center_delayed, stratum, stratum_index)
                            end

                            patients_per_stratum[stratum] += 1

                            i_1 = trunc(Int16, (i - patients_per_stratum[2]) - tot_delayed[1] - length(S.delayed_patients[1]))
                            i_2 = trunc(Int16, (i - patients_per_stratum[1]) - tot_delayed[2] - length(S.delayed_patients[2]))
                            dlgz1s[si, sim, min(max_z1, stratum_index)] = patients_per_stratum[1] - i_1
                            dlgz2s[si, sim, min(max_z2, stratum_index)] = patients_per_stratum[2] - i_2

                            i += 1
                            if i > S.num_patients
                                patients[S.num_patients, 2] = next_resupply
                                break_loop = true
                                break
                            end
                        end

                        S.delayed_patients[1] = S.delayed_patients[1][num_delayed[1]+1:end]
                        S.delayed_patients[2] = S.delayed_patients[2][num_delayed[2]+1:end]
                        next_resupply = next_supply_check + params.resupply_time
                        empty!(sent_supply)
                    end

                    if break_loop
                        break
                    end

                    current_delayed = length(S.delayed_patients[1]) + length(S.delayed_patients[2])
                    patient_index   = i + sum(S.treatments_skipped) - sum(S.patients_sent_home) - sum(tot_delayed) - current_delayed
                    stratum         = strata[patient_index]
                    stratum_index   = trunc(Int16,
                        (i - patients_per_stratum[abs(stratum - 3)]) -
                        tot_delayed[stratum] -
                        length(S.delayed_patients[stratum])
                    )

                    if !fr_allowed
                        allocate_f0a!(S, center, stratum, stratum_index)
                    elseif cap == 0
                        allocate_f0b!(S, center, stratum, stratum_index)
                    elseif !backfill_enabled
                        fr_allowed = allocate_f1a!(S, center, stratum, stratum_index)
                        if stratum == 1
                            unfilled_slots_z1[si, sim, min(max_z1, stratum_index)] = S.treatments_skipped[1]
                        else
                            unfilled_slots_z2[si, sim, min(max_z2, stratum_index)] = S.treatments_skipped[2]
                        end
                    else
                        fr_allowed = allocate_f1b!(S, center, stratum, stratum_index)
                    end

                    patients_per_stratum[stratum] += 1

                    i_1 = trunc(Int16, i - patients_per_stratum[2])
                    i_2 = trunc(Int16, i - patients_per_stratum[1])
                    dlgz1s[si, sim, min(max_z1, stratum_index)] = patients_per_stratum[1] - i_1
                    dlgz2s[si, sim, min(max_z2, stratum_index)] = patients_per_stratum[2] - i_2

                    i += 1
                end  # patient loop

                # Record recruitment time as arrival time of last original patient
                recruitment_times[si, sim] = patients[params.sample_size, 2]

                normalise = (scenario == 7)

                for c in 1:min(max_z1, length(S.treatments_used[1]))
                    cm1 = countmap(S.treatments_used[1][1:c])
                    t1  = get(cm1, 1, 0)
                    t2  = get(cm1, 2, 0)
                    dm1s[si, sim, c] = normalise ? (t1 - t2) / sqrt(c) : (t1 - t2)
                end
                for c in 1:min(max_z2, length(S.treatments_used[2]))
                    cm2 = countmap(S.treatments_used[2][1:c])
                    t1  = get(cm2, 1, 0)
                    t2  = get(cm2, 2, 0)
                    dm2s[si, sim, c] = normalise ? (t1 - t2) / sqrt(c) : (t1 - t2)
                end

                cm1_full = countmap(S.treatments_used[1])
                n1 = length(S.treatments_used[1])
                d1 = get(cm1_full, 1, 0) - get(cm1_full, 2, 0)
                d500z1s[si, sim] = normalise ? d1 / sqrt(n1) : d1

                cm2_full = countmap(S.treatments_used[2])
                n2 = length(S.treatments_used[2])
                d2 = get(cm2_full, 1, 0) - get(cm2_full, 2, 0)
                d500z2s[si, sim] = normalise ? d2 / sqrt(n2) : d2

                characteristics[si, 1, 1, sim] = S.patients_force_allocated[1] / length(S.treatments_used[1])
                characteristics[si, 1, 2, sim] = S.patients_force_allocated[2] / length(S.treatments_used[2])
                characteristics[si, 1, 3, sim] = sum(S.patients_force_allocated) / params.sample_size

            end  # scenario
        end  # sim

        push!(all_results, (
            beta              = BETA,
            dm1s              = dm1s,
            dm2s              = dm2s,
            d500z1s           = d500z1s,
            d500z2s           = d500z2s,
            characteristics   = characteristics,
            dlgz1s            = dlgz1s,
            dlgz2s            = dlgz2s,
            unfilled_slots_z1 = unfilled_slots_z1,
            unfilled_slots_z2 = unfilled_slots_z2,
            recruitment_times = recruitment_times,
            max_z1            = max_z1,
            max_z2            = max_z2,
        ))

    end  # BETA

    return all_results
end  # run_simulation


# ── Thin wrapper: preserves include("stratifiedMain.jl") workflow ────────────
# Skip auto-run when loaded from Pluto (notebook manages execution via run button).
if !isdefined(Main, :PlutoRunner)
params  = default_params()
results = run_simulation(params)

r = results[1]  # first (only) beta

f1a_panels = [
    imbalance_line_panel(r.dm1s[1:1, :, :],   "Var[dm(1)]",    var;        cutoff=r.max_z1, start=4),
    imbalance_line_panel(r.dm2s[1:1, :, :],   "Var[dm(2)]",    var;        cutoff=r.max_z2, start=4),
    imbalance_line_panel(r.dm1s[1:1, :, :],   "Q90[dm(1)]",    quantile_90; cutoff=r.max_z1, start=4),
    imbalance_line_panel(r.dm2s[1:1, :, :],   "Q90[dm(2)]",    quantile_90; cutoff=r.max_z2, start=4),
]
save_line_summary(f1a_panels, 2, 2, joinpath("plots", "f1a_low.png"), "F1a Low Supply")

f1b_panels = [
    imbalance_line_panel(r.dm1s[2:2, :, :],   "Var[Dm(1)]",    var;        cutoff=r.max_z1, start=4),
    imbalance_line_panel(r.dm2s[2:2, :, :],   "Var[Dm(2)]",    var;        cutoff=r.max_z2, start=4),
    imbalance_line_panel(r.dm1s[2:2, :, :],   "Q90[Dm(1)]",    quantile_90; cutoff=r.max_z1, start=4),
    imbalance_line_panel(r.dm2s[2:2, :, :],   "Q90[Dm(2)]",    quantile_90; cutoff=r.max_z2, start=4),
    imbalance_line_panel(r.dlgz1s[2:2, :, :], "Q90[DLG(z=1)]", quantile_90; cutoff=r.max_z1),
    imbalance_line_panel(r.dlgz2s[2:2, :, :], "Q90[DLG(z=2)]", quantile_90; cutoff=r.max_z2),
    imbalance_line_panel(r.dlgz1s[2:2, :, :], "max[DLG(z=1)]", maximum;    cutoff=r.max_z1),
    imbalance_line_panel(r.dlgz2s[2:2, :, :], "max[DLG(z=2)]", maximum;    cutoff=r.max_z2),
]
save_line_summary(f1b_panels, 4, 2, joinpath("plots", "f1b_low.png"), "F1b Low Supply")

plot_imbalance_histograms(r.dm1s[1:1, :, :], r.dm2s[1:1, :, :], r.max_z1 - 4, r.max_z2 - 4,
    joinpath("plots", "dm_hists_f1a.png"), "F1a Low"; num_panels=1, unnormalize=true)
plot_imbalance_histograms(r.dm1s[2:2, :, :], r.dm2s[2:2, :, :], r.max_z1 - 4, r.max_z2 - 4,
    joinpath("plots", "dm_hists_f1b.png"), "F1b Low")

plot_joint_normality_mahalanobis(
    r.d500z1s[1, :], r.d500z2s[1, :],
    joinpath("plots", "joint_normality_f1a.png"),
    "F1a Low Supply — Joint Normality of End-of-Trial Imbalances",
)

for (label, idx) in [("F1a Low", 1), ("F1b Low", 2)]
    A = cov(hcat(r.d500z1s[idx, :], r.d500z2s[idx, :]))
    println("Var-Cov $label: ", A)
    println("\tEigenvalues: ", eigvals(A))
end

@printf("Average FA (F1a low, z=1): %.4f\n", mean(r.characteristics[1, 1, 1, :]))
@printf("Average FA (F1a low, z=2): %.4f\n", mean(r.characteristics[1, 1, 2, :]))
@printf("Average FA (F1a low, tot): %.4f\n", mean(r.characteristics[1, 1, 3, :]))
@printf("Average FA (F1b low, z=1): %.4f\n", mean(r.characteristics[2, 1, 1, :]))
@printf("Average FA (F1b low, z=2): %.4f\n", mean(r.characteristics[2, 1, 2, :]))
@printf("Average FA (F1b low, tot): %.4f\n", mean(r.characteristics[2, 1, 3, :]))
end  # !isdefined PlutoRunner
