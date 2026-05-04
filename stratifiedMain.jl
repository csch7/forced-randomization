using DataFrames
using Dates
using DelimitedFiles
using Distributions
using Random
using Printf
using Statistics
using StatsBase
using LinearAlgebra
using ProgressMeter
include("funcs.jl")
include("plotting.jl")
include("Simulation.jl")

const SCENARIOS_TO_RUN = [7, 10]

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

        ns  = params.num_strata
        nsc = length(SCENARIOS_TO_RUN)

        max_z      = [Int(round(params.strata_probabilities[k] * params.sample_size)) for k in 1:ns]
        max_z[end] = params.sample_size - sum(max_z[1:end-1])

        # characteristics[scenario_offset, characteristic, stratum_or_total, simulation]
        characteristics       = zeros(Float64, nsc, 1, ns + 1, params.number_simulations)

        # dms[stratum]: (nsc, simulations, patients_in_stratum)
        dms                   = [fill(NaN, nsc, params.number_simulations, max_z[k]) for k in 1:ns]
        dlgs                  = [fill(NaN, nsc, params.number_simulations, max_z[k]) for k in 1:ns]

        # d500s[scenario, simulation, stratum]
        d500s                 = zeros(Float64, nsc, params.number_simulations, ns)

        # unfilled_slots[stratum]: (1, simulations, patients_in_stratum) — F1a only
        unfilled_slots        = [zeros(Float64, 1, params.number_simulations, max_z[k]) for k in 1:ns]

        # recruitment_times[scenario_offset, simulation] — arrival time of last patient
        recruitment_times     = zeros(Float64, nsc, params.number_simulations)

        # avg_slots_skipped_f1a[stratum, simulation] — mean block positions skipped per FR event
        avg_slots_skipped_f1a = zeros(Float64, ns, params.number_simulations)

        p = Progress(params.number_simulations)
        Threads.@threads for sim in 1:params.number_simulations

            center_rates = rand(Gamma(params.alpha, 1/BETA), params.centers)
            center_acts  = rand(0:4, params.centers)
            patients     = generate_patient_arrivals(center_rates, center_acts, params.centers, params.sample_size)

            raw_blocks       = generate_treatment_blocks(params.allocation_ratio, params.sample_size, params.treatment_arms, params.block_size, ns)
            treatment_blocks = Matrix(transpose(reshape(raw_blocks, (params.sample_size, ns))))

            strata = generate_strata_assignments(params.sample_size * 2, params.strata_probabilities)

            for (si, scenario) in enumerate(SCENARIOS_TO_RUN)

                resupply_amt, init_supply, critical_pt, fr_allowed, backfill_enabled, cap =
                    scenario_params(scenario, params)

                center_supplies = Dict{Int, Vector{Int}}()
                for i in 1:params.centers
                    center_supplies[i] = copy(init_supply)
                end

                total_drugs = params.centers * sum(init_supply)
                total_cost  = params.centers * sum(init_supply) * params.kit_cost

                S = Simulation(
                    [Int16[] for _ in 1:ns],
                    zeros(Int, ns),
                    Set{Int16}(),
                    [Int16[] for _ in 1:ns],
                    zeros(Int, ns),
                    zeros(Int, ns),
                    [Set{Int16}() for _ in 1:ns],
                    params.sample_size,
                    center_supplies,
                    critical_pt,
                    cap,
                    treatment_blocks,
                )

                next_supply_check    = params.resupply_period
                next_resupply        = next_supply_check + params.resupply_time
                sent_supply          = Dict{Int, Vector{Int}}()
                tot_delayed          = zeros(Int, ns)
                num_waitlisted       = zeros(Int, ns)
                patients_per_stratum = zeros(Int, ns)
                break_loop           = false

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

                        num_delayed   = [length(S.delayed_patients[k]) for k in 1:ns]
                        tot_delayed  .+= num_delayed

                        for z in 1:ns
                            cts = countmap(S.delayed_patients[z])
                            for (_, num) in cts
                                num_waitlisted[z] += min(num, resupply_amt * params.treatment_arms)
                            end
                        end

                        flattened_delayed = vcat([S.delayed_patients[z] for z in 1:ns]...)

                        for j in eachindex(flattened_delayed)
                            new_delayed   = sum(length(S.delayed_patients[z]) for z in 1:ns) - sum(num_delayed)
                            patient_index = i - sum(S.patients_sent_home) - sum(tot_delayed) - new_delayed
                            stratum       = strata[patient_index]
                            center_delayed = Int16(flattened_delayed[j])
                            other_strata_patients = sum(patients_per_stratum) - patients_per_stratum[stratum]
                            stratum_index = trunc(Int16,
                                (i - other_strata_patients) -
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
                                unfilled_slots[stratum][1, sim, min(max_z[stratum], stratum_index)] = S.treatments_skipped[stratum]
                            else
                                fr_allowed = allocate_f1b!(S, center_delayed, stratum, stratum_index)
                            end

                            patients_per_stratum[stratum] += 1

                            dlgs[stratum][si, sim, min(max_z[stratum], stratum_index)] = compute_dlg(S, stratum, stratum_index)

                            i += 1
                            if i > S.num_patients
                                patients[S.num_patients, 2] = next_resupply
                                break_loop = true
                                break
                            end
                        end

                        for z in 1:ns
                            S.delayed_patients[z] = S.delayed_patients[z][num_delayed[z]+1:end]
                        end
                        next_resupply = next_supply_check + params.resupply_time
                        empty!(sent_supply)
                    end

                    if break_loop
                        break
                    end

                    current_delayed = sum(length(S.delayed_patients[z]) for z in 1:ns)
                    patient_index   = i + sum(S.treatments_skipped) - sum(S.patients_sent_home) - sum(tot_delayed) - current_delayed
                    stratum         = strata[patient_index]
                    other_strata_patients = sum(patients_per_stratum) - patients_per_stratum[stratum]
                    stratum_index   = trunc(Int16,
                        (i - other_strata_patients) -
                        tot_delayed[stratum] -
                        length(S.delayed_patients[stratum])
                    )

                    if !fr_allowed
                        allocate_f0a!(S, center, stratum, stratum_index)
                    elseif cap == 0
                        allocate_f0b!(S, center, stratum, stratum_index)
                    elseif !backfill_enabled
                        fr_allowed = allocate_f1a!(S, center, stratum, stratum_index)
                        unfilled_slots[stratum][1, sim, min(max_z[stratum], stratum_index)] = S.treatments_skipped[stratum]
                    else
                        fr_allowed = allocate_f1b!(S, center, stratum, stratum_index)
                    end

                    patients_per_stratum[stratum] += 1

                    dlgs[stratum][si, sim, min(max_z[stratum], stratum_index)] = compute_dlg(S, stratum, stratum_index)

                    i += 1
                end  # patient loop

                # Record recruitment time as arrival time of last original patient
                recruitment_times[si, sim] = patients[params.sample_size, 2]

                normalise = (scenario == 7)

                for k in 1:ns
                    t1 = 0; t2 = 0
                    nk = length(S.treatments_used[k])
                    for c in 1:nk
                        if S.treatments_used[k][c] == 1; t1 += 1 else t2 += 1 end
                        if c <= max_z[k]
                            dms[k][si, sim, c] = normalise ? (t1 - t2) / sqrt(c) : Float64(t1 - t2)
                        end
                    end
                    d500s[si, sim, k] = normalise ? Float64(t1 - t2) / sqrt(nk) : Float64(t1 - t2)
                end

                for k in 1:ns
                    characteristics[si, 1, k, sim] = S.patients_force_allocated[k] / length(S.treatments_used[k])
                end
                characteristics[si, 1, ns + 1, sim] = sum(S.patients_force_allocated) / params.sample_size

                if si == 1  # F1a only
                    for z in 1:ns
                        pfa = S.patients_force_allocated[z]
                        avg_slots_skipped_f1a[z, sim] = pfa > 0 ? S.treatments_skipped[z] / pfa : 0.0
                    end
                end

            end  # scenario
            next!(p)
        end  # sim

        push!(all_results, (
            beta                  = BETA,
            dms                   = dms,
            d500s                 = d500s,
            characteristics       = characteristics,
            dlgs                  = dlgs,
            unfilled_slots        = unfilled_slots,
            recruitment_times     = recruitment_times,
            avg_slots_skipped_f1a = avg_slots_skipped_f1a,
            max_z                 = max_z,
        ))

    end  # BETA

    return all_results
end  # run_simulation


# ── Thin wrapper: preserves include("stratifiedMain.jl") workflow ────────────
# Skip auto-run when loaded from Pluto (notebook manages execution via run button).
if !isdefined(Main, :PlutoRunner)
params  = default_params()
results = run_simulation(params)

r  = results[1]  # first (only) beta
ns = params.num_strata

plot_dir = isempty(params.run_label) ? "plots" : joinpath("plots", params.run_label)
mkpath(plot_dir)

f1a_panels = []
for k in 1:ns
    push!(f1a_panels, imbalance_line_panel(r.dms[k][1:1, :, :], "Var[dm($k)]",    var;         cutoff=r.max_z[k], start=4))
    push!(f1a_panels, imbalance_line_panel(r.dms[k][1:1, :, :], "Q90[dm($k)]",    quantile_90; cutoff=r.max_z[k], start=4))
end
save_line_summary(f1a_panels, ns, 2, joinpath(plot_dir, "f1a_low.png"), "F1a Low Supply")

f1b_panels = []
for k in 1:ns
    push!(f1b_panels, imbalance_line_panel(r.dms[k][2:2, :, :],  "Var[Dm($k)]",     var;         cutoff=r.max_z[k], start=4))
    push!(f1b_panels, imbalance_line_panel(r.dms[k][2:2, :, :],  "Q90[Dm($k)]",     quantile_90; cutoff=r.max_z[k], start=4))
    push!(f1b_panels, imbalance_line_panel(r.dlgs[k][2:2, :, :], "Q90[DLG(z=$k)]",  quantile_90; cutoff=r.max_z[k]))
    push!(f1b_panels, imbalance_line_panel(r.dlgs[k][2:2, :, :], "max[DLG(z=$k)]",  maximum;     cutoff=r.max_z[k]))
end
save_line_summary(f1b_panels, ns, 4, joinpath(plot_dir, "f1b_low.png"), "F1b Low Supply")

plot_imbalance_histograms(
    [r.dms[k][1:1, :, :] for k in 1:ns],
    [r.max_z[k] - 4 for k in 1:ns],
    joinpath(plot_dir, "dm_hists_f1a.png"), "F1a Low"; num_panels=1, unnormalize=true)
plot_imbalance_histograms(
    [r.dms[k][2:2, :, :] for k in 1:ns],
    [r.max_z[k] - 4 for k in 1:ns],
    joinpath(plot_dir, "dm_hists_f1b.png"), "F1b Low")

plot_joint_normality_mahalanobis(
    [r.d500s[1, :, k] for k in 1:ns],
    joinpath(plot_dir, "joint_normality_f1a.png"),
    "F1a Low Supply — Joint Normality of End-of-Trial Imbalances",
)

open(joinpath(plot_dir, "params.txt"), "w") do file
    println(file, "Last run: ", Dates.now(), " US Eastern")
    println(file, "")
    println(file, "--- Scenarios Run ---")
    for s in SCENARIOS_TO_RUN
        println(file, "  Scenario $s: ", SCENARIO_MAP[s])
    end
    println(file, "")
    println(file, "--- Trial Design ---")
    println(file, "  SAMPLE_SIZE                   = ", SAMPLE_SIZE)
    println(file, "  TREATMENT_ARMS                = ", TREATMENT_ARMS)
    println(file, "  ALLOCATION_RATIO              = ", ALLOCATION_RATIO)
    println(file, "  BLOCK_SIZE                    = ", BLOCK_SIZE)
    println(file, "  CENTERS                       = ", CENTERS)
    println(file, "")
    println(file, "--- Resupply Logistics ---")
    println(file, "  RESUPPLY_PERIOD               = ", RESUPPLY_PERIOD, " days")
    println(file, "  RESUPPLY_TIME                 = ", RESUPPLY_TIME, " days")
    println(file, "  KIT_COST                      = ", KIT_COST)
    println(file, "  SHIP_COST                     = ", SHIP_COST)
    println(file, "")
    println(file, "--- Forced Randomization ---")
    println(file, "  INITIAL_CAP                   = ", INITIAL_CAP)
    println(file, "")
    println(file, "--- Recruitment Rate (Gamma(α, 1/β)) ---")
    println(file, "  ALPHA                         = ", ALPHA)
    println(file, "  BETA_OPTIONS                  = ", BETA_OPTIONS)
    println(file, "")
    println(file, "--- Stratification ---")
    println(file, "  NUM_STRATA                    = ", NUM_STRATA)
    println(file, "  STRATA_PROBABILITIES          = ", STRATA_PROBABILITIES)
    println(file, "")
    println(file, "--- Simulation ---")
    println(file, "  NUMBER_SIMULATIONS            = ", NUMBER_SIMULATIONS)
    println(file, "")
    println(file, "--- Supply Strategy Parameters ---")
    println(file, "  Low:  RESUPPLY=", LOW_RESUPPLY,  ", INIT=", LOW_INIT,  ", CRITICAL=", LOW_CRITICAL)
    println(file, "  Med:  RESUPPLY=", MED_RESUPPLY,  ", INIT=", MED_INIT,  ", CRITICAL=", MED_CRITICAL)
    println(file, "  High: RESUPPLY=", HIGH_RESUPPLY, ", INIT=", HIGH_INIT, ", CRITICAL=", HIGH_CRITICAL)
end

open(joinpath(plot_dir, "output.txt"), "w") do file
    for k in 1:ns
        print(file, @sprintf("Average FA (F1a low, z=%d): %.4f\n", k, mean(r.characteristics[1, 1, k, :])))
    end
    print(file, @sprintf("Average FA (F1a low, tot): %.4f\n", mean(r.characteristics[1, 1, ns + 1, :])))
    for k in 1:ns
        print(file, @sprintf("Average FA (F1b low, z=%d): %.4f\n", k, mean(r.characteristics[2, 1, k, :])))
    end
    print(file, @sprintf("Average FA (F1b low, tot): %.4f\n", mean(r.characteristics[2, 1, ns + 1, :])))
    println(file, "")
    for k in 1:ns
        print(file, @sprintf("Avg slots skipped/FR (F1a low, z=%d): %.4f\n", k, mean(r.avg_slots_skipped_f1a[k, :])))
    end
    println(file, "")
    for k in 1:ns
        dm_vals = r.d500s[1, :, k]
        print(file, @sprintf("Mean norm imbalance  (F1a low, z=%d): %.4f\n", k, mean(dm_vals)))
        print(file, @sprintf("Var  norm imbalance  (F1a low, z=%d): %.4f\n", k, var(dm_vals)))
    end
    print(file, @sprintf("Mean norm imbalance  (F1a low, overall): %.4f\n", mean(r.d500s[1, :, :])))
    print(file, @sprintf("Var  norm imbalance  (F1a low, overall): %.4f\n", var(r.d500s[1, :, :])))
    println(file, "")
    for (label, idx) in [("F1a Low", 1), ("F1b Low", 2)]
        A = cov(hcat([r.d500s[idx, :, k] for k in 1:ns]...))
        println(file, "Var-Cov $label: ", A)
        println(file, "\tEigenvalues: ", eigvals(A))
    end
end

end  # !isdefined PlutoRunner
