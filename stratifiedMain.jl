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

plts = Array{Any}(undef, 8, length(BETA_OPTIONS))

for (bi, BETA) in enumerate(BETA_OPTIONS)

    # characteristics[scenario_offset, characteristic, stratum_or_total, simulation]
    # scenario_offset: 1=F1a Low, 2=F1b Low  (scenarios 7 and 10)
    # characteristic:  1=FA_rate
    # stratum_or_total: 1=z1, 2=z2, 3=overall
    characteristics    = zeros(Float64, 2, 1, 3, NUMBER_SIMULATIONS)

    # dm(z, m): standardised treatment imbalance at patient index m within stratum z
    dm1s               = zeros(Float64, 2, NUMBER_SIMULATIONS, Int(STRATA_ASSIGNMENT_PROBABILITY * SAMPLE_SIZE))
    dm2s               = zeros(Float64, 2, NUMBER_SIMULATIONS, Int((1 - STRATA_ASSIGNMENT_PROBABILITY) * SAMPLE_SIZE))

    # dlg(z, m): lag between expected and actual stratum-z allocations at overall index i
    dlgz1s             = zeros(Float64, 2, NUMBER_SIMULATIONS, Int(STRATA_ASSIGNMENT_PROBABILITY * SAMPLE_SIZE))
    dlgz2s             = zeros(Float64, 2, NUMBER_SIMULATIONS, Int((1 - STRATA_ASSIGNMENT_PROBABILITY) * SAMPLE_SIZE))

    # d500(z): total imbalance at end of trial for stratum z
    d500z1s            = zeros(Float64, 2, NUMBER_SIMULATIONS)
    d500z2s            = zeros(Float64, 2, NUMBER_SIMULATIONS)

    # unfilled slots: how many block positions have been skipped at index m (F1a only)
    unfilled_slots_z1  = zeros(Float64, 1, NUMBER_SIMULATIONS, Int(STRATA_ASSIGNMENT_PROBABILITY * SAMPLE_SIZE))
    unfilled_slots_z2  = zeros(Float64, 1, NUMBER_SIMULATIONS, Int((1 - STRATA_ASSIGNMENT_PROBABILITY) * SAMPLE_SIZE))

    for sim in tqdm(1:NUMBER_SIMULATIONS)

        center_rates = rand(Gamma(ALPHA, 1/BETA), CENTERS)
        center_acts  = rand(0:4, CENTERS)
        patients     = generate_patient_arrivals(center_rates, center_acts, CENTERS, SAMPLE_SIZE)

        # Build a 2×N treatment block matrix (one independent PBD sequence per stratum).
        # generate_treatment_blocks returns 1×(2*SAMPLE_SIZE); reshape splits it into two rows.
        raw_blocks       = generate_treatment_blocks(ALLOCATION_RATIO, SAMPLE_SIZE, TREATMENT_ARMS, BLOCK_SIZE)
        treatment_blocks = Matrix(transpose(reshape(raw_blocks, (SAMPLE_SIZE, 2))))

        # Pre-generate stratum assignments for all possible patient slots
        strata = generate_strata_assignments(SAMPLE_SIZE * 2, STRATA_ASSIGNMENT_PROBABILITY)

        for (si, scenario) in enumerate([7, 10])

            resupply_amt, init_supply, critical_pt, fr_allowed, backfill_enabled, cap =
                scenario_params(scenario, INITIAL_CAP)

            center_supplies = Dict{Int, Vector{Int}}()
            for i in 1:CENTERS
                center_supplies[i] = copy(init_supply)
            end

            total_drugs = CENTERS * sum(init_supply)
            total_cost  = CENTERS * sum(init_supply) * KIT_COST

            S = Simulation(
                [Int16[] for _ in 1:2],    # treatments_used
                zeros(Int16, 2),           # treatments_skipped
                Set{Int16}(),              # need_supply
                [Int16[] for _ in 1:2],   # delayed_patients
                zeros(Int16, 2),           # patients_sent_home
                zeros(Int16, 2),           # patients_force_allocated
                [Int16[] for _ in 1:2],   # forward_treated
                SAMPLE_SIZE,               # num_patients
                center_supplies,           # supplies
                critical_pt,               # critical_point
                cap,                       # cap
                treatment_blocks,          # treatment_blocks
            )

            next_supply_check  = RESUPPLY_PERIOD
            next_resupply      = next_supply_check + RESUPPLY_TIME
            sent_supply        = Dict{Int, Vector{Int}}()
            tot_delayed        = zeros(Int16, 2)
            num_waitlisted     = zeros(Int16, 2)
            patients_per_stratum = zeros(Int16, 2)
            break_loop         = false

            i = 1
            while i <= S.num_patients
                center = Int16(patients[i, 3])

                # ── Resupply order ───────────────────────────────────────────
                if patients[i, 2] >= next_supply_check
                    for j in S.need_supply
                        new_supply = zeros(Int, TREATMENT_ARMS)
                        for k in eachindex(S.supplies[j])
                            if S.supplies[j][k] <= critical_pt
                                top_up         = resupply_amt - S.supplies[j][k]
                                total_cost    += top_up * KIT_COST
                                total_drugs   += top_up
                                new_supply[k]  = resupply_amt
                            end
                        end
                        sent_supply[j] = new_supply
                        total_cost    += SHIP_COST
                    end
                    next_supply_check += RESUPPLY_PERIOD
                    empty!(S.need_supply)
                end

                # ── Resupply delivery + delayed patient allocation ───────────
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
                            num_waitlisted[z] += min(num, resupply_amt * TREATMENT_ARMS)
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
                                unfilled_slots_z1[si, sim, min(600, stratum_index)] = S.treatments_skipped[1]
                            else
                                unfilled_slots_z2[si, sim, min(400, stratum_index)] = S.treatments_skipped[2]
                            end
                        else
                            fr_allowed = allocate_f1b!(S, center_delayed, stratum, stratum_index)
                        end

                        patients_per_stratum[stratum] += 1

                        i_1 = trunc(Int16, (i - patients_per_stratum[2]) - tot_delayed[1] - length(S.delayed_patients[1]))
                        i_2 = trunc(Int16, (i - patients_per_stratum[1]) - tot_delayed[2] - length(S.delayed_patients[2]))
                        dlgz1s[si, sim, min(600, stratum_index)] = patients_per_stratum[1] - i_1
                        dlgz2s[si, sim, min(400, stratum_index)] = patients_per_stratum[2] - i_2

                        i += 1
                        if i > S.num_patients
                            patients[S.num_patients, 2] = next_resupply
                            break_loop = true
                            break
                        end
                    end

                    S.delayed_patients[1] = S.delayed_patients[1][num_delayed[1]+1:end]
                    S.delayed_patients[2] = S.delayed_patients[2][num_delayed[2]+1:end]
                    next_resupply = next_supply_check + RESUPPLY_TIME
                    empty!(sent_supply)
                end

                if break_loop
                    break
                end

                # ── Normal patient allocation ────────────────────────────────
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
                        unfilled_slots_z1[si, sim, min(600, stratum_index)] = S.treatments_skipped[1]
                    else
                        unfilled_slots_z2[si, sim, min(400, stratum_index)] = S.treatments_skipped[2]
                    end
                else
                    fr_allowed = allocate_f1b!(S, center, stratum, stratum_index)
                end

                patients_per_stratum[stratum] += 1

                i_1 = trunc(Int16, i - patients_per_stratum[2])
                i_2 = trunc(Int16, i - patients_per_stratum[1])
                dlgz1s[si, sim, min(600, stratum_index)] = patients_per_stratum[1] - i_1
                dlgz2s[si, sim, min(400, stratum_index)] = patients_per_stratum[2] - i_2

                i += 1
            end  # patient loop

            # ── Compute dm(z, m) — standardised imbalance at each index m ───
            # F1a (scenarios 7–9): divide by sqrt(m) to get a normalised statistic
            # F1b (scenarios 10–12): use raw imbalance
            normalise = false
            # scenario in 7:9

            for c in 1:min(600, length(S.treatments_used[1]))
                cm1 = countmap(S.treatments_used[1][1:c])
                t1  = get(cm1, 1, 0)
                t2  = get(cm1, 2, 0)
                dm1s[si, sim, c] = normalise ? (t1 - t2) / sqrt(c) : (t1 - t2)
            end
            for c in 1:min(400, length(S.treatments_used[2]))
                cm2 = countmap(S.treatments_used[2][1:c])
                t1  = get(cm2, 1, 0)
                t2  = get(cm2, 2, 0)
                dm2s[si, sim, c] = normalise ? (t1 - t2) / sqrt(c) : (t1 - t2)
            end

            # ── End-of-trial imbalance d500(z) ───────────────────────────────────
            cm1_full = countmap(S.treatments_used[1])
            n1 = length(S.treatments_used[1])
            d1 = get(cm1_full, 1, 0) - get(cm1_full, 2, 0)
            d500z1s[si, sim] = normalise ? d1 / sqrt(n1) : d1
            cm2_full = countmap(S.treatments_used[2])
            n2 = length(S.treatments_used[2])
            d2 = get(cm2_full, 1, 0) - get(cm2_full, 2, 0)
            d500z2s[si, sim] = normalise ? d2 / sqrt(n2) : d2

            # ── FA rate characteristics ───────────────────────────────────────
            characteristics[si, 1, 1, sim] = S.patients_force_allocated[1] / length(S.treatments_used[1])
            characteristics[si, 1, 2, sim] = S.patients_force_allocated[2] / length(S.treatments_used[2])
            characteristics[si, 1, 3, sim] = sum(S.patients_force_allocated) / SAMPLE_SIZE

        end  # scenario
    end  # sim

    # ── Plots ────────────────────────────────────────────────────────────────
    # Common cutoff: trim both strata to the shorter one (z2) so x-axes match.

    # rm1_data = replace!(dm1s[1:1, :, :] ./ sqrt.(unfilled_slots_z1), NaN => 0)
    # rm2_data = replace!(dm2s[1:1, :, :] ./ sqrt.(unfilled_slots_z2), NaN => 0)

    # # F1a Low: 3×2 grid — (Var, Q90, Var[r]) × (stratum 1, stratum 2)
    # f1a_panels = [
    #     imbalance_line_panel(dm1s[1:1, :, :], "Var[dm(1)]",  var;        max_value=1, cutoff=528, start=4),
    #     imbalance_line_panel(dm2s[1:1, :, :], "Var[dm(2)]",  var;        max_value=1, cutoff=328, start=4),
    #     imbalance_line_panel(dm1s[1:1, :, :], "Q90[dm(1)]",  quantile_90; max_value=1, cutoff=528, start=4),
    #     imbalance_line_panel(dm2s[1:1, :, :], "Q90[dm(2)]",  quantile_90; max_value=1, cutoff=328, start=4),
    #     # imbalance_line_panel(rm1_data,         "Var[rm(1)]",  var;        cutoff=528, start=4),
    #     # imbalance_line_panel(rm2_data,         "Var[rm(2)]",  var;        cutoff=328, start=4),
    # ]
    # save_line_summary(f1a_panels, 2, 2, joinpath("plots", "f1a_low.png"), "F1a Low Supply")

    # # F1b Low: 4×2 grid — (Var, Q90, Q90[DLG], max[DLG]) × (stratum 1, stratum 2)
    # f1b_panels = [
    #     imbalance_line_panel(dm1s[2:2, :, :],   "Var[Dm(1)]",    var;        cutoff=528, start=4),
    #     imbalance_line_panel(dm2s[2:2, :, :],   "Var[Dm(2)]",    var;        cutoff=328, start=4),
    #     imbalance_line_panel(dm1s[2:2, :, :],   "Q90[Dm(1)]",    quantile_90; cutoff=528, start=4),
    #     imbalance_line_panel(dm2s[2:2, :, :],   "Q90[Dm(2)]",    quantile_90; cutoff=328, start=4),
    #     imbalance_line_panel(dlgz1s[2:2, :, :], "Q90[DLG(z=1)]", quantile_90; cutoff=528),
    #     imbalance_line_panel(dlgz2s[2:2, :, :], "Q90[DLG(z=2)]", quantile_90; cutoff=328),
    #     imbalance_line_panel(dlgz1s[2:2, :, :], "max[DLG(z=1)]", maximum;    cutoff=528),
    #     imbalance_line_panel(dlgz2s[2:2, :, :], "max[DLG(z=2)]", maximum;    cutoff=328),
    # ]
    # save_line_summary(f1b_panels, 4, 2, joinpath("plots", "f1b_low.png"), "F1b Low Supply")

    plot_imbalance_histograms(dm1s[1:1, :, :], dm2s[1:1, :, :], 528, 328, joinpath("plots", "dm_hists_f1a.png"), "F1a Low")
    plot_imbalance_histograms(dm1s[2:2, :, :], dm2s[2:2, :, :], 528, 328, joinpath("plots", "dm_hists_f1b.png"), "F1b Low")

    plot_joint_normality_mahalanobis(
        d500z1s[1, :], d500z2s[1, :],
        joinpath("plots", "joint_normality_f1a.png"),
        "F1a Low Supply — Joint Normality of End-of-Trial Imbalances",
    )

    # ── Variance-covariance and eigenvalue analysis ───────────────────────────
    labels_vc = [("F1a Low", 1), ("F1b Low", 2)]
    for (label, idx) in labels_vc
        A = cov(hcat(d500z1s[idx, :], d500z2s[idx, :]))
        println("Var-Cov $label: ", A)
        println("\tEigenvalues: ", eigvals(A))
    end

    # ── FA rate summary ───────────────────────────────────────────────────────
    @printf("Average FA (F1a low, z=1): %.4f\n", mean(characteristics[1, 1, 1, :]))
    @printf("Average FA (F1a low, z=2): %.4f\n", mean(characteristics[1, 1, 2, :]))
    @printf("Average FA (F1a low, tot): %.4f\n", mean(characteristics[1, 1, 3, :]))
    @printf("Average FA (F1b low, z=1): %.4f\n", mean(characteristics[2, 1, 1, :]))
    @printf("Average FA (F1b low, z=2): %.4f\n", mean(characteristics[2, 1, 2, :]))
    @printf("Average FA (F1b low, tot): %.4f\n", mean(characteristics[2, 1, 3, :]))

end  # BETA
