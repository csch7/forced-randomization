using DataFrames
using Distributions
using Random
using Printf
using Statistics
using StatsBase
include("funcs.jl")
include("plotting.jl")
include("Simulation.jl")

mkpath("plots")

# Column labels for CSV exports (first two are scenario identifiers, rest are characteristics)
const HEADINGS = [
    "IRT Approach", "Resupply Strategy",
    "Treatment_Imbalance", "Pct_FAs", "Pct_Patients_Sent_Home",
    "Drug_Overage", "Cost", "Time", "Patients_NA", "Patients_Waitlisted",
]

# Plots are collected across beta values for side-by-side comparison figures
plts = Array{Any}(undef, 8, length(BETA_OPTIONS))

for (bi, BETA) in enumerate(BETA_OPTIONS)

    # characteristics[scenario, characteristic, simulation]
    characteristics = zeros(Float64, 12, 8, NUMBER_SIMULATIONS)

    for sim in 1:NUMBER_SIMULATIONS

        center_rates = rand(Gamma(ALPHA, 1/BETA), CENTERS)
        center_acts  = rand(0:4, CENTERS)
        patients     = generate_patient_arrivals(center_rates, center_acts, CENTERS, SAMPLE_SIZE)

        # 1×N treatment block matrix; accessed as treatment_blocks[1, position]
        treatment_blocks = generate_treatment_blocks(ALLOCATION_RATIO, SAMPLE_SIZE, TREATMENT_ARMS, BLOCK_SIZE)

        for scenario in 1:12

            resupply_amt, init_supply, critical_pt, fr_allowed, backfill_enabled, cap =
                scenario_params(scenario, INITIAL_CAP)

            center_supplies = Dict{Int, Vector{Int}}()
            for i in 1:CENTERS
                center_supplies[i] = copy(init_supply)
            end

            total_drugs = CENTERS * sum(init_supply)
            total_cost  = CENTERS * sum(init_supply) * KIT_COST

            # All state is held in a Simulation struct with a single stratum (index 1).
            S = Simulation(
                [Int64[]],          # treatments_used
                zeros(Int16, 1),    # treatments_skipped
                Set{Int}(),         # need_supply
                [Int[]],            # delayed_patients
                zeros(Int16, 1),    # patients_sent_home
                zeros(Int16, 1),    # patients_force_allocated
                [Int[]],            # forward_treated
                SAMPLE_SIZE,        # num_patients
                center_supplies,    # supplies
                critical_pt,        # critical_point
                cap,                # cap
                treatment_blocks,   # treatment_blocks
            )

            next_supply_check = RESUPPLY_PERIOD
            next_resupply     = next_supply_check + RESUPPLY_TIME
            sent_supply       = Dict{Int, Vector{Int}}()
            tot_delayed       = 0
            num_waitlisted    = 0
            break_loop        = false

            i = 1
            while i <= S.num_patients
                center = Int16(patients[i, 3])

                # ── Resupply order ───────────────────────────────────────────
                if patients[i, 2] >= next_supply_check
                    for j in S.need_supply
                        new_supply = zeros(Int, TREATMENT_ARMS)
                        for k in eachindex(S.supplies[j])
                            if S.supplies[j][k] <= critical_pt
                                top_up          = resupply_amt - S.supplies[j][k]
                                total_cost     += top_up * KIT_COST
                                total_drugs    += top_up
                                new_supply[k]   = resupply_amt
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

                    num_delayed  = length(S.delayed_patients[1])
                    tot_delayed += num_delayed

                    cts = countmap(S.delayed_patients[1])
                    for (_, num) in cts
                        num_waitlisted += min(num, resupply_amt * TREATMENT_ARMS)
                    end

                    for j in 1:num_delayed
                        new_delayed = length(S.delayed_patients[1]) - num_delayed
                        if !fr_allowed
                            ti = Int16(i - S.patients_sent_home[1] - tot_delayed - new_delayed)
                            allocate_f0a!(S, Int16(S.delayed_patients[1][j]), Int8(1), ti)
                        elseif cap == 0
                            ti = Int16(i - S.patients_sent_home[1] - tot_delayed - new_delayed)
                            allocate_f0b!(S, Int16(S.delayed_patients[1][j]), Int8(1), ti)
                        elseif !backfill_enabled
                            ti = Int16(i - tot_delayed - new_delayed)
                            fr_allowed = allocate_f1a!(S, Int16(S.delayed_patients[1][j]), Int8(1), ti)
                        else
                            ti = Int16(i - tot_delayed - new_delayed)
                            fr_allowed = allocate_f1b!(S, Int16(S.delayed_patients[1][j]), Int8(1), ti)
                        end

                        i += 1
                        if i > S.num_patients
                            patients[S.num_patients, 2] = next_resupply
                            break_loop = true
                            break
                        end
                    end

                    S.delayed_patients[1] = S.delayed_patients[1][num_delayed+1:end]
                    next_resupply = next_supply_check + RESUPPLY_TIME
                    empty!(sent_supply)
                end

                if break_loop
                    break
                end

                # ── Normal patient allocation ────────────────────────────────
                if !fr_allowed
                    ti = Int16(i + S.treatments_skipped[1] - S.patients_sent_home[1] - tot_delayed - length(S.delayed_patients[1]))
                    allocate_f0a!(S, center, Int8(1), ti)
                elseif cap == 0
                    ti = Int16(i - S.patients_sent_home[1] - tot_delayed - length(S.delayed_patients[1]))
                    allocate_f0b!(S, center, Int8(1), ti)
                elseif !backfill_enabled
                    ti = Int16(i - tot_delayed - length(S.delayed_patients[1]))
                    fr_allowed = allocate_f1a!(S, center, Int8(1), ti)
                else
                    ti = Int16(i - tot_delayed - length(S.delayed_patients[1]))
                    fr_allowed = allocate_f1b!(S, center, Int8(1), ti)
                end

                i += 1
            end

            # ── Per-scenario summary ─────────────────────────────────────────
            total_time = patients[S.num_patients - tot_delayed, 2]
            tu = S.treatments_used[1]
            treatment_imbalance = abs(count(==(1), tu) - count(==(2), tu))

            @printf("\nSim: %i  Scenario: %i\n", sim, scenario)
            @printf("  Treatment Imbalance:  %i\n", treatment_imbalance)
            @printf("  Patients Sent Home:   %i\n", S.patients_sent_home[1])
            @printf("  Patients FA:          %i\n", S.patients_force_allocated[1])
            @printf("  Unused Drugs:         %.4f\n", Float64(total_drugs - SAMPLE_SIZE) / SAMPLE_SIZE)
            @printf("  Total Time:           %.2f\n", total_time)
            @printf("  Total Cost:           %i\n", total_cost)
            @printf("  Patients left (end):  %i\n", length(S.delayed_patients[1]))
            @printf("  Total waitlisted:     %i\n", num_waitlisted + length(S.delayed_patients[1]))

            characteristics[scenario, 1, sim] = treatment_imbalance
            characteristics[scenario, 2, sim] = S.patients_force_allocated[1] / SAMPLE_SIZE
            characteristics[scenario, 3, sim] = S.patients_sent_home[1] / SAMPLE_SIZE
            characteristics[scenario, 4, sim] = (total_drugs - SAMPLE_SIZE) / SAMPLE_SIZE
            characteristics[scenario, 5, sim] = total_cost
            characteristics[scenario, 6, sim] = total_time
            characteristics[scenario, 7, sim] = length(S.delayed_patients[1])
            characteristics[scenario, 8, sim] = num_waitlisted + length(S.delayed_patients[1])

        end  # scenario
    end  # sim

    # ── CSV exports ──────────────────────────────────────────────────────────
    tag = string("b_", BETA, "-n", SAMPLE_SIZE)
    export_time_statistics(characteristics[:, 6, :],  joinpath("plots", "time-stats-$tag.csv"))
    export_time_statistics(characteristics[:, 2, :],  joinpath("plots", "FA-stats-$tag.csv"))
    export_scenario_statistics(characteristics, HEADINGS, joinpath("plots", "mean-characteristics-$tag.csv"), mean)
    export_scenario_statistics(characteristics, HEADINGS, joinpath("plots", "std-characteristics-$tag.csv"),  std)
    export_scenario_statistics(characteristics, HEADINGS, joinpath("plots", "q90-characteristics-$tag.csv"),  quantile_90)

    # ── Bar charts ───────────────────────────────────────────────────────────
    scenario_labels = repeat(["FR0a", "FR0b", "FR1a", "FR1b"], inner=3)
    recruitment_label = (BETA == 28 ? "Base Case " : "Faster Recruitment ") *
                        "(α=$ALPHA, β=$BETA)"

    function save_bar(idx, fig_label, ylab, char_col)
        p = grouped_bar_with_error(
            fig_label, recruitment_label, scenario_labels, ylab,
            mean(characteristics[:, char_col, :], dims=2),
            1.645 .* std(characteristics[:, char_col, :], dims=2),
        )
        plts[idx, bi] = p
        savefig(p, joinpath("plots", "$(replace(fig_label, r"[^A-Za-z0-9]" => "-"))-$tag.png"))
    end

    save_bar(1, "Fig.4a: Treatment Imbalance",         "Imbalance (num treatments)",     1)
    save_bar(2, "Fig.4b: Forced Allocations",          "Proportion of FAs",              2)
    save_bar(3, "Fig.4c: Patients Sent Home",          "Proportion of Patients",         3)
    save_bar(4, "Fig.4d: Drug Overage",                "Overage (Pct)",                  4)
    save_bar(5, "Fig.4x: Total Cost",                  "Average Cost (USD)",             5)
    save_bar(6, "Fig.4e: Time to Complete Recruitment","Average Time (days)",            6)
    save_bar(7, "Fig.4f: Patients Not Allocated",      "Number of patients",             7)
    save_bar(8, "Fig.4g: Patients Waitlisted",         "Number of patients",             8)

    plot_recruitment_time_histograms(
        characteristics[:, 6, :],
        joinpath("plots", "Time-histograms-$tag.png"),
    )

end  # BETA

# ── Combined side-by-side plots (when multiple beta values tested) ────────────
if length(BETA_OPTIONS) > 1
    for (ri, row) in enumerate(eachrow(plts))
        combined = plot(row..., layout=(1, length(row)))
        savefig(combined, joinpath("plots", "Combined-$(HEADINGS[ri+2]).png"))
    end
end
