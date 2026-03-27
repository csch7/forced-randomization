include("constants.jl")

# Holds all mutable state for a single scenario run.
# treatments_used records the full sequence of treatment arms assigned per stratum,
# which allows computing imbalance at any point during the trial.
mutable struct Simulation
    treatments_used::Vector           # Vector{Vector{Int16}} — treatment arm sequence per stratum
    treatments_skipped::Vector        # Vector{Int16} — block position offset per stratum (F1a/F1b)
    need_supply::Set                  # Set{Int16} — centers that have hit their critical point
    delayed_patients::Vector          # Vector{Vector{Int16}} — centers with queued patients per stratum
    patients_sent_home::Vector        # Vector{Int16} — count of patients excluded per stratum
    patients_force_allocated::Vector  # Vector{Int16} — count of forced allocations per stratum
    forward_treated::Vector           # Vector{Vector{Int16}} — block positions already used (F1b)
    num_patients::Int                 # Running total of patient slots to process
    supplies::Dict                    # Dict{center_id => Vector{supply per arm}}
    critical_point::Int               # Supply level that triggers a resupply order
    cap::Int                          # Remaining forced-allocation budget
    treatment_blocks::Matrix          # Matrix{Int} — [stratum × block_position] treatment assignment
end


# ── Allocation step: F0a ─────────────────────────────────────────────────────
# No forced randomization. Patient is allocated only when ALL treatments are
# available at the center; sent home if any supply is zero.
function allocate_f0a!(S::Simulation, center::Int16, stratum::Int8, treatment_index::Int16)
    center_supplies   = S.supplies[center]
    patient_treatment = S.treatment_blocks[stratum, treatment_index]

    if length(findall(iszero, center_supplies)) == 0
        center_supplies[patient_treatment] -= 1
        push!(S.treatments_used[stratum], patient_treatment)

        if center_supplies[patient_treatment] <= S.critical_point
            push!(S.need_supply, center)
        end

    elseif count(iszero, center_supplies) == TREATMENT_ARMS
        # All supplies exhausted — delay patient until next resupply
        push!(S.delayed_patients[stratum], center)
        S.num_patients += 1
    else
        # Only the patient's assigned treatment is missing — send home
        S.patients_sent_home[stratum] += 1
        S.num_patients += 1
    end
end


# ── Allocation step: F0b ─────────────────────────────────────────────────────
# Forced randomization is enabled in principle but cap = 0, so it never
# triggers. Patient is sent home only if their specific assigned treatment is
# unavailable (other arms being zero does not block allocation).
function allocate_f0b!(S::Simulation, center::Int16, stratum::Int8, treatment_index::Int16)
    center_supplies   = S.supplies[center]
    patient_treatment = S.treatment_blocks[stratum, treatment_index]

    if center_supplies[patient_treatment] != 0
        center_supplies[patient_treatment] -= 1
        push!(S.treatments_used[stratum], patient_treatment)

        if center_supplies[patient_treatment] <= S.critical_point
            push!(S.need_supply, center)
        end

    elseif count(iszero, center_supplies) == TREATMENT_ARMS
        push!(S.delayed_patients[stratum], center)
        S.num_patients += 1
    else
        S.patients_sent_home[stratum] += 1
        S.num_patients += 1
    end
end


# ── Allocation step: F1a ─────────────────────────────────────────────────────
# Forced randomization without backfilling. When the assigned treatment is
# unavailable, advance forward in the block to the next available treatment.
# The skipped position is permanently lost (no backfill).
# Returns fr_enabled — false once the cap is exhausted.
function allocate_f1a!(S::Simulation, center::Int16, stratum::Int8, treatment_index::Int16)::Bool
    center_supplies   = S.supplies[center]
    patient_treatment = S.treatment_blocks[stratum, treatment_index + S.treatments_skipped[stratum]]

    if center_supplies[patient_treatment] != 0
        center_supplies[patient_treatment] -= 1
        push!(S.treatments_used[stratum], patient_treatment)

        if center_supplies[patient_treatment] <= S.critical_point
            push!(S.need_supply, center)
        end

    elseif count(iszero, center_supplies) == TREATMENT_ARMS
        push!(S.delayed_patients[stratum], center)
        S.num_patients += 1

    else
        # Advance through the block until an available treatment is found
        while true
            S.treatments_skipped[stratum] += 1
            patient_treatment = S.treatment_blocks[stratum, treatment_index + S.treatments_skipped[stratum]]

            if center_supplies[patient_treatment] != 0
                center_supplies[patient_treatment] -= 1
                push!(S.treatments_used[stratum], patient_treatment)

                if center_supplies[patient_treatment] <= S.critical_point
                    push!(S.need_supply, center)
                end
                break
            end
        end

        S.patients_force_allocated[stratum] += 1
        S.cap -= 1
    end

    return S.cap > 0
end


# ── Allocation step: F1b ─────────────────────────────────────────────────────
# Forced randomization with backfilling. When the assigned treatment is
# unavailable, skip forward in the block; the skipped position is recorded in
# forward_treated so a later patient can backfill it.
# Returns fr_enabled — false once the cap is exhausted.
function allocate_f1b!(S::Simulation, center::Int16, stratum::Int8, treatment_index::Int16)::Bool
    # Skip any positions that have already been forward-allocated
    while count(x -> x == treatment_index + S.treatments_skipped[stratum], S.forward_treated[stratum]) != 0
        S.treatments_skipped[stratum] += 1
    end

    center_supplies   = S.supplies[center]
    patient_treatment = S.treatment_blocks[stratum, treatment_index + S.treatments_skipped[stratum]]

    if center_supplies[patient_treatment] != 0
        center_supplies[patient_treatment] -= 1
        push!(S.treatments_used[stratum], patient_treatment)

        if center_supplies[patient_treatment] <= S.critical_point
            push!(S.need_supply, center)
        end

    elseif count(iszero, center_supplies) == TREATMENT_ARMS
        push!(S.delayed_patients[stratum], center)
        S.num_patients += 1

    else
        # Back up one position (current slot will be backfilled later) and search forward
        S.treatments_skipped[stratum] -= 1
        j = 1

        while true
            # Skip positions already forward-allocated
            while count(x -> x == treatment_index + S.treatments_skipped[stratum] + j, S.forward_treated[stratum]) != 0
                j += 1
            end

            forward_pos      = treatment_index + S.treatments_skipped[stratum] + j
            forward_treatment = S.treatment_blocks[stratum, forward_pos]

            if center_supplies[forward_treatment] != 0
                center_supplies[forward_treatment] -= 1
                push!(S.treatments_used[stratum], forward_treatment)
                push!(S.forward_treated[stratum], forward_pos)

                if center_supplies[forward_treatment] <= S.critical_point
                    push!(S.need_supply, center)
                end
                break
            end
            j += 1
        end

        S.patients_force_allocated[stratum] += 1
        S.cap -= 1
    end

    return S.cap > 0
end
