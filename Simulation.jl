
# Perform simulation with forced randomization disabled
function F0a(supplies::Dict, center::Int64, treatment_index::Int64, treatments_used::Vector{Int64}, blocks::Array{Int64}, need_supply::Set, 
    delayed_patients::Array, patients_sent_home::Int64, num_patients::Int64, critical_pt::Int64, TREATMENT_ARMS::Int64)

    # When all supplies are available
    if (length(findall(iszero, (supplies[center])))==0)
        supplies[center][blocks[treatment_index]]-=1        # Take one supply out of the center
        treatments_used[blocks[treatment_index]]+=1         # Update treatments used (needed for imbalance)

        # If this supply has reached the critical point needed for resupply
        if (supplies[center][blocks[treatment_index]]<=critical_pt)
            push!(need_supply, center)                      # Add this center to the centers that need supply
        end

    # If supplies is equal to zero (not possible with this scenario)
    elseif (count(iszero, supplies[center]) == TREATMENT_ARMS)
        push!(delayed_patients, center)
        num_patients+=1
    # When a supply is missing 
    else 
        patients_sent_home+=1
        num_patients+=1         # Take one more patient from patients list to accomadate for missing patient
    end

    return supplies, treatments_used, need_supply, patients_sent_home, num_patients # Return all modified parameters
end



# Perform simulation with forced randomization enabled, but initial cap set to 0
function F0b(supplies::Dict, center::Int64, treatment_index::Int64, treatments_used::Vector{Int64}, blocks::Array{Int64}, need_supply::Set, 
    delayed_patients::Array, patients_sent_home::Int64, num_patients::Int64, critical_pt::Int64, TREATMENT_ARMS::Int64)

    
    # When the associated supply is available
    if (supplies[center][blocks[treatment_index]]!=0)       # Only difference between F0a and F0b
        supplies[center][blocks[treatment_index]]-=1        # Take one supply out of the center
        treatments_used[blocks[treatment_index]]+=1
        
        # If this supply has reached the critical point needed for resupply
        if (supplies[center][blocks[treatment_index]]<=critical_pt)
            push!(need_supply, center)
        end
    
    # If supplies is equal to zero, delay this patient
    elseif (count(iszero, supplies[center]) ==TREATMENT_ARMS)
        push!(delayed_patients, center)
        num_patients+=1
    # When the supply is missing 
    else 
        patients_sent_home+=1
        num_patients+=1         # Take one more patient from patients list to accommodate for missing patient
    end

    return supplies, treatments_used, need_supply, patients_sent_home, num_patients
end



# Perform simulation with forced randomization enabled, with no backfilling
function F1a(supplies::Dict, center::Int64, treatment_index::Int64, treatments_used::Vector{Int64}, blocks::Array{Int64}, need_supply::Set, delayed_patients::Array,
    treatments_skipped::Int64, num_patients::Int64, critical_pt::Int64, patients_FA::Int64, cap::Int64, TREATMENT_ARMS::Int64)

    # FR always enabled at first
    fr_enabled = true

    # If the supply is not missing, assign as normal
    if (supplies[center][blocks[treatment_index+treatments_skipped]]!=0)
        supplies[center][blocks[treatment_index+treatments_skipped]]-=1
        treatments_used[blocks[treatment_index+treatments_skipped]]+=1

        # If this supply has reached the critical point
        if (supplies[center][blocks[treatment_index]]<=critical_pt)
            push!(need_supply, center)
        end

    # If all supplies are missing, add this patient to the delayed_patients array. Add another patient to accommodate (for now)
    elseif (count(iszero, supplies[center])==TREATMENT_ARMS)
        push!(delayed_patients, center)
        num_patients+=1

    # If neither of these are true, perform forced randomization
    else
        # Continue looping until an available treatment is hit.
        # Once this happens, assign as normal
        while true

            treatments_skipped+=1   # Increment treatments skipped 

            # Assign as normal once a non-missing supply is found
            if (supplies[center][blocks[treatment_index+treatments_skipped]]!=0)
                supplies[center][blocks[treatment_index+treatments_skipped]]-=1
                treatments_used[blocks[treatment_index+treatments_skipped]]+=1

                if (supplies[center][blocks[treatment_index+treatments_skipped]]<=critical_pt)
                    push!(need_supply, center)
                end

                break
            end

        end
        
        patients_FA+=1

        cap-=1
        # Once the cap hits 0, default to F0a
        if (cap==0) 
            fr_enabled = false 
        end
    end

    return supplies, treatments_used, need_supply, delayed_patients, treatments_skipped, patients_FA, num_patients, fr_enabled
end



# Perform simulation with forced randomization enabled, with backfilling.
function F1b(supplies::Dict, center::Int64, treatment_index::Int64, treatments_used::Vector{Int64}, blocks::Array{Int64}, need_supply::Set, delayed_patients::Array,
    forward_treated::Vector{Int64}, treatments_skipped::Int64, num_patients::Int64, critical_pt::Int64, patients_FA::Int64, cap::Int64, TREATMENT_ARMS::Int64)

    fr_enabled = true

    # If the current index is in the list of already treated patients, skip it
    while (count( x->x==treatment_index+treatments_skipped, forward_treated) != 0)
        treatments_skipped+=1
    end

    # Assign as normal
    if (supplies[center][blocks[treatment_index+treatments_skipped]]!=0)
        supplies[center][blocks[treatment_index+treatments_skipped]]-=1
        treatments_used[blocks[treatment_index+treatments_skipped]]+=1

        if (supplies[center][blocks[treatment_index+treatments_skipped]]<=critical_pt)
            push!(need_supply, center)
        end

    # If all treatments missing, delay the patient
    elseif (count(iszero, supplies[center]) == TREATMENT_ARMS)
        push!(delayed_patients, center)
        num_patients+=1

    # Otherwise, perform FR with backfilling
    else
        # Here, treatments_skipped defines the "offset" due to patients FA'd, it's the difference between the patient index and the treatment index
        treatments_skipped-=1 

        j=1
        while true

            # If the current index has already been treated, skip it
            while (count( x->x==treatment_index+treatments_skipped+j, forward_treated) != 0)
                j+=1
            end

            # Assign as normal if the treatment is available
            if (supplies[center][blocks[treatment_index+treatments_skipped+j]]!=0)
                supplies[center][blocks[treatment_index+treatments_skipped+j]]-=1
                treatments_used[blocks[treatment_index+treatments_skipped+j]]+=1

                # Add this patient to the forward_treated array
                push!(forward_treated,treatment_index+treatments_skipped+j)

                if (supplies[center][blocks[treatment_index+treatments_skipped+j]]<=critical_pt)
                    push!(need_supply, center)
                end

                # Break out of loop once available treatment is found
                break
            end

            j+=1

        end

        patients_FA+=1

        cap-=1
        if (cap==0) 
            fr_enabled = false 
        end
    end


    return supplies, treatments_used, need_supply, delayed_patients, forward_treated, treatments_skipped, patients_FA, num_patients, fr_enabled
end
