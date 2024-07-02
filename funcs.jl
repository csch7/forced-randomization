using Distributions

low_init = [2,2]
low_resup = 2
low_crit = 1
med_init = [3,3]
med_resup = 4
med_crit = 1
high_init = [4,4]
high_resup = 5
high_crit = 2


# Generates patient arrival times
@views function generatePatList(rates::Array{Float64}, patients::Matrix{Float64}, center_activations::Array{Int64}, numCenters::Int64, numPatients::Int64)

    @simd for i in 1:numCenters

        # Generate numbers as stated in BMC paper
        uis = rand(Uniform(0,1),numPatients)
        yis = -1 .*log.(1 .-uis)./rates[i]
        xis = cumsum(yis)

        # Construct matrix of patient ids (A[:,1]), patient arrival times (A[:,2]), and patient centers (A[:,3])
        A = hcat((i-1)*numPatients.+collect(1:numPatients), xis.+(center_activations[i]*30), fill(i,numPatients))

        patients = vcat(patients, A) # Note that the passed in matrix is always empty

    end

    patients = patients[sortperm(patients[:,2]), :] # Sort by arrival times
    return patients
end

# Generate treatment list for PBD
@views function generateTreatList(ratio::Tuple, blocks::Array{Int64}, sample_size::Int64, treat_arms::Int64, block_size::Int64 )
    block_options = [t for t in 1:treat_arms for a in 1:(ratio[t]/sum(ratio)*block_size)]
    for i in 1:convert(Int64, ceil(2*sample_size/block_size)) # Generate double the treatments needed (for patients sent home/FA) (may need more for very fast recruitment rate)
        blocks = vcat(blocks, shuffle!(block_options))
    end

    return blocks
end

# NOTE: Currently only works for 2 treatment arms, 1:1 ratio
# Returns (Resupply, Initial supply, Critical pt, fr_allowed, backfill_enabled, cap), given a scenario (1-12)
function setParams( scen::Int64, cap::Int64 )

    if(scen%3==1) # Low supply strategy
        if(scen<5)
            if(scen==1)  # F0
                return low_resup , low_init , low_crit , false , false , 0
            else         # F1
                return low_resup , low_init , low_crit , true , false , 0
            end
        elseif (scen==7) # F2a
            return low_resup , low_init , low_crit , true , false , cap
        else             # F2b
            return low_resup , low_init , low_crit , true , true , cap
        end

    elseif (scen%3==2) # Medium supply strategy
        if(scen<6)
            if(scen==2)
                return med_resup , med_init , med_crit , false , false , 0
            else
                return med_resup , med_init , med_crit , true , false , 0
            end
        elseif (scen==8)
            return med_resup , med_init , med_crit , true , false , cap
        else
            return med_resup , med_init , med_crit , true , true , cap
        end

    else # High supply strategy
        if(scen<7)
            if(scen==3)
                return high_resup , high_init , high_crit , false , false , 0
            else
                return high_resup , high_init , high_crit , true , false , 0
            end
        elseif (scen==9)
            return high_resup , high_init , high_crit , true , false , cap
        else
            return high_resup , high_init , high_crit , true , true , cap
        end
        

    end
end

function q90(chars; dims::Int64)
    return mean(chars, dims=dims).+1.645.*std(chars, dims=dims)
end




