using DataFrames
using Distributions
using Random
using Printf
using Statistics
using StatsBase
include("funcs.jl")
include("plotting.jl")
include("Simulation.jl")

# Constant Parameters

SAMPLE_SIZE = 500                   # Total number of patients
# SAMPLE_SIZE = 100                   # Smaller sample size
TREATMENT_ARMS = 2                  # Number of treatments
ALLOC_RATIO = (1,1)                 # Ratio of treatments in each block
BLOCK_SIZE = 4                      # Number of patients allocated to each block
CENTERS = 80                        # Number of centers
# CENTERS = 16                        # Smaller number of centers (for smaller sample size)
RESUPPLY_PERIOD = 7                 # Number of days after which resupply is evaluated and requested
RESUPPLY_TIME = 3                   # Time it takes supplies to reach sites
KIT_COST = 10000                    # Price of one treatment kit
SHIP_COST = 100                     # Price of shipping
INIT_CAP = Int64(0.3*SAMPLE_SIZE)   # Set cap equal to 30% of sample size


# Parameters for patient recruitment
ALPHA = 1.2
BETA_OPTS = [28,16,5]                # Options for beta (all in this array will be run)



num_simulations = 5000                        # Number of simulations per tested beta
heading = ["IRT Approach", "Resupply Strategy","Treatment_Imbalance", "Pct_FAs", "Pct_Patients_Sent_Home", "Drug_Overage","Cost","Time", "Patients_NA", "Patients_Waitlisted"] # Heading for csv files, and used for file names
plts = Array{Any}(undef, 8, length(BETA_OPTS))  # Save individual plots for each beta to be used for side-by-side

for (bi, BETA) ∈ enumerate(BETA_OPTS)

    characteristics = zeros(Float64, 12,8,num_simulations)


    for sim ∈ 1:num_simulations
        center_rates = rand(Gamma(ALPHA, 1/BETA), CENTERS)      # Randomly generate center recruitment rates
        center_acts = Array{Int64}(undef, 1, CENTERS)           # Holds center activation times
        center_supplies = Dict()                                # Keep track of supplies in each center (key: center, value: Array(num A, num B, ...))
        blocks = Array{Int64}(undef, 0,1)                       # Treatment list via PBD
        blocks = generateTreatList(ALLOC_RATIO, blocks, SAMPLE_SIZE, TREATMENT_ARMS, BLOCK_SIZE)

        patients = Array{Float64}(undef, 0, 3)                  # first column: ID, second column: arrival time, third column: center
        rand!(center_acts, 0:4)                                 # Randomly generate center activation times between 0 and 4 months
        patients = generatePatList(center_rates, patients, center_acts, CENTERS, SAMPLE_SIZE)

        for scenario ∈ 1:12

            resupply_amt, init_supply, critical_pt, fr_allowed, backfill_enabled, cap = setParams(scenario, INIT_CAP)  # Set supply/FR parameters, dependent on the scenario
            num_patients = SAMPLE_SIZE         # Keeps track of total patients which need to be "seen"
            patients_force_allocated = 0
            patients_sent_home = 0          
            total_drugs = CENTERS*sum(init_supply)          # Initial amount of total drugs equal to total drugs at all centers
            total_cost = CENTERS*sum(init_supply)*KIT_COST  # Iniital total cost equal to total drugs * cost of one kit

            for i ∈ 1:CENTERS
                center_supplies[i]=copy(init_supply)        # Assign initial supply to each center's supplies
            end

            next_supply_check = RESUPPLY_PERIOD             # Keeps track of when the supply needs to be checked           
            next_resupply = next_supply_check+RESUPPLY_TIME # Keeps track of when the supply actually arrives at the centers (constant for all centers)
            treatments_skipped = 0                          # Keeps track of how many treatments have been skipped over in specification F1a
            forward_treated = Vector{Int64}(undef, 0)       # Keeps track of which treatments (forward) have already been allocated (relevant in F1b)
            sent_supply = Dict()                            # Keeps track of the supplies which will be in the center once shipment arrives (keys: centers, values: supply)
            need_supply = Set()                             # Keeps track of which centers need supply (set so that no duplicates happen)
            treatments_used = zeros(Int64, TREATMENT_ARMS)  # List of length TREATMENT_ARMS, where each index is the amount of that supply used (needed for treatment imbalance)

            tot_delayed = 0                                 # Amount of total patients delayed due to missing supply
            num_waitlisted = 0                              # Amount of UNIQUE patients delayed (note difference with tot_delayed)
            delayed_patients = []                           # Centers where a patient was scheduled to arrive but no treatment was available
            break_loop = false                              # Needed in case the sample size is hit while delayed patients are being allocated

            i=1
            while i ≤ num_patients
                center = Int64(patients[i,3])

                if (patients[i,2] ≥ next_supply_check)       # If the next patient arrives after the next supply check
                    for j in need_supply  # Go through each center that needs supply
                        new_supply = []
                        for k in eachindex(center_supplies[j])  # Go through each treatment in that center
                            if (center_supplies[j][k] <= critical_pt)   # If that treatment is below the critical point, set the sent_supply equal to the resupply amount
                                total_cost+=(resupply_amt-center_supplies[j][k])*KIT_COST
                                total_drugs+=(resupply_amt-center_supplies[j][k])
                                push!(new_supply, resupply_amt)
                            else
                                push!(new_supply, 0)    # Otherwise, don't send any supply (keep it the same)
                            end
                        end
                        sent_supply[j] = new_supply
                        total_cost+=SHIP_COST
                    end
                    
                    next_supply_check+=RESUPPLY_PERIOD
                    empty!(need_supply)
                end



                if (patients[i,2]>=next_resupply)            # If the next patient arrives after the resupply arrives
                    for (cent, newsupply) in sent_supply
                        for k in eachindex(newsupply)
                            if (newsupply[k]!=0)    # If more supply is sent, add it to the center
                                center_supplies[cent][k] = newsupply[k]
                            end
                        end
                    end
                    
                    num_delayed = length(delayed_patients)  # Save the current number of delayed patients for later
                    tot_delayed += num_delayed 

                    cts = countmap(delayed_patients)
                    for (cent, num) in cts
                        if (num >= resupply_amt*TREATMENT_ARMS)
                            num_waitlisted+=resupply_amt*TREATMENT_ARMS
                        else
                            num_waitlisted+=num
                        end
                    end

                    # For each patient/center in the delayed list,
                    for j in 1:num_delayed
                        # Perform simulation.
                        if (!fr_allowed)
                            center_supplies, treatments_used, need_supply, patients_sent_home, num_patients = F0a(
                                center_supplies,        # Current supplies
                                delayed_patients[j],    # Relevant center
                                i-patients_sent_home-tot_delayed-(length(delayed_patients)-num_delayed),  # Index of next unused treatment. When patients are waitlisted, treatment index should not increase.
                                treatments_used,        # Passed in to modify
                                blocks,                 # Treatment list
                                need_supply,            # Passed in to modify
                                delayed_patients,       
                                patients_sent_home,     # Passed in to modify/needed for index
                                num_patients,           # Passed in to modify
                                critical_pt,            # Needed to check resupply
                                TREATMENT_ARMS)         # Needed to check if all treatments are missing
                        elseif (cap==0)
                            center_supplies, treatments_used, need_supply, patients_sent_home, num_patients = F0b(
                                center_supplies, 
                                delayed_patients[j], 
                                i-patients_sent_home-tot_delayed-(length(delayed_patients)-num_delayed),
                                treatments_used, 
                                blocks, 
                                need_supply,
                                delayed_patients, 
                                patients_sent_home, 
                                num_patients, 
                                critical_pt,
                                TREATMENT_ARMS)
                        elseif (!backfill_enabled)
                            center_supplies, treatments_used, need_supply, delayed_patients, treatments_skipped, patients_force_allocated, num_patients, fr_allowed = F1a(
                                center_supplies,
                                delayed_patients[j], 
                                i-tot_delayed-(length(delayed_patients)-num_delayed), 
                                treatments_used, 
                                blocks, 
                                need_supply, 
                                delayed_patients,           # Passed in to modify
                                treatments_skipped,         # Needed for treatment index and modified
                                num_patients, 
                                critical_pt,
                                patients_force_allocated,   # Passed in to modify
                                cap,                        # Passed in to modify/check if FR gets disabled
                                TREATMENT_ARMS          )   # Needed to check if all supply is missing
                        else
                            center_supplies, treatments_used, need_supply, delayed_patients, forward_treated, treatments_skipped, patients_force_allocated, num_patients, fr_allowed = F1b(
                                center_supplies, 
                                delayed_patients[j], 
                                i-tot_delayed-(length(delayed_patients)-num_delayed), 
                                treatments_used, 
                                blocks, 
                                need_supply, 
                                delayed_patients,
                                forward_treated,
                                treatments_skipped,
                                num_patients, 
                                critical_pt,
                                patients_force_allocated,
                                cap, 
                                TREATMENT_ARMS)
                        end
                        i+=1 # Increment overall index

                        if (i > num_patients)   # Check if SAMPLE_SIZE has been hit
                            patients[num_patients,2] = next_resupply    # Set end time equal to the time this supply arrived
                            break_loop = true
                            break
                        end
                    end

                    delayed_patients = delayed_patients[num_delayed+1:end]  # Any new delayed patients stay in the waitlist
                    next_resupply=next_supply_check+RESUPPLY_TIME           # Update next supply time
                    empty!(sent_supply)
                end

                if (break_loop) # If SAMPLE_SIZE was hit during delayed patients allocation, break out of outside loop
                    break
                end


                # Perform simulation as normal
                if (!fr_allowed) # If forced randomization isn't allowed
                    center_supplies, treatments_used, need_supply, patients_sent_home, num_patients = F0a(
                        center_supplies, 
                        center, 
                        i+treatments_skipped-patients_sent_home-tot_delayed-length(delayed_patients), # Treatment index, treatments_skipped=0 if this is the initial specification
                        treatments_used, 
                        blocks, 
                        need_supply, 
                        delayed_patients,
                        patients_sent_home,
                        num_patients, 
                        critical_pt,
                        TREATMENT_ARMS)
                else # If forced randomization is allowed, with cap=0
                    if (cap==0)
                        center_supplies, treatments_used, need_supply, patients_sent_home, num_patients = F0b(
                            center_supplies, 
                            center, 
                            i-patients_sent_home-tot_delayed-length(delayed_patients),      # i increases on patient sent home, but treatment index shouldn't.
                            treatments_used, 
                            blocks, 
                            need_supply,
                            delayed_patients, 
                            patients_sent_home, 
                            num_patients, 
                            critical_pt,
                            TREATMENT_ARMS)
                    elseif (!backfill_enabled) # If forced randomization is allowed, with no backfilling
                        center_supplies, treatments_used, need_supply, delayed_patients, treatments_skipped, patients_force_allocated, num_patients, fr_allowed = F1a(
                            center_supplies, 
                            center, 
                            i-tot_delayed-length(delayed_patients), # i increases once when a "delayed patient" is hit, and once again when the delayed patient is allocated.
                            treatments_used, 
                            blocks, 
                            need_supply, 
                            delayed_patients,
                            treatments_skipped,
                            num_patients, 
                            critical_pt,
                            patients_force_allocated,
                            cap, 
                            TREATMENT_ARMS)
                    else # If forced randomization is allowed, with backfilling
                        center_supplies, treatments_used, need_supply, delayed_patients, forward_treated, treatments_skipped, patients_force_allocated, num_patients, fr_allowed = F1b(
                            center_supplies, 
                            center, 
                            i-tot_delayed-length(delayed_patients), 
                            treatments_used, 
                            blocks, 
                            need_supply, 
                            delayed_patients,
                            forward_treated,
                            treatments_skipped,
                            num_patients, 
                            critical_pt,
                            patients_force_allocated,
                            cap, 
                            TREATMENT_ARMS)
                    end

                end

                i+=1
            end
            
            # Print out relevant stats for that scenario/simulation, helpful for debuggging
            total_time = patients[num_patients-tot_delayed,2]
            @printf("\nGen: %i Scenario: %i\n", sim, scenario)
            @printf("Treatment Imbalance: %i\n", abs(treatments_used[1]-treatments_used[2]))
            @printf("Patients Sent Home: %i\n", patients_sent_home)
            @printf("Patients FA: %i\n", patients_force_allocated)
            @printf("Unused Drugs: %f\n", Float64(total_drugs-SAMPLE_SIZE)/SAMPLE_SIZE)
            @printf("Total Time: %f\n",total_time)
            @printf("Total Cost: %i\n", total_cost)
            @printf("Patients left at end of trial: %i\n", length(delayed_patients))
            @printf("Total patients waitlisted: %i\n", num_waitlisted+length(delayed_patients))

            # Record all characteristics of this simulation
            characteristics[scenario,1,sim]=abs(treatments_used[1]-treatments_used[2])
            characteristics[scenario,2,sim]=patients_force_allocated/SAMPLE_SIZE
            characteristics[scenario,3,sim]=patients_sent_home/SAMPLE_SIZE
            characteristics[scenario,4,sim]=(total_drugs-SAMPLE_SIZE)/SAMPLE_SIZE
            characteristics[scenario,5,sim]=total_cost
            characteristics[scenario,6,sim]=total_time
            characteristics[scenario,7,sim]=length(delayed_patients)
            characteristics[scenario,8,sim]=num_waitlisted+length(delayed_patients)

            println(num_patients, " ",treatments_used, " ", delayed_patients)   # Used for debugging number of treatments allocated/treatment imbalance


        end

    end

    # Export stats in csv files
    export_time_stats_csv(characteristics[:,6,:], string("time-stats-b_",BETA,"-n",SAMPLE_SIZE,".csv"))
    export_time_stats_csv(characteristics[:,2,:], string("FA-stats-b_",BETA,"-n",SAMPLE_SIZE,".csv"))
    export_csv_stat(characteristics, heading, string("mean-characteristics-b_",BETA,"-n",SAMPLE_SIZE,".csv"), mean)
    export_csv_stat(characteristics, heading, string("std-characteristics-b_",BETA,"-n",SAMPLE_SIZE,".csv"), std)
    export_csv_stat(characteristics, heading, string("q90-characteristics-b_",BETA,"-n",SAMPLE_SIZE,".csv"), q90)


    # Plot each bar graph individually
    scenario_labels = repeat(["FR0a", "FR0b", "FR1a", "FR1b"], inner=3)

    plt = bar_with_err("Fig.4a: Treatment Imbalance",(BETA==28 ? "Base Case " : "Faster Recruitment ")*string("(α=",ALPHA)*string(", β=",BETA)*")", 
    scenario_labels, "Imbalance (num treatments)", mean(characteristics[:,1,:],dims=2), 1.645.*std(characteristics[:,1,:],dims=2))
    plts[ 1, bi ] = plt
    savefig(plt, string("Treatment-Imbalance-Plot-b_",BETA,"-n",SAMPLE_SIZE,".png"))

    plt = bar_with_err("Fig.4b: Forced Allocations",(BETA==28 ? "Base Case " : "Faster Recruitment ")*string("(α=",ALPHA)*string(", β=",BETA)*")", 
    scenario_labels, "Proportion of FAs", mean(characteristics[:,2,:],dims=2), 1.645.*std(characteristics[:,2,:],dims=2))
    plts[ 2, bi ] = plt
    savefig(plt, string("Forced-Allocations-Plot-b_",BETA,"-n",SAMPLE_SIZE,".png"))

    plt = bar_with_err("Fig.4c: Patients Sent Home",(BETA==28 ? "Base Case " : "Faster Recruitment ")*string("(α=",ALPHA)*string(", β=",BETA)*")", 
    scenario_labels, "Proportion of Patients", mean(characteristics[:,3,:],dims=2), 1.645.*std(characteristics[:,3,:],dims=2))
    plts[ 3, bi ] = plt
    savefig(plt, string("Patients-Sent-Home-Plot-b_",BETA,"-n",SAMPLE_SIZE,".png"))

    plt = bar_with_err("Fig.4d: Drug Overage",(BETA==28 ? "Base Case " : "Faster Recruitment ")*string("(α=",ALPHA)*string(", β=",BETA)*")",
    scenario_labels, "Overage (Pct)", mean(characteristics[:,4,:],dims=2), 1.645.*std(characteristics[:,4,:],dims=2))
    plts[ 4, bi ] = plt
    savefig(plt, string("Drug-Overage-Plot-b_",BETA,"-n",SAMPLE_SIZE,".png"))

    plt = bar_with_err("Fig.4x: Total Cost",(BETA==28 ? "Base Case " : "Faster Recruitment ")*string("(α=",ALPHA)*string("β=",BETA)*")", scenario_labels, "Average Cost (USD)", mean(characteristics[:,5,:],dims=2), 
    1.645.*std(characteristics[:,5,:],dims=2))
    plts[ 5, bi ] = plt
    savefig(plt, string("Total-Cost-Plot-b_",BETA,"-n",SAMPLE_SIZE,".png"))

    plt = bar_with_err("Fig.4e: Time to Complete Recruitment",(BETA==28 ? "Base Case " : "Faster Recruitment ")*string("(α=",ALPHA)*string(", β=",BETA)*")",
    scenario_labels, "Average Time (days)", mean(characteristics[:,6,:],dims=2), 1.645.*std(characteristics[:,6,:],dims=2))
    plts[ 6, bi ] = plt
    savefig(plt, string("Total-Time-Plot-b_",BETA,"-n",SAMPLE_SIZE,".png"))

    plt = bar_with_err("Fig.4f: Patients Not Allocated",(BETA==28 ? "Base Case " : "Faster Recruitment ")*string("(α=",ALPHA)*string(", β=",BETA)*")",
    scenario_labels, "Number of patients", mean(characteristics[:,7,:],dims=2), 1.645.*std(characteristics[:,7,:],dims=2))
    plts[ 7, bi ] = plt
    savefig(plt, string("Patients-NA-Plot-b_",BETA,"-n",SAMPLE_SIZE,".png"))

    plt = bar_with_err("Fig.4g: Patients Waitlisted",(BETA==28 ? "Base Case " : "Faster Recruitment ")*string("(α=",ALPHA)*string(", β=",BETA)*")",
    scenario_labels, "Number of patients", mean(characteristics[:,8,:],dims=2), 1.645.*std(characteristics[:,8,:],dims=2))
    plts[ 8, bi ] = plt
    savefig(plt, string("Patients-Waitlisted-Plot-b_",BETA,"-n",SAMPLE_SIZE,".png"))

    # Plot 3x4 time histograms
    time_hists(characteristics[:,6,:], string("Time-histograms-b_",BETA,"-n",SAMPLE_SIZE,".png"))

end

# Plot combined graphs for each beta tested
if(length(eachcol(plts))>1)
    for (ri,r) ∈ enumerate(eachrow(plts))
        plt_combined = plot(r..., layout = (1,length(r)))
        savefig(plt_combined, string("Combined-Plot-",heading[ri+2],".png"))
    end

end


