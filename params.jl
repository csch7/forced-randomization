include("constants.jl")

Base.@kwdef struct SimParams
    # Trial design
    sample_size::Int
    treatment_arms::Int
    allocation_ratio::Tuple{Int,Int}
    block_size::Int
    centers::Int

    # Resupply logistics
    resupply_period::Int
    resupply_time::Int
    kit_cost::Float64
    ship_cost::Float64

    # Forced randomization
    initial_cap::Int

    # Recruitment rate
    alpha::Float64
    beta_options::Vector{Float64}

    # Stratification
    strata_assignment_probability::Float64

    # Simulation
    number_simulations::Int

    # Supply strategies
    low_resupply::Int
    low_init::Vector{Int}
    low_critical::Int
    med_resupply::Int
    med_init::Vector{Int}
    med_critical::Int
    high_resupply::Int
    high_init::Vector{Int}
    high_critical::Int
end

function default_params()
    return SimParams(
        sample_size                   = SAMPLE_SIZE,
        treatment_arms                = TREATMENT_ARMS,
        allocation_ratio              = ALLOCATION_RATIO,
        block_size                    = BLOCK_SIZE,
        centers                       = CENTERS,
        resupply_period               = RESUPPLY_PERIOD,
        resupply_time                 = RESUPPLY_TIME,
        kit_cost                      = Float64(KIT_COST),
        ship_cost                     = Float64(SHIP_COST),
        initial_cap                   = INITIAL_CAP,
        alpha                         = ALPHA,
        beta_options                  = Vector{Float64}(BETA_OPTIONS),
        strata_assignment_probability = STRATA_ASSIGNMENT_PROBABILITY,
        number_simulations            = NUMBER_SIMULATIONS,
        low_resupply                  = LOW_RESUPPLY,
        low_init                      = copy(LOW_INIT),
        low_critical                  = LOW_CRITICAL,
        med_resupply                  = MED_RESUPPLY,
        med_init                      = copy(MED_INIT),
        med_critical                  = MED_CRITICAL,
        high_resupply                 = HIGH_RESUPPLY,
        high_init                     = copy(HIGH_INIT),
        high_critical                 = HIGH_CRITICAL,
    )
end
