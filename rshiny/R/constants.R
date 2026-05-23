default_sim_params <- function() {
  ss <- 500L
  list(
    # Trial design
    sample_size          = ss,
    treatment_arms       = 2L,
    allocation_ratio     = c(1L, 1L),
    block_size           = 4L,
    centers              = 80L,

    # Resupply logistics
    resupply_period      = 7L,
    resupply_time        = 3L,
    kit_cost             = 10000,
    ship_cost            = 100,

    # Forced randomization
    initial_cap          = as.integer(0.5 * ss),

    # Recruitment rate — Gamma(alpha, 1/beta)
    alpha                = 1.2,
    beta_options         = c(16),

    # Stratification
    num_strata           = 2L,
    strata_probabilities = c(0.6, 0.4),

    # Simulation
    number_simulations   = 2000L,

    # Supply strategies
    low_resupply         = 2L,
    low_init             = c(2L, 2L),
    low_critical         = 1L,
    med_resupply         = 4L,
    med_init             = c(3L, 3L),
    med_critical         = 1L,
    high_resupply        = 5L,
    high_init            = c(4L, 4L),
    high_critical        = 2L
  )
}

SCENARIOS_TO_RUN <- c(7L, 10L)

SCENARIO_MAP <- c(
  "1"  = "F0a Low",  "2"  = "F0a Med",  "3"  = "F0a High",
  "4"  = "F0b Low",  "5"  = "F0b Med",  "6"  = "F0b High",
  "7"  = "F1a Low",  "8"  = "F1a Med",  "9"  = "F1a High",
  "10" = "F1b Low",  "11" = "F1b Med",  "12" = "F1b High"
)
