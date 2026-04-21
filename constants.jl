# ── Trial Design ────────────────────────────────────────────────────────────
const SAMPLE_SIZE                   = 2000
const TREATMENT_ARMS                = 2
const ALLOCATION_RATIO              = (1, 1)
const BLOCK_SIZE                    = 4
const CENTERS                       = 80

# ── Resupply Logistics ───────────────────────────────────────────────────────
const RESUPPLY_PERIOD               = 7     # Days between resupply checks
const RESUPPLY_TIME                 = 3     # Days for shipment to arrive after order
const KIT_COST                      = 10_000
const SHIP_COST                     = 100

# ── Forced Randomization ─────────────────────────────────────────────────────
const INITIAL_CAP                   = Int(0.3 * SAMPLE_SIZE)

# ── Center Recruitment Rate Distribution (Gamma(α, 1/β)) ────────────────────
const ALPHA                         = 1.2
const BETA_OPTIONS                  = [16]

# ── Stratification ───────────────────────────────────────────────────────────
const STRATA_ASSIGNMENT_PROBABILITY = 0.6   # Probability of stratum 1 assignment

# ── Simulation ───────────────────────────────────────────────────────────────
const NUMBER_SIMULATIONS            = 1000

# ── Supply Strategy Parameters ───────────────────────────────────────────────
# Low supply
const LOW_RESUPPLY                  = 2
const LOW_INIT                      = [2, 2]
const LOW_CRITICAL                  = 1

# Medium supply
const MED_RESUPPLY                  = 4
const MED_INIT                      = [3, 3]
const MED_CRITICAL                  = 1

# High supply
const HIGH_RESUPPLY                 = 5
const HIGH_INIT                     = [4, 4]
const HIGH_CRITICAL                 = 2

const SCENARIO_MAP = Dict{Int, String}(
    1 => "F0a Low",
    2 => "F0a Med",
    3 => "F0a High",
    4 => "F0b Low",
    5 => "F0b Med",
    6 => "F0b High",
    7 => "F1a Low",
    8 => "F1a Med",
    9 => "F1a High",
    10 => "F1b Low",
    11 => "F1b Med",
    12 => "F1b High"
)