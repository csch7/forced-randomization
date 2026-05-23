# Validation script: run a quick simulation and print key metrics.
# Compare output against Julia with the same parameters (seeds won't match —
# compare distribution summaries: mean FA rate, mean imbalance, etc.)
#
# Run from RStudio: open this file, then Session > Set Working Directory >
# To Source File Location, then source("validate.R") or press Ctrl+Shift+S.

# Resolve paths relative to this script regardless of working directory
this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

source(file.path(this_dir, "R/constants.R"))
source(file.path(this_dir, "R/funcs.R"))
source(file.path(this_dir, "R/simulation.R"))
source(file.path(this_dir, "R/run_simulation.R"))

set.seed(42)

params <- default_sim_params()

# Small run for fast validation
params$sample_size         <- 500L
params$number_simulations  <- 200L
params$initial_cap         <- as.integer(0.5 * params$sample_size)
params$beta_options        <- c(3)   # match Julia BETA=3 exploratory runs

cat("Running validation: n =", params$sample_size,
    "| sims =", params$number_simulations, "\n")

t0      <- proc.time()
results <- run_simulation(params)
elapsed <- (proc.time() - t0)[["elapsed"]]

cat(sprintf("Elapsed: %.1f sec\n\n", elapsed))

r  <- results[[1]]
ns <- params$num_strata

cat("=== F1a Low (scenario 7) ===\n")
for (k in seq_len(ns)) {
  fa_mean <- mean(r$characteristics[1L, 1L, k, ])
  cat(sprintf("  Mean FA rate  (z=%d): %.4f\n", k, fa_mean))
}
cat(sprintf("  Mean FA rate  (total): %.4f\n", mean(r$characteristics[1L, 1L, ns + 1L, ])))
for (k in seq_len(ns)) {
  dm_end <- r$d500s[1L, , k]
  cat(sprintf("  Mean end imbalance (z=%d): %.4f  |  Var: %.4f\n",
              k, mean(dm_end), var(dm_end)))
}
for (k in seq_len(ns)) {
  cat(sprintf("  Avg slots skipped/FA (z=%d): %.4f\n",
              k, mean(r$avg_slots_skipped_f1a[k, ])))
}

cat("\n=== F1b Low (scenario 10) ===\n")
for (k in seq_len(ns)) {
  fa_mean <- mean(r$characteristics[2L, 1L, k, ])
  cat(sprintf("  Mean FA rate  (z=%d): %.4f\n", k, fa_mean))
}
cat(sprintf("  Mean FA rate  (total): %.4f\n", mean(r$characteristics[2L, 1L, ns + 1L, ])))
for (k in seq_len(ns)) {
  dm_end <- r$d500s[2L, , k]
  cat(sprintf("  Mean end imbalance (z=%d): %.4f  |  Var: %.4f\n",
              k, mean(dm_end), var(dm_end)))
}

cat("\n=== Recruitment times ===\n")
cat(sprintf("  F1a: mean=%.1f  F1b: mean=%.1f\n",
            mean(r$recruitment_times[1L, ]), mean(r$recruitment_times[2L, ])))
