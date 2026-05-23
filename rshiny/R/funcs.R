# Generate patient arrivals via inverse-transform exponential sampling.
# Returns an N×3 matrix: [patient_id, arrival_time, center_id]
generate_patient_arrivals <- function(rates, center_activations, num_centers, num_patients) {
  blocks <- vector("list", num_centers)
  for (i in seq_len(num_centers)) {
    u              <- runif(num_patients)
    inter_arrivals <- -log(1 - u) / rates[i]
    arrival_times  <- cumsum(inter_arrivals) + center_activations[i] * 30
    blocks[[i]]    <- cbind(
      (i - 1L) * num_patients + seq_len(num_patients),
      arrival_times,
      rep(i, num_patients)
    )
  }
  all_arrivals <- do.call(rbind, blocks)
  all_arrivals[order(all_arrivals[, 2L]), ]
}

# Assign each patient to a stratum. Returns an integer vector of stratum indices.
# Matches Julia's findfirst(>=(u), cumprobs) logic.
generate_strata_assignments <- function(num_patients, probs) {
  cum_probs <- cumsum(probs)
  u         <- runif(num_patients)
  # rowSums counts how many cumprobs u exceeds; +1 gives 1-based stratum index
  as.integer(rowSums(outer(u, cum_probs, ">")) + 1L)
}

# Generate a permuted block treatment sequence.
# Returns a flat integer vector of length num_strata * sample_size.
# Caller reshapes: t(matrix(raw, nrow = sample_size, ncol = num_strata))
generate_treatment_blocks <- function(ratio, sample_size, treat_arms, block_size, num_strata) {
  block_template <- rep(seq_len(treat_arms), times = ratio * (block_size %/% sum(ratio)))
  n_blocks       <- ceiling(num_strata * sample_size / block_size)
  unlist(replicate(n_blocks, sample(block_template), simplify = FALSE))
}

# Return supply and FR parameters for a given scenario number (1-12).
scenario_params <- function(scenario, params) {
  supply_configs <- list(
    list(resupply_amt = params$low_resupply,  init_supply = params$low_init,  critical_pt = params$low_critical),
    list(resupply_amt = params$med_resupply,  init_supply = params$med_init,  critical_pt = params$med_critical),
    list(resupply_amt = params$high_resupply, init_supply = params$high_init, critical_pt = params$high_critical)
  )
  fr_configs <- list(
    list(fr_allowed = FALSE, backfill_enabled = FALSE, cap = 0L),
    list(fr_allowed = TRUE,  backfill_enabled = FALSE, cap = 0L),
    list(fr_allowed = TRUE,  backfill_enabled = FALSE, cap = params$initial_cap),
    list(fr_allowed = TRUE,  backfill_enabled = TRUE,  cap = params$initial_cap)
  )
  supply_idx <- (scenario - 1L) %% 3L + 1L
  fr_idx     <- (scenario - 1L) %/% 3L + 1L
  c(supply_configs[[supply_idx]], fr_configs[[fr_idx]])
}

# 90th quantile approximation: mean + 1.645 * sd, over the first dimension of an array.
quantile_90 <- function(x, dims) {
  apply(x, dims, mean) + 1.645 * apply(x, dims, sd)
}
