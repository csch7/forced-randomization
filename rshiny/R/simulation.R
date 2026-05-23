# Simulation state is stored in an R environment for reference semantics,
# mirroring Julia's mutable struct. All allocate_* functions modify S in place.

new_simulation <- function(ns, sample_size, centers, init_supply, critical_pt, cap, treatment_blocks) {
  S <- new.env(parent = emptyenv())

  max_pos <- ncol(treatment_blocks)

  # Treatment sequences: [stratum × position], pre-allocated to max sample_size
  S$treatments_used   <- matrix(0L, nrow = ns, ncol = sample_size)
  S$treatments_used_n <- integer(ns)   # per-stratum fill counters

  # Block skip offsets per stratum (F1a/F1b)
  S$treatments_skipped <- integer(ns)

  # Centers needing resupply — logical vector for O(1) add/check/clear
  S$need_supply <- logical(centers)

  # Delayed patient queues: [stratum × position], pre-allocated
  S$delayed_patients   <- matrix(0L, nrow = ns, ncol = sample_size)
  S$delayed_patients_n <- integer(ns)

  # Outcome counters
  S$patients_sent_home       <- integer(ns)
  S$patients_force_allocated <- integer(ns)

  # Forward-allocated block positions per stratum (F1b): logical for O(1) lookup
  S$forward_treated <- matrix(FALSE, nrow = ns, ncol = max_pos)

  S$num_patients <- as.integer(sample_size)

  # Supply matrix: treatment_arms × centers (column = one center's supply vector)
  S$supplies <- matrix(
    rep(as.integer(init_supply), centers),
    nrow = length(init_supply),
    ncol = centers
  )

  S$critical_point   <- as.integer(critical_pt)
  S$cap              <- as.integer(cap)
  S$treatment_blocks <- treatment_blocks   # ns × sample_size integer matrix

  S
}


# ── F0a: allocate only when all arms available; send home if any arm is zero ──
allocate_f0a <- function(S, center, stratum, treatment_index) {
  pt <- S$treatment_blocks[stratum, treatment_index]

  if (!any(S$supplies[, center] == 0L)) {
    S$supplies[pt, center] <- S$supplies[pt, center] - 1L
    n <- S$treatments_used_n[stratum] + 1L
    S$treatments_used[stratum, n] <- pt
    S$treatments_used_n[stratum]  <- n
    if (S$supplies[pt, center] <= S$critical_point) {
      S$need_supply[center] <- TRUE
    }

  } else if (all(S$supplies[, center] == 0L)) {
    nd <- S$delayed_patients_n[stratum] + 1L
    S$delayed_patients[stratum, nd] <- center
    S$delayed_patients_n[stratum]   <- nd
    S$num_patients <- S$num_patients + 1L

  } else {
    S$patients_sent_home[stratum] <- S$patients_sent_home[stratum] + 1L
    S$num_patients <- S$num_patients + 1L
  }
}


# ── F0b: send home only if patient's specific arm is zero ─────────────────────
allocate_f0b <- function(S, center, stratum, treatment_index) {
  pt <- S$treatment_blocks[stratum, treatment_index]

  if (S$supplies[pt, center] != 0L) {
    S$supplies[pt, center] <- S$supplies[pt, center] - 1L
    n <- S$treatments_used_n[stratum] + 1L
    S$treatments_used[stratum, n] <- pt
    S$treatments_used_n[stratum]  <- n
    if (S$supplies[pt, center] <= S$critical_point) {
      S$need_supply[center] <- TRUE
    }

  } else if (all(S$supplies[, center] == 0L)) {
    nd <- S$delayed_patients_n[stratum] + 1L
    S$delayed_patients[stratum, nd] <- center
    S$delayed_patients_n[stratum]   <- nd
    S$num_patients <- S$num_patients + 1L

  } else {
    S$patients_sent_home[stratum] <- S$patients_sent_home[stratum] + 1L
    S$num_patients <- S$num_patients + 1L
  }
}


# ── F1a: forced randomization without backfill ────────────────────────────────
# Returns TRUE if the FR cap has not yet been exhausted.
allocate_f1a <- function(S, center, stratum, treatment_index) {
  pt <- S$treatment_blocks[stratum, treatment_index + S$treatments_skipped[stratum]]

  if (S$supplies[pt, center] != 0L) {
    S$supplies[pt, center] <- S$supplies[pt, center] - 1L
    n <- S$treatments_used_n[stratum] + 1L
    S$treatments_used[stratum, n] <- pt
    S$treatments_used_n[stratum]  <- n
    if (S$supplies[pt, center] <= S$critical_point) {
      S$need_supply[center] <- TRUE
    }

  } else if (all(S$supplies[, center] == 0L)) {
    nd <- S$delayed_patients_n[stratum] + 1L
    S$delayed_patients[stratum, nd] <- center
    S$delayed_patients_n[stratum]   <- nd
    S$num_patients <- S$num_patients + 1L

  } else {
    repeat {
      S$treatments_skipped[stratum] <- S$treatments_skipped[stratum] + 1L
      pt <- S$treatment_blocks[stratum, treatment_index + S$treatments_skipped[stratum]]
      if (S$supplies[pt, center] != 0L) {
        S$supplies[pt, center] <- S$supplies[pt, center] - 1L
        n <- S$treatments_used_n[stratum] + 1L
        S$treatments_used[stratum, n] <- pt
        S$treatments_used_n[stratum]  <- n
        if (S$supplies[pt, center] <= S$critical_point) {
          S$need_supply[center] <- TRUE
        }
        break
      }
    }
    S$patients_force_allocated[stratum] <- S$patients_force_allocated[stratum] + 1L
    S$cap <- S$cap - 1L
  }

  S$cap > 0L
}


# ── DLG helper ────────────────────────────────────────────────────────────────
# Distance to Last Gap: furthest forward-allocated position minus first unfilled gap.
compute_dlg <- function(S, stratum, stratum_index) {
  if (!any(S$forward_treated[stratum, ])) return(0L)
  ts       <- S$treatments_skipped[stratum]
  next_pos <- stratum_index + 1L
  while (S$forward_treated[stratum, next_pos + ts]) {
    ts <- ts + 1L
  }
  max(0L, max(which(S$forward_treated[stratum, ])) - (next_pos + ts))
}


# ── F1b: forced randomization with backfill ───────────────────────────────────
# Returns TRUE if the FR cap has not yet been exhausted.
allocate_f1b <- function(S, center, stratum, treatment_index) {
  # Skip any positions already forward-allocated
  while (S$forward_treated[stratum, treatment_index + S$treatments_skipped[stratum]]) {
    S$treatments_skipped[stratum] <- S$treatments_skipped[stratum] + 1L
  }

  pt <- S$treatment_blocks[stratum, treatment_index + S$treatments_skipped[stratum]]

  if (S$supplies[pt, center] != 0L) {
    S$supplies[pt, center] <- S$supplies[pt, center] - 1L
    n <- S$treatments_used_n[stratum] + 1L
    S$treatments_used[stratum, n] <- pt
    S$treatments_used_n[stratum]  <- n
    if (S$supplies[pt, center] <= S$critical_point) {
      S$need_supply[center] <- TRUE
    }

  } else if (all(S$supplies[, center] == 0L)) {
    nd <- S$delayed_patients_n[stratum] + 1L
    S$delayed_patients[stratum, nd] <- center
    S$delayed_patients_n[stratum]   <- nd
    S$num_patients <- S$num_patients + 1L

  } else {
    # Back up one — current slot will be backfilled — then search forward
    S$treatments_skipped[stratum] <- S$treatments_skipped[stratum] - 1L
    j <- 1L

    repeat {
      while (S$forward_treated[stratum, treatment_index + S$treatments_skipped[stratum] + j]) {
        j <- j + 1L
      }
      forward_pos <- treatment_index + S$treatments_skipped[stratum] + j
      forward_pt  <- S$treatment_blocks[stratum, forward_pos]

      if (S$supplies[forward_pt, center] != 0L) {
        S$supplies[forward_pt, center] <- S$supplies[forward_pt, center] - 1L
        n <- S$treatments_used_n[stratum] + 1L
        S$treatments_used[stratum, n]      <- forward_pt
        S$treatments_used_n[stratum]       <- n
        S$forward_treated[stratum, forward_pos] <- TRUE
        if (S$supplies[forward_pt, center] <= S$critical_point) {
          S$need_supply[center] <- TRUE
        }
        break
      }
      j <- j + 1L
    }

    S$patients_force_allocated[stratum] <- S$patients_force_allocated[stratum] + 1L
    S$cap <- S$cap - 1L
  }

  S$cap > 0L
}
