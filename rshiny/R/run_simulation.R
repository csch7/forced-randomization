run_simulation <- function(params, progress_fn = NULL) {
  all_results <- list()

  for (BETA in params$beta_options) {
    ns  <- params$num_strata
    nsc <- length(SCENARIOS_TO_RUN)

    max_z     <- round(params$strata_probabilities * params$sample_size)
    max_z[ns] <- params$sample_size - sum(max_z[-ns])
    max_z     <- as.integer(max_z)

    # [scenario_offset, characteristic=1, stratum_or_total, simulation]
    characteristics <- array(0.0, dim = c(nsc, 1L, ns + 1L, params$number_simulations))

    # Per-stratum imbalance/DLG arrays
    dms  <- lapply(seq_len(ns), function(k) array(NaN, dim = c(nsc, params$number_simulations, max_z[k])))
    dlgs <- lapply(seq_len(ns), function(k) array(NaN, dim = c(nsc, params$number_simulations, max_z[k])))

    # [scenario, simulation, stratum]
    d500s <- array(0.0, dim = c(nsc, params$number_simulations, ns))

    # F1a-only skip tracking
    unfilled_slots        <- lapply(seq_len(ns), function(k) array(0.0, dim = c(1L, params$number_simulations, max_z[k])))
    avg_slots_skipped_f1a <- matrix(0.0, nrow = ns, ncol = params$number_simulations)

    # [scenario_offset, simulation]
    recruitment_times <- matrix(0.0, nrow = nsc, ncol = params$number_simulations)

    for (sim in seq_len(params$number_simulations)) {

      center_rates <- rgamma(params$centers, shape = params$alpha, rate = BETA)
      center_acts  <- sample(0L:4L, params$centers, replace = TRUE)
      patients     <- generate_patient_arrivals(center_rates, center_acts, params$centers, params$sample_size)

      raw_blocks       <- generate_treatment_blocks(params$allocation_ratio, params$sample_size,
                                                    params$treatment_arms, params$block_size, ns)
      treatment_blocks <- t(matrix(raw_blocks, nrow = params$sample_size, ncol = ns))
      mode(treatment_blocks) <- "integer"

      strata <- generate_strata_assignments(params$sample_size * 2L, params$strata_probabilities)

      for (si in seq_along(SCENARIOS_TO_RUN)) {
        scenario <- SCENARIOS_TO_RUN[si]
        sp       <- scenario_params(scenario, params)

        total_drugs <- params$centers * sum(sp$init_supply)
        total_cost  <- params$centers * sum(sp$init_supply) * params$kit_cost

        S <- new_simulation(ns, params$sample_size, params$centers,
                            sp$init_supply, sp$critical_pt, sp$cap, treatment_blocks)

        next_supply_check    <- params$resupply_period
        next_resupply        <- next_supply_check + params$resupply_time
        sent_supply          <- list()   # character(center_id) -> integer supply vector
        tot_delayed          <- integer(ns)
        num_waitlisted       <- integer(ns)
        patients_per_stratum <- integer(ns)
        fr_allowed           <- sp$fr_allowed
        cap                  <- sp$cap
        break_loop           <- FALSE

        i <- 1L
        while (i <= S$num_patients) {
          center <- as.integer(patients[i, 3L])

          # ── Supply check ────────────────────────────────────────────────────
          if (patients[i, 2L] >= next_supply_check) {
            for (j in which(S$need_supply)) {
              new_supply <- integer(params$treatment_arms)
              for (k in seq_len(params$treatment_arms)) {
                if (S$supplies[k, j] <= sp$critical_pt) {
                  top_up      <- sp$resupply_amt - S$supplies[k, j]
                  total_cost  <- total_cost + top_up * params$kit_cost
                  total_drugs <- total_drugs + top_up
                  new_supply[k] <- sp$resupply_amt
                }
              }
              sent_supply[[as.character(j)]] <- new_supply
              total_cost <- total_cost + params$ship_cost
            }
            next_supply_check  <- next_supply_check + params$resupply_period
            S$need_supply[]    <- FALSE
          }

          # ── Resupply arrival ─────────────────────────────────────────────
          if (patients[i, 2L] >= next_resupply) {
            for (cent_name in names(sent_supply)) {
              cent   <- as.integer(cent_name)
              ns_vec <- sent_supply[[cent_name]]
              for (k in seq_along(ns_vec)) {
                if (ns_vec[k] != 0L) S$supplies[k, cent] <- ns_vec[k]
              }
            }

            num_delayed  <- S$delayed_patients_n
            tot_delayed  <- tot_delayed + num_delayed

            for (z in seq_len(ns)) {
              if (num_delayed[z] > 0L) {
                dp  <- S$delayed_patients[z, seq_len(S$delayed_patients_n[z])]
                cts <- tabulate(dp, nbins = params$centers)
                for (cnt in cts[cts > 0L]) {
                  num_waitlisted[z] <- num_waitlisted[z] + min(cnt, sp$resupply_amt * params$treatment_arms)
                }
              }
            }

            # Flatten delayed queues across strata
            flattened_delayed <- unlist(lapply(seq_len(ns), function(z) {
              if (S$delayed_patients_n[z] > 0L)
                S$delayed_patients[z, seq_len(S$delayed_patients_n[z])]
              else
                integer(0L)
            }))

            for (j_idx in seq_along(flattened_delayed)) {
              new_delayed   <- sum(S$delayed_patients_n) - sum(num_delayed)
              patient_index <- i - sum(S$patients_sent_home) - sum(tot_delayed) - new_delayed
              stratum       <- strata[patient_index]
              center_delayed <- as.integer(flattened_delayed[j_idx])
              other_strata_patients <- sum(patients_per_stratum) - patients_per_stratum[stratum]
              stratum_index <- as.integer(trunc(
                (i - other_strata_patients) -
                S$patients_sent_home[stratum] -
                tot_delayed[stratum] -
                (S$delayed_patients_n[stratum] - num_delayed[stratum])
              ))

              if (!fr_allowed) {
                allocate_f0a(S, center_delayed, stratum, stratum_index)
              } else if (cap == 0L) {
                allocate_f0b(S, center_delayed, stratum, stratum_index)
              } else if (!sp$backfill_enabled) {
                fr_allowed <- allocate_f1a(S, center_delayed, stratum, stratum_index)
                unfilled_slots[[stratum]][1L, sim, min(max_z[stratum], stratum_index)] <-
                  S$treatments_skipped[stratum]
              } else {
                fr_allowed <- allocate_f1b(S, center_delayed, stratum, stratum_index)
              }

              patients_per_stratum[stratum] <- patients_per_stratum[stratum] + 1L
              dlgs[[stratum]][si, sim, min(max_z[stratum], stratum_index)] <-
                compute_dlg(S, stratum, stratum_index)

              i <- i + 1L
              if (i > S$num_patients) {
                # Clamp to valid row range — S$num_patients can exceed sample_size
                patients[min(S$num_patients, nrow(patients)), 2L] <- next_resupply
                break_loop <- TRUE
                break
              }
            }

            # Trim processed entries from delayed queues
            for (z in seq_len(ns)) {
              nd        <- num_delayed[z]
              remaining <- S$delayed_patients_n[z] - nd
              if (remaining > 0L) {
                S$delayed_patients[z, seq_len(remaining)] <-
                  S$delayed_patients[z, (nd + 1L):S$delayed_patients_n[z]]
              }
              S$delayed_patients_n[z] <- remaining
            }

            next_resupply <- next_supply_check + params$resupply_time
            sent_supply   <- list()
          }

          if (break_loop) break

          # ── Normal patient allocation ────────────────────────────────────
          current_delayed <- sum(S$delayed_patients_n)
          patient_index   <- i + sum(S$treatments_skipped) - sum(S$patients_sent_home) -
                             sum(tot_delayed) - current_delayed
          stratum         <- strata[patient_index]
          other_strata_patients <- sum(patients_per_stratum) - patients_per_stratum[stratum]
          stratum_index   <- as.integer(trunc(
            (i - other_strata_patients) -
            tot_delayed[stratum] -
            S$delayed_patients_n[stratum]
          ))

          if (!fr_allowed) {
            allocate_f0a(S, center, stratum, stratum_index)
          } else if (cap == 0L) {
            allocate_f0b(S, center, stratum, stratum_index)
          } else if (!sp$backfill_enabled) {
            fr_allowed <- allocate_f1a(S, center, stratum, stratum_index)
            unfilled_slots[[stratum]][1L, sim, min(max_z[stratum], stratum_index)] <-
              S$treatments_skipped[stratum]
          } else {
            fr_allowed <- allocate_f1b(S, center, stratum, stratum_index)
          }

          patients_per_stratum[stratum] <- patients_per_stratum[stratum] + 1L
          dlgs[[stratum]][si, sim, min(max_z[stratum], stratum_index)] <-
            compute_dlg(S, stratum, stratum_index)

          i <- i + 1L
        }  # while patients

        recruitment_times[si, sim] <- patients[params$sample_size, 2L]

        # ── Imbalance metrics (vectorized) ───────────────────────────────
        normalise <- (scenario == 7L)

        for (k in seq_len(ns)) {
          nk      <- S$treatments_used_n[k]
          trt     <- S$treatments_used[k, seq_len(nk)]
          t1_cs   <- cumsum(trt == 1L)
          imbal   <- as.double(2L * t1_cs - seq_len(nk))   # t1 - t2 = 2*t1 - n
          n_fill  <- min(nk, max_z[k])

          if (normalise) {
            dms[[k]][si, sim, seq_len(n_fill)] <- imbal[seq_len(n_fill)] / sqrt(seq_len(n_fill))
            d500s[si, sim, k]                  <- imbal[nk] / sqrt(nk)
          } else {
            dms[[k]][si, sim, seq_len(n_fill)] <- imbal[seq_len(n_fill)]
            d500s[si, sim, k]                  <- imbal[nk]
          }
        }

        for (k in seq_len(ns)) {
          characteristics[si, 1L, k, sim] <-
            S$patients_force_allocated[k] / S$treatments_used_n[k]
        }
        characteristics[si, 1L, ns + 1L, sim] <-
          sum(S$patients_force_allocated) / params$sample_size

        if (si == 1L) {   # F1a only
          for (z in seq_len(ns)) {
            pfa <- S$patients_force_allocated[z]
            avg_slots_skipped_f1a[z, sim] <-
              if (pfa > 0L) S$treatments_skipped[z] / pfa else 0.0
          }
        }

      }  # scenarios

      if (!is.null(progress_fn)) progress_fn(sim, params$number_simulations)
    }  # simulations

    all_results[[length(all_results) + 1L]] <- list(
      beta                  = BETA,
      dms                   = dms,
      d500s                 = d500s,
      characteristics       = characteristics,
      dlgs                  = dlgs,
      unfilled_slots        = unfilled_slots,
      recruitment_times     = recruitment_times,
      avg_slots_skipped_f1a = avg_slots_skipped_f1a,
      max_z                 = max_z
    )

  }  # beta

  all_results
}
