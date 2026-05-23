# quantile_90_vec: single-vector version of quantile_90 for use with apply()
quantile_90_vec <- function(x) mean(x, na.rm = TRUE) + 1.645 * sd(x, na.rm = TRUE)

# make_line_panel: compute a summary statistic over simulations at each patient
# position and return a ggplot line chart.
#
# dm_slice          : n_sims x n_patients matrix (one scenario, one stratum)
# stat_fn           : function applied column-wise, e.g. var, quantile_90_vec, max
# cutoff            : last patient position to include
# start             : first patient position to include (skip high-variance opening)
# stored_normalized : TRUE if values were divided by sqrt(n) at storage time (F1a)
# display_normalized: TRUE to show normalized values in the plot
make_line_panel <- function(dm_slice, label, stat_fn, cutoff, start = 1L,
                             stored_normalized = FALSE, display_normalized = TRUE) {
  m <- dm_slice[, seq_len(cutoff), drop = FALSE]
  if (stored_normalized && !display_normalized) {
    m <- sweep(m, 2L, sqrt(seq_len(cutoff)), `*`)
  } else if (!stored_normalized && display_normalized) {
    m <- sweep(m, 2L, sqrt(seq_len(cutoff)), `/`)
  }
  vals <- apply(m, 2L, stat_fn)
  df   <- data.frame(patient = start:cutoff, value = vals[start:cutoff])
  ggplot(df, aes(x = patient, y = value)) +
    geom_line(color = "#2c7bb6", linewidth = 0.8) +
    labs(x = "Patient", y = label) +
    theme_minimal(base_size = 11) +
    theme(plot.margin = margin(4, 8, 4, 8))
}

make_imbalance_histogram <- function(values, title) {
  df <- data.frame(value = values[!is.nan(values)])
  ggplot(df, aes(x = value)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", color = "white", alpha = 0.8) +
    labs(x = "Imbalance (Dm)", y = "Count", title = title) +
    theme_minimal(base_size = 11)
}

make_recruitment_histogram <- function(times, title) {
  df <- data.frame(time = times)
  ggplot(df, aes(x = time)) +
    geom_histogram(bins = 30, fill = "#d7191c", color = "white", alpha = 0.8) +
    labs(x = "Days to Full Enrollment", y = "Count", title = title) +
    theme_minimal(base_size = 11)
}

plot_joint_normality <- function(r, ns) {
  panels <- lapply(seq_len(ns), function(k) {
    y   <- sort(r$d500s[1L, , k])
    n   <- length(y)
    x   <- qnorm(ppoints(n))
    mu  <- mean(y)
    sig <- sd(y)
    df  <- data.frame(x = x, y = y)
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.6, color = "#2c7bb6", size = 1.5) +
      geom_abline(intercept = mu, slope = sig,
                  color = "#d7191c", linetype = "dashed", linewidth = 0.8) +
      labs(
        x     = "Theoretical N(0,1) quantiles",
        y     = "Empirical end-of-trial imbalance",
        title = sprintf("F1a Low — Stratum %d", k)
      ) +
      theme_minimal(base_size = 11)
  })
  wrap_plots(panels, ncol = ns)
}

build_summary_table <- function(r, ns) {
  scenario_labels <- c("F1a Low", "F1b Low")
  rows <- list()
  for (si in seq_along(scenario_labels)) {
    for (k in seq_len(ns)) {
      avg_skip <- if (si == 1L) sprintf("%.4f", mean(r$avg_slots_skipped_f1a[k, ])) else "—"
      rows[[length(rows) + 1L]] <- data.frame(
        Scenario                             = scenario_labels[si],
        Stratum                              = sprintf("z=%d", k),
        `Mean FA Rate`                       = sprintf("%.4f", mean(r$characteristics[si, 1L, k, ])),
        `Var(End Imbalance)`                 = sprintf("%.4f", var(r$d500s[si, , k])),
        `Avg Block Positions Skipped per FA` = avg_skip,
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }
    rows[[length(rows) + 1L]] <- data.frame(
      Scenario                             = scenario_labels[si],
      Stratum                              = "Total",
      `Mean FA Rate`                       = sprintf("%.4f", mean(r$characteristics[si, 1L, ns + 1L, ])),
      `Var(End Imbalance)`                 = "—",
      `Avg Block Positions Skipped per FA` = "—",
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
