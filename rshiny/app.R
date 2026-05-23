library(shiny)
library(bslib)
library(ggplot2)
library(patchwork)
library(shinycssloaders)
library(shinyjs)
library(callr)

# Build a SimParams list from Shiny input values.
params_from_inputs <- function(input) {
  p                      <- default_sim_params()
  p$sample_size          <- as.integer(input$sample_size)
  p$treatment_arms       <- as.integer(input$treatment_arms)
  p$allocation_ratio     <- c(as.integer(input$alloc_ratio_1), as.integer(input$alloc_ratio_2))
  p$block_size           <- as.integer(input$block_size)
  p$centers              <- as.integer(input$centers)
  p$resupply_period      <- as.integer(input$resupply_period)
  p$resupply_time        <- as.integer(input$resupply_time)
  p$kit_cost             <- input$kit_cost
  p$ship_cost            <- input$ship_cost
  p$initial_cap          <- as.integer(input$initial_cap)
  p$alpha                <- input$alpha
  p$beta_options         <- c(input$beta)
  p$strata_probabilities <- c(input$z1_prob, 1 - input$z1_prob)
  p$number_simulations   <- as.integer(input$number_simulations)
  p$run_label            <- ""
  p$low_resupply         <- as.integer(input$low_resupply)
  p$low_init             <- c(as.integer(input$low_init_1), as.integer(input$low_init_2))
  p$low_critical         <- as.integer(input$low_critical)
  p$med_resupply         <- as.integer(input$med_resupply)
  p$med_init             <- c(as.integer(input$med_init_1), as.integer(input$med_init_2))
  p$med_critical         <- as.integer(input$med_critical)
  p$high_resupply        <- as.integer(input$high_resupply)
  p$high_init            <- c(as.integer(input$high_init_1), as.integer(input$high_init_2))
  p$high_critical        <- as.integer(input$high_critical)
  p
}

d <- default_sim_params()   # default values for input initialisation

ui <- page_sidebar(
  title = "Forced Randomization Simulator",
  theme = bs_theme(version = 5),
  useShinyjs(),

  sidebar = sidebar(
    width = 320,

    accordion(
      open = FALSE,

      accordion_panel(
        "Simulation",
        numericInput("number_simulations", "Number of simulations",
                     value = d$number_simulations, min = 1L, step = 100L)
      ),

      accordion_panel(
        "Trial Design",
        numericInput("sample_size",    "Sample size",    value = d$sample_size,    min = 1L),
        numericInput("centers",        "Centers",        value = d$centers,        min = 1L),
        numericInput("treatment_arms", "Treatment arms", value = d$treatment_arms, min = 2L),
        tags$label(class = "form-label", "Allocation ratio (arm 1 : arm 2)"),
        layout_columns(
          col_widths = c(6, 6),
          numericInput("alloc_ratio_1", "Arm 1", value = d$allocation_ratio[1], min = 1L),
          numericInput("alloc_ratio_2", "Arm 2", value = d$allocation_ratio[2], min = 1L)
        ),
        numericInput("block_size", "Block size", value = d$block_size, min = 2L)
      ),

      accordion_panel(
        "Stratification",
        sliderInput("z1_prob", "Z1 probability",
                    min = 0.01, max = 0.99,
                    value = d$strata_probabilities[1], step = 0.01),
        uiOutput("z2_label")
      ),

      accordion_panel(
        "Recruitment",
        numericInput("alpha", "Alpha — Gamma shape",
                     value = d$alpha, min = 0.01, step = 0.1),
        numericInput("beta",  "Beta — Gamma rate",
                     value = d$beta_options[1], min = 0.01, step = 0.1)
      ),

      accordion_panel(
        "Forced Randomization (FR)",
        numericInput("initial_cap", "FR cap (max forced allocations)",
                     value = d$initial_cap, min = 0L)
      ),

      accordion_panel(
        "Resupply Logistics",
        numericInput("resupply_period", "Resupply check period (days)",
                     value = d$resupply_period, min = 1L),
        numericInput("resupply_time",   "Shipment transit time (days)",
                     value = d$resupply_time,   min = 1L),
        numericInput("kit_cost",  "Kit cost ($)",                value = d$kit_cost,  min = 0),
        numericInput("ship_cost", "Shipping cost per order ($)", value = d$ship_cost, min = 0)
      ),

      accordion_panel(
        "Supply Strategies",
        tags$strong("Low"),
        layout_columns(
          col_widths = c(6, 6),
          numericInput("low_init_1", "Initial — arm 1", value = d$low_init[1], min = 0L),
          numericInput("low_init_2", "Initial — arm 2", value = d$low_init[2], min = 0L)
        ),
        numericInput("low_resupply", "Resupply target",    value = d$low_resupply, min = 0L),
        numericInput("low_critical", "Critical threshold", value = d$low_critical, min = 0L),

        tags$hr(),
        tags$strong("Medium"),
        layout_columns(
          col_widths = c(6, 6),
          numericInput("med_init_1", "Initial — arm 1", value = d$med_init[1], min = 0L),
          numericInput("med_init_2", "Initial — arm 2", value = d$med_init[2], min = 0L)
        ),
        numericInput("med_resupply", "Resupply target",    value = d$med_resupply, min = 0L),
        numericInput("med_critical", "Critical threshold", value = d$med_critical, min = 0L),

        tags$hr(),
        tags$strong("High"),
        layout_columns(
          col_widths = c(6, 6),
          numericInput("high_init_1", "Initial — arm 1", value = d$high_init[1], min = 0L),
          numericInput("high_init_2", "Initial — arm 2", value = d$high_init[2], min = 0L)
        ),
        numericInput("high_resupply", "Resupply target",    value = d$high_resupply, min = 0L),
        numericInput("high_critical", "Critical threshold", value = d$high_critical, min = 0L)
      )
    ),

    actionButton("run_btn", "Run Simulation",
                 class = "btn-primary w-100 mt-3", icon = icon("play")),
    actionButton("cancel_btn", "Cancel / Reset",
                 class = "btn-outline-danger w-100 mt-1")
  ),

  withSpinner(uiOutput("main_ui"), type = 6)
)

server <- function(input, output, session) {
  rv      <- reactiveValues(results = NULL, ns = 2L,
                            running = FALSE,
                            plot_f1a = NULL, plot_f1b = NULL,
                            plot_dlg = NULL, plot_hist = NULL,
                            plot_normality = NULL)
  bg_proc <- reactiveVal(NULL)   # callr background process handle

  # Z2 read-only label beneath the Z1 slider
  output$z2_label <- renderUI({
    tags$p(sprintf("Z2 probability: %.2f", 1 - input$z1_prob),
           class = "text-muted small mt-n2")
  })

  # Launch simulation in a background R process so the UI stays responsive
  observeEvent(input$run_btn, {
    if (isTRUE(rv$running)) return()
    rv$running <- TRUE
    shinyjs::disable("run_btn")

    params <- params_from_inputs(input)
    wd     <- getwd()

    p <- callr::r_bg(
      func = function(params, wd) {
        setwd(wd)
        source("R/constants.R")
        source("R/funcs.R")
        source("R/simulation.R")
        source("R/run_simulation.R")
        run_simulation(params)
      },
      args       = list(params = params, wd = wd),
      supervise  = TRUE,
      stderr     = NULL
    )
    bg_proc(p)
  })

  # Poll the background process every 500 ms until it finishes
  observe({
    p <- bg_proc()
    req(!is.null(p))

    if (p$is_alive()) {
      invalidateLater(500)
      return()
    }

    # Process finished — harvest result
    bg_proc(NULL)
    rv$running <- FALSE
    shinyjs::enable("run_btn")

    result <- tryCatch(p$get_result(), error = function(e) NULL)
    if (!is.null(result)) {
      rv$results <- result
      rv$ns      <- length(result[[1L]]$max_z)
    }
  })

  # Cancel: kill the child process and reset all state
  observeEvent(input$cancel_btn, {
    p <- bg_proc()
    if (!is.null(p)) {
      if (p$is_alive()) p$kill()
      bg_proc(NULL)
    }
    rv$running        <- FALSE
    rv$results        <- NULL
    rv$ns             <- 2L
    rv$plot_f1a       <- NULL
    rv$plot_f1b       <- NULL
    rv$plot_dlg       <- NULL
    rv$plot_hist      <- NULL
    rv$plot_normality <- NULL
    shinyjs::enable("run_btn")
  })

  # Main panel: placeholder / running indicator / tabs
  output$main_ui <- renderUI({
    if (isTRUE(rv$running)) {
      div(class = "text-center p-5 text-muted",
          h4(icon("circle-notch", class = "fa-spin me-2"),
             "Simulation running…"),
          p("Press", strong("Cancel / Reset"), "to abort."))
    } else if (is.null(rv$results)) {
      div(class = "text-center p-5 text-muted",
          h4("Configure parameters and press", strong("Run Simulation"),
             "to generate results."))
    } else {
      navset_tab(
        id = "output_tabs",
        nav_panel("Imbalance Over Time",
          div(style = "text-align:right; margin-bottom:4px;",
              downloadButton("dl_f1a", "F1a",  class = "btn-sm btn-outline-secondary me-1"),
              downloadButton("dl_f1b", "F1b",  class = "btn-sm btn-outline-secondary")),
          radioButtons("normalize", NULL,
                       choices = c("Normalized (÷ √n)" = "norm", "Raw" = "raw"),
                       selected = "norm", inline = TRUE),
          plotOutput("plot_imbalance_f1a", height = "400px"),
          plotOutput("plot_imbalance_f1b", height = "400px")
        ),
        nav_panel("Distance to Last Gap (DLG)",
          div(style = "text-align:right; margin-bottom:4px;",
              downloadButton("dl_dlg", "Download", class = "btn-sm btn-outline-secondary")),
          plotOutput("plot_dlg", height = "400px")
        ),
        nav_panel("End-of-Trial Imbalance",
          div(style = "text-align:right; margin-bottom:4px;",
              downloadButton("dl_hist", "Download", class = "btn-sm btn-outline-secondary")),
          plotOutput("plot_hist_imbalance", height = "400px")
        ),
        nav_panel("Summary Statistics",
          div(style = "text-align:right; margin-bottom:4px;",
              downloadButton("dl_table", "Download CSV", class = "btn-sm btn-outline-secondary")),
          tableOutput("table_summary")
        ),
        nav_panel("Joint Normality",
          div(style = "text-align:right; margin-bottom:4px;",
              downloadButton("dl_normality", "Download", class = "btn-sm btn-outline-secondary")),
          plotOutput("plot_normality", height = "400px")
        )
      )
    }
  })

  # Tab 1a: F1a imbalance over time
  output$plot_imbalance_f1a <- renderPlot({
    req(rv$results, input$normalize)
    r    <- rv$results[[1L]]
    ns   <- rv$ns
    norm <- input$normalize == "norm"
    panels <- lapply(seq_len(ns), function(k) {
      dm <- r$dms[[k]][1L, , ]   # F1a: stored normalized
      list(
        make_line_panel(dm, sprintf("Var[Dm(z=%d)]", k),    var,             r$max_z[k],
                        start = 4L, stored_normalized = TRUE, display_normalized = norm),
        make_line_panel(dm, sprintf("Q90[Dm(z=%d)]", k),    quantile_90_vec, r$max_z[k],
                        start = 4L, stored_normalized = TRUE, display_normalized = norm)
      )
    })
    rv$plot_f1a <- wrap_plots(unlist(panels, recursive = FALSE), ncol = 2L) +
      plot_annotation(title = "F1a Low — Imbalance Over Time")
    rv$plot_f1a
  })

  # Tab 1b: F1b imbalance over time
  output$plot_imbalance_f1b <- renderPlot({
    req(rv$results, input$normalize)
    r    <- rv$results[[1L]]
    ns   <- rv$ns
    norm <- input$normalize == "norm"
    panels <- lapply(seq_len(ns), function(k) {
      dm <- r$dms[[k]][2L, , ]   # F1b: stored raw
      list(
        make_line_panel(dm, sprintf("Var[Dm(z=%d)]", k),    var,             r$max_z[k],
                        start = 4L, stored_normalized = FALSE, display_normalized = norm),
        make_line_panel(dm, sprintf("Q90[Dm(z=%d)]", k),    quantile_90_vec, r$max_z[k],
                        start = 4L, stored_normalized = FALSE, display_normalized = norm)
      )
    })
    rv$plot_f1b <- wrap_plots(unlist(panels, recursive = FALSE), ncol = 2L) +
      plot_annotation(title = "F1b Low — Imbalance Over Time")
    rv$plot_f1b
  })

  # Tab 2: DLG (F1b only)
  output$plot_dlg <- renderPlot({
    req(rv$results)
    r  <- rv$results[[1L]]
    ns <- rv$ns
    panels <- lapply(seq_len(ns), function(k) {
      dlg <- r$dlgs[[k]][2L, , ]
      list(
        make_line_panel(dlg, sprintf("Q90[DLG(z=%d)]", k), quantile_90_vec, r$max_z[k]),
        make_line_panel(dlg, sprintf("max[DLG(z=%d)]", k), max,             r$max_z[k])
      )
    })
    rv$plot_dlg <- wrap_plots(unlist(panels, recursive = FALSE), ncol = 2L) +
      plot_annotation(title = "F1b Low — Distance to Last Gap (DLG)")
    rv$plot_dlg
  })

  # Tab 3: end-of-trial imbalance histograms
  output$plot_hist_imbalance <- renderPlot({
    req(rv$results)
    r  <- rv$results[[1L]]
    ns <- rv$ns
    panels <- c(
      lapply(seq_len(ns), function(k)
        make_imbalance_histogram(r$d500s[1L, , k], sprintf("F1a Low — z=%d", k))),
      lapply(seq_len(ns), function(k)
        make_imbalance_histogram(r$d500s[2L, , k], sprintf("F1b Low — z=%d", k)))
    )
    rv$plot_hist <- wrap_plots(panels, ncol = 2L) +
      plot_annotation(title = "End-of-Trial Imbalance Distribution")
    rv$plot_hist
  })

  # Tab 4: summary statistics table
  output$table_summary <- renderTable({
    req(rv$results)
    build_summary_table(rv$results[[1L]], rv$ns)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # Tab 5: joint normality check
  output$plot_normality <- renderPlot({
    req(rv$results)
    rv$plot_normality <- plot_joint_normality(rv$results[[1L]], rv$ns)
    rv$plot_normality
  })
  output$dl_f1a <- downloadHandler(
    filename = function() "imbalance_f1a.png",
    content  = function(file) {
      req(rv$plot_f1a)
      ggsave(file, plot = rv$plot_f1a, width = 12, height = 8, dpi = 150)
    }
  )

  output$dl_f1b <- downloadHandler(
    filename = function() "imbalance_f1b.png",
    content  = function(file) {
      req(rv$plot_f1b)
      ggsave(file, plot = rv$plot_f1b, width = 12, height = 8, dpi = 150)
    }
  )

  output$dl_dlg <- downloadHandler(
    filename = function() "dlg.png",
    content  = function(file) {
      req(rv$plot_dlg)
      ggsave(file, plot = rv$plot_dlg, width = 12, height = 8, dpi = 150)
    }
  )

  output$dl_hist <- downloadHandler(
    filename = function() "imbalance_histograms.png",
    content  = function(file) {
      req(rv$plot_hist)
      ggsave(file, plot = rv$plot_hist, width = 12, height = 8, dpi = 150)
    }
  )

  output$dl_table <- downloadHandler(
    filename = function() "summary_statistics.csv",
    content  = function(file) {
      req(rv$results)
      write.csv(build_summary_table(rv$results[[1L]], rv$ns), file, row.names = FALSE)
    }
  )

  output$dl_normality <- downloadHandler(
    filename = function() "joint_normality.png",
    content  = function(file) {
      req(rv$plot_normality)
      ggsave(file, plot = rv$plot_normality, width = 10, height = 5, dpi = 150)
    }
  )
}

shinyApp(ui, server)
