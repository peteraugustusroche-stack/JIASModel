# Shiny diagnostics for the HIV per-acquisition model.  Run from R/:  shiny::runApp("app.R")
# Live tabs re-run the engine on each change; the PSA tab reads out/Table1_psa_draws_*.csv (make with Rscript run_all.R psa).
library(shiny)
source("model_functions.R")

PAPER_C <- c(dalys_accrued = 5.46, dalys_no_care = 15.52, dalys_averted_care = 10.06, daly_ideal = 1.378, cost_art_drug = 1009,
             cost_art_svc = 2171, cost_hiv_care = 1281, cost_total = 4461, icer = 199, pct_ever_art = 0.758, med_cd4_at_dx = 312.6, med_lag_yr = 2.25)
PAPER_L <- c(dalys_accrued = 9.68, dalys_averted_care = 5.85, cost_art_drug = 519, cost_art_svc = 1207, cost_hiv_care = 1802,
             cost_total = 3528, icer = 182, pct_ever_art = 0.4386)
PAPER_DR <- data.frame(dr = c(0, .03, .05), acc = c(11.29, 5.46, 3.67), nc = c(32.82, 15.52, 10.19), av = c(21.53, 10.06, 6.52))
LODI <- c(7.9, 5.4, 15.8)
MODEL <- "#0f766e"; PAPER <- "#7c3a5e"; GREY <- "#9aa1ab"

ui <- fluidPage(
  tags$style("body{font-family:'Source Sans 3','Segoe UI',sans-serif;max-width:1100px;margin:auto} h4{margin-top:22px} .shiny-input-container{margin-bottom:6px}"),
  titlePanel("HIV per-acquisition model — diagnostics"),
  sidebarLayout(
    sidebarPanel(width = 3,
      sliderInput("n", "Individuals", 250, 4000, 1000, step = 250),
      numericInput("seed", "Seed", 42),
      sliderInput("dr", "Discount rate", 0, 0.07, 0.03, step = 0.01),
      h5("Cascade anchor"),
      sliderInput("p1", "Diagnosed", 0, 1, 0.93, 0.01), sliderInput("p2", "On ART | diagnosed", 0, 1, 0.91, 0.01),
      sliderInput("p3", "Suppressed | on ART", 0.5, 1, 0.95, 0.01), sliderInput("ret", "12-month retention", 0.5, 1, 0.86, 0.01),
      h5("Uncertain structure"),
      numericInput("dx_rate", "Diagnostic hazard after threshold (/yr; 0 = immediate)", 4),
      numericInput("floor", "On-ART SMR floor", 1.20, step = 0.05),
      checkboxInput("l2", "Second-line rescue available", FALSE),
      selectInput("band", "Mortality age band", c("attained", "acquisition")),
      actionButton("go", "Run model", class = "btn-primary")),
    mainPanel(width = 9, tabsetPanel(
      tabPanel("Replication", h4("Model vs paper (contemporary, then lower coverage)"), tableOutput("repC"), tableOutput("repL"),
               p("Lower-coverage arm is run with 70 / 70 / 80 and 0.70 retention regardless of the sliders; the sliders drive the contemporary arm.")),
      tabPanel("Cascade", plotOutput("casc", height = 320), plotOutput("cd4dx", height = 300), plotOutput("ptyear", height = 300)),
      tabPanel("Natural history", tableOutput("nh"), plotOutput("nhplot", height = 340)),
      tabPanel("Discount rate", p("Sweep of 0–7% at the current settings (runs eight times)."), actionButton("sweep", "Run sweep"), plotOutput("drplot", height = 360), tableOutput("drtab")),
      tabPanel("PSA", uiOutput("psaui")),
      tabPanel("Life table", tableOutput("lt"), plotOutput("ltplot", height = 300)))))
)

server <- function(input, output, session) {
  pars <- function(...) default_params(n = input$n, seed = input$seed, dr = input$dr, p_diagnosed = input$p1, p_art_if_dx = input$p2,
                                       p3 = input$p3, retention_12m = input$ret, dx_rate = if (input$dx_rate <= 0) Inf else input$dx_rate,
                                       art_smr_floor = input$floor, l2_rescue = input$l2, age_band_mode = input$band, ...)
  runs <- eventReactive(input$go, {
    withProgress(message = "Running three arms...", {
      p <- pars(); s <- run_acquisition_cohort(p); incProgress(.6)
      sl <- run_acquisition_cohort(pars(p_diagnosed = .70, p_art_if_dx = .70, p3 = .80, retention_12m = .70), arms = c("no_art", "real_world"))
      list(p = p, s = s, sl = sl) })
  }, ignoreNULL = FALSE)

  summ <- function(s) { f <- cost_daly_frame(s); ic <- compute_icer(s$no_art, s$real_world); cd <- cascade_diagnostics(s$real_world)
    c(unlist(f), icer = ic$icer, inc_cost = ic$inc_cost, unlist(cd)) }
  reptab <- function(m, paper) { k <- names(paper); dev <- 100 * (m[k] - paper) / abs(paper)
    data.frame(Output = k, Model = signif(m[k], 4), Paper = paper, `Deviation %` = round(dev, 1), check.names = FALSE) }
  output$repC <- renderTable(reptab(summ(runs()$s), PAPER_C), striped = TRUE)
  output$repL <- renderTable(reptab(summ(runs()$sl), PAPER_L), striped = TRUE)

  output$casc <- renderPlot({ r <- runs(); rw <- r$s$real_world; cd <- cascade_diagnostics(rw)
    m <- rbind(anchor = c(input$p1, input$p2, input$p3), realised = c(cd$pct_ever_dx, cd$art_given_dx, NA),
               person_time = c(cd$pt_dx, cd$pt_art_given_dx, cd$pt_supp_given_art))
    barplot(m, beside = TRUE, names.arg = c("Diagnosed", "On ART | dx", "Suppressed | ART"), col = c(PAPER, MODEL, "#cfe8e4"), ylim = c(0, 1),
            main = "Cascade: assigned anchor, lifetime realised, person-time", border = NA)
    legend("topright", c("Anchor", "Realised", "Person-time"), fill = c(PAPER, MODEL, "#cfe8e4"), bty = "n") })
  output$cd4dx <- renderPlot({ r <- runs(); rw <- r$s$real_world; cd <- cascade_diagnostics(rw)
    h <- hist(rw$cd4_dx[is.finite(rw$cd4_dx)], breaks = c(0, 50, 100, 200, 350, 500, 700), plot = FALSE)
    m <- rbind(anchor = r$p$cd4_schedule_props, model = h$counts / sum(h$counts))
    barplot(m, beside = TRUE, names.arg = c("<50", "50-99", "100-199", "200-349", "350-499", ">=500"), col = c(PAPER, MODEL), border = NA,
            main = sprintf("CD4 at diagnosis: median %.0f (paper 312.6), lag %.2f yr (paper 2.25)", cd$med_cd4_at_dx, cd$med_lag_yr))
    legend("topleft", c("IeDEA schedule", "Model"), fill = c(PAPER, MODEL), bty = "n") })
  output$ptyear <- renderPlot({ rw <- runs()$s$real_world; pt <- rw$pt_year; k <- pt[, 1] > 5
    plot(rw$years[k], pt[k, 2] / pt[k, 1], type = "l", col = PAPER, lwd = 2, ylim = c(0, 1), xlab = "Calendar year", ylab = "Fraction of person-time", main = "Person-time cascade by calendar year")
    lines(rw$years[k], pt[k, 3] / pt[k, 2], col = MODEL, lwd = 2); lines(rw$years[k], pt[k, 4] / pt[k, 3], col = GREY, lwd = 2, lty = 2)
    legend("bottomright", c("Diagnosed", "On ART | dx", "Suppressed | ART"), col = c(PAPER, MODEL, GREY), lwd = 2, lty = c(1, 1, 2), bty = "n") })

  nh <- reactive({ r <- runs(); p <- r$p; crn <- r$s$crn; ind <- r$s$ind; n <- p$n
    cd4 <- ind$cd4_0; alive <- rep(TRUE, n); age <- ind$age0; t200 <- td <- rep(NA_real_, n); surv <- med <- numeric(0)
    for (k in 1:crn$T) { t <- (k - 1) * DT; if (!any(alive)) break
      surv <- c(surv, mean(alive)); med <- c(med, median(cd4[alive]))
      h <- bg_hazard(age, p) + glaubius_rate(cd4, ind$age_band); dies <- alive & (crn$u_death[k, ] < 1 - exp(-h * DT)); td[dies] <- t; alive <- alive & !dies
      lvl <- ifelse(cd4 < 200, ind$vsp + (200 - cd4) / 200 * (7 - ind$vsp), ind$vsp); Dd <- pmax(0, lvl - LOG200)
      cd4 <- ifelse(alive, pmax(cd4 + DT * (p$s_C * pmin(1, cd4 / p$kappa) - (p$d_C + p$phi * Dd) * cd4), 1), cd4)
      new <- alive & cd4 < 200 & is.na(t200); t200[new] <- t + DT; age <- age + DT }
    list(surv = surv, med = med, t200 = t200, td = td, vsp = ind$vsp) })
  output$nh <- renderTable({ x <- nh(); data.frame(Quantity = c("Years to CD4<200, all", "VL >= 1e5", "VL < 1e4", "Mean untreated survival"),
    Model = round(c(mean(x$t200, na.rm = TRUE), mean(x$t200[x$vsp >= 5], na.rm = TRUE), mean(x$t200[x$vsp < 4], na.rm = TRUE), mean(x$td, na.rm = TRUE)), 2), Lodi = c(LODI, NA)) })
  output$nhplot <- renderPlot({ x <- nh(); t <- (seq_along(x$surv) - 1) * DT; par(mar = c(4, 4, 2, 4))
    plot(t, x$surv, type = "l", col = MODEL, lwd = 2, ylim = c(0, 1), xlab = "Years since acquisition", ylab = "Surviving, untreated", main = "Untreated survival and median CD4")
    par(new = TRUE); plot(t, x$med, type = "l", col = PAPER, lwd = 2, axes = FALSE, xlab = "", ylab = ""); axis(4); mtext("Median CD4", 4, 2.5) })

  sw <- eventReactive(input$sweep, { withProgress(message = "Discount sweep", { do.call(rbind, lapply(seq(0, .07, .01), function(d) {
    s <- run_acquisition_cohort(pars(dr = d)); f <- cost_daly_frame(s); incProgress(1 / 8); cbind(dr = d, f, icer = compute_icer(s$no_art, s$real_world)$icer) })) }) })
  output$drplot <- renderPlot({ r <- sw(); par(mfrow = c(1, 2))
    matplot(r$dr * 100, r[, c("dalys_accrued", "dalys_no_care", "dalys_averted_care")], type = "l", lty = c(1, 2, 3), col = MODEL, lwd = 2, xlab = "Discount rate (%)", ylab = "DALYs per acquisition", main = "DALYs")
    points(rep(PAPER_DR$dr * 100, 3), unlist(PAPER_DR[, c("acc", "nc", "av")]), col = PAPER, pch = 19)
    legend("topright", c("Accrued", "No care", "Averted", "Paper"), lty = c(1, 2, 3, NA), pch = c(NA, NA, NA, 19), col = c(MODEL, MODEL, MODEL, PAPER), bty = "n")
    plot(r$dr * 100, r$icer, type = "l", col = MODEL, lwd = 2, xlab = "Discount rate (%)", ylab = "US$ per DALY averted", main = "Treatment ICER") })
  output$drtab <- renderTable({ r <- sw(); data.frame(`dr %` = r$dr * 100, Accrued = round(r$dalys_accrued, 2), `No care` = round(r$dalys_no_care, 2), Averted = round(r$dalys_averted_care, 2),
    Ideal = round(r$daly_ideal, 3), Cost = round(r$cost_total), ICER = round(r$icer), check.names = FALSE) })

  output$psaui <- renderUI({ fC <- "out/Table1_psa_draws_contemporary.csv"; fL <- "out/Table1_psa_draws_lower.csv"
    if (!file.exists(fC)) return(p("No PSA draws found. Run  Rscript run_all.R psa  from this directory first."))
    tagList(tableOutput("psatab"), plotOutput("ceplane", height = 360)) })
  psa <- reactive(list(C = read.csv("out/Table1_psa_draws_contemporary.csv"), L = read.csv("out/Table1_psa_draws_lower.csv")))
  output$psatab <- renderTable({ x <- psa(); keys <- c("dalys_accrued", "dalys_no_care", "dalys_averted_care", "daly_ideal", "cost_total", "icer")
    q <- function(v) sprintf("%.2f (%.2f – %.2f)", mean(v), quantile(v, .025), quantile(v, .975))
    data.frame(Output = c(keys, "P(net-saving $125)", "P(net-saving $55)"),
      Contemporary = c(sapply(keys, function(k) q(x$C[[k]])), sprintf("%.0f%%", 100 * mean(x$C$net_saving_125)), sprintf("%.0f%%", 100 * mean(x$C$net_saving_55))),
      Paper.C = c("5.27 (4.33–6.10)", "15.25 (11.92–17.72)", "9.98 (7.53–11.92)", "1.378 (1.019–1.925)", "4491 (3218–6193)", "197 (48–382)", "34%", "83%"),
      Lower = c(sapply(keys, function(k) q(x$L[[k]])), sprintf("%.0f%%", 100 * mean(x$L$net_saving_125)), sprintf("%.0f%%", 100 * mean(x$L$net_saving_55))),
      Paper.L = c("9.50 (7.59–10.89)", "—", "5.74 (4.36–7.07)", "—", "3549 (2649–4671)", "180 (44–342)", "24%", "68%")) })
  output$ceplane <- renderPlot({ x <- psa()$C
    plot(x$dalys_averted_care, x$inc_cost, pch = 19, col = adjustcolor(MODEL, .6), xlim = c(0, max(x$dalys_averted_care) * 1.05), ylim = c(0, max(x$inc_cost) * 1.05),
         xlab = "DALYs averted by care", ylab = "Incremental cost (US$)", main = "Cost-effectiveness plane, contemporary")
    abline(0, 500, lty = 2, col = GREY); text(1, 900, "US$500 / DALY", col = GREY, adj = 0) })

  output$lt <- renderTable({ p <- default_params(); e0 <- life_expectancy_at(0, p); e28 <- life_expectancy_at(28, p)
    data.frame(Quantity = c("Life expectancy at birth", "Remaining LE at 28", "Expected age at death, acquiring at 28"), Value = round(c(e0, e28, 28 + e28), 1), Paper = c(NA, NA, 72)) })
  output$ltplot <- renderPlot({ plot(0:100, 1 - exp(-BG_HAZ), log = "y", type = "l", col = MODEL, lwd = 2, xlab = "Age", ylab = "qx (non-HIV)", main = "Background mortality, ESA") })
}
shinyApp(ui, server)
