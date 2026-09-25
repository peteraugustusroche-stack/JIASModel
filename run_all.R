# Usage: Rscript run_all.R [checks|psa|fig2|fig3|all]   (run from the directory holding bg_life_table.csv)
source("model_functions.R")
what <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "checks"
dir.create("out", showWarnings = FALSE)
LOWER <- list(p_diagnosed = 0.70, p_art_if_dx = 0.70, p3 = 0.80, retention_12m = 0.70)

# ---------------------------------------------------------------- checks (spec s11)
if (what %in% c("checks", "all")) {
  p <- default_params()
  s <- run_acquisition_cohort(p); f <- cost_daly_frame(s); ic <- compute_icer(s$no_art, s$real_world)
  tgt <- c(dalys_accrued = 5.46, dalys_no_care = 15.52, dalys_averted_care = 10.06, daly_ideal = 1.378,
           cost_art_drug = 1009, cost_art_svc = 2171, cost_hiv_care = 1281, cost_art = 3180, cost_total = 4461)
  cat("Table 1, contemporary\n")
  for (k in names(tgt)) cat(sprintf("  %-20s model %9.2f   paper %s\n", k, f[[k]], tgt[[k]]))
  cat(sprintf("  %-20s model %9.0f   paper 1998\n  %-20s model %9.0f   paper 199\n", "inc_cost", ic$inc_cost, "ICER", ic$icer))
  for (cc in c(125, 55)) cat(sprintf("  prevention ICER $%d  model %9.0f\n", cc, compare_prevention_treatment(s, cc, 56)$prev_icer))
  d <- cascade_diagnostics(s$real_world); d$assigned_p1 <- mean(s$crn$u_p1 < p$p_diagnosed)
  cat("Realised cascade:\n"); for (k in names(d)) cat(sprintf("  %-18s %.3f\n", k, d[[k]]))
  sl <- run_acquisition_cohort(do.call(default_params, LOWER), arms = c("no_art", "real_world"))
  fl <- cost_daly_frame(sl); icl <- compute_icer(sl$no_art, sl$real_world)
  cat("Table 1, lower coverage\n")
  for (k in c("dalys_accrued", "dalys_averted_care", "cost_art_drug", "cost_art_svc", "cost_hiv_care", "cost_total"))
    cat(sprintf("  %-20s model %9.2f\n", k, fl[[k]]))
  cat(sprintf("  ICER %.0f (182); ever_art %.3f (0.4386)\n", icl$icer, cascade_diagnostics(sl$real_world)$pct_ever_art))
  cat("Table 2 (accrued / no care / averted)\n")
  for (kw in list(list(dr = 0), list(dr = 0.05), list(mean_age = 22), list(mean_age = 34))) {
    r <- cost_daly_frame(run_acquisition_cohort(do.call(default_params, kw), arms = c("no_art", "real_world")))
    cat(sprintf("  %s: %.2f / %.2f / %.2f\n", paste(names(kw), kw), r$dalys_accrued, r$dalys_no_care, r$dalys_averted_care))
  }
}

# ---------------------------------------------------------------- PSA (spec s10.1; p1-p3 fixed)
PSA_SPEC <- list(
  mu_VL = c(4.4, 0.15, 1), sigma_VL = c(0.80, 0.10, 1), mean_age = c(28, 3, 1), sd_age = c(8, 2, 1),
  phi = c(0.090, 0.015, 1), r_rec = c(0.40, 0.08, 1), art_smr_floor = c(1.20, 0.10, 1), bg_mort_mult = c(1.0, 0.10, 1),
  art_init_rate = c(8, 1, 1), vf_rate_L1 = c(0.010, 0.003, 2), vf_rate_L2 = c(0.05, 0.015, 2), rescue_rate = c(0.90, 0.05, 2),
  phi_ltfu_mult = c(1.5, 0.15, 1), retention_12m = c(0.86, 0.03, 2),
  art_cost_L1 = c(45, 8, 1), art_cost_L2 = c(300, 60, 1), cost_provision = c(150, 50, 1), art_startup_cost = c(50, 20, 1),
  hc_cost_200_350 = c(250, 60, 1), hc_cost_100_200 = c(450, 100, 1), hc_cost_LT100 = c(650, 150, 1))   # 1 = gamma, 2 = beta
DW_SPEC <- list(dw_noart = rbind(c(0.08, 0.015), c(0.57, 0.06), c(0.008, 0.001)),
                dw_art = rbind(c(0.04, 0.008), c(0.19, 0.03), c(0.010, 0.002)))
draw1 <- function(m, s, dist) {
  if (dist == 1) { k <- (m / s)^2; rgamma(1, k, scale = m / k) }
  else { c <- m * (1 - m) / s^2 - 1; rbeta(1, m * c, (1 - m) * c) }
}

psa <- function(n_draws = 200, n_ind = 1000, seed = 42) {
  set.seed(seed); rows <- list(contemporary = list(), lower = list())
  for (d in 1:n_draws) {
    kw <- lapply(PSA_SPEC, function(v) draw1(v[1], v[2], v[3]))
    for (nm in names(DW_SPEC)) kw[[nm]] <- apply(DW_SPEC[[nm]], 1, function(v) draw1(v[1], v[2], 2))
    nnp <- draw1(56, 32, 1)
    for (a in names(rows)) {
      kwa <- if (a == "lower") modifyList(kw, LOWER) else kw
      p <- do.call(default_params, c(list(n = n_ind, seed = seed + d), kwa))
      s <- run_acquisition_cohort(p); f <- cost_daly_frame(s); ic <- compute_icer(s$no_art, s$real_world); cd <- cascade_diagnostics(s$real_world)
      r <- c(draw = d, nnp = nnp, as.list(f), inc_cost = ic$inc_cost, icer = ic$icer,
             pct_ever_dx = cd$pct_ever_dx, pct_ever_art = cd$pct_ever_art, assigned_p1 = mean(s$crn$u_p1 < p$p_diagnosed),
             treat_ce_500 = ic$icer < 500)
      for (cc in c(125, 55)) { pv <- compare_prevention_treatment(s, cc, nnp); r[[paste0("prev_icer_", cc)]] <- pv$prev_icer; r[[paste0("net_saving_", cc)]] <- pv$net_cost_saving }
      rows[[a]][[d]] <- as.data.frame(r)
    }
    if (d %% 20 == 0) cat(sprintf("  PSA draw %d/%d\n", d, n_draws))
  }
  lapply(rows, function(x) do.call(rbind, x))
}

if (what %in% c("psa", "all")) {
  res <- psa()
  for (a in names(res)) {
    write.csv(res[[a]], sprintf("out/Table1_psa_draws_%s.csv", a), row.names = FALSE)
    cat(sprintf("\n== PSA %s: mean (95%% UI) ==\n", a))
    for (k in c("dalys_accrued", "dalys_no_care", "dalys_averted_care", "daly_ideal", "cost_art_drug", "cost_art_svc", "cost_hiv_care",
                "cost_total", "inc_cost", "icer", "nnp", "prev_icer_125", "prev_icer_55", "pct_ever_art", "assigned_p1")) {
      x <- res[[a]][[k]]; cat(sprintf("  %-18s %9.3f  (%.3f - %.3f)\n", k, mean(x), quantile(x, .025), quantile(x, .975)))
    }
    for (k in c("treat_ce_500", "net_saving_125", "net_saving_55")) cat(sprintf("  P(%s) = %.0f%%\n", k, 100 * mean(res[[a]][[k]])))
  }
}

# ---------------------------------------------------------------- Figure 2 / 3
sweep <- function(years, n, tvp, seed = 42) {
  do.call(rbind, lapply(years, function(y) {
    p <- default_params(n = n, seed = seed, year = y, cascade_tvp = tvp)
    s <- run_acquisition_cohort(p, arms = c("no_art", "real_world")); f <- cost_daly_frame(s)
    ic <- compute_icer(s$no_art, s$real_world); cd <- cascade_diagnostics(s$real_world)
    cbind(year = y, f, inc_cost = ic$inc_cost, pct_ever_dx = 100 * cd$pct_ever_dx, pct_ever_art = 100 * cd$pct_ever_art,
          med_cd4_at_dx = cd$med_cd4_at_dx, med_lag_yr = cd$med_lag_yr,
          pt_dx = cd$pt_dx, pt_art_given_dx = cd$pt_art_given_dx, pt_supp_given_art = cd$pt_supp_given_art)
  }))
}
if (what %in% c("fig2", "all")) { r2 <- sweep(seq(1980, 2030, 2), 800, FALSE); write.csv(r2, "out/Figure2_data.csv", row.names = FALSE); print(r2[, c("year", "dalys_accrued", "dalys_no_care", "yld", "yll", "pct_ever_art")]) }
if (what %in% c("fig3", "all")) { r3 <- sweep(c(2004, 2024), 2000, TRUE); write.csv(r3, "out/Figure3_data.csv", row.names = FALSE); print(t(r3)) }
