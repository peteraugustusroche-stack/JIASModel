# HIV per-acquisition DALY microsimulation -- R port of engine.py (calibrated defaults, Sept 2026)
# Base R only. Three arms on common random numbers: no_art, real_world, ideal.

DT <- 0.25
LOG200 <- log10(200)
AGE_BANDS <- c(25, 35, 45)
CD4_BANDS <- c(50, 100, 200, 250, 350)
GLAUBIUS <- matrix(c(0.712, 0.954, 1.197, 1.439,
                     0.227, 0.305, 0.382, 0.459,
                     0.048, 0.064, 0.080, 0.097,
                     0.007, 0.010, 0.012, 0.015,
                     0.002, 0.002, 0.003, 0.003,
                     0.000, 0.000, 0.000, 0.000), nrow = 6, byrow = TRUE)
CD4_EDGES <- c(0, 50, 100, 200, 350, 500, 700)

CASCADE_TS <- matrix(c(
  1990, 0.01, 0.30, 0.60, .24, .23, .26, .16, .07, .04, 700, 1200, 120,
  2000, 0.06, 0.40, 0.65, .19, .20, .25, .19, .10, .07, 700, 1200, 120,
  2005, 0.20, 0.58, 0.73, .13, .15, .20, .27, .15, .10, 400,  900, 140,
  2010, 0.45, 0.70, 0.77, .11, .13, .18, .29, .17, .12, 200,  600, 150,
  2013, 0.62, 0.78, 0.79, .11, .12, .17, .29, .18, .13, 140,  500, 150,
  2016, 0.78, 0.86, 0.82, .08, .10, .15, .29, .20, .18,  90,  400, 150,
  2020, 0.90, 0.90, 0.85, .06, .09, .13, .28, .21, .23,  64,  350, 150,
  2024, 0.93, 0.92, 0.86, .05, .07, .12, .27, .22, .27,  45,  300, 150), ncol = 13, byrow = TRUE)
colnames(CASCADE_TS) <- c("year", "p1", "link", "ret", "s1", "s2", "s3", "s4", "s5", "s6", "L1", "L2", "prov")

default_params <- function(...) {
  p <- list(
    n = 2000, seed = 42, mean_age = 28, sd_age = 8, mu_VL = 4.4, sigma_VL = 0.80,
    cd4_init_mode = "loglogistic", cd4_ll_medians = c(579, 554, 528, 503), cd4_ll_shape = 2.92, frailty_sd = 0,
    s_C = 38.6, d_C = 0.0428, phi = 0.090, kappa = 500, log10_Vmax = 7, r_rec = 0.40, cd4_ceiling = 800,
    lam_art = 21, log10_Vfloor = log10(50), substeps = 1, death_timing = 0.5,
    age_band_mode = "attained", art_smr_floor = 1.20, smr_at50 = 10, smr_top = 500, bg_mort_mult = 1,
    year = 2024, cascade_tvp = FALSE, policy_mode = TRUE,
    p_diagnosed = 0.93, p_art_if_dx = 0.91, p3 = 0.95, retention_12m = 0.86, p_supp_on_art = 0.95,
    cd4_schedule_props = c(0.05, 0.07, 0.12, 0.27, 0.22, 0.27),
    art_init_rate = 8, dx_rate_bg = 0, dx_rate = 4, phi_ltfu_mult = 1.5, ltfu_phi_years = 2,
    vf_rate_L1 = 0.010, vf_rate_L2 = 0.05, rescue_rate = 0.90, ltfu_rate_override = NA, l2_rescue = FALSE,
    dw_noart = c(0.08, 0.57, 0.008), dw_art = c(0.04, 0.19, 0.010), dr = 0.03, ref_le = 72, yll_hiv_deaths_only = TRUE,
    art_cost_L1 = 45, art_cost_L2 = 300, cost_provision = 150, art_startup_cost = 50,
    hc_cost_200_350 = 250, hc_cost_100_200 = 450, hc_cost_LT100 = 650, max_age = 100)
  modifyList(p, list(...))
}

# ---------------------------------------------------------------- helpers
derive_cascade_params <- function(retention_12m, p3, p_supp_on_art = 0.95) {
  h <- -log(retention_12m); p3_adj <- pmin(0.999, p3 / p_supp_on_art)
  list(h = h, r = h * p3_adj / (1 - p3_adj))
}
q <- function(rate) 1 - exp(-rate * DT)
dw_curve <- function(cd4, pars) pars[1] + (pars[2] - pars[1]) * exp(-pars[3] * cd4)
hc_cost <- function(cd4, p) ifelse(cd4 < 100, p$hc_cost_LT100, ifelse(cd4 < 200, p$hc_cost_100_200, ifelse(cd4 <= 350, p$hc_cost_200_350, 0)))

load_life_table <- function(path = "bg_life_table.csv") {
  lt <- read.csv(path)
  -log(pmax(1 - lt$qx_non_hiv, 1e-9))            # annual hazard by single year of age 0..100
}
BG_HAZ <- if (file.exists("bg_life_table.csv")) load_life_table("bg_life_table.csv") else NULL

bg_hazard <- function(age, p) {
  stopifnot(!is.null(BG_HAZ))
  a <- pmin(pmax(floor(age), 0), length(BG_HAZ) - 1) + 1
  p$bg_mort_mult * BG_HAZ[a]
}
life_expectancy_at <- function(age0, p) {
  ages <- seq(age0, 110, by = DT); S <- exp(-cumsum(bg_hazard(ages, p)) * DT); sum(S) * DT
}

who_art_policy <- function(year) {
  thr <- if (year < 2001) -Inf else if (year < 2010) 200 else if (year < 2013) 350 else if (year < 2016) 500 else Inf
  vf <- if (year < 2016) 0.04 else if (year < 2020) 0.02 else 0.01
  c(thr = thr, vf = vf)
}
cascade_state <- function(years) {
  out <- lapply(colnames(CASCADE_TS)[-1], function(c) approx(CASCADE_TS[, "year"], CASCADE_TS[, c], years, rule = 2)$y)
  names(out) <- colnames(CASCADE_TS)[-1]; out
}
smr_schedule <- function(cd4, p) {
  f <- pmin(pmax((p$smr_top - cd4) / (p$smr_top - 50), 0), 1)
  pmax(p$art_smr_floor, exp(f * log(p$smr_at50)))
}
glaubius_rate <- function(cd4, age_band) GLAUBIUS[cbind(findInterval(cd4, CD4_BANDS) + 1, age_band + 1)]
yll_discounted <- function(t_death, age_death, p) {
  lost <- pmax(0, p$ref_le - age_death)
  if (p$dr == 0) return(lost)
  rho <- log(1 + p$dr); exp(-rho * t_death) * (1 - exp(-rho * lost)) / rho
}
# fixed percentile u -> CD4 target through a six-band piecewise-uniform schedule (vector or n x 6 matrix)
cd4_target_from_percentile <- function(u, props) {
  n <- length(u)
  if (is.null(dim(props))) props <- matrix(props, n, 6, byrow = TRUE)
  props <- props / rowSums(props); cum <- t(apply(props, 1, cumsum))
  band <- pmin(rowSums(u > cum), 5)                  # 0-based
  lo_c <- ifelse(band == 0, 0, cum[cbind(1:n, pmax(band, 1))])
  w <- props[cbind(1:n, band + 1)]
  frac <- pmin(pmax((u - lo_c) / pmax(w, 1e-12), 0), 1)
  CD4_EDGES[band + 1] + frac * (CD4_EDGES[band + 2] - CD4_EDGES[band + 1])
}

# ---------------------------------------------------------------- CRN draws
make_crn <- function(p) {
  set.seed(p$seed); n <- p$n
  Tn <- ceiling((p$max_age - 15) / DT) + 1
  U <- function() matrix(runif(Tn * n), Tn, n)
  list(T = Tn,
       age0 = pmin(pmax(rnorm(n, p$mean_age, p$sd_age), 15), 70),
       vsp = pmin(pmax(rnorm(n, p$mu_VL, p$sigma_VL), 3), 7),
       u_cd4 = runif(n), u_frail = rnorm(n), u_p1 = runif(n), u_p2 = runif(n), u_stratum = runif(n), u_within = runif(n),
       u_death = U(), u_init = U(), u_ltfu = U(), u_reeng = U(), u_vf = U(), u_switch = U(), u_cause = U(), u_dx = U())
}

init_individuals <- function(p, crn) {
  age_band <- findInterval(crn$age0, AGE_BANDS)
  if (p$cd4_init_mode == "loglogistic") {
    u <- crn$u_cd4; cd4 <- pmin(pmax(p$cd4_ll_medians[age_band + 1] * (u / (1 - u))^(1 / p$cd4_ll_shape), 200), 1500)
  } else cd4 <- rep(p$s_C / p$d_C, p$n)
  frail <- if (p$frailty_sd > 0) exp(crn$u_frail * p$frailty_sd - 0.5 * p$frailty_sd^2) else rep(1, p$n)
  props <- p$cd4_schedule_props / sum(p$cd4_schedule_props)
  stratum <- pmin(rowSums(outer(crn$u_stratum, cumsum(props), ">")), 5)
  list(age0 = crn$age0, age_band = age_band, vsp = crn$vsp, cd4_0 = cd4, frail = frail,
       dx_thr = cd4_target_from_percentile(crn$u_stratum, props), stratum = stratum)
}

# ---------------------------------------------------------------- core loop
simulate_arm <- function(p, crn, ind, arm) {
  n <- p$n; Tn <- crn$T
  years <- p$year + (0:(Tn - 1)) * DT
  if (p$cascade_tvp) {
    cs <- cascade_state(years)
    p1_k <- cs$p1; link_k <- cs$link; ret_k <- cs$ret
    sched_k <- cbind(cs$s1, cs$s2, cs$s3, cs$s4, cs$s5, cs$s6)
    L1_k <- cs$L1; L2_k <- cs$L2; prov_k <- cs$prov
  } else {
    p1_k <- rep(p$p_diagnosed, Tn); link_k <- rep(p$p_art_if_dx, Tn); ret_k <- rep(p$retention_12m, Tn)
    sched_k <- matrix(p$cd4_schedule_props, Tn, 6, byrow = TRUE)
    L1_k <- rep(p$art_cost_L1, Tn); L2_k <- rep(p$art_cost_L2, Tn); prov_k <- rep(p$cost_provision, Tn)
  }
  if (p$policy_mode) {
    pol <- t(sapply(years, who_art_policy)); thr_k <- pol[, "thr"]; vf_k <- pol[, "vf"]
    if (!p$cascade_tvp && p$year >= 2020) vf_k[] <- p$vf_rate_L1
  } else { thr_k <- rep(Inf, Tn); vf_k <- rep(p$vf_rate_L1, Tn) }
  if (is.na(p$ltfu_rate_override)) { cp <- derive_cascade_params(ret_k, p$p3, p$p_supp_on_art); h_k <- cp$h; r_k <- cp$r
  } else { p3a <- min(0.999, p$p3 / p$p_supp_on_art); h_k <- rep(p$ltfu_rate_override, Tn); r_k <- h_k * p3a / (1 - p3a) }

  ideal <- arm == "ideal"; noart <- arm == "no_art"
  will_art <- rep(ideal, n)
  alive <- rep(TRUE, n); cd4 <- ind$cd4_0; lvl <- ind$vsp; nadir <- cd4
  diagnosed <- rep(ideal, n); t_dx <- ifelse(diagnosed, 0, NA); cd4_dx <- ifelse(diagnosed, cd4, NA)
  on_art <- ever_art <- line2 <- failing <- ltfu <- rep(FALSE, n); t_ltfu <- rep(-Inf, n)
  age <- ind$age0; phi_i <- p$phi * ind$frail
  yld <- yld_prediag <- yll <- c_drug <- c_svc <- c_start <- c_hc <- py_total <- py_dx <- py_art <- py_supp <- numeric(n)
  ny <- ceiling(Tn * DT) + 1; pt_year <- matrix(0, ny, 4)
  t_death <- age_death <- rep(NA_real_, n); cause_hiv <- rep(FALSE, n)
  rho <- log(1 + p$dr)
  inc_of <- function(nd) ifelse(nd >= 200, 350, ifelse(nd >= 100, 250, 200))

  for (k in 1:Tn) {
    t <- (k - 1) * DT
    if (!any(alive)) break
    idx <- alive; df <- exp(-rho * t)
    suppressed <- on_art & !failing & (lvl <= LOG200)
    ab <- if (p$age_band_mode == "acquisition") ind$age_band else findInterval(age, AGE_BANDS)
    h_bg <- bg_hazard(age, p)
    h_hiv <- ifelse(suppressed, h_bg * (smr_schedule(cd4, p) - 1), glaubius_rate(cd4, ab))
    dies <- idx & (crn$u_death[k, ] < 1 - exp(-(h_bg + h_hiv) * DT))

    w <- ifelse(dies, p$death_timing, 1) * DT * idx
    dw <- ifelse(suppressed, dw_curve(cd4, p$dw_art), dw_curve(cd4, p$dw_noart))
    yld <- yld + w * dw * df; yld_prediag <- yld_prediag + w * dw * df * !diagnosed
    c_hc <- c_hc + w * hc_cost(cd4, p) * df
    c_drug <- c_drug + w * ifelse(on_art, ifelse(line2, L2_k[k], L1_k[k]), 0) * df
    c_svc <- c_svc + w * ifelse(on_art, prov_k[k], 0) * df
    py_total <- py_total + w; py_dx <- py_dx + w * diagnosed; py_art <- py_art + w * on_art; py_supp <- py_supp + w * suppressed
    yi <- floor(t) + 1; pt_year[yi, ] <- pt_year[yi, ] + c(sum(w), sum(w * diagnosed), sum(w * on_art), sum(w * suppressed))

    if (any(dies)) {
      t_death[dies] <- t + DT * p$death_timing; age_death[dies] <- age[dies] + DT * p$death_timing
      hiv_death <- dies & (crn$u_cause[k, ] < h_hiv / pmax(h_bg + h_hiv, 1e-12))
      cause_hiv[hiv_death] <- TRUE
      counted <- if (p$yll_hiv_deaths_only) hiv_death else dies
      yll[counted] <- yll_discounted(t_death[counted], age_death[counted], p)
      alive <- alive & !dies; idx <- alive
    }

    if (!noart && !ideal) {
      target <- cd4_target_from_percentile(crn$u_stratum, sched_k[k, ])
      h_dx <- ifelse(cd4 <= target, p$dx_rate, 0) + p$dx_rate_bg
      newdx <- idx & !diagnosed & (crn$u_p1 < p1_k[k]) & (crn$u_dx[k, ] < q(h_dx))
      diagnosed <- diagnosed | newdx; t_dx[newdx] <- t; cd4_dx[newdx] <- cd4[newdx]
      will_art <- will_art | (newdx & (crn$u_p2 < link_k[k]))
    }
    eligible <- idx & diagnosed & will_art & !ever_art & (cd4 < thr_k[k])
    init <- eligible & (crn$u_init[k, ] < q(p$art_init_rate))
    reeng <- idx & ltfu & (crn$u_reeng[k, ] < q(r_k[k]))
    start <- init | reeng
    if (any(start)) {
      on_art[start] <- TRUE; ever_art[start] <- TRUE; ltfu[start] <- FALSE; failing[start] <- FALSE
      c_start[start] <- c_start[start] + p$art_startup_cost * df
    }
    if (!ideal) {
      drop <- idx & on_art & !start & (crn$u_ltfu[k, ] < q(h_k[k]))
      on_art[drop] <- FALSE; ltfu[drop] <- TRUE; t_ltfu[drop] <- t; failing[drop] <- FALSE; lvl[drop] <- ind$vsp[drop]
      supp_now <- idx & on_art & !failing & (lvl <= LOG200)
      vf <- supp_now & (crn$u_vf[k, ] < ifelse(line2, q(p$vf_rate_L2), q(vf_k[k])))
      failing[vf] <- TRUE; lvl[vf] <- ind$vsp[vf]
      resc <- idx & on_art & failing & !vf & (crn$u_switch[k, ] < q(p$rescue_rate)) & (p$l2_rescue | !line2)
      line2[resc & !line2] <- TRUE; failing[resc] <- FALSE
    }

    decaying <- idx & on_art & !failing
    lvl <- ifelse(decaying, p$log10_Vfloor + (lvl - p$log10_Vfloor) * exp(-p$lam_art * DT), lvl)
    offtx <- idx & !decaying
    lvl <- ifelse(offtx & cd4 < 200, ind$vsp + (200 - cd4) / 200 * (p$log10_Vmax - ind$vsp), lvl)
    lvl <- ifelse(offtx & cd4 >= 200, ind$vsp, lvl)

    supp_now <- idx & on_art & !failing & (lvl <= LOG200)
    cstar <- pmin(p$cd4_ceiling, nadir + inc_of(nadir))
    cd4_art <- cstar - (cstar - cd4) * exp(-p$r_rec * DT)
    phi_eff <- phi_i * ifelse(ltfu & (t - t_ltfu < p$ltfu_phi_years), p$phi_ltfu_mult, 1)
    if (p$substeps == 1) {
      D <- pmax(0, lvl - LOG200)
      cd4_off <- cd4 + DT * (p$s_C * pmin(1, cd4 / p$kappa) - (p$d_C + phi_eff * D) * cd4)
    } else {
      h <- DT / p$substeps; cc <- cd4; v <- lvl
      for (s in 1:p$substeps) {
        Ds <- pmax(0, v - LOG200)
        cc <- pmax(cc + h * (p$s_C * pmin(1, cc / p$kappa) - (p$d_C + phi_eff * Ds) * cc), 1)
        v <- ifelse(offtx & cc < 200, ind$vsp + (200 - cc) / 200 * (p$log10_Vmax - ind$vsp), v)
      }
      cd4_off <- cc
    }
    cd4 <- ifelse(idx, pmax(ifelse(supp_now, cd4_art, cd4_off), 1), cd4)
    nadir <- pmin(nadir, cd4)
    age <- age + DT
  }
  if (any(alive)) { t_death[alive] <- Tn * DT; age_death[alive] <- age[alive]; yll[alive] <- 0 }

  list(yld = yld, yll = yll, dalys = yld + yll, yld_prediag = yld_prediag,
       cost_art_drug = c_drug, cost_art_svc = c_svc + c_start, cost_hiv_care = c_hc,
       cost_art = c_drug + c_svc + c_start, cost_total = c_drug + c_svc + c_start + c_hc,
       ever_dx = diagnosed, ever_art = ever_art, t_dx = t_dx, cd4_dx = cd4_dx,
       t_death = t_death, age_death = age_death, cause_hiv = cause_hiv,
       py_total = py_total, py_dx = py_dx, py_art = py_art, py_supp = py_supp,
       pt_year = pt_year, years = p$year + 0:(ny - 1))
}

run_acquisition_cohort <- function(p, arms = c("no_art", "real_world", "ideal")) {
  crn <- make_crn(p); ind <- init_individuals(p, crn)
  out <- lapply(arms, function(a) simulate_arm(p, crn, ind, a)); names(out) <- arms
  out$ind <- ind; out$crn <- crn; out
}

# ---------------------------------------------------------------- accounting
cost_daly_frame <- function(s) {
  rw <- s$real_world; na <- s$no_art
  data.frame(dalys_accrued = mean(rw$dalys), dalys_no_care = mean(na$dalys), dalys_averted_care = mean(na$dalys - rw$dalys),
             daly_ideal = if (!is.null(s$ideal)) mean(s$ideal$dalys) else NA,
             yld = mean(rw$yld), yll = mean(rw$yll), yld_prediag = mean(rw$yld_prediag),
             cost_art_drug = mean(rw$cost_art_drug), cost_art_svc = mean(rw$cost_art_svc), cost_hiv_care = mean(rw$cost_hiv_care),
             cost_art = mean(rw$cost_art), cost_total = mean(rw$cost_total), cost_no_care = mean(na$cost_total))
}
compute_icer <- function(no_art, real_world) {
  inc_cost <- mean(real_world$cost_total - no_art$cost_total); inc_daly <- mean(no_art$dalys - real_world$dalys)
  list(inc_cost = inc_cost, inc_daly = inc_daly, icer = inc_cost / inc_daly)
}
cascade_diagnostics <- function(rw) {
  dx <- mean(rw$ever_dx); art <- mean(rw$ever_art)
  list(pct_ever_dx = dx, pct_ever_art = art, art_given_dx = art / dx,
       med_cd4_at_dx = median(rw$cd4_dx, na.rm = TRUE), med_lag_yr = median(rw$t_dx, na.rm = TRUE),
       pt_dx = sum(rw$py_dx) / sum(rw$py_total), pt_art_given_dx = sum(rw$py_art) / sum(rw$py_dx),
       pt_supp_given_art = sum(rw$py_supp) / sum(rw$py_art))
}
# prevention averts the whole acquisition: (c_PrEP*NNP - lifetime cost)/DALYs accrued  (Table 1D: 465, -253)
compare_prevention_treatment <- function(s, c_prep, nnp, wtp = 500) {
  rw <- s$real_world; cost <- mean(rw$cost_total); dalys <- mean(rw$dalys); ic <- compute_icer(s$no_art, rw)
  c(list(prev_icer = (c_prep * nnp - cost) / dalys, net_cost_saving = (c_prep * nnp) < cost, treat_ce = ic$icer < wtp), ic)
}
