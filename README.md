# HIV per-acquisition DALY microsimulation — rebuild

Reconstruction of the individual-level Monte Carlo model behind the per-acquisition DALY / cost paper
(eastern & southern Africa), rebuilt from `docs/MODEL_SPEC_rebuild.md` after the original R engine was lost.

- `python/engine.py` — reference implementation (numpy); `run_checks.py` prints model-vs-paper for Table 1/2;
  `run_all.py {psa|fig2|fig3|all}` reproduces the PSA and Figures 2–3; `calibrate.py` is the grid search that chose the four calibrated switches.
- `R/model_functions.R`, `R/run_all.R` — line-for-line R port (base R). `Rscript run_all.R checks` from `R/`.
- `bg_life_table.csv` — UN WPP 2024 ESA × GBD 2021 HIV-deleted single-year qx.
- `out/` — PSA draws (200×1000, seed 42), Figure 2/3 data, calibration grid.

Structural decisions settled during the rebuild: log-logistic CD4 at acquisition (~550); YLL counted for
HIV-attributable deaths only; reference age 72 = expected age at death at acquisition (not e0); percentile
coverage gate ≡ Bernoulli gate when p1 is constant. Calibrated against recorded outputs (not in the param table):
diagnostic hazard 4/yr after threshold crossing, on-ART SMR floor 1.20, second-line failure permanent, attained-age
mortality band. Residual: contemporary-arm costs ~5% low (HIV-care cost; lower-arm ART time).
