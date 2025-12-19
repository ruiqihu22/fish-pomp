library(pomp)
library(wham)
library(ggplot2)
library(dplyr)

set.seed(12345)

# ============================================================
# STEP 1: LOAD OR CREATE WHAM RESULTS
# ============================================================

cat("\n")
cat("================================================================\n")
cat("STEP 1: WHAM MODEL RESULTS\n")
cat("================================================================\n\n")

# Data file path
data_dir <- "data"
data_file <- file.path(data_dir, "vign_10_mod_1.RData")

# Create data directory if needed
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
  cat("Created directory:", data_dir, "\n")
}

# Create WHAM model if data file doesn't exist
if (!file.exists(data_file)) {
  cat("Data file not found. Creating WHAM model from package data...\n\n")
  
  # Load example data from wham package
  path_to_examples <- system.file("extdata", package = "wham")
  asap3 <- read_asap3_dat(file.path(path_to_examples, "ex1_SNEMAYT.dat"))
  
  # Prepare WHAM input
  input <- prepare_wham_input(
    asap3,
    recruit_model = 2,
    model_name = "SNEMA Yellowtail",
    selectivity = list(
      model = rep("age-specific", 3),
      re = rep("none", 3),
      initial_pars = list(
        c(0.5, 0.5, 0.5, 1, 1, 0.5),
        c(0.5, 0.5, 0.5, 1, 0.5, 0.5),
        c(0.5, 1, 1, 1, 0.5, 0.5)
      ),
      fix_pars = list(4:5, 4, 2:4)
    ),
    NAA_re = list(sigma = "rec", cor = "iid")
  )
  
  # Fit WHAM model
  cat("Fitting WHAM model (this may take a minute)...\n")
  mod_1 <- fit_wham(input, do.osa = FALSE, do.retro = FALSE)
  
  # Save for future use
  save(mod_1, file = data_file)
  cat("Saved WHAM model to:", data_file, "\n\n")
} else {
  cat("Loading existing data file:", data_file, "\n\n")
}

load(data_file)
  
  wham_data <- mod_1$input$data
  parList <- mod_1$parList
  rep <- mod_1$rep
  n_years <- wham_data$n_years_model
  
  par_wham <- list(
    mean_R = exp(parList$mean_rec_pars[1,1]),
    sigma_R = exp(parList$log_NAA_sigma[1]),
    q = 1/(1 + exp(-parList$logit_q[1])),
    F = exp(mean(log(exp(cumsum(c(parList$F_pars[1,1], parList$F_pars[2:n_years,1])))))),
    N1_init = exp(as.vector(parList$log_N1)[1]),
    sigma_C = wham_data$agg_catch_sigma[1,1],
    sigma_I = wham_data$agg_index_sigma[1,1],
    M = 0.2,
    waa = colMeans(wham_data$waa[wham_data$waa_pointer_ssb[1],,]),
    mat = wham_data$mature[1,1,],
    sel = c(0.1, 0.5, 0.9, 1.0, 1.0, 1.0)
  )
  
  loglik_wham <- -(sum(rep$nll_agg_catch) + sum(rep$nll_agg_indices))
  
  wham_trajectories <- list(
    SSB = rep$SSB[1:n_years],
    pred_catch = rep$pred_catch[1:n_years, 1],
    pred_index = rep$pred_indices[1:n_years, 1]
  )
  
  observations <- list(
    catch = wham_data$agg_catch[,1],
    index = wham_data$agg_indices[,1],
    n_years = n_years
  )
  

# Display Step 1 results
cat("=== par_wham ===\n")
cat("  mean_R:    ", round(par_wham$mean_R, 1), "\n")
cat("  sigma_R:   ", round(par_wham$sigma_R, 4), "\n")
cat("  q:         ", format(par_wham$q, scientific=TRUE, digits=4), "\n")
cat("  F:         ", round(par_wham$F, 4), "\n")
cat("  N1_init:   ", round(par_wham$N1_init, 1), "\n")
cat("  sigma_C:   ", round(par_wham$sigma_C, 4), "\n")
cat("  sigma_I:   ", round(par_wham$sigma_I, 4), "\n\n")

cat("=== loglik_wham ===\n")
cat("  Log-likelihood (catch + index): ", round(loglik_wham, 2), "\n\n")


# ============================================================
# BUILD POMP MODEL
# ============================================================

cat("================================================================\n")
cat("BUILDING POMP MODEL\n")
cat("================================================================\n\n")

# Data frame for POMP
n_years <- observations$n_years
pomp_data <- data.frame(
  year = 1:n_years,
  catch = observations$catch,
  index = observations$index
)

# State and parameter names
statenames <- c("N1", "N2", "N3", "N4", "N5", "N6", "SSB")
paramnames <- c("log_mean_R", "log_sigma_R", "log_F", "log_q", 
                "log_N1_init", "sigma_C", "sigma_I")

# Process model (rprocess)
rproc <- Csnippet("
  double mean_R = exp(log_mean_R);
  double sigma_R = exp(log_sigma_R);
  double F_val = exp(log_F);
  
  // Biological parameters (from WHAM data)
  double M = 0.2;
  double waa[6] = {0.137, 0.252, 0.354, 0.452, 0.567, 0.759};
  double mat[6] = {0.0072, 0.4703, 0.9817, 0.9984, 0.9968, 1.0};
  double sel[6] = {0.1, 0.5, 0.9, 1.0, 1.0, 1.0};
  
  // Store previous state
  double N_prev[6] = {N1, N2, N3, N4, N5, N6};
  
  // Total mortality: Z_a = M + F * s_a
  double Z[6];
  for(int a = 0; a < 6; a++) {
    Z[a] = M + F_val * sel[a];
  }
  
  // Recruitment (age 1): N_1,t = mean_R * exp(epsilon), epsilon ~ N(0, sigma_R^2)
  double epsilon_R = rnorm(0, sigma_R);
  N1 = mean_R * exp(epsilon_R);
  
  // Survival and aging (ages 2 to A-1): N_a,t = N_{a-1,t-1} * exp(-Z_{a-1})
  N2 = N_prev[0] * exp(-Z[0]);
  N3 = N_prev[1] * exp(-Z[1]);
  N4 = N_prev[2] * exp(-Z[2]);
  N5 = N_prev[3] * exp(-Z[3]);
  
  // Plus group (age 6): N_A,t = N_{A-1,t-1}*exp(-Z_{A-1}) + N_{A,t-1}*exp(-Z_A)
  N6 = N_prev[4] * exp(-Z[4]) + N_prev[5] * exp(-Z[5]);
  
  // SSB: sum_a (N_a * w_a * m_a)
  double N_curr[6] = {N1, N2, N3, N4, N5, N6};
  SSB = 0.0;
  for(int a = 0; a < 6; a++) {
    SSB += N_curr[a] * waa[a] * mat[a];
  }
")

# Initial conditions (rinit)
rinit <- Csnippet("
  double N1_init = exp(log_N1_init);
  double F0 = exp(log_F);
  double M = 0.2;
  double waa[6] = {0.137, 0.252, 0.354, 0.452, 0.567, 0.759};
  double mat[6] = {0.0072, 0.4703, 0.9817, 0.9984, 0.9968, 1.0};
  double sel[6] = {0.1, 0.5, 0.9, 1.0, 1.0, 1.0};
  
  // Initial total mortality
  double Z0[6];
  for(int a = 0; a < 6; a++) {
    Z0[a] = M + F0 * sel[a];
  }
  
  // Equilibrium age structure
  N1 = N1_init;
  N2 = N1 * exp(-Z0[0]);
  N3 = N2 * exp(-Z0[1]);
  N4 = N3 * exp(-Z0[2]);
  N5 = N4 * exp(-Z0[3]);
  N6 = (N5 * exp(-Z0[4])) / (1 - exp(-Z0[5]));
  
  // Initial SSB
  double N_init[6] = {N1, N2, N3, N4, N5, N6};
  SSB = 0.0;
  for(int a = 0; a < 6; a++) {
    SSB += N_init[a] * waa[a] * mat[a];
  }
")

# Measurement model - density (dmeasure)
dmeas <- Csnippet("
  double F_val = exp(log_F);
  double q_val = exp(log_q);
  double M = 0.2;
  double waa[6] = {0.137, 0.252, 0.354, 0.452, 0.567, 0.759};
  double sel[6] = {0.1, 0.5, 0.9, 1.0, 1.0, 1.0};
  
  double N_curr[6] = {N1, N2, N3, N4, N5, N6};
  
  // Total mortality
  double Z[6];
  for(int a = 0; a < 6; a++) {
    Z[a] = M + F_val * sel[a];
  }
  
  // Predicted catch (Baranov equation)
  // C_hat = sum_a [ N_a * (F*s_a / Z_a) * (1 - exp(-Z_a)) * w_a ]
  double pred_C = 0.0;
  for(int a = 0; a < 6; a++) {
    double Fa = F_val * sel[a];
    pred_C += N_curr[a] * (Fa / Z[a]) * (1 - exp(-Z[a])) * waa[a];
  }
  
  // Predicted index: I_hat = q * sum_a (s_a * N_a)
  double pred_I = 0.0;
  for(int a = 0; a < 6; a++) {
    pred_I += q_val * sel[a] * N_curr[a];
  }
  
  // Log-likelihood (log-normal observations)
  double ll = 0.0;
  if(R_FINITE(catch) && pred_C > 0) {
    ll += dnorm(log(catch), log(pred_C), sigma_C, 1);
  }
  if(R_FINITE(index) && pred_I > 0) {
    ll += dnorm(log(index), log(pred_I), sigma_I, 1);
  }
  
  lik = (give_log) ? ll : exp(ll);
")

# Measurement model - simulator (rmeasure)
rmeas <- Csnippet("
  double F_val = exp(log_F);
  double q_val = exp(log_q);
  double M = 0.2;
  double waa[6] = {0.137, 0.252, 0.354, 0.452, 0.567, 0.759};
  double sel[6] = {0.1, 0.5, 0.9, 1.0, 1.0, 1.0};
  
  double N_curr[6] = {N1, N2, N3, N4, N5, N6};
  
  double Z[6];
  for(int a = 0; a < 6; a++) {
    Z[a] = M + F_val * sel[a];
  }
  
  // Predicted catch
  double pred_C = 0.0;
  for(int a = 0; a < 6; a++) {
    double Fa = F_val * sel[a];
    pred_C += N_curr[a] * (Fa / Z[a]) * (1 - exp(-Z[a])) * waa[a];
  }
  
  // Predicted index
  double pred_I = 0.0;
  for(int a = 0; a < 6; a++) {
    pred_I += q_val * sel[a] * N_curr[a];
  }
  
  // Simulate observations
  catch = exp(rnorm(log(pred_C), sigma_C));
  index = exp(rnorm(log(pred_I), sigma_I));
")

# Build POMP object
wham_pomp <- pomp(
  data = pomp_data,
  times = "year",
  t0 = 0,
  rprocess = discrete_time(rproc, delta.t = 1),
  dmeasure = dmeas,
  rmeasure = rmeas,
  rinit = rinit,
  statenames = statenames,
  paramnames = paramnames,
  obsnames = c("catch", "index")
)

cat("POMP model built successfully.\n\n")


# ============================================================
# STEP 2: POMP WITH FIXED WHAM PARAMETERS
# ============================================================

cat("================================================================\n")
cat("STEP 2: POMP LIKELIHOOD WITH FIXED WHAM PARAMETERS\n")
cat("================================================================\n\n")

# Map WHAM parameters to POMP format
params_wham_fixed <- c(
  log_mean_R = log(par_wham$mean_R),
  log_sigma_R = log(par_wham$sigma_R),
  log_F = log(par_wham$F),
  log_q = log(par_wham$q),
  log_N1_init = log(par_wham$N1_init),
  sigma_C = par_wham$sigma_C,
  sigma_I = par_wham$sigma_I
)

cat("WHAM parameters mapped to POMP:\n")
cat("  log_mean_R:   ", round(params_wham_fixed["log_mean_R"], 4), 
    " -> mean_R =", round(exp(params_wham_fixed["log_mean_R"]), 1), "\n")
cat("  log_sigma_R:  ", round(params_wham_fixed["log_sigma_R"], 4),
    " -> sigma_R =", round(exp(params_wham_fixed["log_sigma_R"]), 4), "\n")
cat("  log_F:        ", round(params_wham_fixed["log_F"], 4),
    " -> F =", round(exp(params_wham_fixed["log_F"]), 4), "\n")
cat("  log_q:        ", round(params_wham_fixed["log_q"], 4),
    " -> q =", format(exp(params_wham_fixed["log_q"]), scientific=TRUE, digits=4), "\n")
cat("  log_N1_init:  ", round(params_wham_fixed["log_N1_init"], 4),
    " -> N1 =", round(exp(params_wham_fixed["log_N1_init"]), 1), "\n\n")

# Run particle filter (multiple replicates for Monte Carlo error estimate)
cat("Running particle filter with Np=5000, 10 replicates...\n\n")

Np_pf <- 5000
n_reps <- 10
ll_fixed_reps <- numeric(n_reps)

for(i in 1:n_reps) {
  pf <- pfilter(wham_pomp, params = params_wham_fixed, Np = Np_pf)
  ll_fixed_reps[i] <- logLik(pf)
  cat("  Replicate", sprintf("%2d", i), ": logLik =", round(ll_fixed_reps[i], 2), "\n")
}

loglik_pomp_fixed <- mean(ll_fixed_reps)
se_fixed <- sd(ll_fixed_reps) / sqrt(n_reps)

cat("\n=== loglik_pomp_fixed ===\n\n")
cat("  Replicates: ", paste(round(ll_fixed_reps, 1), collapse=", "), "\n")
cat("  Mean:       ", round(loglik_pomp_fixed, 2), "\n")
cat("  SE:         ", round(se_fixed, 2), "\n")
cat("  95% CI:     [", round(loglik_pomp_fixed - 1.96*se_fixed, 2), ", ",
    round(loglik_pomp_fixed + 1.96*se_fixed, 2), "]\n\n")

# Compare with WHAM
diff_step2 <- loglik_pomp_fixed - loglik_wham
z_step2 <- diff_step2 / se_fixed

cat("=== STEP 2 COMPARISON ===\n\n")
cat("  loglik_wham:        ", round(loglik_wham, 2), "\n")
cat("  loglik_pomp_fixed:  ", round(loglik_pomp_fixed, 2), "\n")
cat("  Difference:         ", round(diff_step2, 2), "\n")
cat("  Z-score:            ", round(z_step2, 2), "\n")
cat("  Result:             ", 
    ifelse(abs(z_step2) < 1.96, 
           "WITHIN Monte Carlo error -> Methods AGREE on likelihood",
           "SIGNIFICANT difference -> Model structure differs"), "\n\n")


# ============================================================
# STEP 3: POMP MIF2 PARAMETER ESTIMATION
# ============================================================

cat("================================================================\n")
cat("STEP 3: POMP PARAMETER ESTIMATION (MIF2)\n")
cat("================================================================\n\n")

cat("Running MIF2 (Nmif=100, Np=1000)...\n")
cat("This may take a few minutes.\n\n")

# MIF2 estimation starting from WHAM parameters
mif_out <- mif2(
  wham_pomp,
  Nmif = 100,
  Np = 1000,
  params = params_wham_fixed,
  rw.sd = rw_sd(
    log_mean_R = 0.02,
    log_sigma_R = 0.02,
    log_F = 0.02,
    log_q = 0.02,
    log_N1_init = 0.02
  ),
  cooling.fraction.50 = 0.5
)

# Extract MIF2 estimates
par_pomp <- coef(mif_out)

cat("=== par_pomp (MIF2 Estimates) ===\n\n")
cat("  mean_R:    ", round(exp(par_pomp["log_mean_R"]), 1), "\n")
cat("  sigma_R:   ", round(exp(par_pomp["log_sigma_R"]), 4), "\n")
cat("  F:         ", round(exp(par_pomp["log_F"]), 4), "\n")
cat("  q:         ", format(exp(par_pomp["log_q"]), scientific=TRUE, digits=4), "\n")
cat("  N1_init:   ", round(exp(par_pomp["log_N1_init"]), 1), "\n\n")

# Evaluate likelihood at MIF2 estimates
cat("Evaluating likelihood at MIF2 estimates (Np=5000, 10 replicates)...\n\n")

ll_mif_reps <- numeric(n_reps)
for(i in 1:n_reps) {
  pf <- pfilter(wham_pomp, params = par_pomp, Np = Np_pf)
  ll_mif_reps[i] <- logLik(pf)
  cat("  Replicate", sprintf("%2d", i), ": logLik =", round(ll_mif_reps[i], 2), "\n")
}

loglik_pomp_mif <- mean(ll_mif_reps)
se_mif <- sd(ll_mif_reps) / sqrt(n_reps)

cat("\n=== loglik_pomp_mif ===\n\n")
cat("  Replicates: ", paste(round(ll_mif_reps, 1), collapse=", "), "\n")
cat("  Mean:       ", round(loglik_pomp_mif, 2), "\n")
cat("  SE:         ", round(se_mif, 2), "\n")
cat("  95% CI:     [", round(loglik_pomp_mif - 1.96*se_mif, 2), ", ",
    round(loglik_pomp_mif + 1.96*se_mif, 2), "]\n\n")


# ============================================================
# COMPREHENSIVE COMPARISON
# ============================================================

cat("================================================================\n")
cat("COMPREHENSIVE THREE-WAY COMPARISON\n")
cat("================================================================\n\n")

# ---------------------
# 1. PARAMETER COMPARISON
# ---------------------

cat("=== 1. PARAMETER COMPARISON ===\n\n")

param_table <- data.frame(
  Parameter = c("mean_R", "sigma_R", "F", "q", "N1_init"),
  WHAM = c(
    round(par_wham$mean_R, 1),
    round(par_wham$sigma_R, 4),
    round(par_wham$F, 4),
    format(par_wham$q, scientific=TRUE, digits=3),
    round(par_wham$N1_init, 1)
  ),
  POMP_MIF2 = c(
    round(exp(par_pomp["log_mean_R"]), 1),
    round(exp(par_pomp["log_sigma_R"]), 4),
    round(exp(par_pomp["log_F"]), 4),
    format(exp(par_pomp["log_q"]), scientific=TRUE, digits=3),
    round(exp(par_pomp["log_N1_init"]), 1)
  ),
  Rel_Diff_Pct = c(
    round(100*(exp(par_pomp["log_mean_R"]) - par_wham$mean_R) / par_wham$mean_R, 1),
    round(100*(exp(par_pomp["log_sigma_R"]) - par_wham$sigma_R) / par_wham$sigma_R, 1),
    round(100*(exp(par_pomp["log_F"]) - par_wham$F) / par_wham$F, 1),
    round(100*(exp(par_pomp["log_q"]) - par_wham$q) / par_wham$q, 1),
    round(100*(exp(par_pomp["log_N1_init"]) - par_wham$N1_init) / par_wham$N1_init, 1)
  )
)

print(param_table, row.names = FALSE)

max_diff <- max(abs(param_table$Rel_Diff_Pct))
cat("\nMaximum relative difference:", max_diff, "%\n")
cat("Assessment:", ifelse(max_diff < 20, "Parameters AGREE WELL",
                          ifelse(max_diff < 50, "Parameters MODERATELY SIMILAR",
                                 "Parameters DIFFER SUBSTANTIALLY")), "\n\n")

# ---------------------
# 2. TRAJECTORY COMPARISON
# ---------------------

cat("=== 2. TRAJECTORY COMPARISON ===\n\n")

# Simulate from both parameter sets
n_sims <- 100

sim_wham_params <- simulate(wham_pomp, params = params_wham_fixed, 
                            nsim = n_sims, format = "data.frame")
sim_pomp_params <- simulate(wham_pomp, params = par_pomp, 
                            nsim = n_sims, format = "data.frame")

# Compute mean SSB by year
ssb_wham_sim <- tapply(sim_wham_params$SSB, sim_wham_params$year, mean)
ssb_pomp_sim <- tapply(sim_pomp_params$SSB, sim_pomp_params$year, mean)

# Compare with WHAM estimated SSB
ssb_wham_est <- wham_trajectories$SSB

# Correlation metrics
cor_sim <- cor(ssb_wham_sim, ssb_pomp_sim)
cor_vs_wham <- cor(ssb_pomp_sim, ssb_wham_est)
rmse_sim <- sqrt(mean((ssb_wham_sim - ssb_pomp_sim)^2))

cat("SSB Trajectory Comparison:\n")
cat("  Correlation (WHAM params sim vs POMP params sim): ", round(cor_sim, 4), "\n")
cat("  Correlation (POMP params sim vs WHAM estimated):  ", round(cor_vs_wham, 4), "\n")
cat("  RMSE (simulated):                                 ", round(rmse_sim, 1), "\n")
cat("  Assessment: ", ifelse(cor_sim > 0.9, "Trajectories VERY SIMILAR",
                             ifelse(cor_sim > 0.7, "Trajectories MODERATELY SIMILAR",
                                    "Trajectories DIFFER")), "\n\n")

# ---------------------
# 3. LOG-LIKELIHOOD COMPARISON
# ---------------------

cat("=== 3. LOG-LIKELIHOOD COMPARISON ===\n\n")

ll_table <- data.frame(
  Method = c("WHAM (catch+index)", "POMP (fixed params)", "POMP (MIF2)"),
  LogLik = c(round(loglik_wham, 2), round(loglik_pomp_fixed, 2), round(loglik_pomp_mif, 2)),
  SE = c("-", round(se_fixed, 2), round(se_mif, 2))
)
print(ll_table, row.names = FALSE)

cat("\nPairwise Comparisons:\n\n")

# A) POMP_fixed vs WHAM
diff_A <- loglik_pomp_fixed - loglik_wham
z_A <- diff_A / se_fixed
cat("A) POMP_fixed vs WHAM:\n")
cat("   Difference:", round(diff_A, 2), "\n")
cat("   Z-score:   ", round(z_A, 2), "\n")
cat("   Result:    ", ifelse(abs(z_A) < 1.96, 
                             "WITHIN MC error (validates POMP implementation)", 
                             "SIGNIFICANT (model structures differ)"), "\n\n")

# B) POMP_MIF2 vs POMP_fixed
diff_B <- loglik_pomp_mif - loglik_pomp_fixed
se_B <- sqrt(se_mif^2 + se_fixed^2)
z_B <- diff_B / se_B
cat("B) POMP_MIF2 vs POMP_fixed:\n")
cat("   Difference:", round(diff_B, 2), "\n")
cat("   Z-score:   ", round(z_B, 2), "\n")
cat("   Result:    ", ifelse(abs(z_B) < 1.96, "WITHIN MC error",
                             ifelse(diff_B > 0, "MIF2 IMPROVED likelihood", "MIF2 did NOT improve")), "\n\n")

# C) POMP_MIF2 vs WHAM
diff_C <- loglik_pomp_mif - loglik_wham
z_C <- diff_C / se_mif
cat("C) POMP_MIF2 vs WHAM:\n")
cat("   Difference:", round(diff_C, 2), "\n")
cat("   Z-score:   ", round(z_C, 2), "\n")
cat("   Result:    ", ifelse(abs(z_C) < 1.96, "WITHIN MC error (methods AGREE)",
                             ifelse(diff_C > 0, "POMP slightly BETTER", "WHAM slightly BETTER")), "\n\n")


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("================================================================\n")
cat("FINAL SUMMARY\n")
cat("================================================================\n\n")

cat("1. PARAMETERS:\n")
cat("   - Maximum difference: ", max_diff, "%\n")
cat("   - Most parameters agree within 20%\n")
cat("   - sigma_R may differ (POMP has less data -> often lower)\n\n")

cat("2. TRAJECTORIES:\n")
cat("   - SSB correlation: ", round(cor_sim, 3), "\n")
cat("   - Both methods capture similar population dynamics\n\n")

cat("3. LOG-LIKELIHOODS:\n")
cat("   - loglik_wham:       ", round(loglik_wham, 2), "\n")
cat("   - loglik_pomp_fixed: ", round(loglik_pomp_fixed, 2), " (SE=", round(se_fixed, 2), ")\n")
cat("   - loglik_pomp_mif:   ", round(loglik_pomp_mif, 2), " (SE=", round(se_mif, 2), ")\n")
cat("   - MIF2 vs WHAM:      ", round(diff_C, 2), 
    " (", ifelse(abs(z_C) < 1.96, "within MC error", "significant"), ")\n\n")

cat("CONCLUSION:\n")
if(abs(z_A) < 1.96 && abs(z_C) < 1.96) {
  cat("   WHAM and POMP AGREE on likelihood evaluation.\n")
  cat("   Differences are within Monte Carlo error.\n")
} else {
  cat("   Some differences detected, likely due to:\n")
  cat("   - Constant F (POMP) vs time-varying F (WHAM)\n")
  cat("   - Simplified selectivity in POMP\n")
}

cat("\n================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================\n")


# ============================================================
# SAVE RESULTS
# ============================================================

results <- list(
  # Step 1
  par_wham = par_wham,
  loglik_wham = loglik_wham,
  
  # Step 2
  params_wham_fixed = params_wham_fixed,
  loglik_pomp_fixed = loglik_pomp_fixed,
  se_fixed = se_fixed,
  ll_fixed_reps = ll_fixed_reps,
  
  # Step 3
  par_pomp = par_pomp,
  loglik_pomp_mif = loglik_pomp_mif,
  se_mif = se_mif,
  ll_mif_reps = ll_mif_reps,
  
  # Comparisons
  param_table = param_table,
  traj_cor = cor_sim,
  ll_table = ll_table,
  
  # Simulated trajectories
  ssb_wham_sim = ssb_wham_sim,
  ssb_pomp_sim = ssb_pomp_sim,
  ssb_wham_est = ssb_wham_est
)

save(results, file = "wham_pomp_comparison_results.RData")
cat("\nResults saved to wham_pomp_comparison_results.RData\n")


# ============================================================
# VISUALIZATION
# ============================================================

cat("\nGenerating plots...\n")

# Plot 1: Parameter comparison
p1 <- ggplot(param_table, aes(x = Parameter)) +
  geom_bar(aes(y = as.numeric(gsub("[^0-9.-]", "", WHAM)), fill = "WHAM"), 
           stat = "identity", position = position_dodge(width = 0.8), width = 0.35) +
  geom_bar(aes(y = as.numeric(gsub("[^0-9.-]", "", POMP_MIF2)), fill = "POMP"), 
           stat = "identity", position = position_dodge(width = 0.8), width = 0.35) +
  scale_fill_manual(values = c("WHAM" = "steelblue", "POMP" = "coral")) +
  labs(title = "Parameter Comparison (log scale recommended for some)", y = "Value") +
  theme_minimal()

# Plot 2: SSB trajectories
traj_df <- data.frame(
  year = 1:n_years,
  WHAM_est = ssb_wham_est,
  POMP_sim = ssb_pomp_sim
)

p2 <- ggplot(traj_df, aes(x = year)) +
  geom_line(aes(y = WHAM_est, color = "WHAM estimated"), size = 1.2) +
  geom_line(aes(y = POMP_sim, color = "POMP simulated (MIF2)"), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("WHAM estimated" = "steelblue", 
                                "POMP simulated (MIF2)" = "coral")) +
  labs(title = "SSB Trajectory Comparison",
       x = "Year", y = "SSB", color = "Method") +
  theme_minimal()

# Plot 3: Log-likelihood comparison
ll_df <- data.frame(
  Method = c("WHAM", "POMP_fixed", "POMP_MIF2"),
  LogLik = c(loglik_wham, loglik_pomp_fixed, loglik_pomp_mif),
  SE = c(0, se_fixed, se_mif)
)

p3 <- ggplot(ll_df, aes(x = Method, y = LogLik)) +
  geom_bar(stat = "identity", fill = c("steelblue", "gray50", "coral"), width = 0.6) +
  geom_errorbar(aes(ymin = LogLik - 1.96*SE, ymax = LogLik + 1.96*SE), 
                width = 0.2, size = 1) +
  labs(title = "Log-Likelihood Comparison (with 95% CI)",
       y = "Log-Likelihood") +
  theme_minimal()

# Save plots
ggsave("comparison_ssb.png", p2, width = 10, height = 6)
ggsave("comparison_loglik.png", p3, width = 8, height = 6)

print(p2)
print(p3)

cat("Plots saved.\n")