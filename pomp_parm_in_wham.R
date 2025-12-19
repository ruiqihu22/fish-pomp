library(wham)

# ============================================================
# Load Data
# ============================================================

cat("\n")
cat("================================================================\n")
cat("EVALUATE POMP PARAMETERS IN WHAM MODEL\n")
cat("================================================================\n\n")

# Load original WHAM model for reference
load("vign_10_mod_1.RData")

# Get WHAM parameter estimates for comparison
parList <- mod_1$parList
rep <- mod_1$rep
n_years <- mod_1$input$data$n_years_model

wham_mean_R <- exp(parList$mean_rec_pars[1,1])
wham_sigma_R <- exp(parList$log_NAA_sigma[1])
wham_q <- 1/(1 + exp(-parList$logit_q[1]))
wham_F <- exp(mean(log(exp(cumsum(c(parList$F_pars[1,1], parList$F_pars[2:n_years,1]))))))
wham_N1 <- exp(as.vector(parList$log_N1)[1])

loglik_wham <- -(sum(rep$nll_agg_catch) + sum(rep$nll_agg_indices))

cat("WHAM Parameter Estimates (reference):\n")
cat("  mean_R:  ", round(wham_mean_R, 1), "\n")
cat("  sigma_R: ", round(wham_sigma_R, 4), "\n")
cat("  F:       ", round(wham_F, 4), "\n")
cat("  q:       ", format(wham_q, scientific=TRUE, digits=4), "\n")
cat("  N1_init: ", round(wham_N1, 1), "\n")
cat("  loglik:  ", round(loglik_wham, 2), "\n\n")

# ============================================================
# POMP Parameters (replace with your MIF2 results)
# ============================================================

# Example POMP MIF2 estimates - REPLACE WITH YOUR ACTUAL VALUES
pomp_mean_R <- 9500
pomp_sigma_R <- 0.7
pomp_F <- 0.55
pomp_q <- 5e-4
pomp_N1 <- 28000

cat("POMP Parameter Estimates (to evaluate):\n")
cat("  mean_R:  ", round(pomp_mean_R, 1), "\n")
cat("  sigma_R: ", round(pomp_sigma_R, 4), "\n")
cat("  F:       ", round(pomp_F, 4), "\n")
cat("  q:       ", format(pomp_q, scientific=TRUE, digits=4), "\n")
cat("  N1_init: ", round(pomp_N1, 1), "\n\n")

# ============================================================
# Create New WHAM Input
# ============================================================

cat("================================================================\n")
cat("SETTING UP WHAM WITH POMP PARAMETERS\n")
cat("================================================================\n\n")

# Read ASAP data
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples, "ex1_SNEMAYT.dat"))

# Create WHAM input (same structure as original)
input <- prepare_wham_input(
  asap3,
  recruit_model = 2,
  model_name = "WHAM with POMP params",
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

# ============================================================
# Check Parameter Dimensions
# ============================================================

cat("Parameter dimensions:\n")
cat("  mean_rec_pars: ", paste(dim(input$par$mean_rec_pars), collapse=" x "), 
    " (length:", length(input$par$mean_rec_pars), ")\n")
cat("  log_NAA_sigma: ", paste(dim(input$par$log_NAA_sigma), collapse=" x "), 
    " (length:", length(input$par$log_NAA_sigma), ")\n")
cat("  logit_q:       ", paste(dim(input$par$logit_q), collapse=" x "), 
    " (length:", length(input$par$logit_q), ")\n")
cat("  F_pars:        ", paste(dim(input$par$F_pars), collapse=" x "), 
    " (length:", length(input$par$F_pars), ")\n")
cat("  log_N1:        ", paste(dim(input$par$log_N1), collapse=" x "), 
    " (length:", length(input$par$log_N1), ")\n")
cat("\n")

# ============================================================
# Set Parameters to POMP Values
# ============================================================

cat("Setting parameters to POMP values...\n\n")

# 1. Mean recruitment (log scale)
input$par$mean_rec_pars[1,1] <- log(pomp_mean_R)

# 2. Recruitment sigma (log scale)
input$par$log_NAA_sigma[1] <- log(pomp_sigma_R)

# 3. Catchability (logit scale)
# logit(q) = log(q / (1-q))
input$par$logit_q[1] <- log(pomp_q / (1 - pomp_q))
input$par$logit_q[2] <- log(pomp_q / (1 - pomp_q))  # Set both indices

# 4. F_pars: first value is log(F_1), rest are deviations
# For constant F: set F_1 = pomp_F, all deviations = 0
input$par$F_pars[1, 1] <- log(pomp_F)
for(i in 2:n_years) {
  input$par$F_pars[i, 1] <- 0  # No deviations = constant F
}

# 5. Initial N at age 1
# log_N1 structure depends on WHAM version
if(length(dim(input$par$log_N1)) == 3) {
  # 3D array: [1, 1, n_ages]
  input$par$log_N1[1, 1, 1] <- log(pomp_N1)
} else if(length(dim(input$par$log_N1)) == 1 || is.vector(input$par$log_N1)) {
  # Vector
  input$par$log_N1[1] <- log(pomp_N1)
} else {
  # 2D matrix
  input$par$log_N1[1, 1] <- log(pomp_N1)
}

# ============================================================
# FIX PARAMETERS USING MAP (Correct Length!)
# ============================================================

cat("Fixing parameters with correct map lengths...\n\n")

# The map must have the SAME LENGTH as the parameter
# Use factor(rep(NA, length)) to fix ALL elements

# 1. Fix mean_rec_pars
n_mean_rec <- length(input$par$mean_rec_pars)
input$map$mean_rec_pars <- factor(rep(NA, n_mean_rec))
cat("  mean_rec_pars map length:", n_mean_rec, "\n")

# 2. Fix log_NAA_sigma
n_sigma <- length(input$par$log_NAA_sigma)
input$map$log_NAA_sigma <- factor(rep(NA, n_sigma))
cat("  log_NAA_sigma map length:", n_sigma, "\n")

# 3. Fix logit_q
n_q <- length(input$par$logit_q)
input$map$logit_q <- factor(rep(NA, n_q))
cat("  logit_q map length:", n_q, "\n")

# 4. Fix F_pars
n_F <- length(input$par$F_pars)
input$map$F_pars <- factor(rep(NA, n_F))
cat("  F_pars map length:", n_F, "\n")

# 5. Fix log_N1
n_N1 <- length(input$par$log_N1)
input$map$log_N1 <- factor(rep(NA, n_N1))
cat("  log_N1 map length:", n_N1, "\n")

cat("\nAll parameters fixed.\n\n")

# ============================================================
# Evaluate Likelihood (Don't Optimize Fixed Params)
# ============================================================

cat("================================================================\n")
cat("EVALUATING WHAM LIKELIHOOD\n")
cat("================================================================\n\n")

cat("Running WHAM model (with fixed parameters)...\n")
cat("Only selectivity parameters will be optimized.\n\n")

# Fit model - with parameters fixed, this mainly evaluates likelihood
# and optimizes only the non-fixed parameters (selectivity)
mod_pomp <- tryCatch({
  fit_wham(input, do.osa = FALSE, do.retro = FALSE, do.sdrep = FALSE)
}, error = function(e) {
  cat("Error in fit_wham:", e$message, "\n")
  return(NULL)
})

if(!is.null(mod_pomp)) {
  # Extract likelihood components
  nll_total <- mod_pomp$opt$objective
  nll_catch <- sum(mod_pomp$rep$nll_agg_catch)
  nll_index <- sum(mod_pomp$rep$nll_agg_indices)
  loglik_wham_at_pomp <- -(nll_catch + nll_index)
  
  cat("=== RESULTS ===\n\n")
  cat("  NLL (catch):           ", round(nll_catch, 2), "\n")
  cat("  NLL (index):           ", round(nll_index, 2), "\n")
  cat("  loglik (catch+index):  ", round(loglik_wham_at_pomp, 2), "\n")
  cat("  NLL (total):           ", round(nll_total, 2), "\n\n")
  
  # ============================================================
  # Comparison
  # ============================================================
  
  cat("================================================================\n")
  cat("COMPARISON\n")
  cat("================================================================\n\n")
  
  cat("  loglik_wham (WHAM params):     ", round(loglik_wham, 2), "\n")
  cat("  loglik_wham_at_pomp (POMP params): ", round(loglik_wham_at_pomp, 2), "\n")
  cat("  Difference:                    ", round(loglik_wham_at_pomp - loglik_wham, 2), "\n\n")
  
  if(loglik_wham_at_pomp < loglik_wham) {
    cat("Interpretation: WHAM parameters achieve BETTER likelihood\n")
    cat("  (Expected, since WHAM optimized for its own likelihood)\n")
  } else {
    cat("Interpretation: POMP parameters achieve BETTER likelihood in WHAM\n")
    cat("  (Unexpected - consider checking parameter mapping)\n")
  }
  
  # Save results
  results <- list(
    pomp_params = list(
      mean_R = pomp_mean_R,
      sigma_R = pomp_sigma_R,
      F = pomp_F,
      q = pomp_q,
      N1_init = pomp_N1
    ),
    loglik_wham = loglik_wham,
    loglik_wham_at_pomp = loglik_wham_at_pomp,
    nll_catch = nll_catch,
    nll_index = nll_index,
    SSB_pomp = mod_pomp$rep$SSB
  )
  
  save(results, file = "wham_at_pomp_results.RData")
  cat("\nResults saved to wham_at_pomp_results.RData\n")
  
} else {
  cat("Model fitting failed. Check parameter settings.\n")
}

cat("\n================================================================\n")
cat("DONE\n")
cat("================================================================\n")