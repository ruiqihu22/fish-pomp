# Comparison of Stock Assessment Model Inference: WHAM in TMB versus POMP Frameworks

This repository contains code and analysis for comparing the Woods Hole Assessment Model (WHAM) implemented in TMB with an equivalent POMP (Partially Observed Markov Process) formulation using iterated filtering (IF2).

## Repository Structure

```
fish-pomp/
├── Makefile                    # Build automation
├── README.md                   # This file
├── thesis.Rmd                  # Main thesis document
├── references.bib              # Bibliography
├── wham_parm_in_pomp.R         # Evaluate WHAM params in POMP
├── pomp_parm_in_wham.R         # Evaluate POMP params in WHAM
```

## Requirements

Install required R packages:

```r
install.packages(c("pomp", "ggplot2", "dplyr", "knitr", "rmarkdown"))

# Install WHAM from GitHub
remotes::install_github("timjmiller/wham", dependencies = TRUE)
```

## Reproducing the Analysis

### Option 1: Using Make (Recommended)

```bash
# Clone the repository
git clone git@github.com:ruiqihu22/fish-pomp.git
cd fish-pomp

# Build everything
make all

# Or step by step:
make wham_in_pomp_results.rds # Direction 1 analysis
make pomp_in_wham_results.rds # Direction 2 analysis
make thesis.pdf               # Render thesis
```

### Option 2: Direct R commands

```r
# Step 1: Evaluate WHAM parameters in POMP
source("wham_parm_in_pomp.R")

# Step 2: Evaluate POMP parameters in WHAM
source("pomp_parm_in_wham.R")

# Render thesis
rmarkdown::render("thesis.Rmd", output_format = "pdf_document")
```

## Cached Results

Results are cached in `.rds` files to avoid re-running time-consuming computations:

- `wham_in_pomp_results.rds`: Direction 1 results (WHAM params → POMP likelihood)
- `pomp_in_wham_results.rds`: Direction 2 results (POMP params → WHAM likelihood)

To force re-computation:

```bash
make clean    # Remove cached results (keeps wham_model.rds)
make cleanall # Remove all cached results
make fresh    # Clean and rebuild everything
```

## Analysis Overview

### Direction 1: WHAM → POMP
Evaluate TMB-optimized parameters in the POMP framework using particle filtering, then run IF2 to estimate parameters directly in POMP.

### Direction 2: POMP → WHAM
Evaluate POMP/IF2-estimated parameters in the WHAM/TMB framework.


## Author

Ruiqi Hu

## References

- Stock, B.C., & Miller, T.J. (2021). The Woods Hole Assessment Model (WHAM). *Fisheries Research*, 240, 105959.
- King, A.A., et al. (2016). Statistical inference for partially observed Markov processes via the R package pomp. *Journal of Statistical Software*, 69(12).
- Ionides, E.L., et al. (2015). Inference for dynamic and latent variable models via iterated, perturbed Bayes maps. *PNAS*, 112(3), 719-724.
