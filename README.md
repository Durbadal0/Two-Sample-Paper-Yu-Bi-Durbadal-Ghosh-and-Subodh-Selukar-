# Replication Code: A Comparative Study of Two-Sample Hypothesis Tests in the Presence of Long-Term Survivors

Bi Y, Ghosh D, Selukar S. *A comparative study of two-sample hypothesis tests in the presence of long-term survivors.*

## R Scripts

- **01_numerical_results.R** — Computes the weighted average hazard difference A(tau) for three scenarios and generates Figure 5 (main text Section 3.2).
- **02_supplement_scenarios.R** — Generates three-panel figures (survival, hazard, A(tau)/Delta H) for five theoretical scenarios (Supplementary Material Section B).
- **03_beacon_analysis.R** — Performs the complete BEACON-Immuno trial analysis: hypothesis tests, mixture cure model fitting with RECeUS-AIC, parameter estimation, and generates Figure 6 and supplement tables (main text Section 4 and Supplementary Material Section C).

## Data

- **ipd_pfs_jco.csv** — Reconstructed individual patient-level progression-free survival data from the BEACON-Immuno trial (NCT02308527), obtained using the WebDigitizer-based tool at https://biostatistics.mdanderson.org/shinyapps/IPDfromKM/.

## Required R Packages

```r
install.packages(c("survival", "flexsurv", "flexsurvcure", "ggplot2", "gridExtra"))
```

## How to Run

Run scripts from the project root directory:

```r
source("01_numerical_results.R")
source("02_supplement_scenarios.R")
source("03_beacon_analysis.R")
```

## Uploading to GitHub

1. Create a new repository on GitHub (https://github.com/new).

2. Initialize the local repo, add files, and push:

```bash
cd "/path/to/this/directory"
git init
git add 01_numerical_results.R 02_supplement_scenarios.R 03_beacon_analysis.R ipd_pfs_jco.csv README.md
git commit -m "Initial commit: replication code for two-sample tests with long-term survivors"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
git push -u origin main
```

Replace `YOUR_USERNAME` and `YOUR_REPO_NAME` with your GitHub username and repository name.
