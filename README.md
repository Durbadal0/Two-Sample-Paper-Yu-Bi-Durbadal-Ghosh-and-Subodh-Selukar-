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



