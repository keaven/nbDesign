# Group Sequential Design and Simulation

``` r
library(gsDesignNB)
library(data.table)
library(ggplot2)
library(gt)
```

This vignette demonstrates how to create a group sequential design for
negative binomial outcomes using
[`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
and simulate the design to confirm design operating characteristicsusing
[`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md).

## Trial Design Parameters

We design a trial with the following characteristics:

- **Enrollment:** 1 year (12 months) with a constant rate
- **Trial Duration:** 2 years (24 months)
- **Analyses:**
  - Interim 1: 10 months
  - Interim 2: 18 months
  - Final: 24 months
- **Event Rates:**
  - Control: 1.5 events per year

  - Experimental: 1.0 events per year (rate ratio = 0.67)
- **Dispersion:** 0.5
- **Power:** 90% (beta = 0.1)
- **Alpha:** 0.025 (one-sided)

### Sample Size Calculation

First, we calculate the required sample size for a fixed design:

``` r
# Sample size calculation
# Enrollment: constant rate over 1 year
# Trial duration: 2 years
nb_ss <- sample_size_nbinom(
  lambda1 = 1.5,           # Control event rate (per year)
  lambda2 = 1.0,           # Experimental event rate (per year)
  dispersion = 0.5,        # Overdispersion parameter
  power = 0.9,             # 90% power
  alpha = 0.025,           # One-sided alpha
  accrual_rate = 100,      # Patients per year (will determine total n)
  accrual_duration = 1,    # 1 year enrollment
  trial_duration = 2,      # 2 year trial
  max_followup = 1,        # 1 year of follow-up per patient
  event_gap = 20 / 365.25, # Minimum gap between events is 20 days
  method = "zhu"           # Zhu and Lakkis sample size method
)

# Print key results
cat("Fixed design\n")
#> Fixed design
nb_ss
#> Sample Size for Negative Binomial Outcome
#> ==========================================
#> 
#> Sample size:     n1 = 178, n2 = 178, total = 356
#> Expected events: 415.5 (n1: 246.7, n2: 168.8)
#> Power: 90%, Alpha: 0.025 (1-sided)
#> Rates: control = 1.5000, treatment = 1.0000 (RR = 0.6667)
#> Dispersion: 0.5000, Avg exposure: 1.00
#> Accrual: 1.0, Trial duration: 2.0
```

### Group Sequential Design

Now we convert this a group sequential design with 3 analyses after 10,
18 and 24 months. Note that the final analysis time must be the same as
for the fixed design. The relative enrollment rates will be increased to
increase the sample size as with standard group sequential design
theory. We specify `usTime = c(0.1, 0.2, 1)` which along with the
`sfLinear()` spending function will spend 10%, 20% and 100% of the
cumulative $\alpha$ at the 3 planned analyses regardless of the observed
statistical information at each analysis.

``` r
# Analysis times (in years, for simulation)
analysis_times <- c(10, 18, 24) / 12  # Convert months to years

# Create group sequential design with integer sample sizes
gs_nb <- gsNBCalendar(
  x = nb_ss,               # Input fixed design for negative binomial
  k = 3,                   # 3 analyses
  test.type = 4,           # Two-sided asymmetric, non-binding futility
  sfu = sfLinear,          # Linear spending function for upper bound
  sfupar = c(.5, .5),      # Identity function
  sfl = sfLinear,          # Same lower spending
  sflpar = c(.5, .5),
  usTime = c(.1, .2, 1),   # Conservative efficacy spending: 10%, 20% at IAs
  lsTime = c(.1, .4, 1),   # Futility spending more aggressive at IA2
  analysis_times = analysis_times * 12  # Calendar times in months
) |> toInteger() # Round to integer sample size
```

Textual group sequential design summary:

``` r
summary(gs_nb)
#> Asymmetric two-sided with non-binding futility bound group sequential design
#> for negative binomial outcomes, 3 analyses, total sample size 378.0, 90 percent
#> power, 2.5 percent (1-sided) Type I error. Control rate 1.5000, treatment rate
#> 1.0000, risk ratio 0.6667, dispersion 0.5000. Accrual duration 1.0, trial
#> duration 2.0, average exposure 1.00. Randomization ratio 1:1.
```

Tabular summary:

``` r
gs_nb |> 
  gsDesign::gsBoundSummary(
    deltaname = "RR", 
    logdelta = TRUE,
    Nname = "Information",
    timename = "Year",
    digits = 4,
    ddigits = 2
  ) |>
  gt() |>
  tab_header(
    title = "Group Sequential Design Bounds for Negative Binomial Outcome",
    subtitle = paste0("N = ", ceiling(gs_nb$n_total[gs_nb$k]), 
                      ", Expected events = ", round(gs_nb$nb_design$total_events, 1))
  )
```

| Group Sequential Design Bounds for Negative Binomial Outcome |                     |          |          |
|--------------------------------------------------------------|---------------------|----------|----------|
| N = 378, Expected events = 415.5                             |                     |          |          |
| Analysis                                                     | Value               | Efficacy | Futility |
| IA 1: 91%                                                    | Z                   | 2.8070   | -0.4001  |
| Information: 160.79                                          | p (1-sided)         | 0.0025   | 0.6554   |
|                                                              | ~RR at bound        | 0.8011   | 1.0321   |
|                                                              | P(Cross) if RR=1    | 0.0025   | 0.3446   |
|                                                              | P(Cross) if RR=0.67 | 0.1892   | 0.0100   |
| IA 2: 98%                                                    | Z                   | 2.7403   | 0.9175   |
| Information: 172.57                                          | p (1-sided)         | 0.0031   | 0.1794   |
|                                                              | ~RR at bound        | 0.8114   | 0.9324   |
|                                                              | P(Cross) if RR=1    | 0.0050   | 0.8259   |
|                                                              | P(Cross) if RR=0.67 | 0.5121   | 0.0400   |
| Final                                                        | Z                   | 1.9964   | 1.9964   |
| Information: 176.48                                          | p (1-sided)         | 0.0229   | 0.0229   |
|                                                              | ~RR at bound        | 0.8603   | 0.8603   |
|                                                              | P(Cross) if RR=1    | 0.0235   | 0.9765   |
|                                                              | P(Cross) if RR=0.67 | 0.9000   | 0.1000   |

## Simulation Study

We now simulate 50 trials to evaluate the operating characteristics of
the group sequential design.

### Simulation Setup

``` r
set.seed(42)
n_sims <- 50

# Enrollment rate (patients per year) to achieve target sample size
n_target <- ceiling(nb_ss$n_total)
enroll_rate_val <- n_target / 1  # All enrolled in 1 year

# Define enrollment
enroll_rate <- data.frame(
  rate = enroll_rate_val,
  duration = 1  # 1 year enrollment
)

# Define failure rates (with dispersion)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(1.5, 1.0),
  dispersion = c(0.5, 0.5)
)

# No dropout for simplicity
dropout_rate <- NULL

# Maximum follow-up (trial duration minus minimum enrollment time)
max_followup <- 2  # 2 years from enrollment start
```

### Run Simulations

``` r
# Storage for results
results <- vector("list", n_sims)

for (sim in 1:n_sims) {
  # Simulate trial data
  sim_data <- nb_sim(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate,
    max_followup = max_followup,
    n = n_target
  )
  
  # Analyze at each interim
  sim_results <- data.frame(
    sim = sim,
    analysis = 1:3,
    analysis_time = analysis_times,
    n_enrolled = NA_integer_,
    events_ctrl = NA_integer_,
    events_exp = NA_integer_,
    events_total = NA_integer_,
    exposure_ctrl = NA_real_,
    exposure_exp = NA_real_,
    z_stat = NA_real_,
    p_value = NA_real_,
    cross_upper = NA,
    cross_lower = NA
  )
  
  stopped <- FALSE
  
  for (k in 1:3) {
    if (stopped) {
      # Trial already stopped at earlier analysis
      sim_results$cross_upper[k] <- FALSE
      sim_results$cross_lower[k] <- FALSE
      next
    }
    
    # Cut data at analysis time
    cut_time <- analysis_times[k]
    cut_data <- cut_data_by_date(sim_data, cut_date = cut_time)
    
    # Count enrolled subjects (those with enroll_time <= cut_time)
    enrolled <- unique(sim_data$id[sim_data$enroll_time <= cut_time])
    cut_data <- cut_data[cut_data$id %in% enrolled, ]
    
    # Summary by treatment
    summary_dt <- as.data.table(cut_data)[
      ,
      .(n = .N, events = sum(events), exposure = sum(tte)),
      by = treatment
    ]
    
    ctrl_row <- summary_dt[treatment == "Control"]
    exp_row <- summary_dt[treatment == "Experimental"]
    
    sim_results$n_enrolled[k] <- nrow(cut_data)
    sim_results$events_ctrl[k] <- if (nrow(ctrl_row) > 0) ctrl_row$events else 0
    sim_results$events_exp[k] <- if (nrow(exp_row) > 0) exp_row$events else 0
    sim_results$events_total[k] <- sim_results$events_ctrl[k] + sim_results$events_exp[k]
    sim_results$exposure_ctrl[k] <- if (nrow(ctrl_row) > 0) ctrl_row$exposure else 0
    sim_results$exposure_exp[k] <- if (nrow(exp_row) > 0) exp_row$exposure else 0
    
    # Run Mütze test
    if (nrow(cut_data) >= 4 && sim_results$events_total[k] >= 2) {
      test_result <- tryCatch(
        mutze_test(cut_data),
        error = function(e) NULL
      )
      
      if (!is.null(test_result)) {
        sim_results$z_stat[k] <- test_result$z
        sim_results$p_value[k] <- test_result$p_value
        
        # Check boundaries (one-sided: reject if z < -upper bound for benefit)
        # For rate ratio < 1 (experimental better), log(RR) < 0, so z < 0
        z_eff <- -test_result$z  # Flip sign for efficacy direction
        
        sim_results$cross_upper[k] <- z_eff > gs_nb$upper$bound[k]
        sim_results$cross_lower[k] <- z_eff < gs_nb$lower$bound[k]
        
        if (sim_results$cross_upper[k] || sim_results$cross_lower[k]) {
          stopped <- TRUE
        }
      }
    }
  }
  
  results[[sim]] <- sim_results
}

# Combine all results
all_results <- do.call(rbind, results)
```

## Simulation Results Summary

### Events and Exposure by Analysis

``` r
# Summarize by analysis
summary_by_analysis <- as.data.table(all_results)[
  ,
  .(
    mean_enrolled = mean(n_enrolled, na.rm = TRUE),
    mean_events_total = mean(events_total, na.rm = TRUE),
    mean_events_ctrl = mean(events_ctrl, na.rm = TRUE),
    mean_events_exp = mean(events_exp, na.rm = TRUE),
    mean_exposure_ctrl = mean(exposure_ctrl, na.rm = TRUE),
    mean_exposure_exp = mean(exposure_exp, na.rm = TRUE),
    mean_z = mean(z_stat, na.rm = TRUE),
    sd_z = sd(z_stat, na.rm = TRUE)
  ),
  by = .(analysis, analysis_time)
]

summary_by_analysis |>
  gt() |>
  tab_header(title = "Summary Statistics by Analysis") |>
  cols_label(
    analysis = "Analysis",
    analysis_time = "Time (yrs)",
    mean_enrolled = "N Enrolled",
    mean_events_total = "Total Events",
    mean_events_ctrl = "Ctrl Events",
    mean_events_exp = "Exp Events",
    mean_exposure_ctrl = "Ctrl Exposure",
    mean_exposure_exp = "Exp Exposure",
    mean_z = "Mean Z",
    sd_z = "SD Z"
  ) |>
  fmt_number(decimals = 2)
```

| Summary Statistics by Analysis |            |            |              |             |            |               |              |        |      |
|--------------------------------|------------|------------|--------------|-------------|------------|---------------|--------------|--------|------|
| Analysis                       | Time (yrs) | N Enrolled | Total Events | Ctrl Events | Exp Events | Ctrl Exposure | Exp Exposure | Mean Z | SD Z |
| 1.00                           | 0.83       | 297.20     | 153.54       | 91.52       | 62.02      | 60.21         | 60.61        | −2.10  | 0.87 |
| 2.00                           | 1.50       | 356.00     | 435.36       | 257.69      | 177.67     | 174.08        | 175.16       | −3.08  | 0.83 |
| 3.00                           | 2.00       | 356.00     | 657.31       | 380.12      | 277.19     | 261.45        | 262.85       | −2.92  | 0.47 |

### Statistical Information

The statistical information at each analysis is proportional to the
precision of the treatment effect estimate. For negative binomial
outcomes, this relates to the total exposure and event counts.

``` r
# Information proxy: inverse variance of log rate ratio
# For large samples: Var(log RR) ≈ 1/events_ctrl + 1/events_exp
info_by_analysis <- as.data.table(all_results)[
  ,
  .(
    mean_info = mean(1 / (1/events_ctrl + 1/events_exp), na.rm = TRUE),
    planned_info_frac = unique(gs_nb$timing[analysis])
  ),
  by = analysis
]

# Normalize to get observed information fractions
info_by_analysis[, observed_info_frac := mean_info / max(mean_info)]

info_by_analysis |>
  gt() |>
  tab_header(title = "Information by Analysis") |>
  cols_label(
    analysis = "Analysis",
    mean_info = "Mean Information",
    planned_info_frac = "Planned Info Frac",
    observed_info_frac = "Observed Info Frac"
  ) |>
  fmt_number(decimals = 3)
```

| Information by Analysis |                  |                   |                    |
|-------------------------|------------------|-------------------|--------------------|
| Analysis                | Mean Information | Planned Info Frac | Observed Info Frac |
| 1.000                   | 36.708           | 0.911             | 0.229              |
| 2.000                   | 104.918          | 0.978             | 0.655              |
| 3.000                   | 160.199          | 1.000             | 1.000              |

### Boundary Crossings and Power

``` r
# Calculate crossing probabilities
crossing_summary <- as.data.table(all_results)[
  ,
  .(
    n_cross_upper = sum(cross_upper, na.rm = TRUE),
    n_cross_lower = sum(cross_lower, na.rm = TRUE),
    n_continue = sum(!cross_upper & !cross_lower, na.rm = TRUE)
  ),
  by = analysis
]

crossing_summary[, prob_cross_upper := n_cross_upper / n_sims]
crossing_summary[, prob_cross_lower := n_cross_lower / n_sims]

crossing_summary[, .(analysis, n_cross_upper, prob_cross_upper, n_cross_lower, prob_cross_lower)] |>
  gt() |>
  tab_header(title = "Boundary Crossing by Analysis") |>
  cols_label(
    analysis = "Analysis",
    n_cross_upper = "N Cross Upper",
    prob_cross_upper = "P(Cross Upper)",
    n_cross_lower = "N Cross Lower",
    prob_cross_lower = "P(Cross Lower)"
  ) |>
  fmt_number(columns = starts_with("prob"), decimals = 3)
```

| Boundary Crossing by Analysis |               |                |               |                |
|-------------------------------|---------------|----------------|---------------|----------------|
| Analysis                      | N Cross Upper | P(Cross Upper) | N Cross Lower | P(Cross Lower) |
| 1                             | 11            | 0.220          | 0             | 0.000          |
| 2                             | 23            | 0.460          | 0             | 0.000          |
| 3                             | 16            | 0.320          | 0             | 0.000          |

### Overall Power

``` r
# Determine if each simulation crossed the efficacy boundary at any analysis
efficacy_by_sim <- as.data.table(all_results)[
  ,
  .(efficacy = any(cross_upper, na.rm = TRUE)),
  by = sim
]

overall_power <- mean(efficacy_by_sim$efficacy, na.rm = TRUE)

# Futility stopping
futility_by_sim <- as.data.table(all_results)[
  ,
  .(futility = any(cross_lower, na.rm = TRUE) & !any(cross_upper, na.rm = TRUE)),
  by = sim
]

overall_futility <- mean(futility_by_sim$futility, na.rm = TRUE)

cat("\n=== Overall Operating Characteristics ===\n")
#> 
#> === Overall Operating Characteristics ===
cat(sprintf("Number of simulations: %d\n", n_sims))
#> Number of simulations: 50
cat(sprintf("Overall Power (P[reject H0]): %.1f%%\n", overall_power * 100))
#> Overall Power (P[reject H0]): 100.0%
cat(sprintf("Futility Stopping Rate: %.1f%%\n", overall_futility * 100))
#> Futility Stopping Rate: 0.0%
cat(sprintf("Design Power (target): %.1f%%\n", (1 - gs_nb$beta) * 100))
#> Design Power (target): 90.0%
```

### Visualization of Z-Statistics

``` r
# Prepare data for plotting
plot_data <- all_results
plot_data$z_flipped <- -plot_data$z_stat  # Flip for efficacy direction

# Boundary data
bounds_df <- data.frame(
  analysis = 1:3,
  upper = gs_nb$upper$bound,
  lower = gs_nb$lower$bound
)

ggplot(plot_data, aes(x = factor(analysis), y = z_flipped)) +
  geom_violin(fill = "steelblue", alpha = 0.5, color = "steelblue") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_hline(data = bounds_df, aes(yintercept = upper), 
             linetype = "dashed", color = "darkgreen", linewidth = 1) +
  geom_hline(data = bounds_df, aes(yintercept = lower), 
             linetype = "dashed", color = "darkred", linewidth = 1) +
  geom_hline(yintercept = 0, color = "gray50") +
  labs(
    title = "Simulated Z-Statistics by Analysis",
    subtitle = "Green dashed = efficacy bound, Red dashed = futility bound",
    x = "Analysis",
    y = "Z-statistic (positive = favors experimental)"
  ) +
  theme_minimal() +
  ylim(c(-4, 6))
#> Warning: Removed 45 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 45 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
```

![Z-statistics across analyses with group sequential
boundaries](group_sequential_simulation_files/figure-html/plot-z-stats-1.png)

## Design Summary

The [`summary()`](https://rdrr.io/r/base/summary.html) function provides
a textual overview of the group sequential design:

``` r
summary(gs_nb)
#> Asymmetric two-sided with non-binding futility bound group sequential design
#> for negative binomial outcomes, 3 analyses, total sample size 378.0, 90 percent
#> power, 2.5 percent (1-sided) Type I error. Control rate 1.5000, treatment rate
#> 1.0000, risk ratio 0.6667, dispersion 0.5000. Accrual duration 1.0, trial
#> duration 2.0, average exposure 1.00. Randomization ratio 1:1.
```

For detailed boundary information, use `gsBoundSummary()`. We set
`logdelta = TRUE` since the test statistic is based on the log rate
ratio. The `~RR at bound` displays the rate ratio
($\lambda_{2}/\lambda_{1}$), where values \< 1 indicate treatment
benefit:

``` r
gsDesign::gsBoundSummary(gs_nb, 
                         deltaname = "RR", 
                         logdelta = TRUE,
                         Nname = "Information",
                         timename = "Month",
                         digits = 4,
                         ddigits = 2) |>
  gt() |>
  tab_header(
    title = "Group Sequential Design Bounds for Negative Binomial Outcome",
    subtitle = paste0("N = ", ceiling(gs_nb$n_total[gs_nb$k]), 
                      ", Expected events = ", round(gs_nb$nb_design$total_events, 1))
  )
```

| Group Sequential Design Bounds for Negative Binomial Outcome |                     |          |          |
|--------------------------------------------------------------|---------------------|----------|----------|
| N = 378, Expected events = 415.5                             |                     |          |          |
| Analysis                                                     | Value               | Efficacy | Futility |
| IA 1: 91%                                                    | Z                   | 2.8070   | -0.4001  |
| Information: 160.79                                          | p (1-sided)         | 0.0025   | 0.6554   |
|                                                              | ~RR at bound        | 0.8011   | 1.0321   |
|                                                              | P(Cross) if RR=1    | 0.0025   | 0.3446   |
|                                                              | P(Cross) if RR=0.67 | 0.1892   | 0.0100   |
| IA 2: 98%                                                    | Z                   | 2.7403   | 0.9175   |
| Information: 172.57                                          | p (1-sided)         | 0.0031   | 0.1794   |
|                                                              | ~RR at bound        | 0.8114   | 0.9324   |
|                                                              | P(Cross) if RR=1    | 0.0050   | 0.8259   |
|                                                              | P(Cross) if RR=0.67 | 0.5121   | 0.0400   |
| Final                                                        | Z                   | 1.9964   | 1.9964   |
| Information: 176.48                                          | p (1-sided)         | 0.0229   | 0.0229   |
|                                                              | ~RR at bound        | 0.8603   | 0.8603   |
|                                                              | P(Cross) if RR=1    | 0.0235   | 0.9765   |
|                                                              | P(Cross) if RR=0.67 | 0.9000   | 0.1000   |

Note that `P(Cross) if RR=0.67` corresponds to the design’s alternate
hypothesis (treatment rate / control rate = 0.67).

Sample sizes at each analysis:

``` r
data.frame(
  Analysis = 1:gs_nb$k,
  n1 = gs_nb$n1,
  n2 = gs_nb$n2,
  n_total = gs_nb$n_total
) |>
  gt() |>
  tab_header(title = "Sample Sizes at Each Analysis") |>
  fmt_number(columns = c(n1, n2, n_total), decimals = 1)
```

| Sample Sizes at Each Analysis |       |       |         |
|-------------------------------|-------|-------|---------|
| Analysis                      | n1    | n2    | n_total |
| 1                             | 63.0  | 63.0  | 126.0   |
| 2                             | 126.0 | 126.0 | 252.0   |
| 3                             | 189.0 | 189.0 | 378.0   |

After rounding to integer sample sizes with
[`toInteger()`](https://keaven.github.io/gsDesignNB/reference/toInteger.md):

``` r
gs_nb_int <- toInteger(gs_nb)
summary(gs_nb_int)
#> Asymmetric two-sided with non-binding futility bound group sequential design
#> for negative binomial outcomes, 3 analyses, total sample size 378.0, 90 percent
#> power, 2.5 percent (1-sided) Type I error. Control rate 1.5000, treatment rate
#> 1.0000, risk ratio 0.6667, dispersion 0.5000. Accrual duration 1.0, trial
#> duration 2.0, average exposure 1.00. Randomization ratio 1:1.
```

## Notes

This simulation demonstrates the basic workflow for group sequential
designs with negative binomial outcomes:

1.  **Sample size calculation** using
    [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
    for a fixed design
2.  **Group sequential design** using
    [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
    to add interim analyses
3.  **Simulation** using
    [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md)
    to generate trial data
4.  **Analysis** using
    [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md)
    and
    [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md)
    at each interim
5.  **Boundary checking** against the group sequential bounds

The `usTime = c(0.1, 0.2, 1)` specification provides conservative alpha
spending at early analyses, preserving most of the Type I error for
later analyses when more information is available.

With only 50 simulations, the estimated power has substantial Monte
Carlo error. For more precise estimates, increase `n_sims` to 1000 or
more.
