# Blinded Sample Size Re-estimation for Recurrent Events

Estimates the blinded dispersion and event rate from aggregated interim
data and calculates the required sample size to maintain power, assuming
the planned treatment effect holds. This function supports constant
rates (Friede & Schmidli 2010) and accommodates future extensions for
time-varying rates (Schneider et al. 2013) by using the
exposure-adjusted rate.

## Usage

``` r
blinded_ssr(
  data,
  ratio = 1,
  lambda1_planning,
  lambda2_planning,
  power = 0.8,
  alpha = 0.025,
  method = "friede",
  accrual_rate,
  accrual_duration,
  trial_duration,
  dropout_rate = 0,
  max_followup = NULL,
  event_gap = NULL
)
```

## Arguments

- data:

  A data frame containing the blinded interim data. Must include columns
  `events` (number of events) and `tte` (total exposure/follow-up time).
  This is typically the output of
  [`cut_data_by_date()`](https://keaven.github.io/nbDesign/reference/cut_data_by_date.md).

- ratio:

  Planned allocation ratio (experimental / control). Default is 1.

- lambda1_planning:

  Planned event rate for the control group used in original calculation.

- lambda2_planning:

  Planned event rate for the experimental group used in original
  calculation.

- power:

  Target power (1 - beta). Default is 0.8.

- alpha:

  One-sided significance level. Default is 0.025.

- method:

  Method for sample size recalculation. Currently "friede" (Friede &
  Schmidli 2010) is implemented, which uses the blinded nuisance
  parameter estimates.

## Value

A list containing:

- n_total_unadjusted:

  Original planned total sample size (based on planning parameters).

- n_total_blinded:

  Re-estimated total sample size using blinded estimates.

- dispersion_blinded:

  Estimated dispersion parameter (k) from blinded data.

- lambda_blinded:

  Estimated overall event rate from blinded data.

- info_fraction:

  Estimated information fraction at interim (blinded information /
  target information).

- blinded_info:

  Estimated statistical information from the blinded interim data.

- target_info:

  Target statistical information required for the planned power.
