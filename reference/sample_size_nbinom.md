# Sample Size Calculation for Negative Binomial Distribution

Computes the sample size for comparing two treatment groups assuming a
negative binomial distribution for the outcome.

## Usage

``` r
sample_size_nbinom(
  lambda1,
  lambda2,
  dispersion,
  power = NULL,
  alpha = 0.025,
  sided = 1,
  ratio = 1,
  accrual_rate,
  accrual_duration,
  trial_duration,
  dropout_rate = 0,
  max_followup = NULL,
  event_gap = NULL,
  method = "zhu"
)
```

## Arguments

- lambda1:

  Rate in group 1 (control).

- lambda2:

  Rate in group 2 (treatment).

- dispersion:

  Dispersion parameter `k` such that \\Var(Y) = \mu + k \mu^2\\. Note
  that this is equivalent to `1/size` in R's `rnbinom` parameterization.

- power:

  Power of the test (1 - beta). Default is 0.9.

- alpha:

  Significance level. Default is 0.025.

- sided:

  One-sided or two-sided test. 1 for one-sided, 2 for two-sided. Default
  is 1.

- ratio:

  Allocation ratio n2/n1. Default is 1.

- accrual_rate:

  Vector of accrual rates (patients per unit time).

- accrual_duration:

  Vector of durations for each accrual rate. Must be same length as
  `accrual_rate`.

- trial_duration:

  Total planned duration of the trial.

- dropout_rate:

  Dropout rate (hazard rate). Default is 0.

- max_followup:

  Maximum follow-up time for any patient. Default is NULL (infinite).

- event_gap:

  Gap duration after each event during which no new events are counted.
  Default is NULL (no gap). If provided, the effective event rate is
  reduced.

- method:

  Method for sample size calculation. "zhu" for Zhu and Lakkis (2014),
  "friede" for Friede and Schmidli (2010) / Mütze et al. (2018).

## Value

A list containing:

- n1:

  Sample size for group 1

- n2:

  Sample size for group 2

- n_total:

  Total sample size

- exposure:

  Average exposure time used in calculation

## References

Zhu, H., & Lakkis, H. (2014). Sample size calculation for comparing two
negative binomial rates in clinical trials. *Statistics in
Biopharmaceutical Research*, 6(1), 107-115.

Friede, T., & Schmidli, H. (2010). Sample size estimation for clinical
trials with negative binomial rates. *Methods of Information in
Medicine*, 49(6), 623-631.

Mütze, T., Glimm, E., Schmidli, H., & Friede, T. (2018). Group
sequential designs for negative binomial outcomes. *Statistical Methods
in Medical Research*, 27(10), 2978-2993.

## Examples

``` r
# Calculate sample size for lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1
# with fixed recruitment of 10/month for 20 months, 24 month trial duration
sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
#> $inputs
#> $inputs$lambda1
#> [1] 0.5
#> 
#> $inputs$lambda2
#> [1] 0.3
#> 
#> $inputs$dispersion
#> [1] 0.1
#> 
#> $inputs$power
#> [1] 0.8
#> 
#> $inputs$alpha
#> [1] 0.025
#> 
#> $inputs$sided
#> [1] 1
#> 
#> $inputs$ratio
#> [1] 1
#> 
#> $inputs$accrual_rate
#> [1] 10
#> 
#> $inputs$accrual_duration
#> [1] 20
#> 
#> $inputs$trial_duration
#> [1] 24
#> 
#> $inputs$dropout_rate
#> [1] 0
#> 
#> $inputs$max_followup
#> NULL
#> 
#> $inputs$event_gap
#> NULL
#> 
#> $inputs$method
#> [1] "zhu"
#> 
#> 
#> $n1
#> [1] 18
#> 
#> $n2
#> [1] 18
#> 
#> $n_total
#> [1] 36
#> 
#> $alpha
#> [1] 0.025
#> 
#> $sided
#> [1] 1
#> 
#> $power
#> [1] 0.8
#> 
#> $exposure
#> [1] 14
#> 
#> $events_n1
#> [1] 126
#> 
#> $events_n2
#> [1] 75.6
#> 
#> $total_events
#> [1] 201.6
#> 
#> $variance
#> [1] 0.03227513
#> 
#> $accrual_rate
#> [1] 1.8
#> 
#> $accrual_duration
#> [1] 20
#> 
#> attr(,"class")
#> [1] "sample_size_nbinom_result" "list"                     

# With piecewise accrual
# 5 patients/month for 3 months, then 10 patients/month for 3 months
# Trial ends at month 12.
sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
  accrual_rate = c(5, 10), accrual_duration = c(3, 3),
  trial_duration = 12
)
#> $inputs
#> $inputs$lambda1
#> [1] 0.5
#> 
#> $inputs$lambda2
#> [1] 0.3
#> 
#> $inputs$dispersion
#> [1] 0.1
#> 
#> $inputs$power
#> [1] 0.8
#> 
#> $inputs$alpha
#> [1] 0.025
#> 
#> $inputs$sided
#> [1] 1
#> 
#> $inputs$ratio
#> [1] 1
#> 
#> $inputs$accrual_rate
#> [1]  5 10
#> 
#> $inputs$accrual_duration
#> [1] 3 3
#> 
#> $inputs$trial_duration
#> [1] 12
#> 
#> $inputs$dropout_rate
#> [1] 0
#> 
#> $inputs$max_followup
#> NULL
#> 
#> $inputs$event_gap
#> NULL
#> 
#> $inputs$method
#> [1] "zhu"
#> 
#> 
#> $n1
#> [1] 25
#> 
#> $n2
#> [1] 25
#> 
#> $n_total
#> [1] 50
#> 
#> $alpha
#> [1] 0.025
#> 
#> $sided
#> [1] 1
#> 
#> $power
#> [1] 0.8
#> 
#> $exposure
#> [1] 8.5
#> 
#> $events_n1
#> [1] 106.25
#> 
#> $events_n2
#> [1] 63.75
#> 
#> $total_events
#> [1] 170
#> 
#> $variance
#> [1] 0.03309804
#> 
#> $accrual_rate
#> [1]  5.555556 11.111111
#> 
#> $accrual_duration
#> [1] 3 3
#> 
#> attr(,"class")
#> [1] "sample_size_nbinom_result" "list"                     
```
