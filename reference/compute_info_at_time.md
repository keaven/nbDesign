# Compute Statistical Information at Analysis Time

Computes the statistical information for the log rate ratio at a given
analysis time, accounting for staggered enrollment and varying exposure
times.

## Usage

``` r
compute_info_at_time(
  analysis_time,
  accrual_rate,
  accrual_duration,
  lambda1,
  lambda2,
  dispersion,
  ratio = 1
)
```

## Arguments

- analysis_time:

  The calendar time of the analysis.

- accrual_rate:

  The enrollment rate (subjects per time unit).

- accrual_duration:

  The duration of the enrollment period.

- lambda1:

  Event rate for group 1 (control).

- lambda2:

  Event rate for group 2 (treatment).

- dispersion:

  The negative binomial dispersion parameter.

- ratio:

  Allocation ratio (n2/n1). Default is 1.

## Value

The statistical information (inverse of variance) at the analysis time.

## Details

For subjects enrolled at time `t`, their exposure at analysis time `T`
is `T - t`. The function integrates over all enrolled subjects to
compute the total expected events and variance of the log rate ratio.

The variance formula for the log rate ratio is: \$\$Var = \sum_i
\frac{1/\mu\_{1i} + k}{n\_{1i}} + \frac{1/\mu\_{2i} + k}{n\_{2i}}\$\$

where \\\mu\_{ji} = \lambda_j \times exposure_i\\ and \\k\\ is the
dispersion.

## Examples

``` r
# Compute information at month 12 for a trial with:
# - 10 subjects/month enrollment rate
# - 20 month enrollment period
# - Control event rate 0.5, treatment rate 0.3
# - Dispersion 0.1
compute_info_at_time(
  analysis_time = 12,
  accrual_rate = 10,
  accrual_duration = 20,
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1
)
#> [1] 55.10204
```
