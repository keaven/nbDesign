# Wald test for treatment effect using Negative Binomial model (Mütze et al.)

Fits a Negative Binomial (or Poisson) log-rate model to the aggregated
subject-level data produced by
[`cut_data_by_date()`](https://keaven.github.io/nbDesign/reference/cut_data_by_date.md).
The method matches the Wald test described by Mütze et al. (2018) for
comparing treatment arms with recurrent event outcomes.

## Usage

``` r
mutze_test(data, method = c("nb", "poisson"), conf_level = 0.95)
```

## Arguments

- data:

  A data frame with at least the columns `treatment`, `events`, and
  `tte` (follow-up time). Typically output from
  [`cut_data_by_date()`](https://keaven.github.io/nbDesign/reference/cut_data_by_date.md).

- method:

  Type of model to fit: "nb" (default) uses a Negative Binomial GLM via
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html),
  "poisson" fits a Poisson GLM.

- conf_level:

  Confidence level for the rate ratio interval. Default 0.95.

## Value

A list containing the fitted model summary with elements:

- `estimate`: log rate ratio (experimental vs control).

- `se`: standard error for the log rate ratio.

- `z`: Wald statistic.

- `p_value`: two-sided p-value.

- `rate_ratio`: estimated rate ratio and its confidence interval.

- `dispersion`: estimated dispersion (theta) when `method = "nb"`.

- `group_summary`: observed subjects/events/exposure per treatment.

## Examples

``` r
enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.1, 0.05), duration = c(100, 100)
)
sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 40)
cut <- cut_data_by_date(sim, cut_date = 1.5)
mutze_test(cut)
#> $method
#> [1] "Negative Binomial Wald"
#> 
#> $estimate
#> [1] -0.5353677
#> 
#> $se
#> [1] 0.4022541
#> 
#> $z
#> [1] -1.330919
#> 
#> $p_value
#> [1] 0.1832157
#> 
#> $rate_ratio
#> [1] 0.585454
#> 
#> $conf_int
#> [1] 0.2661298 1.2879295
#> 
#> $conf_level
#> [1] 0.95
#> 
#> $dispersion
#> [1] 36.08517
#> 
#> $model
#> 
#> Call:  MASS::glm.nb(formula = events ~ treatment + offset(log(tte)), 
#>     data = df, init.theta = 36.08517111, link = log)
#> 
#> Coefficients:
#>           (Intercept)  treatmentExperimental  
#>               -0.1949                -0.5354  
#> 
#> Degrees of Freedom: 39 Total (i.e. Null);  38 Residual
#> Null Deviance:       42.72 
#> Residual Deviance: 40.89     AIC: 89.83
#> 
#> $group_summary
#>      treatment subjects events exposure
#> 1 Experimental       20     10 20.78206
#> 2      Control       20     17 20.64381
#> 
```
