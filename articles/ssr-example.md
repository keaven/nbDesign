# Sample size re-estimation example

``` r
library(gsDesignNB)
library(data.table)
library(MASS)
```

## Introduction

This vignette demonstrates a group sequential trial with a single
interim analysis where information-based sample size re-estimation (SSR)
is performed. We follow the methodology suggested by Friede and Schmidli
(2010) for blinded sample size re-estimation with negative binomial
counts.

We simulate a single trial where the actual dispersion parameter is
higher than the planned dispersion parameter. This scenario typically
leads to a loss of power if the sample size is not adjusted. We will
show how to increase the sample size based on both blinded and unblinded
interim data to maintain the desired power.

## Trial setup and initial design

We plan a trial to compare two treatment groups (Control
vs. Experimental) with respect to recurrent event rates.

**Planned parameters:**

- Control rate (\\\lambda_1\\): 2.0 events/year
- Experimental rate (\\\lambda_2\\): 1.5 events/year (Hazard Ratio =
  0.75)
- Dispersion (\\k\\): 0.5
- Power: 90%
- One-sided Type I error (\\\alpha\\): 0.025
- Enrollment: 20 patients/month for 12 months
- Study duration: 24 months (12 months accrual + 12 months follow-up)

**Actual parameters (simulation truth):**

- Dispersion (\\k\\): 0.8 (Higher than planned, implying higher
  variability)
- Rates are as planned.

### Initial sample size calculation

First, we calculate the required sample size under the *planned*
parameters.

``` r
# Planned parameters
lambda1_plan <- 2.0
lambda2_plan <- 1.5
k_plan <- 0.5
power_plan <- 0.9
alpha_plan <- 0.025
accrual_rate_plan <- 20
accrual_dur_plan <- 12
trial_dur_plan <- 24

# Calculate sample size
design_plan <- sample_size_nbinom(
  lambda1 = lambda1_plan,
  lambda2 = lambda2_plan,
  dispersion = k_plan,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan
)

print(design_plan)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 140, n2 = 140, total = 280
#> Expected events: 8820.0 (n1: 5040.0, n2: 3780.0)
#> Power: 90%, Alpha: 0.025 (1-sided)
#> Rates: control = 2.0000, treatment = 1.5000 (RR = 0.7500)
#> Dispersion: 0.5000, Avg exposure (calendar): 18.00
#> Accrual: 12.0, Trial duration: 24.0
```

The design requires a total of 280 patients.

## Simulation of a single trial

We simulate a single realization of the trial using the *actual*
parameters (higher dispersion).

``` r
set.seed(1234)

# Actual parameters
k_true <- 0.8
lambda1_true <- 2.0
lambda2_true <- 1.5

# Enrollment and rates for simulation
enroll_rate <- data.frame(rate = accrual_rate_plan, duration = accrual_dur_plan)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(lambda1_true, lambda2_true),
  dispersion = k_true
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0, 0),
  duration = c(100, 100)
)

# Simulate trial
# We simulate enough patients to cover the planned sample size and potential increase
# Let's simulate a large pool and then we will cut it.
sim_data <- nb_sim(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = trial_dur_plan, # Max follow-up is effectively trial duration here
  n = 400 # Simulate more than planned ~200 to allow for SSR increase
)

head(sim_data)
#>   id id    treatment enroll_time        tte calendar_time event
#> 1  1  1 Experimental   0.1250879 0.05426115     0.1793491     1
#> 2  1  1 Experimental   0.1250879 0.45662549     0.5817134     1
#> 3  1  1 Experimental   0.1250879 0.86838555     0.9934735     1
#> 4  1  1 Experimental   0.1250879 1.13483540     1.2599233     1
#> 5  1  1 Experimental   0.1250879 1.58643344     1.7115214     1
#> 6  1  1 Experimental   0.1250879 1.58735933     1.7124473     1
```

## Interim analysis

We perform an interim analysis 9 months after the start of the trial. At
this point, enrollment is still ongoing (planned 12 months).

``` r
interim_time <- 9

# Cut data at interim
interim_data <- cut_data_by_date(sim_data, cut_date = interim_time)

# Summary of interim data
table(interim_data$treatment)
#> 
#>      Control Experimental 
#>           90           91
sum(interim_data$events)
#> [1] 1347
```

## Blinded sample size re-estimation

We use the
[`blinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/blinded_ssr.md)
function to estimate the nuisance parameters (dispersion and overall
rate) from the blinded data and recalculate the sample size.

The method assumes the treatment effect (rate ratio) is maintained as
planned, but updates the sample size based on the observed variance
structure.

``` r
# Perform blinded SSR
blinded_result <- blinded_ssr(
  data = interim_data,
  ratio = 1,
  lambda1_planning = lambda1_plan,
  lambda2_planning = lambda2_plan,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan
)

print(blinded_result)
#> $n_total_blinded
#> (Intercept) 
#>         440 
#> 
#> $dispersion_blinded
#> [1] 0.8021666
#> 
#> $lambda_blinded
#> (Intercept) 
#>     1.73056 
#> 
#> $lambda1_adjusted
#> (Intercept) 
#>    1.977783 
#> 
#> $lambda2_adjusted
#> (Intercept) 
#>    1.483337 
#> 
#> $info_fraction
#> [1] 0.3482652
#> 
#> $blinded_info
#> [1] 44.21613
#> 
#> $target_info
#> [1] 126.9611
```

The blinded estimate of dispersion is 0.802, which is higher than the
planned 0.5. Consequently, the re-estimated sample size 440 is larger
than the planned 280.

## Unblinded sample size re-estimation

Alternatively, if the interim analysis is performed by an Independent
Data Monitoring Committee (IDMC) with access to unblinded data, we can
estimate the rates and dispersion separately for each group.

``` r
# Fit Negative Binomial model to unblinded interim data
# We use glm.nb from MASS package
fit <- glm.nb(events ~ treatment + offset(log(tte)), data = interim_data)
summary(fit)
#> 
#> Call:
#> glm.nb(formula = events ~ treatment + offset(log(tte)), data = interim_data, 
#>     init.theta = 1.269845637, link = log)
#> 
#> Coefficients:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)             0.6875     0.1042   6.596 4.21e-11 ***
#> treatmentExperimental  -0.2938     0.1492  -1.969    0.049 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(1.2698) family taken to be 1)
#> 
#>     Null deviance: 200.87  on 180  degrees of freedom
#> Residual deviance: 197.03  on 179  degrees of freedom
#> AIC: 1040.3
#> 
#> Number of Fisher Scoring iterations: 1
#> 
#> 
#>               Theta:  1.270 
#>           Std. Err.:  0.167 
#> 
#>  2 x log-likelihood:  -1034.284

# Extract estimates
# Intercept is log(lambda_control)
# treatmentExperimental is log(RR)
# theta is 1/k

est_lambda1 <- exp(coef(fit)["(Intercept)"])
est_lambda2 <- exp(coef(fit)["(Intercept)"] + coef(fit)["treatmentExperimental"])
est_k <- 1 / fit$theta

cat("Unblinded Estimates:\n")
#> Unblinded Estimates:
cat("Lambda 1:", est_lambda1, "\n")
#> Lambda 1: 1.98877
cat("Lambda 2:", est_lambda2, "\n")
#> Lambda 2: 1.482545
cat("Dispersion (k):", est_k, "\n")
#> Dispersion (k): 0.7874973

# Recalculate sample size with unblinded estimates
# We keep the original power and alpha
unblinded_design <- sample_size_nbinom(
  lambda1 = est_lambda1,
  lambda2 = est_lambda2,
  dispersion = est_k,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan
)

print(unblinded_design)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 207, n2 = 207, total = 414
#> Expected events: 12934.1 (n1: 7410.2, n2: 5524.0)
#> Power: 90%, Alpha: 0.025 (1-sided)
#> Rates: control = 1.9888, treatment = 1.4825 (RR = 0.7455)
#> Dispersion: 0.7875, Avg exposure (calendar): 18.00
#> Accrual: 12.0, Trial duration: 24.0
```

The unblinded re-estimation uses the observed rates and dispersion. In
this specific simulation, the observed dispersion 0.787 is close to the
true value 0.8. The rates might also fluctuate. The resulting sample
size 414 reflects the actual variability and effect size observed so
far.

## Comparison and conclusion

``` r
comparison <- data.frame(
  Method = c("Planned", "Blinded SSR", "Unblinded SSR"),
  Dispersion = c(k_plan, blinded_result$dispersion_blinded, est_k),
  Sample_Size = c(design_plan$n_total, blinded_result$n_total_blinded, unblinded_design$n_total)
)

knitr::kable(comparison, digits = 3, caption = "Comparison of Sample Size Estimates")
```

| Method        | Dispersion | Sample_Size |
|:--------------|-----------:|------------:|
| Planned       |      0.500 |         280 |
| Blinded SSR   |      0.802 |         440 |
| Unblinded SSR |      0.787 |         414 |

Comparison of Sample Size Estimates

In this example, the true dispersion was higher than planned. Both
blinded and unblinded SSR detected the increased variability and
suggested increasing the sample size to maintain power. The blinded
method is often preferred to maintain trial integrity, assuming the
treatment effect is not far from the planning assumption. The unblinded
method adapts to both the variance and the observed effect size, which
can lead to larger or smaller sample sizes depending on the observed
efficacy.

**Note:** In practice, sample size is typically not reduced if the
re-estimated size is smaller than planned. We would take
`max(n_planned, n_reestimated)`.

## References

Friede, T, and H Schmidli. 2010. “Blinded Sample Size Reestimation with
Negative Binomial Counts in Superiority and Non-Inferiority Trials.”
*Methods of Information in Medicine* 49 (06): 618–24.
