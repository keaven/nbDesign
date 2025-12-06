# Sample Size Calculation for Negative Binomial Outcomes

``` r
library(gsDesignNB)
```

This vignette describes the methodology used in the `sample_size_nbinom`
function for calculating sample sizes when comparing two treatment
groups with negative binomial outcomes. It covers two available methods
and how to account for variable accrual rates.

## Methodology

We assume the outcome $Y$ follows a negative binomial distribution with
mean $\mu$ and dispersion parameter $k$, such that the variance is given
by:

$$Var(Y) = \mu + k\mu^{2}$$

Note that in R’s `rnbinom` parameterization, $k = 1/\texttt{𝚜𝚒𝚣𝚎}$.

We wish to test the null hypothesis $H_{0}:\lambda_{1} = \lambda_{2}$
against the alternative $H_{1}:\lambda_{1} \neq \lambda_{2}$, where
$\lambda_{1}$ and $\lambda_{2}$ are the event rates in the control and
treatment groups, respectively.

The function implements two methods:

### Method 1: Zhu and Lakkis (2014)

This method is based on the asymptotic normality of the log rate ratio.
The sample size for the control group ($n_{1}$) is:

$$n_{1} = \frac{\left( z_{\alpha/s} + z_{\beta} \right)^{2} \cdot V}{\left( \log\left( \mu_{1}/\mu_{2} \right) \right)^{2}}$$

where $V$ is the variance component:
$$V = \left( \frac{1}{\mu_{1}} + k \right) + \frac{1}{r}\left( \frac{1}{\mu_{2}} + k \right)$$

where $r = n_{2}/n_{1}$ is the allocation ratio, and
$\mu_{i} = \lambda_{i} \cdot t$ is the expected mean count over exposure
duration $t$.

### Method 2: Friede and Schmidli (2010) / Mütze et al. (2018)

This method uses a Wald test statistic and is commonly used in group
sequential designs (as implemented in the `gscounts` package). The total
sample size $n_{total}$ is calculated as:

$$n_{total} = \frac{\left( z_{\alpha/s} + z_{\beta} \right)^{2} \cdot \bar{V}}{\left( \log\left( \lambda_{1}/\lambda_{2} \right) \right)^{2}}$$

where $\bar{V}$ is the average variance per subject:
$$\bar{V} = \frac{1/\mu_{1} + k}{p_{1}} + \frac{1/\mu_{2} + k}{p_{2}}$$

where $p_{1} = n_{1}/n_{total}$ and $p_{2} = n_{2}/n_{total}$ are the
allocation proportions.

**Note:** For a fixed design with equal allocation, both methods yield
identical sample sizes.

## Average Exposure with Variable Accrual and Dropout

When the accrual rate is not constant or when the trial has a fixed
duration with ongoing recruitment, the exposure time for patients will
vary. The function calculates an *average exposure* time to use in the
sample size formula. Additionally, if a `dropout_rate` is specified, the
exposure is adjusted to account for patients leaving the study early.

Let $T$ be the total trial duration. Suppose recruitment occurs in $J$
segments, where the $j$-th segment has accrual rate $R_{j}$ and duration
$D_{j}$.

If `dropout_rate` is 0:

1.  The number of patients recruited in segment $j$ is
    $N_{j} = R_{j} \cdot D_{j}$.
2.  The start time of segment $j$ is
    $S_{j - 1} = \sum_{i = 1}^{j - 1}D_{i}$ (with $S_{0} = 0$).
3.  The midpoint of recruitment for segment $j$ is
    $M_{j} = S_{j - 1} + D_{j}/2$.
4.  The average follow-up (exposure) time for patients recruited in
    segment $j$ is approximately $E_{j} = T - M_{j}$.

If `dropout_rate` ($\delta$) \> 0, the average exposure is calculated by
integrating the exposure function over the recruitment interval:

$$E_{j} = \frac{1}{\delta} - \frac{1}{\delta^{2}D_{j}}\left( e^{- \delta u_{min}} - e^{- \delta u_{max}} \right)$$
where $u_{max}$ and $u_{min}$ are the maximum and minimum potential
follow-up times for patients in that segment ($T - S_{j - 1}$ and
$T - \left( S_{j - 1} + D_{j} \right)$ respectively).

### Maximum Follow-up

If `max_followup` ($F$) is specified, the follow-up time for any
individual is capped at $F$. This creates three scenarios for a
recruitment segment:

1.  **All truncated:** If $u_{min} \geq F$, all patients in the segment
    have potential follow-up $\geq F$, so their actual follow-up is $F$
    (subject to dropout).
2.  **None truncated:** If $u_{max} \leq F$, no patients reach the cap
    $F$ before the trial ends. The calculation is as above.
3.  **Partial truncation:** If $u_{min} < F < u_{max}$, patients
    recruited earlier in the segment are capped at $F$, while those
    recruited later are followed until the trial end. The segment is
    split into two parts for calculation.

The overall average exposure used for the calculation is the weighted
average:

$$\bar{t} = \frac{\sum\limits_{j = 1}^{J}N_{j}E_{j}}{\sum\limits_{j = 1}^{J}N_{j}}$$

## Examples

### Basic Calculation (Zhu and Lakkis)

Calculate sample size for: \* Control rate $\lambda_{1} = 0.5$ \*
Treatment rate $\lambda_{2} = 0.3$ \* Dispersion $k = 0.1$ \* Power =
80% \* Alpha = 0.025 (one-sided) \* Accrual over 12 months \* Trial
duration 12 months (implying exposure approx 6 months)

``` r
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  alpha = 0.025,
  sided = 1,
  accrual_rate = 10, # arbitrary, just for average exposure
  accrual_duration = 12,
  trial_duration = 12,
  method = "zhu"
)
#> Sample Size for Negative Binomial Outcome
#> ==========================================
#> 
#> Method:          zhu
#> Sample size:     n1 = 33, n2 = 33, total = 66
#> Expected events: 158.4 (n1: 99.0, n2: 59.4)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 6.00
#> Accrual: 12.0, Trial duration: 12.0
```

### Using Friede et al. Method

``` r
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  alpha = 0.025,
  sided = 1,
  accrual_rate = 10,
  accrual_duration = 12,
  trial_duration = 12,
  method = "friede"
)
#> Sample Size for Negative Binomial Outcome
#> ==========================================
#> 
#> Method:          friede
#> Sample size:     n1 = 33, n2 = 33, total = 66
#> Expected events: 158.4 (n1: 99.0, n2: 59.4)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 6.00
#> Accrual: 12.0, Trial duration: 12.0
```

### Piecewise Constant Accrual

Consider a trial where recruitment ramps up: \* 5 patients/month for the
first 3 months \* 10 patients/month for the next 3 months \* Total trial
duration is 12 months

The function automatically calculates the average exposure based on this
accrual pattern.

``` r
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = c(5, 10),
  accrual_duration = c(3, 3),
  trial_duration = 12
)
#> Sample Size for Negative Binomial Outcome
#> ==========================================
#> 
#> Method:          zhu
#> Sample size:     n1 = 25, n2 = 25, total = 50
#> Expected events: 170.0 (n1: 106.2, n2: 63.7)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 8.50
#> Accrual: 6.0, Trial duration: 12.0
```

### Accrual with Dropout and Max Follow-up

Same design as above, but with a 5% dropout rate per unit time and a
maximum follow-up of 6 months.

``` r
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = c(5, 10),
  accrual_duration = c(3, 3),
  trial_duration = 12,
  dropout_rate = 0.05,
  max_followup = 6
)
#> Sample Size for Negative Binomial Outcome
#> ==========================================
#> 
#> Method:          zhu
#> Sample size:     n1 = 37, n2 = 37, total = 74
#> Expected events: 153.4 (n1: 95.9, n2: 57.5)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 5.18
#> Accrual: 6.0, Trial duration: 12.0
```

### Calculating Power for Fixed Design

Using the accrual rates and design from the previous example, suppose we
want to calculate the power if the treatment effect is smaller
($\lambda_{2} = 0.4$ instead of $0.3$). We use the `accrual_rate`
computed in the previous step.

``` r
# Store the result from the previous calculation
design_result <- sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = c(5, 10),
  accrual_duration = c(3, 3),
  trial_duration = 12,
  dropout_rate = 0.05,
  max_followup = 6
)

# Use the computed accrual rates to calculate power for a smaller effect size
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.4, # Smaller effect size
  dispersion = 0.1,
  power = NULL, # Request power calculation
  accrual_rate = design_result$accrual_rate, # Use computed rates
  accrual_duration = c(3, 3),
  trial_duration = 12,
  dropout_rate = 0.05,
  max_followup = 6
)
#> Sample Size for Negative Binomial Outcome
#> ==========================================
#> 
#> Method:          zhu
#> Sample size:     n1 = 37, n2 = 37, total = 74
#> Expected events: 172.6 (n1: 95.9, n2: 76.7)
#> Power: 26%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.4000 (RR = 0.8000)
#> Dispersion: 0.1000, Avg exposure (calendar): 5.18
#> Accrual: 6.0, Trial duration: 12.0
```

### Unequal Allocation

Sample size with a 2:1 allocation ratio ($n_{2} = 2n_{1}$).

``` r
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  ratio = 2,
  accrual_rate = 10,
  accrual_duration = 12,
  trial_duration = 12
)
#> Sample Size for Negative Binomial Outcome
#> ==========================================
#> 
#> Method:          zhu
#> Sample size:     n1 = 40, n2 = 80, total = 120
#> Expected events: 264.0 (n1: 120.0, n2: 144.0)
#> Power: 96%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 6.00
#> Accrual: 12.0, Trial duration: 12.0
```

## Accounting for Event Gaps

In some recurrent event trials, there may be a mandatory “gap” period
after each event during which no new events can be recorded (e.g., a
recovery period or administrative window). This effectively reduces the
time at risk.

If an `event_gap` ($g$) is specified, the function adjusts the
calculation as follows:

1.  **Effective Rates**: The event rates are adjusted to
    $\lambda^{*} = \lambda/(1 + \lambda g)$ for the sample size
    calculation (Zhu and Lakkis method).
2.  **At-Risk Exposure**: The function reports the “average at-risk
    exposure” $E_{risk} = E_{cal}/(1 + \lambda g)$ alongside the
    standard calendar exposure. This provides transparency on the actual
    time subjects are at risk for events.

Since the gap reduction depends on the event rate ($\lambda$), the
at-risk exposure differs between treatment groups if their rates differ.

### Example with Event Gap

Calculate sample size assuming a 5-day gap after each event (approx
0.0137 years).

``` r
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = 10,
  accrual_duration = 12,
  trial_duration = 12,
  event_gap = 5 / 365.25
)
#> Sample Size for Negative Binomial Outcome
#> ==========================================
#> 
#> Method:          zhu
#> Sample size:     n1 = 33, n2 = 33, total = 66
#> Expected events: 157.5 (n1: 98.3, n2: 59.2)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 6.00
#> Avg exposure (at-risk): n1 = 5.96, n2 = 5.98
#> Accrual: 12.0, Trial duration: 12.0
```

The output shows both the “Avg exposure (calendar)” and the “Avg
exposure (at-risk)” for each group.

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
