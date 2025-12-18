# Sample size calculation for negative binomial outcomes

``` r
library(gsDesignNB)
```

This vignette describes the methodology used in the `sample_size_nbinom`
function for calculating sample sizes when comparing two treatment
groups with negative binomial outcomes. It covers two available methods
and how to account for variable accrual rates.

## Methodology

We assume the outcome $Y$ follows a negative binomial distribution with
mean $\mu$ and a common dispersion parameter $k$ for both treatment
groups, such that the variance is given by:

$$Var(Y) = \mu + k\mu^{2}$$

Note that in R’s `rnbinom` parameterization, $k = 1/\texttt{𝚜𝚒𝚣𝚎}$.

We wish to test the null hypothesis $H_{0}:\lambda_{1} = \lambda_{2}$
against the alternative $H_{1}:\lambda_{1} \neq \lambda_{2}$, where
$\lambda_{1}$ and $\lambda_{2}$ are the event rates in the control and
treatment groups, respectively.

### Connection between rates and counts

The negative binomial distribution can be motivated by a Gamma-Poisson
mixture model. Suppose that for each subject $i$, the event rate
$\Lambda_{i}$ is a random variable following a Gamma distribution with
shape $\alpha = 1/k$ and rate $\beta = 1/(k\lambda)$, where $\lambda$ is
the underlying population event rate and $k$ is the dispersion
parameter. The mean of this Gamma distribution is
$E\left\lbrack \Lambda_{i} \right\rbrack = \alpha/\beta = \lambda$ and
the variance is
$Var\left( \Lambda_{i} \right) = \alpha/\beta^{2} = k\lambda^{2}$.

Given the subject-specific rate $\Lambda_{i}$, the number of events
$Y_{i}$ observed over a time period $t_{i}$ follows a Poisson
distribution with rate $\Lambda_{i}t_{i}$:
$$Y_{i}|\Lambda_{i} \sim \text{Poisson}\left( \Lambda_{i}t_{i} \right)$$

The marginal distribution of $Y_{i}$ (integrating out $\Lambda_{i}$) is
then a negative binomial distribution with mean
$\mu_{i} = \lambda t_{i}$ and dispersion parameter $k$. The variance is:
$$Var\left( Y_{i} \right) = E\left\lbrack Var\left( Y_{i}|\Lambda_{i} \right) \right\rbrack + Var\left( E\left\lbrack Y_{i}|\Lambda_{i} \right\rbrack \right)$$$$= E\left\lbrack \Lambda_{i}t_{i} \right\rbrack + Var\left( \Lambda_{i}t_{i} \right)$$$$= \lambda t_{i} + t_{i}^{2}\left( k\lambda^{2} \right)$$$$= \mu_{i} + k\mu_{i}^{2}$$

This formulation connects the event rate $\lambda$ (used in the
hypothesis testing framework) with the expected count $\mu$ (used in the
negative binomial parameterization), showing how heterogeneity in
individual rates leads to the overdispersion characteristic of the
negative binomial distribution.

### Sample size formula

The sample size calculation is based on the asymptotic normality of the
log rate ratio. This approach corresponds to **Method 3** of Zhu and
Lakkis (2014), which uses a Wald statistic for the log rate ratio. It is
also the method described by Friede and Schmidli (2010) and Mütze et al.
(2019) (as implemented in the `gscounts` package).

The total sample size $n_{total}$ is calculated as:

$$n_{total} = \frac{\left( z_{\alpha/s} + z_{\beta} \right)^{2} \cdot \widetilde{V}}{\left( \log\left( \lambda_{1}/\lambda_{2} \right) \right)^{2}}$$

where:

- $z_{\alpha/s}$ and $z_{\beta}$ are the standard normal critical values
  for the significance level $\alpha$ (with $s = 1$ or $2$ sided) and
  power $1 - \beta$.
- $\lambda_{1}$ and $\lambda_{2}$ are the event rates in the control and
  treatment groups.
- $\widetilde{V}$ is the average variance per subject, defined as:

$$\widetilde{V} = \frac{1/\mu_{1} + k}{p_{1}} + \frac{1/\mu_{2} + k}{p_{2}}$$

where:

- $p_{1} = n_{1}/n_{total}$ and $p_{2} = n_{2}/n_{total}$ are the
  allocation proportions.
- $\mu_{i} = \lambda_{i} \cdot \bar{t}$ is the expected mean count for
  group $i$ over the average exposure duration $\bar{t}$.
- $k$ is the dispersion parameter. While the cited literature assumes a
  common dispersion parameter across groups, `sample_size_nbinom` allows
  for different dispersion parameters $k_{1}$ and $k_{2}$ for the
  control and treatment groups, respectively. In this case, the variance
  formula becomes:

$$\widetilde{V} = \frac{1/\mu_{1} + k_{1}}{p_{1}} + \frac{1/\mu_{2} + k_{2}}{p_{2}}$$

This formula assumes that the exposure duration is the same for both
groups (or uses an average exposure $\bar{t}$). However,
`sample_size_nbinom` allows for different dropout rates for the two
groups, which results in different average exposure durations
${\bar{t}}_{1}$ and ${\bar{t}}_{2}$. In this case,
$\mu_{1} = \lambda_{1}{\bar{t}}_{1}$ and
$\mu_{2} = \lambda_{2}{\bar{t}}_{2}$. Also, when the assumed gap between
events is \> 0, the expected exposure is impacted differently for each
group based on their event rates.

## Average exposure with variable accrual and dropout

When the accrual rate is not constant or when the trial has a fixed
duration with ongoing recruitment, the exposure time for patients will
vary. The function calculates an *average exposure* time to use in the
sample size formula. Additionally, if a `dropout_rate` is specified, the
exposure is adjusted to account for patients leaving the study early.

Let $T$ be the total trial duration. Suppose recruitment occurs in $J$
segments, where the $j$-th segment has accrual rate $R_{j}$ and duration
$D_{j}$.

If `dropout_rate` is 0:

1.  The expected number of patients recruited in segment $j$ is
    $N_{j} = R_{j} \cdot D_{j}$.
2.  The start time of segment $j$ is
    $S_{j - 1} = \sum_{i = 1}^{j - 1}D_{i}$ (with $S_{0} = 0$).
3.  The midpoint of recruitment for segment $j$ is
    $M_{j} = S_{j - 1} + D_{j}/2$.
4.  The average follow-up (exposure) time for patients recruited in
    segment $j$ is approximately $E_{j} = T - M_{j}$.

If `dropout_rate` ($\delta$) \> 0, the average exposure is calculated by
integrating the exposure function over the recruitment interval:

$$\begin{aligned}
E_{j} & {= \frac{1}{D_{j}}\int_{u_{min}}^{u_{max}}\frac{1 - e^{- \delta u}}{\delta}du} \\
 & {= \frac{1}{\delta} - \frac{1}{\delta^{2}D_{j}}\left( e^{- \delta u_{min}} - e^{- \delta u_{max}} \right)}
\end{aligned}$$ where $u_{max}$ and $u_{min}$ are the maximum and
minimum potential follow-up times for patients in that segment
($T - S_{j - 1}$ and $T - \left( S_{j - 1} + D_{j} \right)$
respectively).

### Group-specific parameters

The function supports specifying `dropout_rate` as a vector of length 2,
corresponding to the control and treatment groups respectively. This
allows for scenarios where the dropout rate differs between arms.

If these parameters differ, the average exposure is calculated
separately for each group (${\bar{t}}_{1}$ and ${\bar{t}}_{2}$). The
variance inflation factor $Q$ (see below) is also calculated separately
for each group ($Q_{1}$ and $Q_{2}$).

### Maximum follow-up

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

### Variance inflation for variable follow-up

When follow-up times are variable (due to accrual, dropout, or
administrative censoring), simply using the average follow-up time
$\bar{t}$ in the variance formula underestimates the true variance of
the rate estimator. This is because the variance of the negative
binomial distribution depends on the exposure time in a non-linear way
($Var(Y) = \mu + k\mu^{2} = \lambda t + k(\lambda t)^{2}$).

To account for this, we apply a variance inflation factor $Q$ to the
dispersion parameter $k$, as derived by Zhu and Lakkis (2014):

$$Q = \frac{E\left\lbrack t^{2} \right\rbrack}{\left( E\lbrack t\rbrack \right)^{2}}$$

The adjusted dispersion parameter used in the sample size calculation is
$k_{adj} = k \cdot Q$. The function automatically calculates
$E\lbrack t\rbrack$ and $E\left\lbrack t^{2} \right\rbrack$ based on the
accrual, dropout, and trial duration parameters. If exposure differs
between groups, $Q$ is calculated separately for each group.

### Event gaps

In some clinical trials, there is a mandatory “dead time” or gap after
an event during which no new events can occur (e.g., a recovery period).
If an `event_gap` is specified, the effective exposure time for a
subject is reduced by the time spent in these gaps.

The function approximates the effective event rate as:
$$\lambda_{eff} \approx \frac{\lambda}{1 + \lambda \cdot \text{gap}}$$
This adjusted rate is then used in the sample size calculations. The
effective exposure time reported is also adjusted similarly.

## Examples

### Basic calculation (Zhu and Lakkis 2014)

Calculate sample size for:

- Control rate $\lambda_{1} = 0.5$
- Treatment rate $\lambda_{2} = 0.3$
- Dispersion $k = 0.1$
- Power = 80%
- Alpha = 0.025 (one-sided)
- Accrual over 12 months
- Trial duration 12 months (implying exposure approx 6 months)

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
  trial_duration = 12
)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 35, n2 = 35, total = 70
#> Expected events: 168.0 (n1: 105.0, n2: 63.0)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 6.00
#> Accrual: 12.0, Trial duration: 12.0
```

### Piecewise constant accrual

Consider a trial where recruitment ramps up:

- 5 patients/month for the first 3 months
- 10 patients/month for the next 3 months
- Total trial duration is 12 months

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
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 26, n2 = 26, total = 52
#> Expected events: 176.8 (n1: 110.5, n2: 66.3)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 8.50
#> Accrual: 6.0, Trial duration: 12.0
```

### Accrual with dropout and max follow-up

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
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 38, n2 = 38, total = 76
#> Expected events: 157.6 (n1: 98.5, n2: 59.1)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 5.18
#> Dropout rate: 0.0500
#> Accrual: 6.0, Trial duration: 12.0
#> Max follow-up: 6.0
```

### Group-specific dropout rates

Suppose the control group has a higher dropout rate (10%) than the
treatment group (5%).

``` r
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = c(5, 10),
  accrual_duration = c(3, 3),
  trial_duration = 12,
  dropout_rate = c(0.10, 0.05), # Control, Treatment
  max_followup = 6
)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 40, n2 = 40, total = 80
#> Expected events: 152.4 (n1: 90.2, n2: 62.2)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 4.51 (n1), 5.18 (n2)
#> Dropout rate: 0.1000 (n1), 0.0500 (n2)
#> Accrual: 6.0, Trial duration: 12.0
#> Max follow-up: 6.0
```

### Calculating power for fixed design

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
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 38, n2 = 38, total = 76
#> Expected events: 177.3 (n1: 98.5, n2: 78.8)
#> Power: 26%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.4000 (RR = 0.8000)
#> Dispersion: 0.1000, Avg exposure (calendar): 5.18
#> Dropout rate: 0.0500
#> Accrual: 6.0, Trial duration: 12.0
#> Max follow-up: 6.0
```

### Unequal allocation

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
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 40, n2 = 80, total = 120
#> Expected events: 264.0 (n1: 120.0, n2: 144.0)
#> Power: 95%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 6.00
#> Accrual: 12.0, Trial duration: 12.0
```

## Accounting for event gaps

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
at-risk exposure differs between treatment groups if their rates differ,
even if the calendar exposure is the same. This is a key reason why
average exposure may differ between groups in the output.

### Example with event gap

Calculate sample size assuming a 30-day gap after each event (approx
0.082 years). Note how the `exposure_at_risk` differs between groups
because the group with the higher event rate ($\lambda_{1} = 2.0$)
spends more time in the “gap” period than the group with the lower rate
($\lambda_{2} = 1.0$).

``` r
sample_size_nbinom(
  lambda1 = 2.0,
  lambda2 = 1.0,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = 10,
  accrual_duration = 12,
  trial_duration = 12,
  event_gap = 30 / 365.25
)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 9, n2 = 9, total = 18
#> Expected events: 142.7 (n1: 92.8, n2: 49.9)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 2.0000, treatment = 1.0000 (RR = 0.5000)
#> Dispersion: 0.1000, Avg exposure (calendar): 6.00
#> Avg exposure (at-risk): n1 = 5.15, n2 = 5.54
#> Event gap: 0.08
#> Accrual: 12.0, Trial duration: 12.0
```

In this example, even though the calendar time is the same for both
groups, the effective at-risk exposure is lower for Group 1 because they
have more events and thus more “gap” time.

## References

Friede, T, and H Schmidli. 2010. “Blinded Sample Size Reestimation with
Negative Binomial Counts in Superiority and Non-Inferiority Trials.”
*Methods of Information in Medicine* 49 (06): 618–24.

Mütze, Tobias, Ekkehard Glimm, Heinz Schmidli, and Tim Friede. 2019.
“Group Sequential Designs for Negative Binomial Outcomes.” *Statistical
Methods in Medical Research* 28 (8): 2326–47.

Zhu, Haiyuan, and Hassan Lakkis. 2014. “Sample Size Calculation for
Comparing Two Negative Binomial Rates.” *Statistics in Medicine* 33 (3):
376–87.
