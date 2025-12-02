# Simulate Recurrent Events with Fixed Follow-up

Simulates recurrent events for a clinical trial with piecewise constant
enrollment, exponential failure rates (Poisson process), and piecewise
exponential dropout.

## Usage

``` r
nb_sim(
  enroll_rate,
  fail_rate,
  dropout_rate = NULL,
  max_followup = NULL,
  n = NULL,
  block = c(rep("Control", 2), rep("Experimental", 2))
)
```

## Arguments

- enroll_rate:

  A data frame with columns `rate` and `duration` defining the piecewise
  constant enrollment rates.

- fail_rate:

  A data frame with columns `treatment` and `rate` defining the
  exponential failure rate for each treatment group. Optionally, a
  `dispersion` column can be provided to generate data from a Negative
  Binomial distribution. The dispersion parameter `k` is such that
  \\Var(Y) = \mu + k \mu^2\\.

- dropout_rate:

  A data frame with columns `treatment`, `rate`, and `duration` defining
  the piecewise constant dropout rates.

- max_followup:

  Numeric. Maximum duration of follow-up for each individual (relative
  to their randomization time).

- n:

  Total sample size. If NULL, it is estimated from `enroll_rate`. If
  provided, enrollment stops when `n` subjects are recruited.

- block:

  Block vector for treatment allocation. Default is
  `c(rep("Control", 2), rep("Experimental", 2))`. If NULL, simple
  randomization is used (treatments are assigned with equal
  probability). If provided, it specifies the block structure, e.g.,
  `c(rep("A", 2), rep("B", 2))` assigns 2 to group A and 2 to group B in
  each block.

## Value

A data frame (tibble) with columns:

- id:

  Subject identifier

- treatment:

  Treatment group

- enroll_time:

  Time of enrollment relative to trial start

- tte:

  Time to event or censoring relative to randomization

- calendar_time:

  Calendar time of event or censoring (enroll_time + tte)

- event:

  Binary indicator: 1 for event, 0 for censoring

Multiple rows per subject are returned (one for each event, plus one for
the final censoring time).
