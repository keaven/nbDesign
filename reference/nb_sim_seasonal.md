# Simulate recurrent events with seasonal rates

Simulates recurrent events where event rates depend on the season.

## Usage

``` r
nb_sim_seasonal(
  enroll_rate,
  fail_rate,
  dropout_rate = NULL,
  max_followup = NULL,
  randomization_start_date = NULL,
  n = NULL,
  block = c(rep("Control", 2), rep("Experimental", 2))
)
```

## Arguments

- enroll_rate:

  A data frame with columns `rate` and `duration`.

- fail_rate:

  A data frame with columns `treatment`, `season`, `rate`, and
  optionally `dispersion`. Seasons should be "Spring", "Summer", "Fall",
  "Winter".

- dropout_rate:

  A data frame with columns `treatment`, `rate`, `duration`.

- max_followup:

  Numeric. Max follow-up duration (years).

- randomization_start_date:

  Date. Start of randomization.

- n:

  Integer. Total sample size.

- block:

  Character vector for block randomization.

## Value

A data frame of class `nb_sim_seasonal` with columns: `id`, `treatment`,
`season`, `enroll_time`, `start`, `end`, `event`, `calendar_start`,
`calendar_end`. Rows represent intervals of risk or events. `event=1`
indicates an event at `end`. `event=0` indicates censoring or end of a
seasonal interval at `end`.

## Examples

``` r
enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
fail_rate <- data.frame(
  treatment = rep(c("Control", "Experimental"), each = 4),
  season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
  rate = c(0.6, 0.5, 0.4, 0.5, 0.4, 0.3, 0.2, 0.3)
)
sim <- nb_sim_seasonal(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  max_followup = 1,
  randomization_start_date = as.Date("2020-01-01"),
  n = 20
)
head(sim)
#>   id    treatment enroll_time season     start       end event  cal_start
#> 1  1      Control  0.01896022 Winter 0.0000000 0.1453108     0 2020-01-07
#> 2  1      Control  0.01896022 Spring 0.1453108 0.3971931     0 2020-03-01
#> 3  1      Control  0.01896022 Summer 0.3971931 0.6490754     0 2020-06-01
#> 4  1      Control  0.01896022   Fall 0.6490754 0.8982198     0 2020-09-01
#> 5  1      Control  0.01896022 Winter 0.8982198 1.0000000     0 2020-12-01
#> 6  2 Experimental  0.03829783 Winter 0.0000000 0.1259732     0 2020-01-14
#>      cal_end
#> 1 2020-03-01
#> 2 2020-06-01
#> 3 2020-09-01
#> 4 2020-12-01
#> 5 2021-01-07
#> 6 2020-03-01
```
