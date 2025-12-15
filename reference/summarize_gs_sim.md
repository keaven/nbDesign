# Summarize Group Sequential Simulation Results

Provides a summary of the operating characteristics of the group
sequential design based on simulation results.

## Usage

``` r
summarize_gs_sim(x)
```

## Arguments

- x:

  A data frame returned by
  [`check_gs_bound()`](https://keaven.github.io/gsDesignNB/reference/check_gs_bound.md)
  (or
  [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)
  if bounds are manually checked). Must contain columns `cross_upper`,
  `cross_lower`.

## Value

A list containing:

- n_sim:

  Number of simulations

- power:

  Overall power (probability of crossing upper bound)

- futility:

  Overall futility rate (probability of crossing lower bound and not
  upper)

- analysis_summary:

  Data frame with per-analysis statistics (events, crossings)
