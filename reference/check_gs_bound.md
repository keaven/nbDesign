# Check Group Sequential Bounds

Updates the group sequential design boundaries based on observed
information and checks if boundaries have been crossed.

## Usage

``` r
check_gs_bound(sim_results, design, info_scale = c("blinded", "unblinded"))
```

## Arguments

- sim_results:

  Data frame of simulation results (from
  [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)).

- design:

  The planning `gsNB` object.

- info_scale:

  Character. "blinded" (default) or "unblinded" information to use for
  bounds.

## Value

A data frame with added columns:

- cross_upper:

  Logical, true if upper bound crossed (efficacy)

- cross_lower:

  Logical, true if lower bound crossed (futility)
