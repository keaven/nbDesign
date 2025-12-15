# Determine Analysis Date based on Criteria

Finds the earliest calendar date at which all specified criteria are
met. Criteria can include a specific calendar date, a target number of
events, a target number of completers, or a target amount of blinded
information.

## Usage

``` r
get_cut_date(
  data,
  planned_calendar = NULL,
  target_events = NULL,
  target_completers = NULL,
  target_info = NULL,
  event_gap = 0,
  ratio = 1,
  lambda1 = NULL,
  lambda2 = NULL,
  min_date = 0,
  max_date = Inf
)
```

## Arguments

- data:

  A data frame of simulated data (from
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md)).

- planned_calendar:

  Numeric. Target calendar time.

- target_events:

  Integer. Target number of observed events.

- target_completers:

  Integer. Target number of subjects with complete follow-up.

- target_info:

  Numeric. Target blinded information.

- event_gap:

  Numeric. Gap duration for event counting and info calculation.

- ratio:

  Numeric. Randomization ratio (experimental/control) for info
  calculation.

- lambda1:

  Numeric. Planned control rate for info calculation.

- lambda2:

  Numeric. Planned experimental rate for info calculation.

- min_date:

  Numeric. Minimum possible date (e.g., 0 or previous analysis time).

- max_date:

  Numeric. Maximum possible date (e.g., trial duration).

## Value

Numeric. The calendar date satisfying the criteria. If criteria cannot
be met within `max_date` (or data limits), returns `max_date` (or max
data time).
