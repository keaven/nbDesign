# GitHub Copilot Instructions for gsDesignNB

## Project Overview
`gsDesignNB` is an R package for designing and simulating clinical trials with negative binomial outcomes (recurrent events). It integrates with the `gsDesign` package to support group sequential designs.

## Architecture & Key Components
- **Core Logic (`R/`)**:
  - `sample_size_nbinom.R`: Fixed design sample size calculations.
  - `sim_gs_nbinom.R`: Group sequential design simulations.
  - `nb_sim.R`: Underlying data generation for simulations.
- **Data Handling**: Uses `data.table` for high-performance data manipulation in simulations.
- **Integration**: Re-exports `gsDesign` functions (`reexport-gsDesign.R`) to extend its functionality seamlessly.

## Development Workflow
- **Package Management**: Standard R package structure.
- **Build & Reload**: Use `devtools::load_all()` for interactive development.
- **Testing**: Run tests with `devtools::test()` or `testthat::test_dir("tests/testthat")`.
- **Documentation**: Generate documentation with `devtools::document()`. The site is built with `pkgdown`.

## Coding Standards & Patterns
- **Style**: Follows the [tidyverse style guide](https://style.tidyverse.org/).
- **Data Manipulation**:
  - Prefer `data.table` syntax (`dt[i, j, by]`) inside internal functions for performance, especially in simulation loops.
  - Use `utils::globalVariables()` in `R/globals.R` to handle `data.table` non-standard evaluation warnings.
- **Documentation**: Use `roxygen2` with Markdown support.
  - Include `@examples` that are runnable.
  - Use `@importFrom` to manage dependencies explicitly.
- **Error Handling**: Validate inputs early in exported functions (e.g., check for positive rates, valid probabilities).

## Testing Guidelines
- **Framework**: Uses `testthat` (Edition 3).
- **Structure**: Tests are located in `tests/testthat/` and prefixed with `test-`.
- **Best Practices**:
  - Test against known theoretical results (e.g., Poisson limit when dispersion is 0).
  - Verify S3 class structures and return types.
  - Use `expect_error` to verify input validation.
  - Example: See `tests/testthat/test-sample_size_nbinom.R` for testing calculation logic.

## Common Tasks
- **Adding a new simulation feature**:
  1. Update `nb_sim.R` or `nb_sim_seasonal.R` to generate the new data structure.
  2. Update `sim_gs_nbinom.R` to handle the new data in the simulation loop.
  3. Add a test case in `tests/testthat/`.
- **Modifying sample size logic**:
  1. Edit `sample_size_nbinom.R`.
  2. Verify against published literature examples if possible (referenced in documentation).
