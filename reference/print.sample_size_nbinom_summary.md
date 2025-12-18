# Print method for sample_size_nbinom_summary objects

Print method for sample_size_nbinom_summary objects

## Usage

``` r
# S3 method for class 'sample_size_nbinom_summary'
print(x, ...)
```

## Arguments

- x:

  An object of class `sample_size_nbinom_summary`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the input object.

## Examples

``` r
x <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
s <- summary(x)
print(s)
#> 
#> 
```
