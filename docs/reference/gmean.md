# Geometric Mean

Generic function for the (trimmed) geometric mean.

## Usage

``` r
gmean(x, zero.rm = FALSE, trim = 0, na.rm = FALSE)
```

## Arguments

- x:

  A nonnegative vector.

- zero.rm:

  a logical value indicating whether zero values should be stripped
  before the computation proceeds.

- trim:

  the fraction (0 to 0.5) of observations to be trimmed from each end of
  x before the mean is computed. Values of trim outside that range are
  taken as the nearest endpoint.

- na.rm:

  a logical value indicating whether NA values should be stripped before
  the computation proceeds.

## See also

[`center`](https://mcomas.net/coda.base/reference/center.md)
