# Pairwise log-ratio generating system

Construct the system of all pairwise log-ratios between parts.

## Usage

``` r
pairwise_basis(dim)
```

## Arguments

- dim:

  Number of parts. It can be a single integer, a matrix or data frame,
  or a character vector of part names.

## Value

A matrix, or a sparse matrix for large dimensions, whose columns
represent all pairwise log-ratio generators.
