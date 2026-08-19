# Closure operation for compositional data

Applies the closure operation to a numeric vector, matrix or data frame
so that each composition sums to a prescribed constant `k`.

## Usage

``` r
closure(X, k = 1)
```

## Arguments

- X:

  A numeric vector, matrix, data frame, or an object coercible to one of
  these. For matrices and data frames, rows are interpreted as
  compositions.

- k:

  A numeric vector of length 1 or length `nrow(X)`. Must contain
  strictly positive values.

## Value

If `X` is a vector, a numeric vector of the same length.

If `X` is a matrix, a numeric matrix with the same dimensions, dimnames,
and row-wise sums equal to `k`.

If `X` is a data frame, a data frame with the same row and column names,
and row-wise sums equal to `k`.

## Details

If `X` is:

- a vector, the returned vector sums to `k`;

- a matrix or data frame, closure is applied row-wise, and each row sums
  to the corresponding value of `k`.

The argument `k` may be:

- a single positive number, recycled to all rows;

- a numeric vector of length `nrow(X)`, specifying a different closure
  constant for each row.

For a composition \\x = (x_1, \dots, x_D)\\ with positive sum, the
closure to constant \\k\\ is \$\$C(x) = k \frac{x}{\sum\_{j=1}^D
x_j}.\$\$

This function requires all entries of `X` to be finite and non-negative,
and every row sum (or the vector sum) must be strictly positive.

## Examples

``` r
closure(c(2, 3, 5))
#> [1] 0.2 0.3 0.5
closure(c(2, 3, 5), k = 100)
#> [1] 20 30 50

X <- matrix(c(1, 1, 2,
              2, 3, 5), nrow = 2, byrow = TRUE)
closure(X)
#>      [,1] [,2] [,3]
#> [1,] 0.25 0.25  0.5
#> [2,] 0.20 0.30  0.5
closure(X, k = c(1, 100))
#>       [,1]  [,2] [,3]
#> [1,]  0.25  0.25  0.5
#> [2,] 20.00 30.00 50.0

df <- data.frame(a = c(1, 2), b = c(1, 3), c = c(2, 5))
closure(df, k = 10)
#>     a   b c
#> 1 2.5 2.5 5
#> 2 2.0 3.0 5
```
