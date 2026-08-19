# Perturbation of compositional data

The perturbation operation combines two compositions by component-wise
multiplication and then applies closure to ensure the result remains a
valid composition.

## Usage

``` r
perturbation(X, Y)
```

## Arguments

- X:

  A numeric vector, matrix or data.frame containing compositions.

- Y:

  A numeric vector, matrix or data.frame with the same number of parts
  as `X`. If one input is a matrix or data.frame and the other is a
  vector, the vector is applied to every row. If both inputs are
  matrices or data.frames, they must have the same dimensions, or `Y`
  may contain a single composition to be applied to all rows of `X`.

## Value

An object with the same format as `X` containing the perturbed
compositions, except that vector `X` with matrix or data.frame `Y`
returns the same rectangular format as `Y`.

## Details

Perturbation is the analogue of addition in the simplex. Each part of
`X` is multiplied by the corresponding part of `Y`, and the result is
closed with
[`closure`](https://mcomas.net/coda.base/reference/closure.md) so that
each composition has constant sum.

## Examples

``` r
x <- c(a = 1, b = 2, c = 3)
y <- c(a = 1, b = 1, c = 2)
perturbation(x, y)
#> [1] 0.1111111 0.2222222 0.6666667

X <- rbind(
  c(1, 2, 3),
  c(4, 5, 6)
)
perturbation(X, c(1, 1, 2))
#>           [,1]      [,2]      [,3]
#> [1,] 0.1111111 0.2222222 0.6666667
#> [2,] 0.1904762 0.2380952 0.5714286
perturbation(c(1, 1, 2), X)
#>           [,1]      [,2]      [,3]
#> [1,] 0.1111111 0.2222222 0.6666667
#> [2,] 0.1904762 0.2380952 0.5714286
```
