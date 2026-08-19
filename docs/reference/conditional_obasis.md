# Conditional orthonormal basis

Compute orthonormal ilr bases adapted to row-wise conditioning patterns.

## Usage

``` r
conditional_obasis(X, scheme = c("zero", "zero_na"))
```

## Arguments

- X:

  A numeric matrix or data frame with one observation or conditioning
  pattern per row and one part per column.

- scheme:

  Character string indicating the conditioning scheme. Possible values
  are \`"zero"\` and \`"zero_na"\`. Default is \`"zero"\`.

## Value

A three-dimensional array of dimension \`(D - 1, D, nrow(X))\`, where
\`D\` is the number of parts. Each slice contains one orthonormal ilr
basis.

## Details

Each row of \`X\` defines one conditioning pattern on the parts of a
composition. According to \`scheme\`, the parts are split into ordered
blocks:

- \`"zero"\`: parts equal to \`0\` and parts with strictly positive
  values,

- \`"zero_na"\`: missing values (\`NA\`), zeros, and strictly positive
  values.

For each row, the function constructs an orthonormal basis of the
clr-plane preserving the block structure induced by the selected scheme.

Under \`scheme = "zero"\`, if a row contains \`nz\` zeros, then:

- the first \`nz - 1\` coordinates describe the internal log-ratio
  structure of the zero block,

- the coordinate \`nz\` describes the balance between the zero block and
  the positive block,

- the remaining coordinates describe the internal log-ratio structure of
  the positive block.

Under \`scheme = "zero_na"\`, the blocks are ordered as:

- missing values (\`NA\`),

- zeros,

- strictly positive values.

In this case:

- the first coordinates describe the internal structure of the \`NA\`
  block,

- the next coordinate contrasts the \`NA\` block with the positive
  block,

- the following coordinates describe the internal structure of the zero
  block,

- the next coordinate contrasts the zero block with the positive block,

- the remaining coordinates describe the internal structure of the
  positive block.

## Examples

``` r
C <- rbind(
  c(0, 0, 1, 1, 0),
  c(0, 1, 0, 1, 0)
)

conditional_obasis(C)
#> , , 1
#> 
#>           [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] 0.7071068 -0.7071068  0.0000000  0.0000000  0.0000000
#> [2,] 0.4082483  0.4082483  0.0000000  0.0000000 -0.8164966
#> [3,] 0.3651484  0.3651484 -0.5477226 -0.5477226  0.3651484
#> [4,] 0.0000000  0.0000000  0.7071068 -0.7071068  0.0000000
#> 
#> , , 2
#> 
#>           [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] 0.7071068  0.0000000 -0.7071068  0.0000000  0.0000000
#> [2,] 0.4082483  0.0000000  0.4082483  0.0000000 -0.8164966
#> [3,] 0.3651484 -0.5477226  0.3651484 -0.5477226  0.3651484
#> [4,] 0.0000000  0.7071068  0.0000000 -0.7071068  0.0000000
#> 

X <- rbind(
  c(1, NA, 0, 2),
  c(NA, 3, 0, 4),
  c(1, 2, 3, 4)
)

conditional_obasis(X, scheme = "zero_na")
#> , , 1
#> 
#>            [,1]      [,2]      [,3]       [,4]
#> [1,] -0.4082483 0.8164966 0.0000000 -0.4082483
#> [2,] -0.4082483 0.0000000 0.8164966 -0.4082483
#> [3,]  0.7071068 0.0000000 0.0000000 -0.7071068
#> 
#> , , 2
#> 
#>           [,1]       [,2]      [,3]       [,4]
#> [1,] 0.8164966 -0.4082483 0.0000000 -0.4082483
#> [2,] 0.0000000 -0.4082483 0.8164966 -0.4082483
#> [3,] 0.0000000  0.7071068 0.0000000 -0.7071068
#> 
#> , , 3
#> 
#>           [,1]       [,2]       [,3]       [,4]
#> [1,] 0.7071068 -0.7071068  0.0000000  0.0000000
#> [2,] 0.4082483  0.4082483 -0.8164966  0.0000000
#> [3,] 0.2886751  0.2886751  0.2886751 -0.8660254
#> 
```
