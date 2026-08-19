# Plot a balance with node labels under horizontal branches

Plot a balance with node labels under horizontal branches

## Usage

``` r
plot_balance(
  B,
  data = NULL,
  main = "Balance dendrogram",
  summary_fun = NULL,
  cex_node = 0.9,
  offset_node = 0.05,
  ...
)
```

## Arguments

- B:

  Balance basis matrix

- data:

  Optional compositional data used to compute balance summaries

- main:

  Plot title

- summary_fun:

  Optional function applied to each balance coordinate vector. It must
  take a numeric vector and return a character string.

- cex_node:

  Character expansion for node labels

- offset_node:

  Vertical offset below the horizontal branch, relative to max height

- ...:

  Further arguments passed to plot

## Value

Invisibly returns a data.frame with node coordinates and labels

## Examples

``` r
X = waste[,5:9]
B = pb_basis(X, method = 'exact')

plot_balance(B)


plot_balance(B, data = X,
             summary_fun = function(x){
               q = quantile(x, probs = c(0.25, 0.5, 0.75))
               sprintf("%0.2f [%0.2f-%0.2f]", q[2], q[1], q[3])
             })

```
