# Mammals' milk

The \`mammals_milk\` data set contains the percentages of five
constituents of the milk of 26 mammals: \\\[W, P, F, L, A\]\\, where
\`W\` is water, \`P\` is protein, \`F\` is fat, \`L\` is lactose, and
\`A\` is ash.

The blank protein entry for hippopotamus and the blank lactose entry for
seal are recorded as zeros. Before log-ratio
analysis, these values can be imputed using
[`coda_replacement()`](https://mcomas.net/coda.base/reference/coda_replacement.md).

## Usage

``` r
mammals_milk
```

## Format

An object of class `data.frame` with 26 rows and 6 columns.

## Source

Spector, W. S. (ed.) (1956). *Handbook of Biological Data*.
W. B. Saunders, Philadelphia, Table 35, *Milk, chemical composition:
various animals*, p. 50. Values are grams per 100 grams of whole milk,
transcribed from the table.
