# coda.base

# Log-Ratio Coordinates for Compositional Data

This R package provides tools for analyzing compositional data using log-ratio coordinates. It enables users to define coordinate systems tailored to compositional datasets and to generate coordinates based on these systems.

The package focuses on the construction and application of orthonormal and non-orthonormal coordinate systems for representing compositions in real space, facilitating advanced statistical modeling and interpretation.

## Key Functions

- `ilr_basis()`, `alr_basis()`, `clr_basis()` – classical log-ratio bases:
  - **ILR/OLR (Isometric/Orthonormal Log-Ratio)**: orthonormal basis.
  - **ALR (Additive Log-Ratio)**: basis with respect to a reference part.
  - **CLR (Centered Log-Ratio)**: non-orthonormal, but symmetrically treats all parts.

- `pc_basis()`, `pb_basis()`, `pw_basis()` – domain-specific and data-driven bases:
  - **PC basis**: based on principal component analysis of log-ratio coordinates.
  - **PB basis (Principal balances)**: variability explain by balances.
  - **PW (Pairwise Log-Ratios)**: pairwise comparisons.

- `coordinates(x, basis)`: expresses a composition `x` in coordinates with respect to a given `basis`.

- `composition(z, basis)`: reconstructs a composition from coordinates `z` and the associated `basis`.

## Example

```r
library(coda.base)

# Define a simple 3-part composition
x <- c('a' = 0.2, 'b' = 0.3, 'c' = 0.5)

# Create an ILR basis and express x in ILR coordinates
B <- ilr_basis(x)
h <- coordinates(x, B)  
h

# Recover the original composition
composition(h, B)
```


## Installation

You can install the development version from GitHub:

```r
# Install development version from GitHub
remotes::install_github("mcomas/coda.base")
```

and the cran version with:


``` r
# Install release version from CRAN
install.packages("coda.base")
```
