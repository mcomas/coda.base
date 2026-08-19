# Tabu search for a partial principal balance on grouped parts

Finds a single grouped balance by tabu search over a partition of
selected parts. The search is carried out on groups of parts defined by
`lI`, using configurable neighbourhood moves.

## Usage

``` r
partial_pb_tabu_search(
  X,
  lI = NULL,
  min_parts = 2,
  max_parts = NULL,
  iter = 100,
  tabu_size = length(lI),
  ini = NULL,
  remove_active = TRUE,
  add_left = TRUE,
  add_right = TRUE,
  flip_side = FALSE,
  swap_zero = FALSE,
  swap_sides = FALSE,
  debug = FALSE,
  constrained.criterion = "variance"
)
```

## Arguments

- X:

  A numeric matrix with strictly positive finite entries. Rows are
  observations and columns are compositional parts.

- lI:

  A list defining a partition of a subset of the columns of `X`. If
  `NULL`, each column of `X` is used as a singleton group.

- min_parts:

  Integer. Minimum number of active groups.

- max_parts:

  Integer or `NULL`. Maximum number of groups from `lI` allowed to be
  active in the balance. If `NULL`, all groups may be active.

- iter:

  Integer. Maximum number of tabu search iterations.

- tabu_size:

  Integer. Maximum size of the tabu list. If `0`, no tabu memory is used
  and the algorithm performs a greedy local search until no neighbouring
  balance improves the criterion.

- ini:

  Initial grouped split. If `NULL`, the constrained principal balance of
  the grouped subcomposition is used.

- remove_active:

  Logical. Allow moves from `-1` or `+1` to `0`.

- add_left:

  Logical. Allow moves from `0` to `-1`.

- add_right:

  Logical. Allow moves from `0` to `+1`.

- flip_side:

  Logical. Allow direct moves from `-1` to `+1` and from `+1` to `-1`.

- swap_zero:

  Logical. Allow swaps between one active group and one inactive group,
  preserving the active side.

- swap_sides:

  Logical. Allow swaps between one left group and one right group.

- debug:

  Logical. If `TRUE`, progress information is printed during the search.

- constrained.criterion:

  Criterion used to initialise the constrained balance when
  `ini = NULL`. Either `"variance"` (default) or `"angle"`.

## Value

A list with the selected balance, its variance criterion, the search
path, and a `neighbourhoods` element recording the active neighbourhood
types.

## Details

When `ini = NULL`, the constrained grouped balance is adjusted greedily
so that the initial solution has exactly `max_parts` active groups.
