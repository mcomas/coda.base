test_that("formula groups are parsed literally and count each part once", {
  B <- sbp_basis(list(g = a ~ b, h = g + a ~ (c + c)), silent = TRUE)
  expected <- sbp_basis(cbind(c(1, -1, 0), c(1, 1, -1)), silent = TRUE)
  expect_equal(unname(B), unname(expected))
  expect_equal(unname(colSums(B)), c(0, 0), tolerance = 1e-14)
  expect_equal(unname(colSums(B^2)), c(1, 1), tolerance = 1e-14)

  B <- sbp_basis(list(`part one` ~ (`part two` + c)), silent = TRUE)
  expect_equal(rownames(B), c("part one", "part two", "c"))
  expect_equal(unname(B[, 1]), c(sqrt(2 / 3), -sqrt(1 / 6), -sqrt(1 / 6)))

  for (formula in list(a ~ a + b, a + b ~ b + c)) {
    expect_error(sbp_basis(list(formula), silent = TRUE), "disjoint")
  }
  expect_error(sbp_basis(list(g = a ~ b, h = g ~ a + c)), "disjoint")
  expect_error(sbp_basis(list(~ a + b)), "numerator and a denominator")
  for (formula in list(a ~ 1, a ~ b * c, a ~ log(b), a ~ b - c, a ~ b ~ c)) {
    expect_error(sbp_basis(list(formula)), "part names joined")
  }
})

test_that("balance references resolve forward and detect cycles", {
  formulas <- list(g = a ~ b, h = g ~ c, i = h ~ d)
  forward <- sbp_basis(formulas)
  reverse <- sbp_basis(rev(formulas))
  expect_equal(reverse[rownames(forward), colnames(forward)], forward)
  expect_error(sbp_basis(list(g = h ~ a, h = g ~ b)), "Circular balance references")
  expect_error(sbp_basis(list(g = g ~ a)), "Circular balance references")
  expect_error(sbp_basis(list(g = a ~ b), data = c("a", "c")), "Unknown.*b")

  # Long dependency chains must not consume one R call frame per balance.
  D <- 300L
  formulas <- lapply(seq_len(D - 1L), function(i) {
    stats::as.formula(if (i == 1L) "x1 ~ x2" else paste0("b", i - 1L, " ~ x", i + 1L))
  })
  names(formulas) <- paste0("b", seq_along(formulas))
  B <- sbp_basis(rev(formulas), silent = TRUE)
  expect_equal(dim(B), c(D, D - 1L))
  expect_equal(unname(colSums(B)), rep(0, D - 1L), tolerance = 1e-12)
  expect_equal(unname(colSums(B^2)), rep(1, D - 1L), tolerance = 1e-12)
})

test_that("names are validated and generated without ambiguity", {
  expect_error(sbp_basis(list(g = a ~ b, g = c ~ d)), "Balance names")
  formulas <- list(a ~ b)
  names(formulas) <- NA_character_
  expect_error(sbp_basis(formulas), "Balance names")
  expect_error(sbp_basis(list(a = a ~ b), data = c("a", "b")), "also be part names")
  for (parts in list(c("a", "a"), c("a", ""), c("a", NA_character_))) {
    expect_error(sbp_basis(list(a ~ b), data = parts), "Part names")
  }
  expect_error(sbp_basis(list(a ~ b), data = matrix(1, 2, 2)), "Part names")
  S <- matrix(c(1, -1), 2)
  rownames(S) <- c("a", "a")
  expect_error(sbp_basis(S), "Part names")
  rownames(S) <- c("a", "b")
  colnames(S) <- ""
  expect_error(sbp_basis(S), "Balance names")

  B <- sbp_basis(list(a ~ b, b1 = c ~ d), silent = TRUE)
  expect_equal(colnames(B), c("b1.1", "b1"))
  B <- sbp_basis(list(.b2 = a ~ b), data = c("a", "b", "c"), fill = TRUE)
  expect_equal(colnames(B), c(".b2", ".b2.1"))
})

test_that("completion preserves supplied columns and fills both branches", {
  S <- cbind(c(1, 1, 1, -1, -1, -1),
             c(1, 1, -1, 0, 0, 0),
             c(0, 0, 0, 1, 1, -1))
  rownames(S) <- letters[1:6]
  colnames(S) <- c("root", "left", "right")
  B <- sbp_basis(S, fill = TRUE)
  expect_equal(B[, 1:3], sbp_basis(S, silent = TRUE))
  expect_equal(rownames(B), rownames(S))
  expect_equal(colnames(B), c(colnames(S), ".b4", ".b5"))
  expect_equal(unname(crossprod(B)), diag(5), tolerance = 1e-14)

  B <- sbp_basis(list(a ~ b), data = c("a", "b", "c"), fill = TRUE)
  expect_equal(colnames(B), c("b1", ".b2"))
  expect_equal(unname(crossprod(B)), diag(2), tolerance = 1e-14)
  expect_equal(B, sbp_basis(list(a ~ b), data = c("a", "b", "c"),
                            fill = TRUE, silent = TRUE))
})

test_that("completion handles partial trees regardless of ordering and signs", {
  set.seed(131)
  for (D in 2:10) {
    full_bases <- list(sign(ilr_basis(D)))
    if (D > 2L) full_bases[[2]] <- sign(cdp_basis(D))
    for (full in full_bases) {
      for (iteration in 1:10) {
        selected <- sample.int(D - 1L, sample.int(D - 1L, 1))
        S <- sweep(full[, selected, drop = FALSE], 2,
                   sample(c(-1, 1), length(selected), replace = TRUE), "*")
        B <- sbp_basis(S, fill = TRUE)
        expect_equal(dim(B), c(D, D - 1L))
        expect_equal(B[, seq_along(selected), drop = FALSE], sbp_basis(S, silent = TRUE))
        expect_equal(unname(crossprod(B)), diag(D - 1L), tolerance = 1e-12)
        expect_equal(unname(colSums(B)), rep(0, D - 1L), tolerance = 1e-12)
      }
    }
  }
})

test_that("completion rejects incompatible balances even in silent mode", {
  for (S in list(cbind(c(1, -1, 0), c(0, 1, -1)),
                 cbind(c(1, -1, 0), c(1, -1, 0)),
                 cbind(c(1, 1, -1, -1), c(1, -1, 1, -1)))) {
    expect_error(sbp_basis(S, fill = TRUE, silent = TRUE), "hierarchical SBP")
  }
  # Also validate complete subtrees; they can be incompatible without needing
  # any new columns in that branch.
  expect_error(sbp_basis(cbind(c(1, 1, -1, -1), c(1, 0, -1, 0), c(1, 0, 0, -1)),
                         fill = TRUE, silent = TRUE), "hierarchical SBP")
})

test_that("matrix validation rejects invalid entries and preserves sign semantics", {
  for (x in list(NA_real_, NaN, Inf, -Inf)) {
    expect_error(sbp_basis(matrix(c(1, -1, x), 3), silent = TRUE), "finite")
  }
  expect_error(sbp_basis(matrix(c("1", "-1"), 2)), "numeric matrix")
  expect_error(sbp_basis(matrix(c(1i, -1i), 2)), "numeric matrix")
  expect_error(sbp_basis(matrix(c(1, 0), 2)), "positive and one negative")
  expect_error(sbp_basis(matrix(1, 1, 1)), "two parts")
  expect_error(sbp_basis(c(1, -1)), "numeric matrix.*list of formulas")
  expect_error(sbp_basis(list()), "non-empty list")
  expect_error(sbp_basis(list(a ~ b), fill = NA), "TRUE or FALSE")
  expect_error(sbp_basis(list(a ~ b), silent = c(TRUE, FALSE)), "TRUE or FALSE")
  S <- matrix(c(3, -5, 0), 3)
  expect_equal(sbp_basis(S, silent = TRUE), sbp_basis(sign(S), silent = TRUE))
})

test_that("empty sign matrices can be completed without dropping dimensions", {
  S <- matrix(numeric(), 4, 0, dimnames = list(letters[1:4], NULL))
  B <- sbp_basis(S, silent = TRUE)
  expect_equal(dim(B), c(4L, 0L))
  expect_equal(rownames(B), letters[1:4])
  B <- sbp_basis(S, fill = TRUE)
  expect_equal(colnames(B), paste0(".b", 1:3))
  expect_equal(unname(B), unname(ilr_basis(4)))
})

test_that("diagnostics distinguish partial, redundant and non-orthogonal systems", {
  S <- cbind(c(1, -1, 0), c(0, 1, -1))
  expect_warning(sbp_basis(S), "not orthogonal")
  expect_warning(sbp_basis(cbind(S, S[, 1])), "not a basis")
  expect_warning(sbp_basis(cbind(S[, 1], S[, 1])), "not a basis")
  expect_warning(sbp_basis(S[, 1, drop = FALSE]), "not a basis")
  expect_warning(sbp_basis(S, silent = TRUE), NA)
  expect_warning(sbp_basis(cbind(c(1, -1, 0)), fill = TRUE, silent = TRUE), NA)

  # Non-hierarchical orthonormal collections are still valid without fill.
  H <- cbind(c(1, 1, -1, -1), c(1, -1, 1, -1), c(1, -1, -1, 1))
  expect_warning(B <- sbp_basis(H), NA)
  expect_equal(unname(crossprod(B)), diag(3))
})
