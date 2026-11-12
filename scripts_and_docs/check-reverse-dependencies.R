# Run with network access:
# Rscript scripts_and_docs/check-reverse-dependencies.R package.tar.gz output-dir
# Uses the same current CRAN reverse-dependency sources for both versions.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript check-reverse-dependencies.R package.tar.gz output-dir")
}
candidate <- normalizePath(args[1L], mustWork = TRUE)
if (!grepl("^coda[.]base_[0-9].*[.]tar[.]gz$", basename(candidate))) {
  stop("Expected a coda.base source tarball")
}
dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
root <- normalizePath(args[2L], mustWork = TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
db <- available.packages(type = "source")
if (!"coda.base" %in% rownames(db)) {
  stop("Cannot retrieve the current CRAN package index; checks were not run")
}
reverse <- sort(unique(unlist(tools::package_dependencies(
  "coda.base", db = db, reverse = TRUE, recursive = FALSE,
  which = c("Depends", "Imports", "LinkingTo", "Suggests")
), use.names = FALSE)))
writeLines(reverse, file.path(root, "reverse-dependencies.txt"))
write.csv(db[c("coda.base", reverse), , drop = FALSE],
          file.path(root, "package-versions.csv"), row.names = FALSE)

# Download once so baseline and candidate use exactly the same sources.
sources <- file.path(root, "sources")
dir.create(sources, showWarnings = FALSE)
downloads <- download.packages(c("coda.base", reverse), sources,
                              type = "source", available = db)
if (!setequal(downloads[, 1L], c("coda.base", reverse))) {
  stop("Some source packages could not be downloaded")
}
reverse_files <- downloads[downloads[, 1L] %in% reverse, 2L]
baseline <- downloads[downloads[, 1L] == "coda.base", 2L]
original_libraries <- .libPaths()
for (mode in c("baseline", "candidate")) {
  directory <- file.path(root, mode)
  dir.create(directory, showWarnings = FALSE)
  files <- c(if (mode == "baseline") baseline else candidate, reverse_files)
  stopifnot(all(file.copy(files, directory, overwrite = TRUE)))
  library <- file.path(directory, "Library")
  dir.create(library, showWarnings = FALSE)
  .libPaths(c(library, original_libraries))
  status <- system2(file.path(R.home("bin"), "R"), c(
    "CMD", "INSTALL", shQuote(paste0("--library=", library)), shQuote(files[1L])
  ), stdout = file.path(root, paste0(mode, "-install.txt")), stderr = TRUE)
  if (status != 0L) stop("Failed to install ", mode, " coda.base")
  stopifnot(find.package("coda.base") == file.path(library, "coda.base"))
  tools::check_packages_in_dir(
    directory, pfiles = basename(files), check_args = "--no-manual",
    clean = FALSE
  )
  capture.output(tools::summarize_check_packages_in_dir_results(
    directory, full = TRUE
  ), file = file.path(root, paste0(mode, "-summary.txt")))
}

diagnostics <- function(mode, package) {
  log <- file.path(root, mode, paste0(package, ".Rcheck"), "00check.log")
  if (!file.exists(log)) stop("Missing check log: ", log)
  lines <- readLines(log, warn = FALSE)
  if (!any(grepl("^\\* DONE", lines))) stop("Incomplete check: ", log)
  starts <- which(grepl("^\\* checking .* \\.\\.\\. (ERROR|WARNING|NOTE)", lines))
  unlist(lapply(starts, function(i) {
    next_check <- which(seq_along(lines) > i & grepl("^\\* ", lines))
    end <- if (length(next_check)) next_check[1L] - 1L else length(lines)
    lines[i:end][nzchar(lines[i:end])]
  }), use.names = FALSE)
}
regressions <- character()
for (package in reverse) {
  added <- setdiff(diagnostics("candidate", package),
                   diagnostics("baseline", package))
  if (length(added)) {
    regressions <- c(regressions, paste0(package, ":"), added, "")
  }
}
writeLines(regressions, file.path(root, "new-diagnostics.txt"))
if (length(regressions)) {
  stop("New reverse-dependency diagnostics; review new-diagnostics.txt")
}
message("No new diagnostics in ", length(reverse), " reverse dependencies.")
