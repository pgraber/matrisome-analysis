# Continuous validation helper (base R, no dependencies).
#
# Every analysis step whose expected shape is known in advance asserts it here, immediately after the
# step runs. The expected value must come from outside the code that produced the object (an
# annotation file, a sample sheet, a published n), otherwise the check passes by construction.
#
# Usage:
#   source("scripts/validate.R")
#   vcheck_init("results/checks", step = "01_load_counts")
#   vcheck("annotated transcripts", nrow(df), 50000, source = "GENCODE v44 GTF")
#   vcheck("samples retained", length(unique(df$sample_id)), 24, source = "samplesheet.csv")
#   vcheck("IDs preserved after join", setequal(df$id, ref$id), TRUE)
#   vcheck_summary()
#
# A failed check stops the run. Nothing downstream should execute on data of the wrong shape.

.vstate <- new.env(parent = emptyenv())
.vstate$log <- NULL
.vstate$tsv <- NULL
.vstate$step <- NA_character_
.vstate$n_pass <- 0L
.vstate$n_fail <- 0L

#' Start a validation log.
#'
#' @param prefix Path prefix without extension. Writes `<prefix>.log` (human readable) and
#'   `<prefix>.tsv` (machine readable, read by the write-up and review steps).
#' @param step Name of the analysis step, recorded on every line.
#' @param append Keep an existing log (TRUE) or start a fresh one (FALSE, the default).
vcheck_init <- function(prefix = "results/checks", step = NA_character_, append = FALSE) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  .vstate$log <- paste0(prefix, ".log")
  .vstate$tsv <- paste0(prefix, ".tsv")
  .vstate$step <- step
  .vstate$n_pass <- 0L
  .vstate$n_fail <- 0L
  if (!append || !file.exists(.vstate$tsv)) {
    cat("", file = .vstate$log)
    cat("status\ttimestamp\tstep\tlabel\tactual\texpected\tsource\n", file = .vstate$tsv)
  }
  header <- sprintf("=== validation: %s (%s) ===",
                    if (is.na(step)) "analysis" else step,
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  .vwrite(header)
  invisible(prefix)
}

.vwrite <- function(line) {
  cat(line, "\n", sep = "")
  if (!is.null(.vstate$log)) cat(line, "\n", sep = "", file = .vstate$log, append = TRUE)
}

.vfmt <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.logical(x) && length(x) == 1) return(if (isTRUE(x)) "TRUE" else "FALSE")
  if (is.numeric(x) && length(x) == 1) {
    if (is.na(x)) return("NA")
    if (x == round(x) && abs(x) < 1e15) return(format(x, big.mark = ",", scientific = FALSE, trim = TRUE))
    return(format(signif(x, 6), scientific = FALSE, trim = TRUE))
  }
  paste(utils::head(as.character(x), 5), collapse = ", ")
}

#' Assert one expected value and record the result.
#'
#' @param label Plain-language description of what is being checked.
#' @param actual The observed value: a count, a proportion, or a logical.
#' @param expected The value known in advance. Use TRUE for logical assertions.
#' @param tol Absolute tolerance for numeric comparison (default 0, exact).
#' @param source Where the expected value came from. Recorded in the log; strongly recommended.
#' @param warn_only If TRUE, a mismatch is recorded as WARN and the run continues. Use sparingly and
#'   only for informational counts, never for a shape the downstream code depends on.
vcheck <- function(label, actual, expected, tol = 0, source = NA_character_, warn_only = FALSE) {
  if (is.null(.vstate$log)) vcheck_init()

  ok <- if (is.logical(expected) && length(expected) == 1) {
    isTRUE(all.equal(isTRUE(actual), isTRUE(expected)))
  } else if (is.numeric(expected) && is.numeric(actual) &&
             length(expected) == 1 && length(actual) == 1) {
    # Guard the comparison against binary floating-point representation error, so a value exactly
    # on the tolerance boundary (0.51 against 0.50 +/- 0.01) passes.
    slack <- tol * 1e-9 + abs(expected) * 1e-12
    !is.na(actual) && !is.na(expected) && abs(actual - expected) <= tol + slack
  } else {
    identical(actual, expected)
  }

  status <- if (ok) "PASS" else if (warn_only) "WARN" else "FAIL"
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  detail <- sprintf("%s: %s (expected %s%s)", label, .vfmt(actual), .vfmt(expected),
                    if (tol > 0) sprintf(" +/- %s", .vfmt(tol)) else "")
  if (!is.na(source)) detail <- sprintf("%s  [%s]", detail, source)
  .vwrite(sprintf("%-4s %s  %s", status, ts, detail))

  if (!is.null(.vstate$tsv)) {
    cat(paste(status, ts, .vstate$step, label, .vfmt(actual), .vfmt(expected),
              if (is.na(source)) "" else source, sep = "\t"), "\n",
        sep = "", file = .vstate$tsv, append = TRUE)
  }

  if (ok) {
    .vstate$n_pass <- .vstate$n_pass + 1L
  } else {
    .vstate$n_fail <- .vstate$n_fail + 1L
    if (!warn_only) {
      stop(sprintf("Validation failed | %s", detail), call. = FALSE)
    }
  }
  invisible(ok)
}

#' Close the log with a one-line tally. Returns the counts invisibly.
vcheck_summary <- function() {
  .vwrite(sprintf("--- %d checks: %d passed, %d failed ---",
                  .vstate$n_pass + .vstate$n_fail, .vstate$n_pass, .vstate$n_fail))
  invisible(list(pass = .vstate$n_pass, fail = .vstate$n_fail))
}
