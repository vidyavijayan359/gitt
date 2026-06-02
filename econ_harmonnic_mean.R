

print.harmonic_mean <- function(x, ...) {
  cat("Harmonic mean:", round(unclass(x), 4), "\n")
  invisible(x)
}

harmonic_mean <- function(x,
                          na.rm    = TRUE,
                          inf.rm   = FALSE,
                          zero.rm  = FALSE,
                          weights  = NULL,
                          negative = c("warn", "allow"),
                          trim     = 0,
                          margin   = NULL,
                          verbose  = FALSE) {

  negative <- match.arg(negative)

  if (is.null(x) || length(x) == 0)
    stop("Input 'x' is NULL or empty.")

  if (!is.numeric(trim) || length(trim) != 1 || trim < 0 || trim >= 0.5)
    stop("`trim` must be a single number in [0, 0.5).")

  if (!is.null(weights) && trim > 0)
    stop("`weights` and `trim` cannot be used together.\n",
         "  Use `weights` for unequal population shares (surveys, panels).\n",
         "  Use `trim` to discard outliers or measurement error.\n",
         "  Choose one and leave the other at its default.")

  # ── 3D array dispatch ──────────────────────────────────────────────────────
  if (is.array(x) && length(dim(x)) == 3) {

    if (!is.null(weights))
      warning("`weights` not supported for array input and will be ignored.")

    if (is.character(x)) {
      if (verbose) message("Array is character type — returning NA.")
      return(NA_real_)
    }

    if (length(dim(x)) > 3)
      stop("Arrays with more than 3 dimensions are not supported.\n",
           "  Flatten first with: as.vector(x)")

    if (is.null(margin)) {
      if (verbose) message("3D array detected — flattening all values into one result.")
      return(harmonic_mean(as.vector(x), na.rm = na.rm, inf.rm = inf.rm,
                           zero.rm = zero.rm, negative = negative,
                           trim = trim, verbose = verbose))
    }

    if (!margin %in% 1:3)
      stop("`margin` for a 3D array must be 1 (rows), 2 (columns), or 3 (layers).")

    if (verbose) message("3D array detected — computing across margin ", margin, ".")

    result <- apply(x, margin, function(slice) {
      if (is.character(slice) || is.factor(slice)) return(NA_real_)
      tryCatch(
        harmonic_mean(as.numeric(slice), na.rm = na.rm, inf.rm = inf.rm,
                      zero.rm = zero.rm, negative = negative,
                      trim = trim, verbose = FALSE),
        error = function(e) NA_real_
      )
    })

    result <- as.numeric(result)

    dn <- dimnames(x)
    if (!is.null(dn) && !is.null(dn[[margin]]))
      names(result) <- dn[[margin]]
    else
      names(result) <- paste0("slice_", seq_along(result))

    return(result)
  }

  # ── matrix / data frame dispatch ───────────────────────────────────────────
  if (!is.null(nrow(x))) {
    if (!is.null(weights))
      warning("`weights` not supported for column-wise computation and will be ignored.")
    if (verbose) message("Matrix / data frame detected — computing column-wise.")
    result <- sapply(seq_len(ncol(x)), function(j) {
      col <- if (is.data.frame(x)) x[[j]] else x[, j]
      if (is.character(col) || is.factor(col)) {
        if (verbose) message("Column ", j, " is non-numeric — returning NA.")
        return(NA_real_)
      }
      harmonic_mean(col, na.rm = na.rm, inf.rm = inf.rm, zero.rm = zero.rm,
                    weights = NULL, negative = negative, trim = trim, verbose = verbose)
    })
    names(result) <- if (!is.null(colnames(x))) colnames(x) else paste0("V", seq_along(result))
    return(result)
  }

  # ── type coercion ──────────────────────────────────────────────────────────
  if (is.ts(x))      x <- as.numeric(x)
  if (is.logical(x)) x <- as.numeric(x)

  if (is.factor(x))
    stop("Factor input is not supported. Convert first: as.numeric(as.character(x))")
  if (is.character(x))
    stop("Character input is not supported. Fix upstream with: as.numeric(x)")
  if (is.complex(x))
    stop("Complex input is not supported. Use Re(x) if you need only the real part.")
  if (is.list(x))
    stop("List input is not supported. Try: as.numeric(unlist(x))")
  if (inherits(x, c("Date", "POSIXct", "POSIXlt")))
    stop("Date/time input is not meaningful for harmonic mean.")
  if (!is.numeric(x))
    stop("Unsupported type '", class(x)[1], "'.")

  # ── NA handling ────────────────────────────────────────────────────────────
  x[is.nan(x)] <- NA_real_
  n_na <- sum(is.na(x))
  if (!na.rm && n_na > 0) return(NA_real_)
  if (na.rm && n_na > 0) {
    if (verbose) message("Removing ", n_na, " NA/NaN value(s).")
    x <- x[!is.na(x)]
  }
  if (length(x) == 0) stop("No valid observations after NA removal.")

  # ── Inf handling ───────────────────────────────────────────────────────────
  n_inf <- sum(is.infinite(x))
  if (n_inf > 0) {
    if (!inf.rm)
      stop(n_inf, " infinite value(s) detected in 'x'.\n",
           "  1. Clean the source data and investigate why Inf appeared.\n",
           "  2. Set inf.rm = TRUE to remove them and proceed.")
    if (verbose) message("Removing ", n_inf, " Inf/-Inf value(s).")
    x <- x[!is.infinite(x)]
  }
  if (length(x) == 0) stop("No valid observations after Inf removal.")

  # ── zero handling ──────────────────────────────────────────────────────────
  n_zero <- sum(x == 0)
  if (n_zero > 0) {
    if (!zero.rm)
      stop(n_zero, " zero(s) found — harmonic mean is undefined when zeros are present.\n",
           "  Set zero.rm = TRUE to exclude them.")
    if (verbose) message("Removing ", n_zero, " zero(s).")
    x <- x[x != 0]
  }
  if (length(x) == 0) stop("No valid observations after zero removal.")

  # ── negative handling ──────────────────────────────────────────────────────
  n_neg <- sum(x < 0)
  if (n_neg > 0) {
    msg <- sprintf("%d negative value(s) found. Harmonic mean may not be meaningful.", n_neg)
    switch(negative,
           "warn"  = warning(msg, call. = FALSE),
           "allow" = if (verbose) message(msg))
  }

  # ── trimming ───────────────────────────────────────────────────────────────
  n_trimmed <- 0L
  if (trim > 0) {
    n_each <- floor(length(x) * trim)
    if (n_each > 0) {
      x         <- sort(x)[(n_each + 1):(length(x) - n_each)]
      n_trimmed <- 2L * n_each
      if (verbose) message("Trimmed ", n_trimmed, " observations (", n_each, " from each tail).")
    }
  }
  if (length(x) == 0) stop("No observations remaining after trimming.")

  if (length(x) == 1) {
    warning("Only 1 observation remains — result equals that value. Check your data.")
    return(x[[1L]])
  }

  # ── weight validation ──────────────────────────────────────────────────────
  if (!is.null(weights)) {
    if (anyNA(weights))        stop("Weights contain NA.")
    if (length(weights) != length(x))
      stop("Weights length (", length(weights), ") must equal observations after filtering (", length(x), ").\n",
           "  Apply the same filters to your weight vector before passing it in.")
    if (any(weights < 0))      stop("Weights must be non-negative.")
    if (sum(weights) == 0)     stop("Sum of weights is zero — cannot normalise.")
    if (verbose && abs(sum(weights) - 1) > sqrt(.Machine$double.eps))
      message("Normalising weights (original sum = ", round(sum(weights), 6), ").")
    weights <- weights / sum(weights)
  }

  # ── computation ────────────────────────────────────────────────────────────
  result <- if (!is.null(weights))
    1 / sum(weights / x)
  else
    length(x) / sum(1 / x)

  if (verbose) message("Harmonic mean = ", round(result, 6), " (n = ", length(x), ")")

  # ── metadata ───────────────────────────────────────────────────────────────
  attr(result, "n")        <- length(x)
  attr(result, "n_na")     <- n_na
  attr(result, "n_inf")    <- n_inf
  attr(result, "n_zero")   <- n_zero
  attr(result, "n_neg")    <- n_neg
  attr(result, "n_trim")   <- n_trimmed
  attr(result, "weighted") <- !is.null(weights)
  class(result)            <- "harmonic_mean"

  result
}
