
decile <- function(x,
                   weights = NULL,
                   ties = c("mean", "low", "high", "first", "last", "random"),
                   labels = TRUE,
                   summary = FALSE,
                   right = TRUE) {
  
  x_name <- deparse(substitute(x))
  
  x <- .decile_validate_type(x, x_name)
  
  n <- length(x)
  if (n < 10L) {
    stop(
      sprintf("`%s` has only %d element(s). At least 10 observations are required ", x_name, n),
      "to form meaningful deciles.",
      call. = FALSE
    )
  }
  
  .decile_check_bad_values(x, x_name)
  
  u <- length(unique(x))
  if (u < 10L) {
    stop(
      sprintf("`%s` has only %d unique value(s). At least 10 unique values are required ", x_name, u),
      "to form 10 distinct decile intervals.",
      call. = FALSE
    )
  }
  
  if (!is.null(weights)) {
    weights <- .decile_validate_weights(weights, n, x_name)
  }
  
  ties <- match.arg(ties)
  
  decile_int <- .decile_compute(x, weights, ties, right)
  
  result <- .decile_label(decile_int, labels, ties)
  
  if (isTRUE(summary)) {
    .decile_print_summary(x, result, weights, right)
    return(invisible(result))
  }
  
  result
}


.decile_validate_type <- function(x, x_name) {
  
  if (is.numeric(x) && !is.matrix(x) && !is.data.frame(x)) {
    return(x)
  }
  
  if (is.logical(x)) {
    warning(
      sprintf("`%s` is logical; coercing to numeric (FALSE = 0, TRUE = 1). ", x_name),
      "Deciles on a binary variable may not be meaningful.",
      call. = FALSE
    )
    return(as.numeric(x))
  }
  
  if (is.factor(x)) {
    if (!is.ordered(x)) {
      stop(
        sprintf("`%s` is an *unordered* factor. Only ordered factors are supported ", x_name),
        "because unordered factors have no meaningful numeric ranking.\n",
        " Hint: convert with `factor(..., ordered = TRUE, levels = <your order>)`.",
        call. = FALSE
      )
    }
    warning(
      sprintf("`%s` is an ordered factor; coercing to integer via `as.integer()`. ", x_name),
      "Decile ranks reflect the ordering of factor levels.",
      call. = FALSE
    )
    return(as.integer(x))
  }
  
  if (inherits(x, c("Date", "POSIXct", "POSIXlt", "difftime"))) {
    warning(
      sprintf("`%s` is of class '%s'; coercing to numeric. ", x_name, paste(class(x), collapse = " / ")),
      "Deciles reflect chronological order.",
      call. = FALSE
    )
    return(as.numeric(x))
  }
  
  if (inherits(x, "ts")) {
    warning(
      sprintf("`%s` is a `ts` object; extracting numeric values.", x_name),
      call. = FALSE
    )
    return(as.numeric(x))
  }
  
  detected_class <- paste(class(x), collapse = " / ")
  
  hint <- if (is.character(x)) {
    "\n Hint: if the strings represent numbers, use `as.numeric(x)` first."
  } else if (is.matrix(x) || is.data.frame(x)) {
    "\n Hint: pass a single column, e.g. `decile(df$income)`."
  } else if (is.complex(x)) {
    "\n Hint: complex numbers have no natural ordering; use `Re(x)` or `Mod(x)`."
  } else if (is.list(x)) {
    "\n Hint: extract the numeric element, e.g. `decile(x[[1]])`."
  } else {
    ""
  }
  
  stop(
    sprintf("`%s` has unsupported type '%s'.\n", x_name, detected_class),
    " Supported types: numeric, integer, double, logical, ordered factor,\n",
    " Date, POSIXct, POSIXlt, difftime, ts.",
    hint,
    call. = FALSE
  )
}


.decile_check_bad_values <- function(x, x_name) {
  
  is_na <- is.na(x)
  is_nan <- is.nan(x)
  is_inf <- is.infinite(x)
  is_na_proper <- is_na & !is_nan
  
  problems <- list(
    "NA" = which(is_na_proper),
    "NaN" = which(is_nan),
    "Inf" = which(is_inf & x > 0),
    "-Inf" = which(is_inf & x < 0)
  )
  
  bad_counts <- vapply(problems, length, integer(1L))
  
  if (any(bad_counts > 0L)) {
    lines <- character(0L)
    for (kind in names(problems)) {
      idx <- problems[[kind]]
      if (length(idx) == 0L) next
      pos_str <- if (length(idx) <= 10L) {
        paste(idx, collapse = ", ")
      } else {
        paste0(paste(idx[1:10], collapse = ", "), " ... [", length(idx) - 10L, " more]")
      }
      lines <- c(lines, sprintf(" %s : %d value(s) at position(s): %s", kind, length(idx), pos_str))
    }
    
    stop(
      sprintf("`%s` contains bad value(s) that prevent decile computation:\n", x_name),
      paste(lines, collapse = "\n"), "\n\n",
      " Please remove or impute these values before calling decile().\n",
      " Tip: use `na.omit(x)`, `x[!is.na(x)]`, or an imputation package.",
      call. = FALSE
    )
  }
  
  invisible(NULL)
}


.decile_validate_weights <- function(w, n, x_name) {
  
  if (!is.numeric(w)) {
    stop("`weights` must be a numeric vector.", call. = FALSE)
  }
  if (length(w) != n) {
    stop(
      sprintf("`weights` has length %d but `%s` has length %d. They must match.", length(w), x_name, n),
      call. = FALSE
    )
  }
  
  .decile_check_bad_values(w, "weights")
  
  if (any(w < 0)) {
    stop(
      sprintf("All weights must be >= 0. Found %d negative weight(s).", sum(w < 0)),
      call. = FALSE
    )
  }
  if (sum(w) == 0) {
    stop("The sum of `weights` is zero; weights cannot all be zero.", call. = FALSE)
  }
  
  as.numeric(w)
}


.decile_compute <- function(x, weights, ties, right) {
  
  probs <- seq(0, 1, by = 0.1)
  
  if (is.null(weights)) {
    breaks <- stats::quantile(x, probs = probs, na.rm = FALSE, names = FALSE)
  } else {
    breaks <- .decile_weighted_quantile(x, weights, probs)
  }
  
  breaks[1L] <- breaks[1L] - .Machine$double.eps * abs(breaks[1L])
  breaks[length(breaks)] <- breaks[length(breaks)] + .Machine$double.eps * abs(breaks[length(breaks)])
  
  if (ties == "mean") {
    return(base::cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE, right = right))
  }
  
  r <- base::rank(x, ties.method = ties, na.last = "keep")
  n <- length(x)
  dec <- ceiling(r / n * 10)
  dec[dec == 0L] <- 1L
  as.integer(dec)
}


.decile_weighted_quantile <- function(x, w, probs) {
  ord <- order(x)
  xs <- x[ord]
  ws <- w[ord]
  cum_w <- cumsum(ws)
  total_w <- cum_w[length(cum_w)]
  cdf <- cum_w / total_w
  
  vapply(probs, function(p) {
    if (p <= cdf[1L]) return(xs[1L])
    if (p >= cdf[length(cdf)]) return(xs[length(xs)])
    idx <- which(cdf >= p)[1L]
    lo <- cdf[idx - 1L]
    hi <- cdf[idx]
    xs[idx - 1L] + (xs[idx] - xs[idx - 1L]) * (p - lo) / (hi - lo)
  }, numeric(1L))
}


.decile_label <- function(dec, labels, ties) {
  if (!isTRUE(labels)) return(dec)
  lev <- paste0("D", 1:10)
  factor(paste0("D", dec), levels = lev, ordered = TRUE)
}


.decile_print_summary <- function(x, result, weights, right) {
  
  probs <- seq(0, 1, by = 0.1)
  breaks <- if (is.null(weights)) {
    stats::quantile(x, probs = probs, names = TRUE)
  } else {
    q <- .decile_weighted_quantile(x, weights, probs)
    stats::setNames(q, paste0(probs * 100, "%"))
  }
  
  print(sprintf("n = %d | %s | intervals: %s",
                length(x),
                if (!is.null(weights)) "weighted" else "unweighted",
                if (right) "(a, b]" else "[a, b)"))
  
  print(table(result))
  
  bp <- data.frame(
    Decile = paste0("D", 1:10),
    Lower = round(breaks[1:10], 6),
    Upper = round(breaks[2:11], 6)
  )
  rownames(bp) <- NULL
  print(bp, row.names = FALSE)
  
  invisible(NULL)
}

