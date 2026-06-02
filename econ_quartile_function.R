
.wtd_quantile <- function(x, weights = NULL, probs, method = "type2") {
  if (is.null(weights)) {
    type_int <- as.integer(gsub("type", "", tolower(method)))
    return(quantile(x, probs = probs, type = type_int, na.rm = TRUE))
  }
  weights <- weights / sum(weights, na.rm = TRUE)
  ord     <- order(x)
  x       <- x[ord]
  w       <- weights[ord]
  cumw    <- cumsum(w)
  sapply(probs, function(p) {
    if (p <= 0) return(x[1L])
    if (p >= 1) return(x[length(x)])
    lo <- max(c(1L, which(cumw <= p)))
    hi <- min(c(length(x), which(cumw >= p)))
    if (lo == hi && cumw[lo] == p && lo < length(x)) hi <- lo + 1L
    if (lo == hi) x[lo] else mean(x[c(lo, hi)])
  })
}


.bootstrap_se <- function(x, weights = NULL, probs, B = 500, ci_level = 0.95) {
  n      <- length(x)
  boot_q <- replicate(B, {
    idx <- sample(n, replace = TRUE)
    .wtd_quantile(x[idx],
                  weights = if (!is.null(weights)) weights[idx] else NULL,
                  probs   = probs)
  })
  if (length(probs) == 1L) boot_q <- matrix(boot_q, nrow = 1L)
  alpha <- 1 - ci_level
  list(
    se    = apply(boot_q, 1L, sd,       na.rm = TRUE),
    lower = apply(boot_q, 1L, quantile, probs = alpha / 2,     na.rm = TRUE),
    upper = apply(boot_q, 1L, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  )
}


.lorenz_gini <- function(x, weights = NULL) {
  ord        <- order(x)
  x          <- x[ord]
  w          <- if (!is.null(weights)) weights[ord] / sum(weights)
                else rep(1 / length(x), length(x))
  cum_pop    <- cumsum(w)
  cum_income <- cumsum(x * w) / sum(x * w)
  gini       <- round(1 - 2 * sum(diff(c(0, cum_pop)) *
                (c(0, cum_income[-length(cum_income)]) + cum_income) / 2), 4)
  at_pts     <- seq(0, 1, 0.1)
  lorenz_y   <- approx(cum_pop, cum_income, xout = at_pts, rule = 2)$y
  lorenz_y   <- pmax(0, pmin(1, lorenz_y))
  lorenz_y[1] <- 0
  list(gini     = gini,
       lorenz_x = at_pts,
       lorenz_y = round(lorenz_y, 4))
}


.prepare_x <- function(x) {
  if (is.logical(x) || is.complex(x) ||
      inherits(x, c("Date", "POSIXct", "POSIXlt"))) {
    stop(
      "Unsupported data type: '", class(x)[1L], "'.\n",
      "econ_quartile() supports only:\n",
      "  1. numeric        — e.g. income, wages, prices\n",
      "  2. integer        — e.g. years of schooling, age, household size\n",
      "  3. ordered factor — e.g. education level, satisfaction rating\n",
      "For binary/logical data use mean() or prop.test().\n",
      "For datetime data use it as a grouping variable via 'by ='."
    )
  }
  if (is.factor(x)) {
    if (!is.ordered(x))
      stop("Factor must be ordered (use factor(..., ordered = TRUE)).\n",
           "Unordered factors have no natural ranking for quartiles.")
    return(list(x = as.numeric(x), dtype = "ordinal"))
  }
  dtype <- if (is.integer(x)) "integer" else if (is.numeric(x)) "numeric" else
           stop("Unsupported data type: '", class(x)[1L], "'.")
  list(x = as.numeric(x), dtype = dtype)
}


#' @export
print.econ_quartile <- function(x, ...) {
  cat("\n── econ_quartile ────────────────────────────────────────────────\n")
  cat(sprintf("  Data type : %s\n", x$data_type))
  cat(sprintf("  Method    : %s  |  n = %d  |  Weighted: %s\n",
              x$method, x$n, ifelse(x$weighted, "Yes", "No")))
  cat("\n  Quartiles:\n")
  print(round(x$quartiles, 4))
  if (!is.null(x$labels))
    cat("  Labels    :",
        paste(names(x$labels), x$labels, sep = "=", collapse = "  "), "\n")
  if (!is.null(x$ci)) {
    cat(sprintf("\n  Confidence Intervals (%.0f%%):\n", x$ci_level * 100))
    print(x$ci)
  }
  if (!is.null(x$inequality)) {
    iq <- x$inequality
    cat("\n  Inequality Measures:\n")
    cat(sprintf("    P90/P10 ratio  : %s\n",
        ifelse(is.na(iq$P90_P10_ratio), "NA (P10 = 0)", sprintf("%.4f", iq$P90_P10_ratio))))
    cat(sprintf("    P75/P25 ratio  : %s\n",
        ifelse(is.na(iq$P75_P25_ratio), "NA (P25 = 0)", sprintf("%.4f", iq$P75_P25_ratio))))
    cat(sprintf("    IQR            : %.4f\n", iq$IQR))
    if (!anyNA(c(iq$Q1_share, iq$Q2_share, iq$Q3_share, iq$Q4_share)))
      cat(sprintf(
        "    Income shares  : Q1=%.1f%%  Q2=%.1f%%  Q3=%.1f%%  Q4=%.1f%%\n",
        iq$Q1_share * 100, iq$Q2_share * 100,
        iq$Q3_share * 100, iq$Q4_share * 100))
    else
      cat("    Income shares  : NA (duplicate quartile boundaries)\n")
  }
  if (!is.null(x$lorenz)) {
    cat(sprintf("\n  Gini Coefficient: %.4f\n", x$lorenz$gini))
    cat("  Lorenz ordinates:\n")
    print(data.frame(population   = x$lorenz$lorenz_x,
                     income_share = x$lorenz$lorenz_y))
  }
  cat("─────────────────────────────────────────────────────────────────\n")
  invisible(x)
}


#' @export
print.econ_quartile_grouped <- function(x, ...) {
  cat("\n── econ_quartile: Grouped Results ──────────────────────────────\n")
  for (g in names(x)) {
    cat("\nGroup:", g, "\n")
    print(round(x[[g]]$quartiles, 4))
    if (!is.null(x[[g]]$inequality)) {
      iq      <- x[[g]]$inequality
      p_ratio <- if (is.na(iq$P75_P25_ratio)) "NA"
                 else sprintf("%.4f", iq$P75_P25_ratio)
      cat(sprintf("  IQR: %.4f  |  P75/P25: %s\n", iq$IQR, p_ratio))
    }
  }
  invisible(x)
}


econ_quartile <- function(x,
                          weights    = NULL,
                          by         = NULL,
                          method     = "type2",
                          se         = TRUE,
                          B          = 500,
                          ci_level   = 0.95,
                          inequality = TRUE,
                          lorenz     = FALSE,
                          probs      = c(0.25, 0.50, 0.75)) {

  if (length(x) < 10L)
    stop("Need at least 10 observations for reliable quartile estimation.\n",
         "  Provided: ", length(x), " observations.")

  if (!is.null(weights)) {
    stopifnot(
      "'weights' must be the same length as 'x'." = length(weights) == length(x),
      "All weights must be non-negative."          = !any(weights < 0, na.rm = TRUE)
    )
  }

  if (!is.null(by) && length(by) != length(x))
    stop("'by' must be the same length as 'x'.")

  method <- tolower(trimws(method))
  if (!method %in% paste0("type", 1:9))
    stop("'method' must be one of: ",
         paste(paste0("type", 1:9), collapse = ", "),
         " (case-insensitive).")

  ordinal_levels <- if (is.factor(x) && is.ordered(x)) levels(x) else NULL
  prep  <- .prepare_x(x)
  xnum  <- prep$x
  dtype <- prep$dtype

  nas <- is.na(xnum)
  if (any(nas))
    stop(
      "NA values detected. econ_quartile() does not remove NAs automatically.\n",
      "  Total NAs : ", sum(nas), "\n",
      "  Positions : ", paste(which(nas), collapse = ", "), "\n",
      "  Options   :\n",
      "    1. x <- x[!is.na(x)]\n",
      "    2. x[is.na(x)] <- median(x, na.rm = TRUE)"
    )

  if (!is.null(by)) {
    by_fac  <- factor(by, levels = unique(as.character(by)))
    groups  <- split(seq_along(xnum), by_fac)
    results <- lapply(names(groups), function(g) {
      idx <- groups[[g]]
      econ_quartile(
        x          = xnum[idx],
        weights    = if (!is.null(weights)) weights[idx] else NULL,
        by         = NULL,
        method     = method,
        se         = se,
        B          = B,
        ci_level   = ci_level,
        inequality = inequality,
        lorenz     = lorenz,
        probs      = probs
      )
    })
    names(results) <- names(groups)
    class(results) <- c("econ_quartile_grouped", "list")
    return(invisible(results))
  }

  prob_names    <- paste0(probs * 100, "%")
  q_vals        <- .wtd_quantile(xnum, weights, probs, method)
  names(q_vals) <- prob_names

  if (dtype == "ordinal" && !is.null(ordinal_levels)) {
    q_label <- sapply(round(q_vals), function(i)
      if (i >= 1L && i <= length(ordinal_levels)) ordinal_levels[i]
      else as.character(i))
    names(q_label) <- prob_names
  } else {
    q_label <- NULL
  }

  ci_out <- NULL
  if (se) {
    if (dtype == "ordinal") {
      warning("Bootstrap SE not computed for ordered factors: numeric codes ",
              "lack metric meaning. Set se = FALSE to suppress this warning.",
              call. = FALSE)
    } else {
      bs     <- .bootstrap_se(xnum, weights, probs, B, ci_level)
      ci_out <- data.frame(
        prob     = prob_names,
        estimate = round(q_vals,   4),
        SE       = round(bs$se,    4),
        lower_CI = round(bs$lower, 4),
        upper_CI = round(bs$upper, 4)
      )
      rownames(ci_out) <- NULL
    }
  }

  ineq_out <- NULL
  if (inequality && dtype != "ordinal") {
    q10 <- .wtd_quantile(xnum, weights, 0.10, method)
    q25 <- .wtd_quantile(xnum, weights, 0.25, method)
    q50 <- .wtd_quantile(xnum, weights, 0.50, method)
    q75 <- .wtd_quantile(xnum, weights, 0.75, method)
    q90 <- .wtd_quantile(xnum, weights, 0.90, method)

    p90_p10 <- if (isTRUE(q10 == 0)) {
      warning("P10 = 0; P90/P10 ratio set to NA to avoid Inf.", call. = FALSE)
      NA_real_
    } else round(q90 / q10, 4)

    p75_p25 <- if (isTRUE(q25 == 0)) {
      warning("P25 = 0; P75/P25 ratio set to NA to avoid Inf.", call. = FALSE)
      NA_real_
    } else round(q75 / q25, 4)

    breaks <- c(-Inf, q25, q50, q75, Inf)
    if (length(unique(breaks)) < length(breaks)) {
      warning("Duplicate quartile boundaries detected. Income shares set to NA.",
              call. = FALSE)
      sh <- rep(NA_real_, 4)
    } else {
      grp   <- cut(xnum, breaks = breaks, include.lowest = TRUE)
      total <- if (!is.null(weights)) sum(xnum * weights) else sum(xnum)
      sh    <- as.numeric(tapply(
        if (!is.null(weights)) xnum * weights else xnum,
        grp, sum) / total)
    }

    ineq_out <- list(
      P90_P10_ratio = p90_p10,
      P75_P25_ratio = p75_p25,
      IQR           = round(q75 - q25, 4),
      Q1_share      = round(sh[1L], 4),
      Q2_share      = round(sh[2L], 4),
      Q3_share      = round(sh[3L], 4),
      Q4_share      = round(sh[4L], 4)
    )
  }

  lorenz_out <- NULL
  if (lorenz && dtype != "ordinal")
    lorenz_out <- .lorenz_gini(xnum, weights)

  result <- list(
    quartiles  = q_vals,
    labels     = q_label,
    data_type  = dtype,
    method     = method,
    n          = length(xnum),
    weighted   = !is.null(weights),
    ci         = ci_out,
    ci_level   = ci_level,
    inequality = ineq_out,
    lorenz     = lorenz_out
  )
  class(result) <- c("econ_quartile", "list")
  invisible(result)
}
