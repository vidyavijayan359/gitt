
# Functions
.check_na <- function(x, var_name = "x") {
  na_idx <- which(is.na(x))
  if (length(na_idx) > 0L)
    stop(sprintf(
      "NA values detected in '%s' at position(s): %s\n  Handle missing values before calling percentile functions.\n  Options: na.omit(), tidyr::drop_na(), or imputation.",
      var_name, paste(na_idx, collapse = ", ")),
      call. = FALSE)
  invisible(NULL)
}

.check_weights <- function(w, n, var_name = "weights") {
  if (!is.null(w)) {
    if (length(w) != n)
      stop(sprintf("'%s' must have the same length as the data (%d).", var_name, n), call. = FALSE)
    .check_na(w, var_name)
    if (any(w < 0))
      stop(sprintf("'%s' must be non-negative.", var_name), call. = FALSE)
    if (sum(w) == 0)
      stop(sprintf("Sum of '%s' is zero.", var_name), call. = FALSE)
  }
  invisible(NULL)
}

.check_probs <- function(probs) {
  if (any(probs < 0 | probs > 1))
    stop("'probs' must be in [0, 1].", call. = FALSE)
  invisible(NULL)
}

.check_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || digits < 0 || digits != as.integer(digits))
    stop("'digits' must be a single non-negative integer.", call. = FALSE)
  invisible(NULL)
}

.check_cols <- function(df, cols, supported_cls = c("numeric", "integer", "Date", "POSIXct", "currency")) {
  cols[!sapply(cols, function(cn) {
    if (!cn %in% names(df)) return(FALSE)
    cl <- class(df[[cn]])[1L]
    cl %in% supported_cls || is.ordered(df[[cn]])
  })]
}

.to_numeric <- function(x) {
  cls <- class(x)[1L]
  if (cls %in% c("numeric", "integer", "currency")) return(list(num = as.numeric(x), cls = cls))
  if (cls == "Date")     return(list(num = as.numeric(x),  cls = "Date"))
  if (cls == "POSIXct")  return(list(num = as.numeric(x),  cls = "POSIXct"))
  if (is.ordered(x))     return(list(num = as.integer(x),  cls = "ordered"))
  stop(sprintf(
    "Unsupported class '%s'. Supported: numeric, integer, Date, POSIXct, currency, ordered factor.", cls),
    call. = FALSE)
}

.from_numeric <- function(vals, cls, tz = "UTC") {
  switch(cls,
    Date     = structure(round(vals), class = "Date"),
    POSIXct  = structure(vals, class = c("POSIXct", "POSIXt"), tzone = tz),
    currency = structure(vals, class = "currency"),
    vals
  )
}

.wquantile <- function(x_num, probs, w) {
  ord      <- order(x_num, method = "radix")
  xs       <- x_num[ord]
  ws       <- w[ord]
  rle_vals <- rle(xs)
  xs_u     <- rle_vals$values
  ws_u     <- as.numeric(tapply(ws,
                rep(seq_along(rle_vals$lengths), rle_vals$lengths), sum))
  cw       <- cumsum(ws_u) / sum(ws_u)
  n        <- length(xs_u)
  sapply(probs, function(p) {
    if (p <= cw[1L]) return(xs_u[1L])
    if (p >= cw[n])  return(xs_u[n])
    idx <- findInterval(p, cw)
    lo  <- cw[idx]; hi <- cw[idx + 1L]
    if (hi == lo) return(xs_u[idx])
    xs_u[idx] + (xs_u[idx + 1L] - xs_u[idx]) * (p - lo) / (hi - lo)
  })
}

.wprank <- function(x_num, weights) {
  w_total  <- sum(weights)
  ord      <- order(x_num, method = "radix")
  xs       <- x_num[ord]
  ws       <- weights[ord]
  rle_vals <- rle(xs)
  grp_ids  <- rep(seq_along(rle_vals$lengths), rle_vals$lengths)
  grp_w    <- as.numeric(tapply(ws, grp_ids, sum))
  cum_below <- c(0, cumsum(grp_w)[-length(grp_w)])
  pct_grp  <- (cum_below + grp_w / 2) / w_total * 100
  pct_ord  <- rep(pct_grp, rle_vals$lengths)
  pct      <- numeric(length(x_num))
  pct[ord] <- pct_ord
  pct
}

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

print.currency <- function(x, ...) {
  vals <- paste0("$", formatC(as.numeric(x), format = "f", digits = 2, big.mark = ","))
  if (!is.null(names(x))) vals <- paste(names(x), vals, sep = " = ")
  cat(vals, sep = "\n")
  cat("\n")
  invisible(x)
}

as_currency <- function(x) {
  if (!is.numeric(x) && !is.integer(x))
    stop("'x' must be numeric or integer to convert to currency.", call. = FALSE)
  structure(as.numeric(x), class = "currency")
}

eco_percentile <- function(x,
                           probs    = c(0, .10, .25, .50, .75, .90, 1),
                           weights  = NULL,
                           var_name = "x") {
  .check_na(x, var_name)
  .check_probs(probs)
  .check_weights(weights, length(x), "weights")
  parsed <- .to_numeric(x)
  x_num  <- parsed$num
  cls    <- parsed$cls
  vals <- if (is.null(weights))
    stats::quantile(x_num, probs = probs, type = 7, names = FALSE)
  else
    .wquantile(x_num, probs, weights)
  pct_names <- paste0("P", round(probs * 100, 1))
  pct_names[probs == 0]   <- "Min"
  pct_names[probs == 1]   <- "Max"
  pct_names[probs == 0.5] <- "Median"
  result <- .from_numeric(vals, cls,
    tz = if (cls == "POSIXct") attr(x, "tzone") %||% "UTC" else "UTC")
  names(result) <- pct_names
  result
}

eco_prank <- function(x,
                      weights  = NULL,
                      ties     = c("average", "min", "max", "first", "last", "random"),
                      var_name = "x") {
  ties <- match.arg(ties)
  .check_na(x, var_name)
  .check_weights(weights, length(x), "weights")
  x_num <- .to_numeric(x)$num
  n     <- length(x_num)
  pct <- if (is.null(weights)) {
    r <- rank(x_num, ties.method = ties)
    if (n == 1L) 50 else (r - 1) / (n - 1) * 100
  } else {
    .wprank(x_num, weights)
  }
  round(pct, 4L)
}

eco_percentile_df <- function(df,
                               probs      = c(0, .25, .50, .75, 1),
                               group_vars = NULL,
                               weight_var = NULL,
                               cols       = NULL) {
  stopifnot(is.data.frame(df))
  if (!is.null(weight_var) && !weight_var %in% names(df))
    stop(sprintf("weight_var '%s' not found in data frame.", weight_var), call. = FALSE)
  exclude       <- c(group_vars, weight_var)
  supported_cls <- c("numeric", "integer", "Date", "POSIXct", "currency")
  if (is.null(cols)) {
    cols <- names(df)[sapply(names(df), function(cn) {
      cl <- class(df[[cn]])[1L]
      (cl %in% supported_cls || is.ordered(df[[cn]])) && !(cn %in% exclude)
    })]
  } else {
    bad <- .check_cols(df, cols, supported_cls)
    if (length(bad) > 0L)
      stop(sprintf("Unsupported or missing column(s): %s", paste(bad, collapse = ", ")), call. = FALSE)
  }
  if (length(cols) == 0L)
    stop("No supported columns found to analyse.", call. = FALSE)
  .one_group <- function(sub) {
    w    <- if (!is.null(weight_var)) sub[[weight_var]] else NULL
    rows <- lapply(cols, function(cn) {
      vals <- eco_percentile(sub[[cn]], probs = probs, weights = w, var_name = cn)
      as.data.frame(t(vals), stringsAsFactors = FALSE)
    })
    tbl          <- do.call(rbind, rows)
    tbl$variable <- cols
    tbl[, c("variable", names(tbl)[names(tbl) != "variable"]), drop = FALSE]
  }
  if (is.null(group_vars)) return(.one_group(df))
  grp_list <- split(df, df[, group_vars, drop = FALSE], drop = TRUE)
  results  <- lapply(names(grp_list), function(g) {
    sub   <- grp_list[[g]]
    tbl   <- .one_group(sub)
    gvals <- sub[1L, group_vars, drop = FALSE]
    row.names(gvals) <- NULL
    cbind(gvals[rep(1L, nrow(tbl)), , drop = FALSE], tbl, row.names = NULL)
  })
  do.call(rbind, c(results, list(make.row.names = FALSE)))
}

eco_prank_df <- function(df,
                          cols       = NULL,
                          group_vars = NULL,
                          weight_var = NULL,
                          ties       = "average") {
  stopifnot(is.data.frame(df))
  if (!is.null(weight_var) && !weight_var %in% names(df))
    stop(sprintf("weight_var '%s' not found in data frame.", weight_var), call. = FALSE)
  exclude       <- c(group_vars, weight_var)
  supported_cls <- c("numeric", "integer", "Date", "POSIXct", "currency")
  if (is.null(cols)) {
    cols <- names(df)[sapply(names(df), function(cn) {
      cl <- class(df[[cn]])[1L]
      (cl %in% supported_cls || is.ordered(df[[cn]])) && !(cn %in% exclude)
    })]
  } else {
    bad <- .check_cols(df, cols, supported_cls)
    if (length(bad) > 0L)
      stop(sprintf("Unsupported or missing column(s): %s", paste(bad, collapse = ", ")), call. = FALSE)
  }
  if (length(cols) == 0L)
    stop("No supported columns found to rank.", call. = FALSE)
  .rank_one_group <- function(sub) {
    w <- if (!is.null(weight_var)) sub[[weight_var]] else NULL
    for (cn in cols)
      sub[[paste0(cn, "_prank")]] <- eco_prank(sub[[cn]], weights = w,
                                               ties = ties, var_name = cn)
    sub
  }
  if (is.null(group_vars)) return(.rank_one_group(df))
  orig_order     <- seq_len(nrow(df))
  df$.orig_order <- orig_order
  grp_keys       <- do.call(paste, c(df[, group_vars, drop = FALSE], sep = "\r"))
  df_out         <- do.call(rbind, c(lapply(split(df, grp_keys), .rank_one_group),
                                     list(make.row.names = FALSE)))
  df_out         <- df_out[order(df_out$.orig_order), ]
  df_out$.orig_order <- NULL
  df_out
}

eco_percentile_summary <- function(df,
                                    cols       = NULL,
                                    weight_var = NULL,
                                    digits     = 3L) {
  stopifnot(is.data.frame(df))
  if (!is.null(weight_var) && !weight_var %in% names(df))
    stop(sprintf("weight_var '%s' not found in data frame.", weight_var), call. = FALSE)
  .check_digits(digits)
  supported_cls <- c("numeric", "integer", "currency")
  if (is.null(cols)) {
    cols <- names(df)[sapply(names(df), function(cn) {
      class(df[[cn]])[1L] %in% supported_cls && cn != weight_var
    })]
  } else {
    bad <- .check_cols(df, cols, supported_cls)
    if (length(bad) > 0L)
      stop(sprintf("Unsupported or missing column(s): %s", paste(bad, collapse = ", ")), call. = FALSE)
  }
  if (length(cols) == 0L)
    stop("No supported columns found to summarise.", call. = FALSE)
  w     <- if (!is.null(weight_var)) df[[weight_var]] else NULL
  probs <- c(0, .10, .25, .50, .75, .90, 1)
  rows <- lapply(cols, function(cn) {
    v   <- df[[cn]]
    pct <- eco_percentile(v, probs = probs, weights = w, var_name = cn)
    xn  <- as.numeric(v)
    mn  <- if (is.null(w)) mean(xn) else sum(xn * w) / sum(w)
    sd_ <- if (is.null(w)) {
      stats::sd(xn)
    } else {
      wm    <- sum(xn * w) / sum(w)
      n_eff <- sum(w)^2 / sum(w^2)
      sqrt(sum(w * (xn - wm)^2) / (sum(w) * (1 - 1 / n_eff)))
    }
    data.frame(
      variable = cn,
      Min      = round(pct["Min"],              digits),
      P10      = round(pct["P10"],              digits),
      P25      = round(pct["P25"],              digits),
      Median   = round(pct["Median"],           digits),
      P75      = round(pct["P75"],              digits),
      P90      = round(pct["P90"],              digits),
      Max      = round(pct["Max"],              digits),
      Mean     = round(mn,                      digits),
      SD       = round(sd_,                     digits),
      IQR      = round(pct["P75"] - pct["P25"], digits),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}

eco_gini <- function(x, weights = NULL, var_name = "x") {
  .check_na(x, var_name)
  x_num <- .to_numeric(x)$num
  .check_na(x_num, var_name)
  n <- length(x_num)
  if (n < 2L)
    stop(sprintf("'%s' must have at least 2 observations to compute Gini.", var_name), call. = FALSE)
  if (is.null(weights)) weights <- rep(1, n)
  .check_weights(weights, n, "weights")
  ord    <- order(x_num)
  x_s    <- x_num[ord];  w_s <- weights[ord]
  w_cum  <- cumsum(w_s)  / sum(w_s)
  xw_cum <- cumsum(x_s * w_s) / sum(x_s * w_s)
  B      <- sum((w_cum[-1] - w_cum[-n]) * (xw_cum[-1] + xw_cum[-n]) / 2)
  round(1 - 2 * B, 4)
}
