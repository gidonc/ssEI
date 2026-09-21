#' @export
build_lambda_rotation <- function(fit, row_margins, col_margins, V_ilr_model,
                                  stage = c("per-area", "full"),
                                  draw = 1, delta = 1e-3,
                                  neutral_logit_mode = 2) {

  stage <- match.arg(stage)
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("numDeriv required")

  rm_mat <- as.matrix(row_margins); storage.mode(rm_mat) <- "double"
  cm_mat <- as.matrix(col_margins); storage.mode(cm_mat) <- "double"
  R <- ncol(rm_mat); C <- ncol(cm_mat); n_areas <- nrow(rm_mat)
  free_R <- rowSums(rm_mat > 0); free_C <- rowSums(cm_mat > 0)
  full_areas <- which(free_R == R & free_C == C)
  if (!length(full_areas)) stop("no areas with all rows and columns non-zero")

  ## --- index map, mirroring Stan's transformed data ----------------------
  pcf <- integer(n_areas); n_param <- 0
  for (j in seq_len(n_areas)) {
    pcf[j] <- n_param
    n_param <- n_param + max(0, (free_R[j] - 1) * (free_C[j] - 1))
  }
  map <- do.call(rbind, lapply(seq_len(n_areas), function(j) {
    if (free_R[j] < 2 || free_C[j] < 2) return(NULL)
    d <- expand.grid(c = 1:(free_C[j] - 1), r = 1:(free_R[j] - 1))
    data.frame(idx = pcf[j] + (d$r - 1) * (free_C[j] - 1) + d$c,
               area = j, r = d$r, c = d$c)
  }))
  map <- map[order(map$idx), ]

  ## --- draws --------------------------------------------------------------
  lam_draws  <- as.matrix(fit, pars = "lambda")
  comp_draws <- as.matrix(fit, pars = "composition_arr")
  cv_draws   <- as.matrix(fit, pars = "cell_values")
  lr_draws   <- as.matrix(fit, pars = "lambda_raw")
  if (ncol(lr_draws) != n_param)
    stop("n_param mismatch: reconstructed ", n_param, ", fit has ", ncol(lr_draws))
  if (draw > nrow(lam_draws)) stop("draw index exceeds available draws")

  pick <- function(draws, j, d, nr, nc, nm) {
    idx <- paste0(nm, "[", j, ",", rep(1:nr, each = nc), ",", rep(1:nc, times = nr), "]")
    matrix(draws[d, idx], nrow = nr, ncol = nc, byrow = TRUE)
  }

  ## --- allocation, mirroring ss_assign_cvals_cpanchor_lp ------------------
  hf <- function(x, d) { z <- x / d
  d * ifelse(z > 30, z, ifelse(z < -30, exp(z), log1p(exp(z)))) }
  hm <- function(x, d) { z <- -x / d; m <- max(z); -d * (m + log(sum(exp(z - m)))) }

  alloc_area <- function(rm_j, cm_j, lam, comp_j, dfloor, dmin, nl_mode) {
    rm_j <- as.numeric(rm_j); cm_j <- as.numeric(cm_j)
    ar <- which(rm_j > 0); ac <- which(cm_j > 0)
    fR <- length(ar); fC <- length(ac)
    sr <- rm_j[ar]; sc <- cm_j[ac]; rt <- sum(sr)
    tmp <- matrix(0, fR, fC)
    for (r in seq_len(fR - 1)) {
      comp_row <- comp_j[ar[r], ac]
      for (cc in seq_len(fC - 1)) {
        lower <- hf(sr[r] - sum(sc[(cc + 1):fC]), dfloor)
        upper <- hm(c(sc[cc], sr[r]), dmin)
        width <- hf(upper - lower, dmin)
        if (nl_mode == 2) {
          share <- comp_row[cc] / max(sum(comp_row[cc:fC]), 1e-10)
          neutral <- log(share / (1 - share))
        } else neutral <- -log(fC - cc)
        p <- 1 / (1 + exp(-(neutral + lam[r, cc])))
        tmp[r, cc] <- lower + p * width
        sc[cc] <- max(sc[cc] - tmp[r, cc], 0)
        sr[r]  <- max(sr[r]  - tmp[r, cc], 0)
        rt     <- max(rt     - tmp[r, cc], 0)
      }
      tmp[r, fC] <- sr[r]
      rt <- max(rt - tmp[r, fC], 0)
      sc[fC] <- max(sc[fC] - tmp[r, fC], 0)
      sr[r] <- 0
    }
    for (cc in seq_len(fC - 1)) { tmp[fR, cc] <- sc[cc]; rt <- max(rt - sc[cc], 0) }
    tmp[fR, fC] <- rt
    out <- matrix(0, length(rm_j), length(cm_j)); out[ar, ac] <- tmp
    out
  }

  ## --- per-area directions ------------------------------------------------
  rm_idx <- as.vector(t(matrix(seq_len(R * C), nrow = R, ncol = C)))

  area_dirs <- function(jj) {
    comp_d <- pick(comp_draws, jj, draw, R, C, "composition_arr")
    cv_d   <- pick(cv_draws,   jj, draw, R, C, "cell_values")
    lam_d  <- pick(lam_draws,  jj, draw, R - 1, C - 1, "lambda")
    f <- function(lv) as.vector(alloc_area(rm_mat[jj, ], cm_mat[jj, ],
                                           matrix(lv, R - 1, C - 1), comp_d,
                                           delta, delta, neutral_logit_mode))
    J  <- numDeriv::jacobian(f, as.vector(lam_d), method = "Richardson")
    Jr <- J[rm_idx, , drop = FALSE]
    p <- as.vector(t(comp_d)); vol <- sum(cv_d); qrf <- qr(Jr)
    out <- sapply(3:6, function(s) {
      dcells <- vol * (p * V_ilr_model[, s] - p * sum(p * V_ilr_model[, s]))
      v <- qr.coef(qrf, dcells)
      c(v / sqrt(sum(v^2)),
        resid = sqrt(sum((dcells - Jr %*% v)^2)) / sqrt(sum(dcells^2)))
    })
    colnames(out) <- paste0("s", 3:6)
    out
  }

  res_all <- lapply(full_areas, area_dirs); names(res_all) <- full_areas

  ## --- assemble -----------------------------------------------------------
  ordmap <- c(1, 3, 2, 4)
  Q <- diag(n_param)
  for (jj in full_areas) {
    x  <- res_all[[as.character(jj)]]
    d4 <- x[ordmap, "s4"]; d5 <- x[ordmap, "s5"]
    d4 <- d4 / sqrt(sum(d4^2))
    d5 <- d5 - sum(d5 * d4) * d4
    n5 <- sqrt(sum(d5^2))
    if (n5 < 1e-8) stop("s4 and s5 nearly parallel in area ", jj)
    d5 <- d5 / n5
    B <- qr.Q(qr(cbind(d4, d5, diag(4))))[, 1:4]
    if (sum(B[, 1] * d4) < 0) B[, 1] <- -B[, 1]
    if (sum(B[, 2] * d5) < 0) B[, 2] <- -B[, 2]
    idx <- map$idx[map$area == jj]
    stopifnot(length(idx) == 4, !is.unsorted(idx))
    Q[idx, idx] <- B
  }

  if (stage == "full") {
    nA <- length(full_areas)
    H  <- cbind(rep(1, nA) / sqrt(nA), make_helmert_basis(nA))
    S  <- diag(n_param)
    for (k in 1:2) {
      idx_k <- vapply(full_areas, function(jj) map$idx[map$area == jj][k], numeric(1))
      S[idx_k, idx_k] <- H
    }
    Q <- Q %*% S
  }

  e <- max(abs(crossprod(Q) - diag(n_param)))
  if (e > 1e-9) stop("rotation not orthonormal, max dev ", signif(e, 3))

  attr(Q, "stage") <- stage
  attr(Q, "neutral_logit_mode") <- neutral_logit_mode
  attr(Q, "draw") <- draw
  attr(Q, "resid") <- t(sapply(res_all, function(x) x["resid", ]))
  attr(Q, "full_areas") <- full_areas
  Q
}

check_orth <- function(Q) {
  e <- max(abs(crossprod(Q) - diag(ncol(Q))))
  if (e > 1e-9) stop("rotation not orthonormal, max dev ", signif(e, 3))
  invisible(TRUE)
}
