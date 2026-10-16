## Goodman-level check of log-linear families for the average (row-conditional) table, before building any into ssEI.
## Every family: log P(c | r) = col_c + sum_k theta_k L_k[r, c] (normalised within rows), or a nonlinear version (estimated scores).
## For each family:  margins   = fitted from the areas' margins (weighted least squares of area column shares on row shares, as Goodman)
##                   truth     = best member fitted to the true table (the representation error: the price of the restriction)
## Errors: national table (row-weighted), and area error after raking the fitted table to each area's margins.
## Families:
##   QI          loyalty per row + one column popularity (quasi-independence)
##   dist-k      QI + a distance term -beta * (x_r - x_c)^2 with stated party positions (one beta per dimension). Rows/columns without a
##               position (spoilt, Other) get no distance term. Positions are rough placements written here, not survey estimates.
##   unfold-k    as dist-k with the k-dimensional positions estimated (unfolding; bilinear, so not a linear family)
##   QS          quasi-symmetry: column popularity + own cell per row + one parameter per unordered pair of parties (lambda_rc = lambda_cr)
##   RE          census: row effects, lambda_rc = mu_r v_c with stated ordered column scores
##   RE+oth      as RE with each row's share of the unordered 'Oth' column free
##   RC1         census: lambda_rc = mu_r v_c with the column scores estimated
## usage: Rscript experiments/structure/families.R [scotland | nz | census | all]
suppressMessages({ source("experiments/data_prep.R"); library(quadprog) })
which_data <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "all"
if (which_data == "none") which_data <- "skip"

setup <- function(kc) {
  J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]; rm <- apply(kc, 1:2, sum); cm <- apply(kc, c(1, 3), sum); w <- rowSums(rm)
  list(kc = kc, J = J, R = R, C = C, rm = rm, cm = cm, w = w, xs = rm / w, ys = cm / w, N = sum(kc), nr = colSums(rm), Bt = apply(kc, 2:3, sum) / colSums(rm))
}
softrow <- function(eta) { eta <- eta - apply(eta, 1, max); P <- exp(eta); P / rowSums(P) }
## a family = list(np, eta = function(theta) R x C matrix of log-linear predictors); column effects included in theta
fit_family <- function(d, fam, starts = 8, seed = 1) {
  Pf <- function(th) softrow(fam$eta(th))
  ob_m <- function(th) sum(d$w * (d$ys - d$xs %*% Pf(th))^2) / d$N
  ob_t <- function(th) sum(d$nr * (d$Bt - Pf(th))^2) / d$N
  ## gradients: d objective / d P in closed form, through the softmax, then a finite-difference Jacobian of eta (no data involved, so cheap)
  gP_m <- function(P) -2 * crossprod(d$xs, d$w * (d$ys - d$xs %*% P)) / d$N
  gP_t <- function(P) -2 * d$nr * (d$Bt - P) / d$N
  grad <- function(th, gP) { P <- Pf(th); g <- gP(P); ge <- P * (g - rowSums(g * P)); h <- 1e-6
    vapply(seq_along(th), function(k) { e <- th; e[k] <- e[k] + h; f <- th; f[k] <- f[k] - h; sum(ge * (fam$eta(e) - fam$eta(f))) / (2 * h) }, 0) }
  best <- function(ob) { gP <- if (identical(ob, ob_m)) gP_m else gP_t; set.seed(seed); b <- NULL
    for (s in seq_len(starts)) { th0 <- if (s == 1 && !is.null(fam$start)) fam$start else rnorm(fam$np, 0, if (s == 1) 0.1 else 1)
      o <- tryCatch(optim(th0, ob, function(th) grad(th, gP), method = "BFGS", control = list(maxit = 3000, reltol = 1e-11)), error = function(e) NULL)
      if (!is.null(o) && (is.null(b) || o$value < b$value)) b <- o }
    structure(Pf(b$par), par = b$par) }
  list(margins = best(ob_m), truth = best(ob_t))
}
goodman <- function(d) {   # unrestricted, constrained to the simplex
  R <- d$R; C <- d$C; Dq <- matrix(0, R * C, R * C); dv <- numeric(R * C)
  for (j in 1:d$J) { Dq <- Dq + d$w[j] * kronecker(tcrossprod(d$xs[j, ]), diag(C)); dv <- dv + d$w[j] * kronecker(d$xs[j, ], d$ys[j, ]) }
  A <- cbind(t(kronecker(diag(R), matrix(1, 1, C))), diag(R * C))
  matrix(solve.QP(2 * Dq / d$N + diag(1e-9, R * C), 2 * dv / d$N, A, c(rep(1, R), rep(0, R * C)), meq = R)$solution, R, C, byrow = TRUE)
}
nat <- function(d, P) 50 * sum(d$nr * abs(P - d$Bt)) / d$N
rake <- function(d, P) { e <- array(0, dim(d$kc)); for (j in 1:d$J) { X <- P * d$rm[j, ] + 1e-9; X[d$rm[j, ] == 0, ] <- 0
  for (i in 1:1500) { rs <- rowSums(X); X <- X * ifelse(rs > 0, d$rm[j, ] / pmax(rs, 1e-300), 0); cs <- colSums(X); X <- sweep(X, 2, ifelse(cs > 0, d$cm[j, ] / pmax(cs, 1e-300), 0), "*")
    if (i %% 25 == 0 && max(abs(rowSums(X) - d$rm[j, ])) < 1e-6) break }; e[j, , ] <- X }; 50 * sum(abs(e - d$kc)) / d$N }

## family builders (square tables, rows = columns = the same categories)
colfx <- function(th, R, C) matrix(c(0, th[seq_len(C - 1)]), R, C, byrow = TRUE)
fam_qi <- function(R, C, extra = NULL) {   # extra = list(np, eta_add(th_extra)) on top of QI
  ne <- if (is.null(extra)) 0 else extra$np
  list(np = C - 1 + R + ne, eta = function(th) { e <- colfx(th, R, C); diag(e) <- diag(e) + th[C - 1 + seq_len(R)]
    if (ne) e <- e + extra$eta_add(th[C - 1 + R + seq_len(ne)]); e }, start = c(rep(0, C - 1), rep(2, R), rep(0.1, ne)))
}
dist_term <- function(X) {   # X = rows x dims positions (NA = no position); one beta per dimension, beta = exp(.) > 0
  K <- ncol(X); D <- lapply(seq_len(K), function(k) { M <- outer(X[, k], X[, k], function(a, b) (a - b)^2); M[is.na(M)] <- 0; M / max(M) })
  list(np = K, eta_add = function(b) Reduce(`+`, lapply(seq_len(K), function(k) -exp(b[k]) * D[[k]])))
}
unfold_term <- function(has, K) {   # positions estimated for the rows with has = TRUE; scale carried by the positions
  idx <- which(has); n <- length(idx)
  list(np = n * K, eta_add = function(p) { X <- matrix(p, n, K); M <- matrix(0, length(has), length(has))
    M[idx, idx] <- -as.matrix(dist(X))^2; M })
}
fam_qs <- function(R) {
  pairs <- which(upper.tri(diag(R)), arr.ind = TRUE)
  list(np = R - 1 + R + nrow(pairs), eta = function(th) { e <- colfx(th, R, R); diag(e) <- diag(e) + th[R - 1 + seq_len(R)]
    l <- th[2 * R - 1 + seq_len(nrow(pairs))]; e[pairs] <- e[pairs] + l; e[pairs[, 2:1]] <- e[pairs[, 2:1]] + l; e },
    start = c(rep(0, R - 1), rep(2, R), rep(0, nrow(pairs))))
}
## census: rows unordered, columns (partly) ordered
fam_re <- function(R, C, v, oth = NULL) {   # lambda_rc = mu_r v_c, mu_1 = 0; oth = column index whose row-specific level is free
  no <- if (is.null(oth)) 0 else R - 1
  list(np = C - 1 + R - 1 + no, eta = function(th) { e <- colfx(th, R, C) + outer(c(0, th[C - 1 + seq_len(R - 1)]), v)
    if (no) e[, oth] <- e[, oth] + c(0, th[C - 1 + R - 1 + seq_len(no)]); e })
}
fam_rc1 <- function(R, C) list(np = C - 1 + R - 1 + C, eta = function(th) colfx(th, R, C) + outer(c(0, th[C - 1 + seq_len(R - 1)]), th[C + R - 2 + seq_len(C)]))
fam_indep <- function(R, C) list(np = C - 1, eta = function(th) colfx(th, R, C))

run <- function(label, d, fams) {
  Pg <- goodman(d)
  out <- data.frame(family = "Goodman (unrestricted)", params = d$R * (d$C - 1), margins_nat = nat(d, Pg), margins_area = rake(d, Pg), truth_nat = 0, truth_area = rake(d, d$Bt))
  for (nm in names(fams)) { f <- fit_family(d, fams[[nm]])
    out <- rbind(out, data.frame(family = nm, params = fams[[nm]]$np, margins_nat = nat(d, f$margins), margins_area = rake(d, f$margins), truth_nat = nat(d, f$truth), truth_area = rake(d, f$truth)))
    cat(sprintf("  %-28s done\n", nm)) }
  cat("\n", label, "\n"); out[, 3:6] <- round(out[, 3:6], 2); print(out, row.names = FALSE); invisible(out)
}

res <- list()
if (which_data %in% c("scotland", "all")) {
  s <- prep_scotland(7, 73); d <- setup(s$kc); R <- 7
  ## Con, Lab, LD, SNP, spoilt, Other, Green. Left-right on 0-10, independence (support) on 0-10; rough placements
  lr <- c(7.5, 4.5, 5.0, 4.0, NA, NA, 2.5); ind <- c(0, 1, 2, 10, NA, NA, 7)
  has <- !is.na(lr)
  res$scotland <- run("Scotland 2007, 7 x 7 (national error | raked area error; fitted from margins, and best member on truth)", d, list(
    "QI" = fam_qi(R, R),
    "QI + left-right (stated)" = fam_qi(R, R, dist_term(cbind(lr))),
    "QI + independence (stated)" = fam_qi(R, R, dist_term(cbind(ind))),
    "QI + both (stated)" = fam_qi(R, R, dist_term(cbind(lr, ind))),
    "QI + 1D unfolding" = fam_qi(R, R, unfold_term(has, 1)),
    "QI + 2D unfolding" = fam_qi(R, R, unfold_term(has, 2)),
    "QS" = fam_qs(R)))
}
if (which_data %in% c("nz", "all")) {
  z <- prep_nz(2020, 6, "candidate", 1e6); d <- setup(z$kc); R <- 6
  lr <- c(4, 7, 2, 5.5, 8.5, NA); has <- !is.na(lr)   # Labour, National, Green, NZFirst, ACT, Other; rough placements
  res$nz <- run("NZ 2020, 6 x 6", d, list(
    "QI" = fam_qi(R, R),
    "QI + left-right (stated)" = fam_qi(R, R, dist_term(cbind(lr))),
    "QI + 1D unfolding" = fam_qi(R, R, unfold_term(has, 1)),
    "QI + 2D unfolding" = fam_qi(R, R, unfold_term(has, 2)),
    "QS" = fam_qs(R)))
}
if (which_data %in% c("census", "all")) {
  k <- readRDS("/home/claude/data_census/eth6_ed7.rds"); k <- array(as.numeric(k), dim(k)); d <- setup(k); R <- d$R; C <- d$C
  ## rows Asian, Black, Mixed, White, Other; columns None, L1, L2, L3, L4+, Oth (apprenticeship / other: unordered, score 0 = the middle)
  v <- c(-2, -1, 0, 1, 2, 0)
  res$census <- run("Census 2021, ethnicity x qualification, 309 LTLAs", d, list(
    "independence" = fam_indep(R, C),
    "row effects (stated scores)" = fam_re(R, C, v),
    "row effects + Oth free" = fam_re(R, C, v, oth = 6),
    "RC(1), column scores estimated" = fam_rc1(R, C)))
}
saveRDS(res, file.path(tempdir(), "families.rds"))
