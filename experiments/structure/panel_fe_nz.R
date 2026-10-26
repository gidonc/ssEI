## Panel (fixed-effects) ecological regression on NZ election pairs with unchanged electorate boundaries (2002-05, 2008-11, 2014-17).
## Idea: an area's column propensity (the sorting term) is roughly stable between two elections, so differencing removes it:
##   y_jt = x_jt B_t + a_j,   B_t = B + 1 d_t'  (a national swing common to all rows)   =>   dy_j = dx_j B + (d_2 - d_1).
## Differences identify B only up to a common column shift (dx sums to 0), so the shift for each year is set by the national margins:
## the differenced estimate (clipped to positive, rows renormalised) is used as a seed and raked to that year's NATIONAL row and column totals.
## Compared with: single-year constrained Goodman; pooled Goodman (both years, one table); Goodman raked to the national totals (same last step).
## Errors: national table (row-weighted), and area error after raking the national table to each area. Truth known for both years.
## Diagnostic: stability of the true area column propensities, a_jt = y_jt - x_jt B_t (B_t the true table), between the two years.
## usage: Rscript experiments/structure/panel_fe_nz.R
suppressMessages({ source("experiments/data_prep.R"); library(quadprog) })
nz <- function(y) { d <- prep_nz(y, 6, "candidate", 1e6); env <- new.env(); utils::data(list = sprintf("ei_NZ_%d", y), package = "ei.Datasets", envir = env)
  nm <- as.character(get(sprintf("ei_NZ_%d", y), env)$District); dimnames(d$kc) <- list(nm, d$row_names, d$col_names); d$kc }
ipf <- function(X, rt, ct, it = 3000) { for (i in 1:it) { rs <- rowSums(X); X <- X * ifelse(rs > 0, rt / pmax(rs, 1e-300), 0); cs <- colSums(X)
  X <- sweep(X, 2, ifelse(cs > 0, ct / pmax(cs, 1e-300), 0), "*"); if (i %% 25 == 0 && max(abs(rowSums(X) - rt)) < 1e-7) break }; X }
nat_rake <- function(P, kc) { X <- ipf(pmax(P, 1e-4) * colSums(apply(kc, 1:2, sum)), colSums(apply(kc, 1:2, sum)), colSums(apply(kc, c(1, 3), sum))); X / rowSums(X) }
err_nat <- function(P, kc) { kt <- apply(kc, 2:3, sum); 50 * sum(abs(P * rowSums(kt) - kt)) / sum(kt) }
err_area <- function(P, kc) { e <- 0; for (j in 1:dim(kc)[1]) { rt <- rowSums(kc[j, , ]); e <- e + sum(abs(ipf(pmax(P, 1e-6) * rt, rt, colSums(kc[j, , ])) - kc[j, , ])) }; 50 * e / sum(kc) }
goodman_qp <- function(X, Y, w) {   # constrained (rows on the simplex), weighted least squares of Y on X
  R <- ncol(X); C <- ncol(Y); Dq <- matrix(0, R * C, R * C); dv <- numeric(R * C)
  for (j in seq_len(nrow(X))) { Dq <- Dq + w[j] * kronecker(tcrossprod(X[j, ]), diag(C)); dv <- dv + w[j] * kronecker(X[j, ], Y[j, ]) }
  A <- cbind(t(kronecker(diag(R), matrix(1, 1, C))), diag(R * C))
  matrix(solve.QP(2 * Dq / sum(w) + diag(1e-9, R * C), 2 * dv / sum(w), A, c(rep(1, R), rep(0, R * C)), meq = R)$solution, R, C, byrow = TRUE)
}
fe_diff <- function(X1, Y1, X2, Y2, w) {   # weighted LS of dY on dX with an intercept; B only up to a column shift, fixed later by the national totals
  dX <- X2 - X1; dY <- Y2 - Y1; Z <- cbind(1, dX[, -ncol(dX)])   # dX rows sum to 0: drop one row category (its coefficient row absorbed by the shift)
  coef <- solve(crossprod(Z * sqrt(w)) + diag(1e-9, ncol(Z)), crossprod(Z * w, dY))
  rbind(coef[-1, , drop = FALSE], 0) + matrix(1 / ncol(dY), ncol(dX), ncol(dY))   # rows relative to the dropped category; any common shift is set by raking
}
shift_to_nat <- function(Bh, k) {   # B_t = Bh + 1 d', d set by the national identity xbar_t B_t = ybar_t; then onto the simplex
  rt <- colSums(apply(k, 1:2, sum)); ct <- colSums(apply(k, c(1, 3), sum)); xb <- rt / sum(rt); yb <- ct / sum(ct)
  B <- Bh + matrix(yb - as.vector(xb %*% Bh), nrow(Bh), ncol(Bh), byrow = TRUE); nat_rake(B, k) }
fe_qp <- function(X1, Y1, X2, Y2, w, k1, k2) {   # constrained: B1 on the simplex, B2 = B1 + 1 d', B2 >= 0, both national identities exact
  R <- ncol(X1); C <- ncol(Y1); n <- R * C + C; Dq <- matrix(0, n, n); dv <- numeric(n); dX <- X2 - X1; dY <- Y2 - Y1
  for (j in seq_len(nrow(dX))) { A <- cbind(kronecker(t(dX[j, ]), diag(C)), diag(C)); Dq <- Dq + w[j] * crossprod(A); dv <- dv + w[j] * as.vector(crossprod(A, dY[j, ])) }
  xb <- function(k) { r <- colSums(apply(k, 1:2, sum)); r / sum(r) }; yb <- function(k) { r <- colSums(apply(k, c(1, 3), sum)); r / sum(r) }
  Eq <- rbind(cbind(kronecker(diag(R), matrix(1, 1, C)), matrix(0, R, C)),          # rows of B1 sum to 1
              c(rep(0, R * C), rep(1, C)),                                              # d sums to 0
              cbind(kronecker(t(xb(k1)), diag(C)), matrix(0, C, C))[-C, ],             # national identity, year 1
              cbind(kronecker(t(xb(k2)), diag(C)), diag(C))[-C, ])                      # national identity, year 2
  eqb <- c(rep(1, R), 0, yb(k1)[-C], yb(k2)[-C])
  Ge <- rbind(cbind(diag(R * C), matrix(0, R * C, C)), cbind(diag(R * C), kronecker(rep(1, R), diag(C))))   # B1 >= 0, B1 + 1 d' >= 0
  sol <- solve.QP(2 * Dq / sum(w) + diag(1e-7, n), 2 * dv / sum(w), t(rbind(Eq, Ge)), c(eqb, rep(0, 2 * R * C)), meq = nrow(Eq))$solution
  B1 <- matrix(sol[1:(R * C)], R, C, byrow = TRUE); d <- sol[R * C + 1:C]; list(B1, B1 + matrix(d, R, C, byrow = TRUE))
}
pairs <- list(c(2002, 2005), c(2008, 2011), c(2014, 2017))
out <- NULL
for (p in pairs) {
  k1 <- nz(p[1]); k2 <- nz(p[2]); common <- intersect(dimnames(k1)[[1]], dimnames(k2)[[1]])
  k1 <- k1[common, , ]; k2 <- k2[common, , ]
  cat(sprintf("\n%d-%d: %d electorates matched (of %d / %d)\n", p[1], p[2], length(common), dim(nz(p[1]))[1], dim(nz(p[2]))[1]))
  sh <- function(k) { rm <- apply(k, 1:2, sum); cm <- apply(k, c(1, 3), sum); n <- rowSums(rm); list(X = rm / n, Y = cm / n, n = n) }
  s1 <- sh(k1); s2 <- sh(k2); w <- (s1$n + s2$n) / 2
  B1 <- apply(k1, 2:3, sum); B1 <- B1 / rowSums(B1); B2 <- apply(k2, 2:3, sum); B2 <- B2 / rowSums(B2)
  ## stability of the true area propensities
  a1 <- s1$Y - s1$X %*% B1; a2 <- s2$Y - s2$X %*% B2
  cat(sprintf("  true area column propensities: correlation between years %.2f (by party: %s)\n", cor(as.vector(a1), as.vector(a2)),
              paste(sprintf("%s %.2f", dimnames(k1)[[3]], sapply(1:6, function(c) cor(a1[, c], a2[, c]))), collapse = ", ")))
  cat(sprintf("  composition change: sd of within-area change %.3f against between-area sd %.3f (mean over rows)\n", mean(apply(s2$X - s1$X, 2, sd)), mean(apply(s1$X, 2, sd))))
  G1 <- goodman_qp(s1$X, s1$Y, s1$n); G2 <- goodman_qp(s2$X, s2$Y, s2$n)
  Gp <- goodman_qp(rbind(s1$X, s2$X), rbind(s1$Y, s2$Y), c(s1$n, s2$n))
  Bfe <- fe_diff(s1$X, s1$Y, s2$X, s2$Y, w)
  est <- list(`Goodman (own year)` = list(G1, G2), `Goodman, raked to national` = list(nat_rake(G1, k1), nat_rake(G2, k2)),
              `pooled Goodman, raked to national` = list(nat_rake(Gp, k1), nat_rake(Gp, k2)),
              `fixed effects, linear (shift from national identity)` = list(shift_to_nat(Bfe, k1), shift_to_nat(Bfe, k2)),
              `fixed effects, constrained QP` = fe_qp(s1$X, s1$Y, s2$X, s2$Y, w, k1, k2),
              `other year's TRUE table, raked to national` = list(nat_rake(B2, k1), nat_rake(B1, k2)))
  for (nm in names(est)) for (i in 1:2) { k <- list(k1, k2)[[i]]
    out <- rbind(out, data.frame(pair = paste(p, collapse = "-"), year = p[i], method = nm, nat = err_nat(est[[nm]][[i]], k), area = err_area(est[[nm]][[i]], k))) }
  cat("  FE (QP) diagonal (loyalty) vs truth, second year:", paste(sprintf("%.2f/%.2f", diag(est[[5]][[2]]), diag(B2)), collapse = " "), "\n")
}
out$method <- factor(out$method, unique(out$method))
cat("\nNational error by election:\n"); print(round(xtabs(nat ~ method + year, out), 2))
cat("\nArea error (national table raked to each area):\n"); print(round(xtabs(area ~ method + year, out), 2))
cat("\nMeans:\n"); mn <- aggregate(cbind(nat, area) ~ method, out, mean); mn[, -1] <- round(mn[, -1], 2); print(mn, row.names = FALSE)
