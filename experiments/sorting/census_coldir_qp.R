## Same diagnostics with a stable base: Goodman constrained to the simplex (quadratic programme), then the column-side correction.
source_lines <- readLines("experiments/sorting/census_kappa.R"); eval(parse(text = source_lines[1:grep("^## 1\\.", source_lines) - 1]))
suppressMessages(library(quadprog))
Dq <- matrix(0, R * C, R * C); dv <- numeric(R * C)
for (j in 1:J) { Dq <- Dq + w[j] * kronecker(tcrossprod(xs[j, ]), diag(C)); dv <- dv + w[j] * kronecker(xs[j, ], ys[j, ]) }
Dq <- 2 * Dq / N + diag(1e-9, R * C); dv <- 2 * dv / N; A <- cbind(t(kronecker(diag(R), matrix(1, 1, C))), diag(R * C))
Bc <- matrix(solve.QP(Dq, dv, A, c(rep(1, R), rep(0, R * C)), meq = R)$solution, R, C, byrow = TRUE); dimnames(Bc) <- dimnames(Bt)
cat(sprintf("Constrained Goodman national error %.2f (unconstrained %.2f, independence %.2f)\n", nat_err(Bc), nat_err(Bh), nat_err(ind)))
print(round(100 * Bc, 1))
D <- Bc - Bt; sv <- svd(sqrt(p) * D); cat(sprintf("bias singular values %s | rank-1 share %.2f\n", paste(round(sv$d, 3), collapse = " "), sv$d[1]^2 / sum(sv$d^2)))
a1 <- sv$u[, 1] / sqrt(p); d1 <- sv$v[, 1]; if (d1[5] < 0) { a1 <- -a1; d1 <- -d1 }
cat("rank-1 row pattern:", paste(rn, round(a1 / max(abs(a1)), 2)), "\nrank-1 column pattern:", paste(cn, round(d1, 2)), "\n")
D1 <- (sv$d[1] * sv$u[, 1] %o% sv$v[, 1]) / sqrt(p); cat(sprintf("error after removing the oracle rank-1 part: %.2f\n", nat_err(Bc - D1)))
cat("bias (Goodman constrained - truth), percentage points:\n"); print(round(100 * D, 1))

## Column-side correction from the constrained base: B = Bc - a delta0', a free per row (optionally sign-restricted), bounds + totals.
lo <- hi <- matrix(0, R, C)
for (j in 1:J) for (r in 1:R) for (c in 1:C) { lo[r, c] <- lo[r, c] + max(0, cm[j, c] - (w[j] - rm[j, r])); hi[r, c] <- hi[r, c] + min(rm[j, r], cm[j, c]) }
lo <- lo / nr; hi <- hi / nr
res <- ys - xs %*% Bc; d_res <- eigen(cov.wt(res, wt = w / N)$cov)$vectors[, 1]
norm_d <- function(d) { d <- d - mean(d); d <- d / sqrt(sum(d^2)); if (d[5] < 0) -d else d }
dirs <- list(residual_pc1 = norm_d(d_res), substantive = norm_d(c(-2, -1, -0.5, 0.5, 3, 0)), oracle = d1)
cat("\nColumn directions (cosine with the oracle bias pattern):", paste(names(dirs), round(sapply(dirs, function(d) sum(d * d1)), 2)), "\n")
lp_range <- function(d, r, sgn = NULL) {   # min and max of a_r subject to bounds on every cell, p'a = 0, optional signs on a
  rng <- c(NA, NA)
  for (s in c(1, -1)) { obj <- numeric(R); obj[r] <- s                       # minimise s * a_r via a QP with a tiny quadratic term
    A <- matrix(p, R, 1); b <- 0
    for (rr in 1:R) for (c in 1:C) if (abs(d[c]) > 1e-12) { e <- numeric(R); e[rr] <- -d[c]; A <- cbind(A, e, -e); b <- c(b, lo[rr, c] - Bc[rr, c], Bc[rr, c] - hi[rr, c]) }
    if (!is.null(sgn)) for (rr in 1:R) if (sgn[rr] != 0) { e <- numeric(R); e[rr] <- sgn[rr]; A <- cbind(A, e); b <- c(b, 0) }
    sol <- tryCatch(solve.QP(diag(1e-6, R), -obj, A, b, meq = 1)$solution, error = function(e) NULL)
    if (is.null(sol)) return(c(NA, NA)); rng[if (s == 1) 1 else 2] <- sol[r] }
  rng }
for (dn in names(dirs)) { d <- dirs[[dn]]; a_or <- as.vector((Bc - Bt) %*% d)
  cat(sprintf("\n%s: error after the oracle a along this direction %.2f\n", dn, nat_err(Bc - a_or %o% d)))
  for (lab in c("no sign restriction", "signs: Asian -, Black +, Mixed +, Other +")) {
    sg <- if (lab == "no sign restriction") NULL else c(-1, 1, 1, 0, 1)
    rg <- t(sapply(1:R, function(r) lp_range(d, r, sg)))
    if (anyNA(rg)) { cat("  ", lab, ": infeasible\n"); next }
    cat("  ", lab, "| L4+ %: truth, constrained Goodman, interval\n")
    for (r in 1:R) cat(sprintf("     %-6s %5.1f %6.1f  [%5.1f, %5.1f]%s\n", rn[r], 100 * Bt[r, 5], 100 * Bc[r, 5], 100 * min(Bc[r, 5] - rg[r, ] * d[5]), 100 * max(Bc[r, 5] - rg[r, ] * d[5]),
      if (Bt[r, 5] >= min(Bc[r, 5] - rg[r, ] * d[5]) - 1e-9 & Bt[r, 5] <= max(Bc[r, 5] - rg[r, ] * d[5]) + 1e-9) "" else "  truth outside")) } }
