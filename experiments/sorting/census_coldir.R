## Column-side-only correction: B = B_hat - a delta0', with delta0 the column direction of the sorting shift (from the Goodman residuals
## or substantive) and a free per row. Constraints: national column totals (p' a = 0) and the aggregated Duncan-Davis bounds.
## Gives each group an interval along delta0; column contrasts orthogonal to delta0 are taken from Goodman.
source_lines <- readLines("experiments/sorting/census_kappa.R"); eval(parse(text = source_lines[1:grep("^## 1\\.", source_lines) - 1]))
suppressMessages(library(quadprog))
res <- ys - xs %*% Bh; d_res <- eigen(cov.wt(res, wt = w / N)$cov)$vectors[, 1]
norm_d <- function(d) { d <- d - mean(d); d <- d / sqrt(sum(d^2)); if (d[5] < 0) -d else d }
dirs <- list(residual_pc1 = norm_d(d_res), substantive = norm_d(c(-2, -1, -0.5, 0.5, 3, 0)))
lo <- hi <- matrix(0, R, C)
for (j in 1:J) for (r in 1:R) for (c in 1:C) { lo[r, c] <- lo[r, c] + max(0, cm[j, c] - (w[j] - rm[j, r])); hi[r, c] <- hi[r, c] + min(rm[j, r], cm[j, c]) }
lo <- lo / nr; hi <- hi / nr
## what is left of Goodman's error once the best a (oracle) is removed along delta0: how much of the bias lies in that column direction
for (dn in names(dirs)) { d <- dirs[[dn]]
  a_or <- as.vector((Bh - Bt) %*% d)                                             # least squares per row (d has unit norm)
  cat(sprintf("%-13s: Goodman error %.2f -> %.2f after removing the oracle a along delta0 (independence %.2f)\n", dn, nat_err(Bh), nat_err(Bh - a_or %o% d), nat_err(ind)))
  ## per-row interval for a_r from the bounds, then the coupling p'a = 0 (each row's range given the others' ranges)
  I <- t(sapply(1:R, function(r) { ks <- c(-Inf, Inf); for (c in 1:C) if (abs(d[c]) > 1e-12) { b1 <- (Bh[r, c] - hi[r, c]) / d[c]; b2 <- (Bh[r, c] - lo[r, c]) / d[c]
    ks <- c(max(ks[1], min(b1, b2)), min(ks[2], max(b1, b2))) }; ks }))
  ok <- all(I[, 1] <= I[, 2]); cat(sprintf("   bounds feasible for every row: %s\n", ok))
  if (!ok) next
  J2 <- I; for (r in 1:R) { o <- setdiff(1:R, r); J2[r, ] <- c(max(I[r, 1], -sum(p[o] * I[o, 2]) / p[r]), min(I[r, 2], -sum(p[o] * I[o, 1]) / p[r])) }
  cat("   group  | a true | a interval      | L4+ % truth | Goodman | interval            | None % truth | interval\n")
  for (r in 1:R) cat(sprintf("   %-6s | %6.3f | [%6.3f, %6.3f] | %5.1f | %6.1f | %6.1f to %6.1f | %5.1f | %6.1f to %6.1f\n", rn[r], a_or[r], J2[r, 1], J2[r, 2],
    100 * Bt[r, 5], 100 * Bh[r, 5], 100 * (Bh[r, 5] - J2[r, 2] * d[5]), 100 * (Bh[r, 5] - J2[r, 1] * d[5]),
    100 * Bt[r, 1], 100 * (Bh[r, 1] - J2[r, 1] * d[1]), 100 * (Bh[r, 1] - J2[r, 2] * d[1])))
  ## a single point: the a closest to zero correction that satisfies bounds and totals (minimum-norm feasible), and the bounds' centre
  Dm <- diag(p); A <- cbind(p, diag(R), -diag(R)); b0 <- c(0, J2[, 1], -J2[, 2])
  amn <- solve.QP(Dm + diag(1e-9, R), rep(0, R), A, b0, meq = 1)$solution
  amid <- rowMeans(J2); amid <- amid - sum(p * amid)                              # centre of each interval, re-centred for the totals
  cat(sprintf("   point choices: minimum correction %.2f | interval centres %.2f | oracle %.2f\n", nat_err(Bh - amn %o% d), nat_err(Bh - amid %o% d), nat_err(Bh - a_or %o% d)))
}
