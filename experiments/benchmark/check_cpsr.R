## Check of claims in an outside report (October 2026): (1) raking a naive 75%-loyalty seed beats ssEI on Scotland;
## (2) 'covariance-projected saddlepoint raking' (CPSR) from the oracle national table beats IPF raking of the same seed.
## CPSR: per area, the minimum of (b - E)' Omega^-1 (b - E) over within-row ILR coordinates b, subject to the area's column margins,
## with Omega(rho) = (1 - rho) I + (rho / R) Wm Wm' (Wm = the shared column-shift directions); rho = 1 is IPF, rho = 0 an isotropic
## Gaussian projection. Implemented here independently (SQP / Newton on the margin constraint), then the result is raked to the exact margins.
## usage: Rscript experiments/benchmark/check_cpsr.R
suppressMessages({ source("experiments/data_prep.R"); source("experiments/benchmark/bm_data.R") })
H <- function(n) { V <- matrix(0, n, n - 1); for (k in 1:(n - 1)) { V[1:k, k] <- 1 / sqrt(k * (k + 1)); V[k + 1, k] <- -k / sqrt(k * (k + 1)) }; V }
ipf <- function(X, rt, ct, it = 2000) { for (i in 1:it) { rs <- rowSums(X); X <- X * ifelse(rs > 0, rt / pmax(rs, 1e-300), 0); cs <- colSums(X)
  X <- sweep(X, 2, ifelse(cs > 0, ct / pmax(cs, 1e-300), 0), "*"); if (i %% 25 == 0 && max(abs(rowSums(X) - rt)) < 1e-7) break }; X }
cpsr <- function(P, kc, rho, n_iter = 50) {
  J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]; Hc <- H(C)
  V <- kronecker(diag(R), Hc)                                # within-row ILR basis; rows of V indexed (r - 1) C + c (row-major)
  Wm <- crossprod(V, kronecker(rep(1, R), Hc))                # shared column shift
  Om <- (1 - rho) * diag(R * (C - 1)) + (rho / R) * tcrossprod(Wm)
  E <- as.vector(crossprod(V, as.vector(t(log(P)))))          # ILR of the seed
  out <- array(0, dim(kc))
  for (j in 1:J) {
    rt <- rowSums(kc[j, , ]); ct <- colSums(kc[j, , ]); N <- sum(rt); if (N == 0) next
    w <- rt / N; m <- ct / N; b <- E
    for (it in 1:n_iter) {
      q <- matrix(V %*% b, R, C, byrow = TRUE); q <- exp(q - apply(q, 1, max)); q <- q / rowSums(q)
      err <- (m - colSums(w * q))[-C]; if (max(abs(err)) < 1e-10) break
      ## d m_c / d logq_rc = w_r q_rc (1[c = c'] - q_rc'), arranged row-major to match V
      Jl <- matrix(0, C - 1, R * C)
      for (cc in 1:(C - 1)) { G <- -w * q[, cc] * q; G[, cc] <- G[, cc] + w * q[, cc]; Jl[cc, ] <- as.vector(t(G)) }
      Jb <- Jl %*% V
      lam <- solve(Jb %*% Om %*% t(Jb) + diag(1e-10, C - 1), err + Jb %*% (b - E))
      bt <- E + as.vector(Om %*% t(Jb) %*% lam)
      step <- 1
      for (ls in 1:12) { bn <- (1 - step) * b + step * bt
        qn <- matrix(V %*% bn, R, C, byrow = TRUE); qn <- exp(qn - apply(qn, 1, max)); qn <- qn / rowSums(qn)
        if (max(abs((m - colSums(w * qn))[-C])) < max(abs(err))) break; step <- step / 2 }
      b <- bn
    }
    q <- matrix(V %*% b, R, C, byrow = TRUE); q <- exp(q - apply(q, 1, max)); q <- q / rowSums(q)
    out[j, , ] <- ipf(q * rt + 1e-12, rt, ct)                 # exact margins (small final correction)
  }
  out
}
rake_seed <- function(P, kc) { out <- array(0, dim(kc)); for (j in 1:dim(kc)[1]) out[j, , ] <- ipf(P * rowSums(kc[j, , ]) + 1e-12, rowSums(kc[j, , ]), colSums(kc[j, , ])); out }
err <- function(est, kc) c(area = 50 * sum(abs(est - kc)) / sum(kc), nat = 50 * sum(abs(apply(est, 2:3, sum) - apply(kc, 2:3, sum))) / sum(kc))
qi_seed <- function(kc, loy = 0.75) { R <- dim(kc)[2]; C <- dim(kc)[3]; w <- colSums(apply(kc, c(1, 3), sum)); w <- w / sum(w)
  P <- matrix(0, R, C); for (r in 1:R) { ww <- w; if (r <= C) ww[r] <- 0; P[r, ] <- (1 - loy) * ww / sum(ww); if (r <= C) P[r, r] <- loy }; P }

dsets <- list(senc = bm_data("senc")$kc, scotland = prep_scotland(7, 73)$kc, nz2020 = prep_nz(2020, 6, "candidate", 1e6)$kc)
for (nm in names(dsets)) {
  kc <- dsets[[nm]]; Bt <- apply(kc, 2:3, sum); Bt <- Bt / rowSums(Bt); Bt <- pmax(Bt, 1e-6); Bt <- Bt / rowSums(Bt)
  cat(sprintf("\n%s (%d x %d x %d): area | national error\n", nm, dim(kc)[1], dim(kc)[2], dim(kc)[3]))
  r <- function(lab, est) { e <- err(est, kc); cat(sprintf("  %-38s %6.2f | %5.2f\n", lab, e[1], e[2])) }
  if (nm != "senc") r("75%-loyalty seed + IPF (claimed: Sco 6.60 | 5.53)", rake_seed(qi_seed(kc), kc))
  r("oracle national seed + IPF", rake_seed(Bt, kc))
  for (rho in c(0, 0.5, 0.85, 0.95)) r(sprintf("oracle national seed + CPSR rho = %.2f", rho), cpsr(Bt, kc, rho))
}
