## Simulated R x C area tables with known truth, for sweeps over the regimes that matter for EI:
##   assoc     strength of the average association (0 = independence; 1 = strong, election-like diagonal)
##   sort      sorting: a shared column shift in each area, correlated with the area's row composition (the census mechanism;
##             0 = none). It moves every row's rates in the same direction, so the cross-area pattern overstates the association.
##   sigma     area-specific variation of each row's rates (ILR coordinates within a row, independent of composition)
##   seg       segregation: spread of the areas' row compositions (log-ratio sd)
##   J, R, C, N_mean: areas, rows, columns, mean area size (sizes lognormal, sd 0.5 on the log scale)
## Returns list(kc = J x R x C counts, E = the average row-conditional table used, rates = J x R x C true expected rates, settings).
##   assim     assimilation of small groups (row-specific contextual effect): row r's rates in area j become (1 - w) x its own rates
##             + w x the other rows' rates weighted by their shares, w = inv_logit(-1.2 (log s_jr - log h)), h = the share at which a row
##             is half assimilated. assim = 0: none; assim > 0: h = 0.05 * assim (so 1 gives half assimilation at a 5% share).
sim_tables <- function(J = 60, R = 4, C = 4, assoc = 1, sort = 0, sigma = 0.3, seg = 0.8, N_mean = 5000, seed = 1, assim = 0) {
  set.seed(seed)
  H <- function(n) { V <- matrix(0, n, n - 1); for (k in 1:(n - 1)) { V[1:k, k] <- 1 / sqrt(k * (k + 1)); V[k + 1, k] <- -k / sqrt(k * (k + 1)) }; V }
  ## average table: column popularity (unequal) + diagonal association (row r prefers column r when R = C; otherwise a monotone pattern)
  col_pop <- log(seq(1.6, 0.6, length.out = C))
  aff <- outer(seq_len(R), seq_len(C), function(r, c) -abs((r - 1) / max(R - 1, 1) - (c - 1) / max(C - 1, 1)))   # closeness of row and column positions
  L0 <- matrix(col_pop, R, C, byrow = TRUE) + 3 * assoc * (aff - rowMeans(aff))
  E <- exp(L0) / rowSums(exp(L0))
  ## compositions: logistic normal around unequal mean shares
  mu_r <- log(seq(2, 0.5, length.out = R)); Zr <- matrix(rnorm(J * (R - 1)), J, R - 1) %*% t(H(R)) * seg
  rs <- exp(sweep(Zr, 2, mu_r, "+")); rs <- rs / rowSums(rs)
  ## sorting: column shift = Phi %*% (centred log composition); Phi random but fixed for the seed, scaled so the shift sd is about `sort`
  lc <- log(rs); lc <- lc - rowMeans(lc); lc <- sweep(lc, 2, colMeans(lc))
  Phi <- matrix(rnorm(C * R), C, R); shift <- lc %*% t(Phi); shift <- shift - rowMeans(shift)
  shift <- if (sort > 0) shift * sort / sd(as.vector(shift)) else shift * 0
  N <- round(exp(rnorm(J, log(N_mean), 0.5))); kc <- array(0, c(J, R, C)); rates <- array(0, c(J, R, C))
  Hc <- H(C)
  for (j in 1:J) {
    nr <- as.vector(rmultinom(1, N[j], rs[j, ]))
    own <- matrix(0, R, C)
    for (r in 1:R) {
      l <- log(E[r, ]) + shift[j, ] + as.vector(Hc %*% rnorm(C - 1, 0, sigma))
      p <- exp(l - max(l)); own[r, ] <- p / sum(p)
    }
    for (r in 1:R) {
      p <- own[r, ]
      if (assim > 0) {
        w <- plogis(-1.2 * (log(rs[j, r]) - log(0.05 * assim))); o <- setdiff(1:R, r)
        p <- (1 - w) * p + w * colSums(rs[j, o] * own[o, , drop = FALSE]) / sum(rs[j, o])
      }
      rates[j, r, ] <- p
      kc[j, r, ] <- if (nr[r] > 0) as.vector(rmultinom(1, nr[r], p)) else 0
    }
  }
  list(kc = kc, E = E, rates = rates, settings = list(J = J, R = R, C = C, assoc = assoc, sort = sort, sigma = sigma, seg = seg, N_mean = N_mean, seed = seed, assim = assim))
}
