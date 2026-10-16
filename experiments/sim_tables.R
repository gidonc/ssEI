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
##   sort_rank rank of the sorting map Phi: 0 = full (as drawn); k > 0 = its best rank-k approximation (after centring), from the same
##             random draw, so the data otherwise match the full-rank dataset with the same seed. Rank 1 = people sort on one latent
##             dimension s_j = v' lc_j (returned as sort_dir = v, the row loading, and sort_score = s).
##   sort_form "log" (default): the shift is linear in the centred log composition; "share": linear in the centred shares themselves
##             (same draws otherwise). Tests whether a method's gain comes from matching the log functional form.
sim_tables <- function(J = 60, R = 4, C = 4, assoc = 1, sort = 0, sigma = 0.3, seg = 0.8, N_mean = 5000, seed = 1, assim = 0, sort_rank = 0, sort_form = "log") {
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
  lc <- if (sort_form == "share") rs else log(rs); lc <- lc - rowMeans(lc); lc <- sweep(lc, 2, colMeans(lc))
  Phi <- matrix(rnorm(C * R), C, R)
  Pc <- diag(C) - 1 / C; Pr <- diag(R) - 1 / R; sv <- svd(Pc %*% Phi %*% Pr)                      # centred map: what reaches the shift
  k <- if (sort_rank > 0) min(sort_rank, length(sv$d)) else length(sv$d)
  Phi_used <- sv$u[, 1:k, drop = FALSE] %*% diag(sv$d[1:k], k) %*% t(sv$v[, 1:k, drop = FALSE])
  if (sort_rank == 0) Phi_used <- Phi                                                            # full rank: exactly as before
  shift <- lc %*% t(Phi_used); shift <- shift - rowMeans(shift)
  sort_dir <- sv$v[, 1:k, drop = FALSE]; sort_score <- lc %*% sort_dir
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
  list(kc = kc, E = E, rates = rates, sort_dir = sort_dir, sort_score = sort_score,
       settings = list(J = J, R = R, C = C, assoc = assoc, sort = sort, sigma = sigma, seg = seg, N_mean = N_mean, seed = seed, assim = assim, sort_rank = sort_rank, sort_form = sort_form))
}
