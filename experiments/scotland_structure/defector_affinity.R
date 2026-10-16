## Where defectors go ON AVERAGE (national true table): affinities after removing row size and column popularity.
## Quasi-independence fitted to the off-diagonal cells (Poisson log-linear, row + column effects, diagonal excluded); residual = log(observed / fitted).
## Scotland 2007 7x7, and New Zealand 2020 (candidate vote x party vote) for contrast.
suppressMessages(source("experiments/data_prep.R"))
aff <- function(kt, lab, show) {
  R <- nrow(kt); d <- expand.grid(r = 1:R, c = 1:R); d$n <- as.vector(kt); d <- d[d$r != d$c, ]
  f <- glm(n ~ factor(r) + factor(c), family = poisson, data = d)
  A <- matrix(NA, R, R, dimnames = list(lab, lab)); A[cbind(d$r, d$c)] <- log(d$n / fitted(f)); print(round(A[show, show], 2))
  cat(sprintf("  share of off-diagonal deviance left after row + column effects: %.2f (1 = no structure from popularity)\n", deviance(f) / glm(n ~ 1, family = poisson, data = d)$deviance)) }
s <- prep_scotland(7, 73); kt <- apply(s$kc, 2:3, sum)
cat("Scotland 2007: rows constituency vote, columns list vote\n"); aff(kt, c("Con", "Lab", "LD", "SNP", "Spoilt", "Other", "Green"), 1:7)
z <- prep_nz(2020, 6, "candidate", 1e6); ktz <- apply(z$kc, 2:3, sum)
cat("\nNew Zealand 2020: rows candidate vote, columns party vote\n"); print(z$row_names); aff(ktz, z$row_names, seq_len(nrow(ktz)))
