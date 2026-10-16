## Measured urbanity from the LTLA code prefix (E09 London borough, E08 metropolitan district, E06 unitary, E07 shire district):
## constrained Goodman with a common column shift D' U_j (all groups shifted alike), national table against the truth.
source_lines <- readLines("experiments/sorting/census_kappa.R"); eval(parse(text = source_lines[1:grep("^## 1\\.", source_lines) - 1]))
suppressMessages(library(quadprog)); codes <- dimnames(readRDS(f))[[1]]; typ <- substr(codes, 1, 3); print(table(typ))
fitU <- function(U) {                                       # variables: vec(B) row-major (R*C), then vec(D) row-major (k*C); B rows in simplex, D rows sum 0
  U <- as.matrix(U); k <- ncol(U); Z <- cbind(xs, U); nv <- (R + k) * C
  Dq <- matrix(0, nv, nv); dv <- numeric(nv)
  for (j in 1:J) { Dq <- Dq + w[j] * kronecker(tcrossprod(Z[j, ]), diag(C)); dv <- dv + w[j] * kronecker(Z[j, ], ys[j, ]) }
  Dq <- 2 * Dq / N + diag(1e-9, nv); dv <- 2 * dv / N
  Aeq <- cbind(t(kronecker(diag(R + k), matrix(1, 1, C)))); beq <- c(rep(1, R), rep(0, k))
  Ain <- rbind(diag(R * C), matrix(0, k * C, R * C)); A <- cbind(Aeq, Ain)
  s <- solve.QP(Dq, dv, A, c(beq, rep(0, R * C)), meq = R + k)$solution
  B <- matrix(s[1:(R * C)], R, C, byrow = TRUE); dimnames(B) <- dimnames(Bt); B }
cat("\nR-squared of each group's share on the urbanity dummies (identification needs these well below 1):\n")
Ud <- cbind(london = typ == "E09", metro = typ == "E08", unitary = typ == "E06") * 1
print(round(sapply(1:R, function(r) summary(lm(xs[, r] ~ Ud, weights = w))$r.squared), 2))
specs <- list(none = NULL, london = Ud[, 1, drop = FALSE], london_metro = Ud[, 1:2], all_three = Ud)
Bc <- fitU(matrix(0, J, 1))
cat(sprintf("\nNational error: independence %.2f | constrained Goodman %.2f\n", nat_err(ind), nat_err(Bc)))
for (nm in names(specs)[-1]) { B <- fitU(specs[[nm]]); cat(sprintf("  + common shift on %-12s %.2f | L4+ %%: %s\n", nm, nat_err(B), paste(rn, round(100 * B[, 5], 1), collapse = " "))) }
cat("  truth L4+ %:", paste(rn, round(100 * Bt[, 5], 1), collapse = " "), "\n")
cat("  constrained Goodman L4+ %:", paste(rn, round(100 * Bc[, 5], 1), collapse = " "), "\n")
