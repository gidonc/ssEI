## Scotland 2007 truth: how much of the area variation is regional, regional list tactics, and Goodman with a regional column shift.
src <- readLines("experiments/scotland_structure/truth_structure.R"); eval(parse(text = src[1:grep("^cat\\(\"\\\\n== 4", src) - 1]))
rg <- read.csv("experiments/scotland_structure/sco2007_regions.csv"); reg <- factor(rg$region)
r2 <- function(y, ...) summary(lm(y ~ ...))$r.squared
cat("\n\n== A. Share of the common column shift explained by region, by composition, and by both (R-squared) ==\n")
for (c in 1:C) cat(sprintf("  %-6s region %.2f | composition %.2f | both %.2f\n", lab[c], r2(S[, c], reg), r2(S[, c], rs[, 1:4]), summary(lm(S[, c] ~ reg + rs[, 1:4]))$r.squared))
Y <- cbind(sapply(maj, function(r) lg(rate[, r, r]) - mean(lg(rate[, r, r]))), do.call(cbind, lapply(maj, function(r) D[[r]][, -1])))
e <- eigen(cov(Y)); pcs <- Y %*% e$vectors[, 1:3]
cat(sprintf("  full deviations (24 dims): PC1 R-squared on region %.2f, PC2 %.2f, PC3 %.2f; all 24 dims together %.2f\n",
    r2(pcs[, 1], reg), r2(pcs[, 2], reg), r2(pcs[, 3], reg), 1 - sum(resid(lm(Y ~ reg))^2) / sum(scale(Y, scale = FALSE)^2)))
cat(sprintf("  loyalty logits: R-squared on region %s\n", paste(sprintf("%s %.2f", lab[maj], sapply(maj, function(r) r2(lg(rate[, r, r]), reg))), collapse = " ")))
cat("\n== B. By region: seats won (constituency winners), list shares, and major-party loyalty (row view, %) ==\n")
win <- maj[apply(rs[, maj], 1, which.max)]
tab <- t(sapply(levels(reg), function(g) { ix <- reg == g
  c(seats = sapply(maj, function(k) sum(win[ix] == k)), list_other = 100 * sum(cmg[ix, 6] * rowSums(N[ix, ])) / sum(N[ix, ]), list_green = 100 * sum(cmg[ix, 7] * rowSums(N[ix, ])) / sum(N[ix, ]),
    loy = sapply(maj, function(k) 100 * sum(kc[ix, k, k]) / sum(N[ix, k]))) }))
colnames(tab) <- c(paste0("won_", lab[maj]), "list_Other%", "list_Green%", paste0("loy_", lab[maj])); print(round(tab, 1))
cat(sprintf("  cor across regions: Labour seats won vs Labour loyalty %.2f; SNP seats won vs SNP loyalty %.2f\n", cor(tab[, "won_Lab"], tab[, "loy_Lab"]), cor(tab[, "won_SNP"], tab[, "loy_SNP"])))
cat("  where Labour constituency voters' list votes go when not Labour, by region (%):\n")
for (g in levels(reg)) { ix <- reg == g; d <- colSums(kc[ix, 2, -2]); cat(sprintf("    %-22s %s\n", g, paste(sprintf("%s %.0f", lab[-2], 100 * d / sum(d)), collapse = " "))) }
cat("\n== C. Goodman (constrained to valid rates) with a common column shift on region dummies: national error ==\n")
suppressMessages(library(quadprog)); xs <- rs; ys <- cmg; w <- rowSums(N); Ntot <- sum(w); Bt <- apply(kc, 2:3, sum) / colSums(N)
fitU <- function(U) { U <- as.matrix(U); k <- ncol(U); Z <- cbind(xs, U); nv <- (R + k) * C; Dq <- matrix(0, nv, nv); dv <- numeric(nv)
  for (j in 1:J) { Dq <- Dq + w[j] * kronecker(tcrossprod(Z[j, ]), diag(C)); dv <- dv + w[j] * kronecker(Z[j, ], ys[j, ]) }
  Dq <- 2 * Dq / Ntot + diag(1e-9, nv); dv <- 2 * dv / Ntot; A <- cbind(t(kronecker(diag(R + k), matrix(1, 1, C))), rbind(diag(R * C), matrix(0, k * C, R * C)))
  s <- solve.QP(Dq, dv, A, c(rep(1, R), rep(0, k), rep(0, R * C)), meq = R + k)$solution; matrix(s[1:(R * C)], R, C, byrow = TRUE) }
ne <- function(B) 50 * sum(colSums(N) * abs(B - Bt)) / Ntot
Ureg <- model.matrix(~ reg)[, -1]; Ureg <- sweep(Ureg, 2, colMeans(Ureg))
B0 <- fitU(matrix(0, J, 1)); B1 <- fitU(Ureg)
cat(sprintf("  Goodman %.2f | + regional column shift %.2f\n", ne(B0), ne(B1)))
dimnames(B1) <- list(lab, lab); cat("  bias with the regional shift (percentage points):\n"); print(round(100 * (B1 - Bt), 1))
