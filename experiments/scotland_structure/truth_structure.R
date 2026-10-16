## Scotland 2007, 7x7, 73 constituencies: how the TRUE tables vary across areas. Is it "loyalty roughly constant + one latent dimension"?
## Uses the four major rows (Con, Lab, LD, SNP: large in every area). Rates smoothed by +0.5.
a <- commandArgs(TRUE); x <- readRDS(if (length(a)) a[1] else "sco7.rds"); kc <- x$kc; J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]
lab <- c("Con", "Lab", "LD", "SNP", "Spoilt", "Other", "Green"); maj <- 1:4
N <- apply(kc, 1:2, sum); rs <- N / rowSums(N); cmg <- apply(kc, c(1, 3), sum) / rowSums(N)
rate <- (kc + 0.5) / array(N + 0.5 * C, dim(kc))
lg <- function(p) log(p / (1 - p))
cat("== 1. Loyalty (own list vote) by area: logit sd across areas against sampling sd ==\n")
for (r in maj) { l <- lg(rate[, r, r]); p <- mean(rate[, r, r]); ssd <- sqrt(mean(1 / (N[, r] * rate[, r, r] * (1 - rate[, r, r]))))
  cat(sprintf("  %-4s mean %.2f range %.2f-%.2f | logit sd %.2f (sampling %.2f) | cor with own row share %.2f, with SNP share %.2f\n", lab[r], p,
              min(rate[, r, r]), max(rate[, r, r]), sd(l), ssd, cor(l, rs[, r]), cor(l, rs[, 4]))) }
cat("\n== 2. Where defectors go: log-ratio of each other column against the row's defector total, sd across areas ==\n")
D <- list()                                                         # per row: J x 6 centred log shares among defectors (clr over the 6)
for (r in maj) { o <- setdiff(1:C, r); q <- rate[, r, o] / rowSums(rate[, r, o]); l <- log(q); l <- l - rowMeans(l)
  D[[r]] <- sweep(l, 2, colMeans(l)); colnames(D[[r]]) <- lab[o]
  noise <- sqrt(colMeans(1 / ((N[, r] * (1 - rate[, r, r])) %o% rep(1, C - 1) * q)))
  cat(sprintf("  %-4s defect to: %s\n", lab[r], paste(sprintf("%s %.2f(%.2f)", lab[o], apply(l, 2, sd), noise), collapse = " "))) }
cat("  (in brackets: rough sampling sd)\n")
cat("\n== 3. A common column shift: does each area move all rows' defectors towards the same parties? ==\n")
## for each area and column c, the shift = mean over major rows r != c of D[[r]][, c]; share of defector-split variance it explains
S <- sapply(1:C, function(c) { rr <- setdiff(maj, c); rowMeans(sapply(rr, function(r) D[[r]][, lab[c]])) }); colnames(S) <- lab
tot <- 0; expl <- 0
for (r in maj) { o <- setdiff(1:C, r); fit <- S[, o] - rowMeans(S[, o]); tot <- tot + sum(D[[r]]^2); expl <- expl + sum(D[[r]]^2) - sum((D[[r]] - fit)^2) }
cat(sprintf("  share of defector-split variation (4 rows) explained by one common shift per area and column: %.2f\n", expl / tot))
cat("  sd of the common shift by column:", paste(sprintf("%s %.2f", lab, apply(S, 2, sd)), collapse = " "), "\n")
cat("\n== 4. Factor structure of all area deviations (loyalty logits + defector splits, 4 rows; 4 + 4 x 5 = 24 dims) ==\n")
Y <- cbind(sapply(maj, function(r) lg(rate[, r, r]) - mean(lg(rate[, r, r]))), do.call(cbind, lapply(maj, function(r) D[[r]][, -1])))
colnames(Y) <- c(paste0("loy_", lab[maj]), unlist(lapply(maj, function(r) paste0(lab[r], ">", colnames(D[[r]])[-1]))))
e <- eigen(cov(Y)); cat("  variance shares of the first 6 components:", round(e$values[1:6] / sum(e$values), 2), "\n")
pc <- Y %*% e$vectors[, 1:3]
for (k in 1:2) { ld <- e$vectors[, k]; o <- order(-abs(ld))[1:8]
  cat(sprintf("  PC%d largest loadings: %s\n", k, paste(sprintf("%s %+.2f", colnames(Y)[o], ld[o]), collapse = ", ")))
  cat(sprintf("       cor of scores with: common shift to Con %.2f, Lab %.2f, LD %.2f, SNP %.2f | row shares Con %.2f Lab %.2f LD %.2f SNP %.2f\n",
    cor(pc[, k], S[, 1]), cor(pc[, k], S[, 2]), cor(pc[, k], S[, 3]), cor(pc[, k], S[, 4]), cor(pc[, k], rs[, 1]), cor(pc[, k], rs[, 2]), cor(pc[, k], rs[, 3]), cor(pc[, k], rs[, 4]))) }
cat("\n== 5. Is the common shift related to composition (sorting)? R-squared of each column's shift on the four row shares ==\n")
for (c in 1:C) cat(sprintf("  %-6s %.2f\n", lab[c], summary(lm(S[, c] ~ rs[, 1:4]))$r.squared))

cat("\n== 6. Local party strength: in areas where party k's constituency share is higher, do OTHER rows' voters also give k more list votes? ==\n")
for (k in maj) { o <- setdiff(maj, k); v <- rowMeans(sapply(o, function(r) lg(rate[, r, k])))
  cat(sprintf("  %-4s: cor(row share of %s, mean logit of other major rows' list vote for %s) = %.2f\n", lab[k], lab[k], lab[k], cor(rs[, k], v))) }
cat("\n== 7. Goodman (constrained to valid rates) against the truth, national row rates in percentage points: bias on the diagonal? ==\n")
suppressMessages(library(quadprog)); xs <- rs; ys <- cmg; w <- rowSums(N); Ntot <- sum(w)
Dq <- matrix(0, R * C, R * C); dv <- numeric(R * C)
for (j in 1:J) { Dq <- Dq + w[j] * kronecker(tcrossprod(xs[j, ]), diag(C)); dv <- dv + w[j] * kronecker(xs[j, ], ys[j, ]) }
Dq <- 2 * Dq / Ntot + diag(1e-9, R * C); dv <- 2 * dv / Ntot; A <- cbind(t(kronecker(diag(R), matrix(1, 1, C))), diag(R * C))
B <- matrix(solve.QP(Dq, dv, A, c(rep(1, R), rep(0, R * C)), meq = R)$solution, R, C, byrow = TRUE)
Bt <- apply(kc, 2:3, sum) / colSums(N); dimnames(B) <- dimnames(Bt) <- list(lab, lab)
cat("  bias (Goodman - truth):\n"); print(round(100 * (B - Bt), 1))
cat(sprintf("  national error: Goodman %.2f\n", 50 * sum(colSums(N) * abs(B - Bt)) / Ntot))
