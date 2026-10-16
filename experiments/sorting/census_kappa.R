## One-scalar sorting analysis on census ethnicity x qualification (309 LTLAs, England 2021), national table.
## Linear Goodman: B_hat = B + a delta' if sorting is one latent dimension U_j = a' x_j + xi_j with a common column shift delta U_j.
## Estimate the row direction a0 from the composition, the column direction delta0 from the residuals (or substantively), then bracket
## kappa in B(kappa) = B_hat - kappa a0 delta0' by aggregated Duncan-Davis bounds and the sign. Truth is known, so every step is checked.
f <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "/home/claude/data_census/eth6_ed7.rds"
k <- readRDS(f); k <- array(as.numeric(k), dim(k)); J <- dim(k)[1]; R <- dim(k)[2]; C <- dim(k)[3]
rn <- c("Asian", "Black", "Mixed", "White", "Other"); cn <- c("None", "L1", "L2", "L3", "L4+", "Oth")
rm <- apply(k, 1:2, sum); cm <- apply(k, c(1, 3), sum); w <- rowSums(rm); N <- sum(w); xs <- rm / w; ys <- cm / w
nr <- colSums(rm); p <- nr / N; Bt <- apply(k, 2:3, sum) / nr; dimnames(Bt) <- list(rn, cn)
nat_err <- function(B) 50 * sum(nr * abs(B - Bt)) / N
gls <- function(X, Y) solve(crossprod(X * sqrt(w)), crossprod(X * w, Y))
Bh <- gls(xs, ys); dimnames(Bh) <- dimnames(Bt)
ind <- matrix(colSums(cm) / N, R, C, byrow = TRUE)
cat(sprintf("National error: Goodman (unconstrained) %.2f | independence %.2f\n", nat_err(Bh), nat_err(ind)))
cat("\nTruth (row %):\n"); print(round(100 * Bt, 1)); cat("\nGoodman (row %):\n"); print(round(100 * Bh, 1))

## 1. Is the Goodman bias close to rank 1? (rows weighted by sqrt(p) so the metric matches the national error)
D <- Bh - Bt; sv <- svd(sqrt(p) * D)
cat(sprintf("\n1. Bias D = Goodman - truth: singular values %s | share of squared bias in rank 1: %.2f\n",
            paste(round(sv$d, 3), collapse = " "), sv$d[1]^2 / sum(sv$d^2)))
a_true <- sv$u[, 1] / sqrt(p); a_true <- a_true / sqrt(sum(p * a_true^2)); d_true <- sv$v[, 1]
if (sum(a_true * (rn != "White")) < 0) { a_true <- -a_true; d_true <- -d_true }
cat("   oracle row direction (minorities +):", paste(rn, round(a_true, 2)), "\n   oracle column direction:", paste(cn, round(d_true, 2)), "\n")
## national-table error left after removing the best rank-1 piece (what a perfect one-scalar correction could reach at best)
D1 <- (sv$d[1] * sv$u[, 1] %o% sv$v[, 1]) / sqrt(p); cat(sprintf("   error after removing the oracle rank-1 part: %.2f\n", nat_err(Bh - D1)))

## 2. Row direction from the composition alone
cw <- cov.wt(xs, wt = w / N)$cov; e <- eigen(cw); v1 <- e$vectors[, 1]
cat(sprintf("\n2. Composition PC1 explains %.2f of compositional variance\n", e$values[1] / sum(e$values)))
lx <- log(xs + 1e-4); lx <- lx - rowMeans(lx); el <- eigen(cov.wt(lx, wt = w / N)$cov); vl <- el$vectors[, 1]
cat(sprintf("   log-composition PC1 explains %.2f\n", el$values[1] / sum(el$values)))
cand_a <- list(share_pc1 = v1, logshare_pc1 = vl / p)                            # log version linearised at the national shares
norm_a <- function(a) { a <- a - sum(p * a); a <- a / sqrt(sum(p * a^2)); if (sum(a * (rn != "White")) < 0) -a else a }
cand_a <- lapply(cand_a, norm_a)
cosw <- function(x, y) sum(p * x * y) / sqrt(sum(p * x^2) * sum(p * y^2))
for (nm in names(cand_a)) cat(sprintf("   %-13s %s | cosine with oracle %.2f\n", nm, paste(rn, round(cand_a[[nm]], 2), collapse = " "), cosw(cand_a[[nm]], a_true)))

## 3. Column direction: from the Goodman residuals, or substantive (urbanity moves people up the scale; 'Other' left at 0)
res <- ys - xs %*% Bh; er <- eigen(cov.wt(res, wt = w / N)$cov); d_res <- er$vectors[, 1]
cat(sprintf("\n3. Residual PC1 explains %.2f of residual column-margin variance\n", er$values[1] / sum(er$values)))
d_sub <- c(-2, -1, -0.5, 0.5, 3, 0)
norm_d <- function(d) { d <- d - mean(d); d <- d / sqrt(sum(d^2)); if (d[5] < 0) -d else d }
cand_d <- lapply(list(residual_pc1 = d_res, substantive = d_sub), norm_d)
for (nm in names(cand_d)) cat(sprintf("   %-13s %s | cosine with oracle %.2f\n", nm, paste(cn, round(cand_d[[nm]], 2), collapse = " "), sum(cand_d[[nm]] * d_true)))

## 4. Aggregated Duncan-Davis bounds for the national table; kappa interval for each pair of directions
lo <- hi <- matrix(0, R, C)
for (j in 1:J) for (r in 1:R) for (c in 1:C) { lo[r, c] <- lo[r, c] + max(0, cm[j, c] - (w[j] - rm[j, r])); hi[r, c] <- hi[r, c] + min(rm[j, r], cm[j, c]) }
lo <- lo / nr; hi <- hi / nr
cat(sprintf("\n4. Truth inside the aggregated bounds: %s | mean bound width %.3f\n", all(Bt >= lo - 1e-9 & Bt <= hi + 1e-9), mean(hi - lo)))
kint <- function(a, d) { M <- a %o% d; ks <- c(-Inf, Inf)        # B_hat - kappa M within [lo, hi]
  for (i in seq_along(M)) if (abs(M[i]) > 1e-12) { b1 <- (Bh[i] - hi[i]) / M[i]; b2 <- (Bh[i] - lo[i]) / M[i]; ks <- c(max(ks[1], min(b1, b2)), min(ks[2], max(b1, b2))) }
  ks }
cat("\n   directions                      | kappa interval (bounds) | with kappa>=0 | error at: 0, ends, mid | best kappa (oracle) and its error\n")
for (an in names(cand_a)) for (dn in names(cand_d)) { a <- cand_a[[an]]; d <- cand_d[[dn]]; M <- a %o% d; ki <- kint(a, d); kp <- c(max(ki[1], 0), ki[2])
  ke <- function(kk) nat_err(Bh - kk * M); kb <- optimize(ke, c(min(ki[1], -1), max(ki[2], 1)))$minimum
  cat(sprintf("   %-13s x %-13s | [%6.3f, %6.3f] | [%6.3f, %6.3f] | %.2f, %.2f / %.2f, mid %.2f | %.3f -> %.2f\n", an, dn, ki[1], ki[2], kp[1], kp[2],
              ke(0), ke(kp[1]), ke(kp[2]), ke(mean(kp)), kb, ke(kb))) }

## 5. What a known row direction identifies: contrasts orthogonal to a0 (p-weighted) against the contrast along a0
a <- cand_a$share_pc1; Q <- qr.Q(qr(cbind(sqrt(p) * a, diag(R))))[, 2:R]          # orthonormal complement in the sqrt(p) metric
Kc <- t(Q) * 0; Kc <- sweep(t(Q), 2, sqrt(p), "*")                                   # contrast rows: sum_r K[r] * p_r^(1/2)...
along <- sum(p * a * 0) ; err_orth <- sqrt(sum((Kc %*% D)^2)); err_along <- sqrt(sum((sum(p * a) * 0 + crossprod(p * a, D))^2))
cat(sprintf("\n5. Goodman bias in contrasts orthogonal to the composition direction: %.3f | along it: %.3f\n", err_orth, err_along))
cat("\n   L4+ (degree-level) share by group, %: truth | Goodman | corrected range (share PC1 x substantive, kappa>=0 within bounds)\n")
a <- cand_a$share_pc1; d <- cand_d$substantive; M <- a %o% d; ki <- kint(a, d); kp <- c(max(ki[1], 0), ki[2])
for (r in 1:R) cat(sprintf("   %-6s %5.1f | %5.1f | %5.1f to %5.1f\n", rn[r], 100 * Bt[r, 5], 100 * Bh[r, 5], 100 * (Bh[r, 5] - kp[2] * M[r, 5]), 100 * (Bh[r, 5] - kp[1] * M[r, 5])))
