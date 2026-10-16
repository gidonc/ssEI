## Does restricting the average table (not just centring its prior) help locate it? Goodman-type fits from the margins, with the row
## conditional table restricted to a low-dimensional family, against unrestricted constrained Goodman. National table error, and area error
## after raking the fitted table to each area's margins.
## Families:  QI: P(c | r) = l_r if c = r, else (1 - l_r) w_c / sum_{c' != r} w_c'    (loyalty per row + one popularity vector for defectors)
##            QI + bloc: as QI with defector weights w_c * exp(a * [c in r's bloc])                (one extra parameter)
## Also: how well each family can represent the TRUE national table (best member fitted to the truth) = the price of the restriction.
suppressMessages({ source("experiments/data_prep.R"); library(quadprog) })
fitfam <- function(kc, bloc = NULL, own_free = NULL) {
  J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]; rm <- apply(kc, 1:2, sum); cm <- apply(kc, c(1, 3), sum); w <- rowSums(rm)
  xs <- rm / w; ys <- cm / w; N <- sum(kc); kt <- apply(kc, 2:3, sum); nr <- colSums(rm); Bt <- kt / nr
  same <- if (is.null(bloc)) matrix(0, R, C) else outer(bloc[seq_len(R)], bloc[seq_len(C)], "==") * 1
  Pfun <- function(p) { l <- plogis(p[1:R]); lw <- c(0, p[R + 1:(C - 1)]); a <- if (is.null(bloc)) 0 else p[R + C]
    P <- matrix(0, R, C); for (r in 1:R) { ww <- exp(lw + a * same[r, ]); if (r <= C) ww[r] <- 0; ww <- ww / sum(ww); P[r, ] <- (1 - l[r]) * ww; if (r <= C) P[r, r] <- l[r] }; P }
  np <- R + C - 1 + !is.null(bloc)
  obj_margin <- function(p) sum(w * (ys - xs %*% Pfun(p))^2)
  obj_truth <- function(p) sum(nr * (Bt - Pfun(p))^2)
  best <- function(obj) { b <- NULL; set.seed(1); for (s in 1:6) { o <- optim(c(rep(qlogis(0.7), R), rnorm(C - 1 + !is.null(bloc), 0, 0.5)), obj, method = "BFGS", control = list(maxit = 2000)); if (is.null(b) || o$value < b$value) b <- o }; Pfun(b$par) }
  Pm <- best(obj_margin); Ptr <- best(obj_truth)
  ## unrestricted constrained Goodman
  Dq <- matrix(0, R * C, R * C); dv <- numeric(R * C)
  for (j in 1:J) { Dq <- Dq + w[j] * kronecker(tcrossprod(xs[j, ]), diag(C)); dv <- dv + w[j] * kronecker(xs[j, ], ys[j, ]) }
  A <- cbind(t(kronecker(diag(R), matrix(1, 1, C))), diag(R * C))
  Pg <- matrix(solve.QP(2 * Dq / N + diag(1e-9, R * C), 2 * dv / N, A, c(rep(1, R), rep(0, R * C)), meq = R)$solution, R, C, byrow = TRUE)
  nat <- function(P) 50 * sum(nr * abs(P - Bt)) / N
  rake <- function(P) { e <- array(0, dim(kc)); for (j in 1:J) { X <- P * rm[j, ] + 1e-9; X[rm[j, ] == 0, ] <- 0
    for (i in 1:1500) { rs <- rowSums(X); X <- X * ifelse(rs > 0, rm[j, ] / pmax(rs, 1e-300), 0); cs <- colSums(X); X <- sweep(X, 2, ifelse(cs > 0, cm[j, ] / pmax(cs, 1e-300), 0), "*")
      if (i %% 25 == 0 && max(abs(rowSums(X) - rm[j, ])) < 1e-6) break }; e[j, , ] <- X }; 50 * sum(abs(e - kc)) / N }
  c(goodman_nat = nat(Pg), goodman_area = rake(Pg), family_nat = nat(Pm), family_area = rake(Pm), family_best_on_truth_nat = nat(Ptr), truth_raked_area = rake(Bt), params = np)
}
s <- prep_scotland(7, 73)
cat("Scotland 2007 (7 x 7): loyalty + quasi-independent defection (13 parameters against 42)\n"); print(round(fitfam(s$kc), 2))
z <- prep_nz(2020, 6, "candidate", 1e6); cat("\nNZ 2020 rows:", z$row_names, "\n")
bl <- c(1, 2, 1, 3, 2, 3)   # Labour, National, Green, NZFirst, ACT, Other: left bloc (Lab, Green), right (National, ACT), other
cat("NZ 2020: loyalty + quasi-independent defection\n"); print(round(fitfam(z$kc), 2))
cat("NZ 2020: + one bloc affinity parameter\n"); print(round(fitfam(z$kc, bloc = bl), 2))
