## Cheap check: replace the spoilt (and Other) constituency rows of existing Scotland 7x7 estimates with structured versions, re-rake the
## remaining rows to each area's margins, and compare errors with the truth. Area error 50*sum|est-true|/N; national error on the summed table.
## Spoilt structure: a share pi of the row is rejected on both halves (spoilt list column, capped by the column margin); the rest splits like
##   the area's valid list vote (independence over the valid columns). pi is not identified from margins: a grid, with 0.15 near the truth (0.156).
## Other structure: (a) independence over the area's list vote; (b) ceiling: the true national Other-row profile as the seed for every area.
S <- "/tmp/claude-0/-home-claude/e71a6b8d-b589-5b64-8f29-846a7938884f/scratchpad"
x <- readRDS(file.path(S, "sco7.rds")); kc <- x$kc; J <- dim(kc)[1]; R <- 7; C <- 7; N <- sum(kc)
rm <- apply(kc, 1:2, sum); cm <- apply(kc, c(1, 3), sum); kt <- apply(kc, 2:3, sum)
err <- function(e) c(area = 50 * sum(abs(e - kc)) / N, national = 50 * sum(abs(apply(e, 2:3, sum) - kt)) / N,
                     spoilt_row = 50 * sum(abs(e[, 5, ] - kc[, 5, ])) / sum(kc[, 5, ]), other_row = 50 * sum(abs(e[, 6, ] - kc[, 6, ])) / sum(kc[, 6, ]))
rake <- function(seed, r, c) { X <- seed + 1e-9; X[r == 0, ] <- 0; X[, c == 0] <- 0
  for (i in 1:2000) { rs <- rowSums(X); X <- X * ifelse(rs > 0, r / pmax(rs, 1e-300), 0); cs <- colSums(X); X <- sweep(X, 2, ifelse(cs > 0, c / pmax(cs, 1e-300), 0), "*")
    if (i %% 25 == 0 && max(abs(rowSums(X) - r)) < 1e-6) break }; X }
## structured rows for area j: returns a list of fixed rows (row index -> counts)
spoilt_row <- function(j, pi) { s <- rm[j, 5]; v <- cm[j, ]; dbl <- min(pi * s, v[5]); rest <- s - dbl; vv <- v; vv[5] <- 0
  out <- rest * vv / sum(vv); out[5] <- out[5] + dbl; out }
other_ind <- function(j) { s <- rm[j, 6]; s * cm[j, ] / sum(cm[j, ]) }
apply_struct <- function(est, fixed_fun, seed_fun = NULL) { o <- est
  for (j in 1:J) { fx <- fixed_fun(j); seed <- est[j, , ]; if (!is.null(seed_fun)) seed <- seed_fun(j, seed)
    free <- setdiff(1:R, as.integer(names(fx))); X <- matrix(0, R, C)
    for (nm in names(fx)) X[as.integer(nm), ] <- fx[[nm]]
    cfree <- pmax(cm[j, ] - colSums(X), 0)
    if (length(free)) X[free, ] <- rake(seed[free, , drop = FALSE], rm[j, free], cfree)
    o[j, , ] <- X }
  o }
## bases
kg <- readRDS(file.path(S, "kgL_scotland7.rds")); gb <- kg$goodman
base <- list(goodman_raked = { e <- array(0, dim(kc)); for (j in 1:J) e[j, , ] <- rake(gb * rm[j, ], rm[j, ], cm[j, ]); e },
             king = kg$king, ssei_default = readRDS(file.path(S, "n2_sco_base.rds"))$est, ssei_sig08 = readRDS(file.path(S, "p_scotland7_sig08.rds"))$est)
stopifnot(all(sapply(base, function(e) max(abs(apply(e, 1:2, sum) - rm))) < 10))   # ssEI posterior means fit the margins to within a few votes
natprof <- kt[6, ] / sum(kt[6, ])
res <- list()
for (b in names(base)) { e <- base[[b]]
  res[[paste(b, "as estimated")]] <- err(e)
  for (pi in c(0, 0.1, 0.15)) res[[sprintf("%s spoilt pi=%.2f", b, pi)]] <- err(apply_struct(e, function(j) list(`5` = spoilt_row(j, pi))))
  res[[paste(b, "spoilt 0.15 + Other independence")]] <- err(apply_struct(e, function(j) list(`5` = spoilt_row(j, 0.15), `6` = other_ind(j))))
  res[[paste(b, "spoilt 0.15 + Other true national profile (ceiling)")]] <- err(apply_struct(e, function(j) list(`5` = spoilt_row(j, 0.15)),
                                                                              function(j, seed) { seed[6, ] <- rm[j, 6] * natprof; seed })) }
out <- do.call(rbind, res); print(round(out, 2))
cat("\nreference: true national table raked to each area (oracle) area error 3.21\n")
