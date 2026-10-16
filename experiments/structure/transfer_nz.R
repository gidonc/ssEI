## Transfer of structure between NZ elections (truth known for 2002-2020), Goodman level.
## The margins identify the association parameters of a family badly (families.R), so: learn them from the PREVIOUS election's true table,
## hold them, and fit only the rest (column popularity, and possibly loyalty) from this election's margins.
## Variants (national error of the fitted table | national error of the raked area tables summed | area error after raking the fitted table to each area's margins):
##   Goodman            unrestricted, from margins
##   QI                 loyalty + popularity, from margins
##   prev table         the previous election's true table, unchanged
##   QI+lr: beta held   stated left-right positions; beta from previous truth; loyalty + popularity from margins
##   unfold1: pos held  positions (1D) and their scale from previous truth; loyalty + popularity from margins
##   unfold1: pos+loy   as above, loyalty held too; only popularity from margins
##   QS: pairs held     quasi-symmetric pair terms from previous truth; loyalty + popularity from margins
##   QS: pairs+loy      as above, loyalty held too
## usage: Rscript experiments/structure/transfer_nz.R
commandArgs <- function(...) "none"
suppressMessages(source("experiments/structure/families.R"))
rake_nat <- function(d, P) { e <- array(0, dim(d$kc)); for (j in 1:d$J) { X <- P * d$rm[j, ] + 1e-9; X[d$rm[j, ] == 0, ] <- 0
  for (i in 1:1500) { rs <- rowSums(X); X <- X * ifelse(rs > 0, d$rm[j, ] / pmax(rs, 1e-300), 0); cs <- colSums(X); X <- sweep(X, 2, ifelse(cs > 0, d$cm[j, ] / pmax(cs, 1e-300), 0), "*")
    if (i %% 25 == 0 && max(abs(rowSums(X) - d$rm[j, ])) < 1e-6) break }; e[j, , ] <- X }; 50 * sum(abs(apply(e, 2:3, sum) - apply(d$kc, 2:3, sum))) / d$N }
hold <- function(fam, idx, val) {   # the family with theta[idx] held at val
  free <- setdiff(seq_len(fam$np), idx)
  full <- function(th) { x <- numeric(fam$np); x[idx] <- val; x[free] <- th; x }
  list(np = length(free), eta = function(th) fam$eta(full(th)), start = if (!is.null(fam$start)) fam$start[free])
}
R <- 6; C <- 6; lr <- c(4, 7, 2, 5.5, 8.5, NA); has <- !is.na(lr)
iloy <- C - 1 + seq_len(R); n_c <- C - 1
years <- c(2002, 2005, 2008, 2011, 2014, 2017, 2020)
dat <- lapply(years, function(y) setup(suppressMessages(prep_nz(y, 6, "candidate", 1e6))$kc)); names(dat) <- years
out <- NULL
for (i in 2:length(years)) {
  d <- dat[[i]]; dp <- dat[[i - 1]]
  f_lr <- fam_qi(R, C, dist_term(cbind(lr))); f_un <- fam_qi(R, C, unfold_term(has, 1)); f_qs <- fam_qs(R)
  p_lr <- attr(fit_family(dp, f_lr)$truth, "par"); p_un <- attr(fit_family(dp, f_un)$truth, "par"); p_qs <- attr(fit_family(dp, f_qs)$truth, "par")
  i_lr <- n_c + R + 1; i_un <- n_c + R + seq_len(sum(has)); i_qs <- 2 * R - 1 + seq_len(R * (R - 1) / 2)
  m <- function(f) fit_family(d, f)$margins
  P <- list(Goodman = goodman(d), QI = m(fam_qi(R, C)), `prev table` = dp$Bt,
            `QI+lr: beta held` = m(hold(f_lr, i_lr, p_lr[i_lr])),
            `unfold1: pos held` = m(hold(f_un, i_un, p_un[i_un])),
            `unfold1: pos+loy held` = m(hold(f_un, c(iloy, i_un), p_un[c(iloy, i_un)])),
            `QS: pairs held` = m(hold(f_qs, i_qs, p_qs[i_qs])),
            `QS: pairs+loy held` = m(hold(f_qs, c(iloy, i_qs), p_qs[c(iloy, i_qs)])))
  for (nm in names(P)) out <- rbind(out, data.frame(year = years[i], from = years[i - 1], variant = nm, nat = nat(d, P[[nm]]), nat_raked = rake_nat(d, P[[nm]]), area = rake(d, P[[nm]])))
  cat(sprintf("%d done\n", years[i]))
}
out$variant <- factor(out$variant, names(P))
cat("\nNational error by election (structure learned from the previous election's truth):\n")
print(round(xtabs(nat ~ variant + year, out), 2))
cat("\nNational error of the raked area tables, summed:\n"); print(round(xtabs(nat_raked ~ variant + year, out), 2))
cat("\nArea error after raking:\n"); print(round(xtabs(area ~ variant + year, out), 2))
cat("\nMeans over 2005-2020 (national | area):\n"); mn <- aggregate(cbind(nat, nat_raked, area) ~ variant, out, mean); mn[, -1] <- round(mn[, -1], 2); print(mn, row.names = FALSE)
saveRDS(out, file.path(tempdir(), "transfer_nz.rds"))
