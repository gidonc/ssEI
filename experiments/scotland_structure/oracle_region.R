## Ceiling on what regional structure could add: rake a TRUE seed to each area's margins, with the national table as the seed
## against each area's own region's table (both oracles: they use the true cells). Area error index 50 * sum|est - true| / N.
x <- readRDS(commandArgs(TRUE)[1]); kc <- x$kc; J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]
rg <- read.csv("experiments/scotland_structure/sco2007_regions.csv")$region; N <- sum(kc)
rm <- apply(kc, 1:2, sum); cm <- apply(kc, c(1, 3), sum)
rake <- function(seed, j) { X <- seed / sum(seed) + 1e-12; X[rm[j, ] == 0, ] <- 0
  for (i in 1:2000) { X <- X * ifelse(rowSums(X) > 0, rm[j, ] / pmax(rowSums(X), 1e-300), 0); X <- sweep(X, 2, ifelse(colSums(X) > 0, cm[j, ] / pmax(colSums(X), 1e-300), 0), "*") }; X }
nat <- apply(kc, 2:3, sum); e_nat <- e_reg <- array(0, dim(kc))
for (j in 1:J) { e_nat[j, , ] <- rake(nat, j); e_reg[j, , ] <- rake(apply(kc[rg == rg[j], , , drop = FALSE], 2:3, sum), j) }
err <- function(e, ix = 1:J) 50 * sum(abs(e[ix, , ] - kc[ix, , ])) / sum(kc[ix, , ])
cat(sprintf("area error: true national table raked %.2f | true regional table raked %.2f\n", err(e_nat), err(e_reg)))
for (g in unique(rg)) cat(sprintf("  %-22s national seed %.2f | regional seed %.2f\n", g, err(e_nat, which(rg == g)), err(e_reg, which(rg == g))))
