## Diagnostic for panel_fe_nz.R: does an area's change in a party's candidate share go with a change in that row's loyalty (true tables)?
## usage: Rscript experiments/structure/panel_fe_nz_loyalty_check.R
src <- readLines("experiments/structure/panel_fe_nz.R"); eval(parse(text = src[grep("^suppressMessages", src):(grep("^ipf <-", src) - 1)]))
for (p in list(c(2002, 2005), c(2008, 2011), c(2014, 2017))) {
  k1 <- nz(p[1]); k2 <- nz(p[2]); cm <- intersect(dimnames(k1)[[1]], dimnames(k2)[[1]]); k1 <- k1[cm, , ]; k2 <- k2[cm, , ]
  n1 <- apply(k1, 1:2, sum); n2 <- apply(k2, 1:2, sum); x1 <- n1 / rowSums(n1); x2 <- n2 / rowSums(n2)
  cat(sprintf("%d-%d: correlation of an area's change in a party's candidate share with the change in that row's loyalty, and loyalty of the area's 'new' candidate voters vs all\n", p[1], p[2]))
  for (r in 1:4) { l1 <- k1[, r, r] / pmax(n1[, r], 1); l2 <- k2[, r, r] / pmax(n2[, r], 1); ok <- n1[, r] > 200 & n2[, r] > 200
    ## marginal loyalty: among areas where the row grew, extra candidate votes and extra loyal votes
    g <- ok & n2[, r] > n1[, r] * 1.1
    marg <- sum(k2[g, r, r] - k1[g, r, r]) / sum(n2[g, r] - n1[g, r])
    cat(sprintf("   %-9s corr %.2f | loyalty all %.2f, of the extra candidate votes where the row grew >10%%: %.2f (%d areas)\n", dimnames(k1)[[2]][r],
                cor(x2[ok, r] - x1[ok, r], l2[ok] - l1[ok]), sum(k2[, r, r]) / sum(n2[, r]), marg, sum(g))) }
}
