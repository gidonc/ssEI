## slopes implied by each fit's posterior-mean area cells: within-row coordinate on standardised own-row log share (weighted by row size)
suppressMessages(source("experiments/rowwise_helpers.R")); data(senc, package = "eiPack")
rown <- c("white", "black", "natam"); coln <- c("dem", "rep", "non"); kc <- array(0, c(nrow(senc), 3, 3), dimnames = list(NULL, rown, coln))
cellnm <- outer(c("wh", "bl", "natam"), coln, paste0); for (r in 1:3) for (c in 1:3) kc[, r, c] <- senc[[cellnm[r, c]]]
Vr <- rw_sbp_to_v(cbind(c(1, 1, -1), c(1, -1, 0))); rm <- apply(kc, 1:2, sum)
sl <- function(tab, minrow = 20) unlist(lapply(1:3, function(r) { ok <- rm[, r] >= minrow; l <- log(rm[ok, r] / rowSums(rm)[ok]); x <- (l - mean(l)) / sd(l)
  sapply(1:2, function(k) coef(lm(as.vector(log(tab[ok, r, ] + 0.5) %*% Vr[, k]) ~ x, weights = rm[ok, r]))[2]) }))
vs <- commandArgs(TRUE); out <- rbind(truth = sl(kc))
for (v in vs) { f <- sprintf("experiments/structure/results/sign_context_senc_%s.rds", v); if (file.exists(f)) out <- rbind(out, setNames(sl(readRDS(f)$est), NULL)); rownames(out)[nrow(out)] <- v }
colnames(out) <- c("wh:reg", "wh:dem-rep", "bl:reg", "bl:dem-rep", "na:reg", "na:dem-rep"); print(round(out, 2))
## rates of each row by its share band: truth and fits
band <- function(tab, r) { s <- rm[, r] / rowSums(rm); b <- cut(s, c(0, .05, .2, .5, 1.01), right = FALSE); sapply(levels(b), function(g) { i <- b == g & rm[, r] > 0; sum(tab[i, r, 2]) / sum(tab[i, r, 1:2]) }) }
for (r in 1:3) { cat("\nRepublican share of", rown[r], "party registrants, by the row's share of the precinct\n"); o <- rbind(truth = band(kc, r))
  for (v in vs) { f <- sprintf("experiments/structure/results/sign_context_senc_%s.rds", v); if (file.exists(f)) { o <- rbind(o, band(readRDS(f)$est, r)); rownames(o)[nrow(o)] <- v } }
  print(round(o, 3)) }
