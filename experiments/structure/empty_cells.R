## Empty cells in the true senc tables against the number each model expects.
## The small-cell term is, to a close approximation, a prior that down-weights each empty cell of the realised table (see
## claude/ssei_small_cell_justification.md), so its implicit claim is "empty cells are rarer than the hierarchy implies". This checks it.
## Expected empties for a row with n voters and rates q: sum_c (1 - q_c)^n (multinomial). Rates: each fit's posterior-mean area table,
## or the true pooled rates of the row applied everywhere (no contextual variation, right on average).
suppressMessages(source("experiments/rowwise_helpers.R")); data(senc, package = "eiPack")
rown <- c("white", "black", "natam"); coln <- c("dem", "rep", "non"); kc <- array(0, c(nrow(senc), 3, 3), dimnames = list(NULL, rown, coln))
cellnm <- outer(c("wh", "bl", "natam"), coln, paste0); for (r in 1:3) for (c in 1:3) kc[, r, c] <- senc[[cellnm[r, c]]]
rm <- apply(kc, 1:2, sum); pooled <- apply(kc, 2:3, sum) / colSums(rm)
exp0 <- function(q, r) { n <- rm[, r]; ok <- n > 0; colSums((1 - q[ok, , drop = FALSE])^n[ok]) }
rates <- function(tab, r) { x <- tab[, r, ] + 1e-9; x / rowSums(x) }
vs <- commandArgs(TRUE)
for (r in 2:3) { ok <- rm[, r] > 0
  o <- rbind("observed (truth)" = colSums(kc[ok, r, ] == 0), "true pooled rates everywhere" = exp0(matrix(pooled[r, ], nrow(kc), 3, byrow = TRUE), r))
  for (v in vs) { f <- sprintf("experiments/structure/results/sign_context_senc_%s.rds", v); if (file.exists(f)) { o <- rbind(o, exp0(rates(readRDS(f)$est, r), r)); rownames(o)[nrow(o)] <- v } }
  cat(sprintf("\n%s row: empty cells among the %d precincts where the row is present\n", rown[r], sum(ok))); print(round(o, 1)) }
