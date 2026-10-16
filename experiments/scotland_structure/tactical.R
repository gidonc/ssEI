## Tactical voting in Scotland 2007 (true tables): is the constituency vote pulled to the locally competitive parties?
## Competitive = the two parties with the largest constituency vote in the area (observable from the row margins).
suppressMessages(source("experiments/data_prep.R")); s <- prep_scotland(7, 73); kc <- s$kc; J <- 73
lab <- c("Con", "Lab", "LD", "SNP", "Spoilt", "Other", "Green"); maj <- 1:4
N <- apply(kc, 1:2, sum); M <- apply(kc, c(1, 3), sum); rs <- N / rowSums(N)
top2 <- t(apply(rs[, maj], 1, function(x) rank(-x) <= 2))                      # J x 4: is party k one of the top two locally
cat("times each party is in the local top two:", paste(lab[maj], colSums(top2)), "\n")
cat("\nA. List voters of party k: share whose constituency vote is also k (by whether k is locally top-two)\n")
for (k in maj) { t2 <- top2[, k]
  f <- function(ix) 100 * sum(kc[ix, k, k]) / sum(M[ix, k])
  cat(sprintf("  %-4s top-two %5.1f%% (%2d areas) | not top-two %5.1f%% (%2d areas)\n", lab[k], f(t2), sum(t2), f(!t2), sum(!t2))) }
cat("\nB. Where list voters go in the constituency vote when their own party is not top-two (or has no candidate):\n")
cat("   share going to the two local front-runners, against the front-runners' share of all constituency votes there\n")
for (k in c(maj, 6, 7)) { pr <- numeric(0); ba <- numeric(0); wt <- numeric(0)
  for (j in 1:J) { if (k %in% maj && top2[j, k]) next
    fr <- maj[top2[j, ]]; tot <- M[j, k]; if (tot < 20) next
    pr <- c(pr, sum(kc[j, fr, k]) / (tot - kc[j, k, k])); ba <- c(ba, sum(N[j, fr]) / (sum(N[j, ]) - N[j, k])); wt <- c(wt, tot) }
  cat(sprintf("  list %-6s (%2d areas): to front-runners %5.1f%% of those not voting for their own candidate | front-runners' share of the other constituency votes %5.1f%%\n",
              lab[k], length(pr), 100 * weighted.mean(pr, wt), 100 * weighted.mean(ba, wt))) }
cat("\nC. What EI estimates (row view): loyalty of constituency row k (share listing k), by whether k is locally top-two\n")
for (k in maj) { t2 <- top2[, k]; f <- function(ix) 100 * sum(kc[ix, k, k]) / sum(N[ix, k])
  cat(sprintf("  %-4s top-two %5.1f%% | not top-two %5.1f%%\n", lab[k], f(t2), f(!t2))) }

cat("\nD. Local pull, fairer comparison: of list party k's voters who do NOT vote for k's candidate, the share going to major party m,\n")
cat("   in areas where m is a local front-runner against areas where it is not (pooled counts; ratio > 1 = pulled to front-runners)\n")
rat <- c()
for (k in c(maj, 6, 7)) { out <- c()
  for (m in setdiff(maj, k)) { d1 <- 0; n1 <- 0; d0 <- 0; n0 <- 0
    for (j in 1:J) { defe <- M[j, k] - (if (k <= 7) kc[j, k, k] else 0); if (defe <= 0) next
      if (top2[j, m]) { d1 <- d1 + kc[j, m, k]; n1 <- n1 + defe } else { d0 <- d0 + kc[j, m, k]; n0 <- n0 + defe } }
    if (n1 > 0 && n0 > 0) { r <- (d1 / n1) / (d0 / n0); rat <- c(rat, log(r)); out <- c(out, sprintf("%s %.1f%% vs %.1f%% (x%.2f)", lab[m], 100 * d1 / n1, 100 * d0 / n0, r)) } }
  cat(sprintf("  list %-6s -> %s\n", lab[k], paste(out, collapse = " | "))) }
cat(sprintf("  geometric mean ratio over all pairs: %.2f\n", exp(mean(rat))))
