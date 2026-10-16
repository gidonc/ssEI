## Scotland 2007 truth: (1) does the spoilt (constituency) row split like the area's list vote? (2) who the Other candidates were.
x <- readRDS(commandArgs(TRUE)[1]); kc <- x$kc; J <- dim(kc)[1]; lab <- c("Con", "Lab", "LD", "SNP", "Spoilt", "Other", "Green")
N <- apply(kc, 1:2, sum); M <- apply(kc, c(1, 3), sum)
cat("== 1. Spoilt constituency row: true national split against 'split like the area's list vote' (independence within each area) ==\n")
tr <- colSums(kc[, 5, ]) / sum(N[, 5])
ind <- colSums(t(sapply(1:J, function(j) N[j, 5] * M[j, ] / sum(M[j, ])))) / sum(N[, 5])
## variant: independence for the valid list columns, with the spoilt-spoilt cell taken at its national rate (a 'fully spoilt' block)
print(round(100 * rbind(truth = tr, independence = ind), 1))
cat(sprintf("  error index for the spoilt row (50 x sum|est - true| / row total): independence %.1f | Goodman (earlier) was about 71 points on Other alone\n",
            50 * sum(abs(sapply(1:J, function(j) N[j, 5] * M[j, ] / sum(M[j, ])) - t(kc[, 5, ]))) / sum(N[, 5])))
cat(sprintf("  area by area: correlation of the spoilt row's Lab share with the area's Lab list share %.2f; SNP %.2f\n",
  cor(kc[N[, 5] > 200, 5, 2] / N[N[, 5] > 200, 5], M[N[, 5] > 200, 2] / rowSums(M[N[, 5] > 200, ])), cor(kc[N[, 5] > 200, 5, 4] / N[N[, 5] > 200, 5], M[N[, 5] > 200, 4] / rowSums(M[N[, 5] > 200, ]))))
cat("\n== 2. The Other constituency row: true national split, and against independence ==\n")
tro <- colSums(kc[, 6, ]) / sum(N[, 6]); indo <- colSums(t(sapply(1:J, function(j) N[j, 6] * M[j, ] / sum(M[j, ])))) / sum(N[, 6])
print(round(100 * rbind(truth = tro, independence = indo), 1))
