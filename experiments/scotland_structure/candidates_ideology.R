## Scotland 2007 truth: (1) candidate presence (Green, Other) and the minor-party columns; (2) theory-specified ideology directions for
## the common column shift. Run after truth_structure.R's set-up (sourced lines 1-24).
src <- readLines("experiments/scotland_structure/truth_structure.R"); eval(parse(text = src[1:grep("^cat\\(\"\\\\n== 4", src) - 1]))
gs <- N[, 7] > 0; os <- N[, 6] > 0
cat("\n\n== A. Candidate presence ==\n")
cat(sprintf("  Green candidate stood in %d of %d constituencies; Other candidate(s) in %d\n", sum(gs), J, sum(os)))
cat("  list vote for Green from each major row (%), areas without / with a Green candidate:\n")
for (r in maj) cat(sprintf("    %-4s %5.1f / %5.1f\n", lab[r], 100 * sum(kc[!gs, r, 7]) / sum(N[!gs, r]), 100 * sum(kc[gs, r, 7]) / sum(N[gs, r])))
cat("  list vote for Other from each major row (%), areas without / with an Other candidate:\n")
for (r in maj) cat(sprintf("    %-4s %5.1f / %5.1f\n", lab[r], 100 * sum(kc[!os, r, 6]) / sum(N[!os, r]), 100 * sum(kc[os, r, 6]) / sum(N[os, r])))
cat("  overall list share Green / Other, areas without vs with that candidate:",
    sprintf("Green %.1f vs %.1f; Other %.1f vs %.1f\n", 100 * mean(cmg[!gs, 7]), 100 * mean(cmg[gs, 7]), 100 * mean(cmg[!os, 6]), 100 * mean(cmg[os, 6])))
Sc <- S - rowMeans(S); Sc <- sweep(Sc, 2, colMeans(Sc))
cat(sprintf("  R-squared of the common shift towards Green on 'Green stood': %.2f; towards Other on 'Other stood': %.2f\n",
            summary(lm(Sc[, 7] ~ gs))$r.squared, summary(lm(Sc[, 6] ~ os))$r.squared))
cat("  where do Green list voters' constituency votes go when no Green stood? (column Green, row shares %):",
    paste(sprintf("%s %.0f", lab, 100 * colSums(kc[!gs, , 7]) / sum(kc[!gs, , 7])), collapse = " "), "\n")
cat("  ... and Other list voters when no Other stood:", paste(sprintf("%s %.0f", lab, 100 * colSums(kc[!os, , 6]) / sum(kc[!os, , 6])), collapse = " "), "\n")
cat("\n== B. Ideology directions for the common column shift (7 list columns, centred) ==\n")
v <- function(pos, neg) { x <- numeric(C); x[match(pos, lab)] <- 1 / length(pos); x[match(neg, lab)] <- -1 / length(neg); x - mean(x) }
dirs <- list(independence = v(c("SNP", "Green"), c("Lab", "Con", "LD")), left_right = v(c("Lab", "SNP", "Green"), c("Con", "LD")),
             lab_con = v("Lab", "Con"), ld = v("LD", c("Con", "Lab", "SNP")))
tot <- sum(Sc^2); proj <- function(M) { Q <- qr.Q(qr(M)); sum((Sc %*% Q)^2) / tot }
ev <- eigen(cov(Sc)); cat("  PCs of the common shift: variance shares", round(ev$values[1:5] / sum(ev$values), 2), "\n")
for (k in 1:2) cat(sprintf("   PC%d loadings: %s\n", k, paste(sprintf("%s %+.2f", lab, ev$vectors[, k]), collapse = " ")))
cat("  share of common-shift variation captured by specified directions (best k PCs for comparison):\n")
cat(sprintf("    independence alone %.2f | left-right alone %.2f | Lab-Con alone %.2f | LD alone %.2f\n", proj(cbind(dirs$independence)), proj(cbind(dirs$left_right)), proj(cbind(dirs$lab_con)), proj(cbind(dirs$ld))))
cat(sprintf("    independence + Lab-Con %.2f | + LD %.2f | + Green-own and Other-own column directions %.2f | best 2 PCs %.2f, best 4 PCs %.2f\n",
    proj(cbind(dirs$independence, dirs$lab_con)), proj(cbind(dirs$independence, dirs$lab_con, dirs$ld)),
    proj(cbind(dirs$independence, dirs$lab_con, dirs$ld, v("Green", lab[-7]), v("Other", lab[-6]))), proj(ev$vectors[, 1:2]), proj(ev$vectors[, 1:4])))
sc_ind <- Sc %*% dirs$independence / sum(dirs$independence^2); sc_lc <- Sc %*% dirs$lab_con / sum(dirs$lab_con^2)
cat(sprintf("  area scores: independence vs SNP row share r = %.2f, vs Green stood r = %.2f; Lab-Con vs Lab share r = %.2f, vs Con share r = %.2f\n",
            cor(sc_ind, rs[, 4]), cor(sc_ind, gs), cor(sc_lc, rs[, 2]), cor(sc_lc, rs[, 1])))
