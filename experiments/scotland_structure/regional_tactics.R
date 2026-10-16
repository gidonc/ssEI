## Regional list tactics in Scotland 2007, with the d'Hondt position of each party.
## Data: sco_list_2003_2007.csv (Wikipedia region pages; parties below about 2% of the regional list vote omitted, which cannot affect
## the 7 list seats). Check: d'Hondt on list votes with divisors (constituency seats + list seats so far + 1) reproduces the published list seats.
## Viability of party k's list vote in a region = log( v_k / (c_k + l_k + 1) / T ), where T is the quotient that won the 7th list seat:
## the party's next quotient against the seat threshold (0 = exactly at the margin; strongly negative = a list vote for k elects nobody more).
## Expectations from 2003 (what voters could know); 2007 outcome as a comparison (endogenous).
d <- read.csv("experiments/scotland_structure/sco_list_2003_2007.csv")
dhondt <- function(v, c, seats = 7) { l <- rep(0, length(v)); q7 <- NA
  for (s in 1:seats) { q <- v / (c + l + 1); k <- which.max(q); q7 <- q[k]; l[k] <- l[k] + 1 }; list(l = l, T = q7) }
cat("== d'Hondt check: computed list seats against published ==\n"); ok <- TRUE
for (y in c(2003, 2007)) for (g in unique(d$region)) { x <- d[d$year == y & d$region == g, ]; r <- dhondt(x$list_votes, x$const_seats)
  if (!all(r$l == x$list_seats)) { ok <- FALSE; cat("  MISMATCH", y, g, "\n"); print(cbind(x[, c("party", "list_votes", "const_seats", "list_seats")], computed = r$l)) } }
cat("  all 16 region-years reproduce the published allocation:", ok, "\n")
via <- do.call(rbind, lapply(c(2003, 2007), function(y) do.call(rbind, lapply(unique(d$region), function(g) { x <- d[d$year == y & d$region == g, ]; r <- dhondt(x$list_votes, x$const_seats)
  data.frame(year = y, region = g, party = x$party, viab = log(x$list_votes / (x$const_seats + r$l + 1) / r$T), list_seats = r$l, const_seats = x$const_seats) }))))
maj <- c("Con", "Lab", "LD", "SNP")
cat("\n== Viability of each major party's list vote, 2003 (expectation) / 2007 ==\n")
w <- reshape(via[via$party %in% maj, c("year", "region", "party", "viab")], idvar = c("region", "party"), timevar = "year", direction = "wide")
print(xtabs(round(viab.2003, 2) ~ region + party, w)); cat("\n"); print(xtabs(round(viab.2007, 2) ~ region + party, w))
## 2007 loyalty by region (row view) from the true tables
src <- readLines("experiments/scotland_structure/truth_structure.R"); eval(parse(text = src[1:grep("^cat\\(\"== 1", src) - 1]))
rg <- read.csv("experiments/scotland_structure/sco2007_regions.csv"); labm <- c("Con", "Lab", "LD", "SNP")
loy <- do.call(rbind, lapply(unique(rg$region), function(g) { ix <- rg$region == g
  data.frame(region = g, party = labm, loyalty = sapply(1:4, function(k) sum(kc[ix, k, k]) / sum(N[ix, k]))) }))
m <- merge(loy, w); m$lgt <- log(m$loyalty / (1 - m$loyalty))
cat("\n== Does lower list viability go with lower loyalty (constituency voters of k also giving k their list vote)? 32 region-party points ==\n")
for (v in c("viab.2003", "viab.2007")) { f <- lm(as.formula(paste("lgt ~ party +", v)), m); s <- summary(f)$coefficients[v, ]
  cat(sprintf("  logit(loyalty) ~ party + %s: slope %+.3f (se %.3f, t %.1f) | within-party correlation %.2f\n", v, s[1], s[2], s[3],
              cor(resid(lm(lgt ~ party, m)), resid(lm(as.formula(paste(v, "~ party")), m))))) }
print(m[order(m$party, m$viab.2003), c("party", "region", "viab.2003", "viab.2007", "loyalty")], row.names = FALSE, digits = 2)

cat("\n== Robustness: without Glasgow; with region fixed effects (comparing parties within a region) ==\n")
for (v in c("viab.2003", "viab.2007")) {
  f1 <- summary(lm(as.formula(paste("lgt ~ party +", v)), m[m$region != "Glasgow", ]))$coefficients[v, ]
  f2 <- summary(lm(as.formula(paste("lgt ~ party + region +", v)), m))$coefficients[v, ]
  cat(sprintf("  %s: without Glasgow slope %+.3f (t %.1f) | party + region fixed effects slope %+.3f (t %.1f)\n", v, f1[1], f1[3], f2[1], f2[3])) }
cat(sprintf("  size: a viability change of -0.5 (e.g. Labour, South to Glasgow in 2003) moves loyalty by about %.1f points at 80%% (party-only model, 2003)\n",
            100 * (plogis(qlogis(0.8) + coef(lm(lgt ~ party + viab.2003, m))["viab.2003"] * -0.5) - 0.8)))
