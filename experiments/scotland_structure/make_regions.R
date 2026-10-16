## Constituency -> region map for Scotland 2007, from House of Commons Library RP07-46, Table 4 (constituency results by region).
## The pasted table text is parsed; constituencies are matched to the ei_SCO_2007 areas by their constituency vote counts (Con, Lab, LD, SNP).
## Writes experiments/scotland_structure/sco2007_regions.csv (area index in prep_scotland(7, 73) order, district, region, majority).
a <- commandArgs(TRUE); txt <- paste(readLines(a[1], warn = FALSE), collapse = " ")
txt <- gsub("18 RESEARCH PAPER 07/46 CON LAB LD SNP Others Total Electorate Turnout Majority", " ", txt)
regions <- c("Central Scotland", "Glasgow", "Highlands and Islands", "Lothians", "Mid Scotland and Fife", "North East Scotland", "South of Scotland", "West of Scotland")
num <- "([0-9,]+)"; pat <- paste0("([A-Z][A-Za-z&,'. -]+?) ", paste(rep(num, 7), collapse = " "), " ([0-9.]+)%( [0-9,]+)?")
m <- regmatches(txt, gregexpr(pat, txt, perl = TRUE))[[1]]
rows <- lapply(m, function(s) { p <- regmatches(s, regexec(pat, s, perl = TRUE))[[1]]; n <- as.numeric(gsub(",", "", p[3:9]))
  data.frame(name = trimws(p[2]), con = n[1], lab = n[2], ld = n[3], snp = n[4], oth = n[5], total = n[6], majority = as.numeric(gsub(",", "", trimws(p[11])))) })
d <- do.call(rbind, rows); d$name <- sub("^.*Turnout Majority ", "", d$name)
d$is_region <- d$name %in% regions; reg <- NA; d$region <- NA
for (i in seq_len(nrow(d))) { if (d$is_region[i]) reg <- d$name[i]; d$region[i] <- reg }
cons <- d[!d$is_region, ]; cat(nrow(cons), "constituencies parsed in", length(unique(cons$region)), "regions\n")
## check: the region rows equal the sums of their constituencies
chk <- aggregate(cbind(con, lab, ld, snp) ~ region, cons, sum); rr <- d[d$is_region, c("name", "con", "lab", "ld", "snp")]
cat("region totals match the constituency sums:", isTRUE(all.equal(as.matrix(chk[match(rr$name, chk$region), -1]), as.matrix(rr[, -1]), check.attributes = FALSE)), "\n")
suppressMessages(source("experiments/data_prep.R")); s <- prep_scotland(7, 73); rm <- apply(s$kc, 1:2, sum)
key_d <- paste(rm[, 1], rm[, 2], rm[, 3], rm[, 4]); key_t <- paste(cons$con, cons$lab, cons$ld, cons$snp)
mt <- match(key_d, key_t); cat("areas matched exactly on Con/Lab/LD/SNP constituency votes:", sum(!is.na(mt)), "of", nrow(rm), "\n")
if (anyNA(mt)) { for (j in which(is.na(mt))) { dd <- rowSums(abs(sweep(as.matrix(cons[, c("con", "lab", "ld", "snp")]), 2, rm[j, 1:4]))); o <- order(dd)
  cat(sprintf("  area %d nearest: %s (abs diff %d; next nearest %s, diff %d)\n", j, cons$name[o[1]], dd[o[1]], cons$name[o[2]], dd[o[2]]))
  if (dd[o[1]] <= 50 && dd[o[2]] > 20 * max(dd[o[1]], 1)) mt[j] <- o[1] } }
stopifnot(!anyNA(mt), !any(duplicated(mt)))
out <- data.frame(area = seq_len(nrow(rm)), constituency = cons$name[mt], region = cons$region[mt], majority = cons$majority[mt])
write.csv(out, "experiments/scotland_structure/sco2007_regions.csv", row.names = FALSE); print(table(out$region))
