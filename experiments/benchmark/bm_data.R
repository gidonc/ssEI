## The elections of the benchmark, each as the true cells (areas x rows x columns) and the margins derived from them.
##   "senc"          eiPack senc, 212 precincts, race (white, black, natam) by vote (dem, rep, non)
##   "scotland7"     Scotland 2007, 73 constituencies, constituency vote (rows) by list vote (columns), all seven categories
##   "nz<year>_<k>"  New Zealand, e.g. "nz2020_6": candidate vote (rows) by party vote (columns), k = 6 or 4 categories, all electorates
## Margins are summed from the true cells, so rows and columns always have the same total in every area.

BM_DIR <- if (exists("BM_DIR")) BM_DIR else "experiments/benchmark"
if (!exists("prep_scotland")) source(file.path(dirname(BM_DIR), "data_prep.R"))

bm_data <- function(name) {
  if (name == "senc") {
    env <- new.env(); utils::data("senc", package = "eiPack", envir = env); senc <- env$senc
    rown <- c("white", "black", "natam"); coln <- c("dem", "rep", "non")
    cellnm <- outer(c("wh", "bl", "natam"), coln, paste0)
    kc <- array(0, c(nrow(senc), 3, 3), dimnames = list(NULL, rown, coln))
    for (r in 1:3) for (c in 1:3) kc[, r, c] <- senc[[cellnm[r, c]]]
    stopifnot(all(apply(kc, c(1, 2), sum) == as.matrix(senc[rown])), all(apply(kc, c(1, 3), sum) == as.matrix(senc[coln])))
    bloc <- NULL
  } else if (name == "scotland7") {
    d <- prep_scotland(7, 73); kc <- d$kc; dimnames(kc) <- list(NULL, d$row_names, d$col_names); bloc <- d$bloc
  } else if (grepl("^nz[0-9]{4}_[46]$", name)) {
    year <- as.integer(substr(name, 3, 6)); size <- as.integer(substring(name, 8))
    d <- prep_nz(year, size, "candidate", n_areas = 1e6); kc <- d$kc; dimnames(kc) <- list(NULL, d$row_names, d$col_names); bloc <- d$bloc
  } else stop("unknown election: ", name)
  stopifnot(all(kc >= 0), all(kc == round(kc)))
  rm <- apply(kc, c(1, 2), sum); cm <- apply(kc, c(1, 3), sum)
  list(name = name, kc = kc, rm = rm, cm = cm, row_names = dimnames(kc)[[2]], col_names = dimnames(kc)[[3]], bloc = bloc,
       n_areas = dim(kc)[1], n_empty_rows = sum(rm == 0), n_empty_cols = sum(cm == 0))
}
