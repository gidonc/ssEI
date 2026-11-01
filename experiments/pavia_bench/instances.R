## Per-district EI instances in the ei.Datasets / Pavia convention: one problem per NZ electorate or Scottish constituency, units = its
## polling places, rows = party (list) votes, columns = candidate (constituency) votes, informal categories included; parties and
## candidates with under 3% of the district's votes merged into one "Others" row / column. Truth = the district's cross table.
## pavia_instances(name) returns a list of list(name, rm = units x R, cm = units x C, truth = R x C).
pavia_instances <- function(dataset = "ei_NZ_2020", thresh = 0.03) {
  env <- new.env(); utils::data(list = dataset, package = "ei.Datasets", envir = env); x <- get(dataset, env)
  lapply(seq_len(nrow(x)), function(j) {
    P <- as.matrix(x$Votes_to_parties[[j]][, -(1:2)]); Cn <- as.matrix(x$Votes_to_candidates[[j]][, -(1:2)])
    P[is.na(P)] <- 0; Cn[is.na(Cn)] <- 0
    cv <- as.data.frame(x$District_cross_votes[[j]]); Tr <- as.matrix(cv[, -1]); Tr[is.na(Tr)] <- 0
    stopifnot(nrow(Tr) == ncol(P), ncol(Tr) == ncol(Cn), max(abs(colSums(P) - rowSums(Tr))) < 1, max(abs(colSums(Cn) - colSums(Tr))) < 1)
    tot <- sum(Tr); kr <- colSums(P) / tot >= thresh; kcl <- colSums(Cn) / tot >= thresh
    agg <- function(M, keep) { M2 <- M[, keep, drop = FALSE]; if (any(!keep)) M2 <- cbind(M2, Others = rowSums(M[, !keep, drop = FALSE])); M2 }
    rm <- agg(P, kr); cm <- agg(Cn, kcl)
    Tr2 <- t(agg(t(agg(Tr, kcl)), kr))
    ok <- rowSums(rm) > 0 & abs(rowSums(rm) - rowSums(cm)) < 0.5
    list(name = as.character(x$District[j]), rm = rm[ok, , drop = FALSE], cm = cm[ok, , drop = FALSE], truth = Tr2)
  })
}
ei_index <- function(est, truth) 50 * sum(abs(est - truth)) / sum(truth)
