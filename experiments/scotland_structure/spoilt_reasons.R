## Scotland 2007: rejected ballots by reason, constituency section (rows) x list section (columns), from the full ei_SCO_2007 cross votes.
utils::data("ei_SCO_2007", package = "ei.Datasets"); x <- ei_SCO_2007
rej <- c("Voting for too many candidates", "Uncertain or Blank", "Writing a mark by which the voter could be identified", "Lack of official mark")
short <- c("too many", "blank/uncertain", "id mark", "no official mark")
cls <- function(n) ifelse(n %in% rej, short[match(n, rej)], "valid")
T <- matrix(0, 5, 5, dimnames = list(const = c(short, "valid"), list = c(short, "valid")))
for (k in seq_len(nrow(x))) { cv <- as.data.frame(x$District_cross_votes[[k]]); lists <- cv[[1]]; cands <- names(cv)[-1]
  for (ci in seq_along(cands)) { v <- as.numeric(cv[[ci + 1]]); a <- cls(cands[ci]); for (li in seq_along(lists)) { b <- cls(lists[li]); T[a, b] <- T[a, b] + v[li] } } }
cat("All ballots, constituency section (rows) x list section (columns):\n"); print(T)
tot <- sum(T); cat(sprintf("\nrejected constituency: %d (%.2f%%) | rejected list: %d (%.2f%%) | both rejected: %d\n",
  sum(T[1:4, ]), 100 * sum(T[1:4, ]) / tot, sum(T[, 1:4]), 100 * sum(T[, 1:4]) / tot, sum(T[1:4, 1:4])))
cat(sprintf("of rejected constituency votes, %.0f%% had a valid list vote; of rejected list votes, %.0f%% had a valid constituency vote\n",
  100 * sum(T[1:4, 5]) / sum(T[1:4, ]), 100 * sum(T[5, 1:4]) / sum(T[, 1:4])))
