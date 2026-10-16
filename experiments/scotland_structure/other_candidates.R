## Who the 'Other' constituency candidates were (2007), with each candidate's votes and where their voters' list votes went.
utils::data("ei_SCO_2007", package = "ei.Datasets"); x <- ei_SCO_2007
major <- c("Conservative and Unionist Party [The]", "Labour Party [The]", "Liberal Democrats", "Scottish National Party", "Scottish Green Party")
rej <- c("Voting for too many candidates", "Uncertain or Blank", "Writing a mark by which the voter could be identified", "Lack of official mark")
out <- list()
for (k in seq_len(nrow(x))) { cv <- as.data.frame(x$District_cross_votes[[k]]); lists <- cv[[1]]; cands <- setdiff(names(cv)[-1], rej)
  for (cn in cands) { v <- as.numeric(cv[[cn]]); valid <- !lists %in% rej; top <- lists[valid][which.max(v[valid])]
    if (!top %in% major) out[[length(out) + 1]] <- data.frame(constituency = x$District[k], candidate = cn, votes = sum(v), top_list = substr(top, 1, 30),
      top_share = round(100 * max(v[valid]) / sum(v[valid])), to_major = round(100 * sum(v[lists %in% major]) / sum(v[valid]))) } }
o <- do.call(rbind, out); o <- o[order(-o$votes), ]
cat(nrow(o), "Other candidates (top list destination not a major party or Green), total votes", sum(o$votes), "\n")
print(head(o, 15), row.names = FALSE)
cat("\nby top list destination:\n"); print(aggregate(votes ~ top_list, o, function(v) c(n = length(v), votes = sum(v))))
