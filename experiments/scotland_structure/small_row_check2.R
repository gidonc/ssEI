## Cheap check, part 2: quasi-independence for the small rows, with no candidate typing.
## Row r in area j: a share pi_r to its own list column ("loyalty"); the rest split like the area's list vote over the other valid columns.
## Applied to spoilt (own = spoilt column), Other (own = Other list column) and Green (own = Green list). The pi values are a grid;
## true national own-column shares: spoilt 0.156, Other 0.447, Green 0.399.
src <- readLines("experiments/scotland_structure/small_row_check.R"); eval(parse(text = src[1:grep("^natprof <-", src)]))
qi_row <- function(j, r, own, pi) { s <- rm[j, r]; v <- cm[j, ]; o <- min(pi * s, v[own]); vv <- v; vv[own] <- 0; if (r != 5) vv[5] <- v[5]
  out <- (s - o) * vv / sum(vv); out[own] <- out[own] + o; out }
specs <- list(
  "spoilt 0.15" = function(j) list(`5` = qi_row(j, 5, 5, 0.15)),
  "spoilt 0.15 + Other 0.30" = function(j) list(`5` = qi_row(j, 5, 5, 0.15), `6` = qi_row(j, 6, 6, 0.30)),
  "spoilt 0.15 + Other 0.45" = function(j) list(`5` = qi_row(j, 5, 5, 0.15), `6` = qi_row(j, 6, 6, 0.45)),
  "spoilt 0.15 + Other 0.45 + Green 0.40" = function(j) list(`5` = qi_row(j, 5, 5, 0.15), `6` = qi_row(j, 6, 6, 0.45), `7` = qi_row(j, 7, 7, 0.40)),
  "spoilt 0.10 + Other 0.30 + Green 0.30 (rougher guesses)" = function(j) list(`5` = qi_row(j, 5, 5, 0.10), `6` = qi_row(j, 6, 6, 0.30), `7` = qi_row(j, 7, 7, 0.30)))
res <- list()
for (b in c("goodman_raked", "king", "ssei_default")) { e <- base[[b]]; res[[paste(b, "| as estimated")]] <- err(e)
  for (nm in names(specs)) res[[paste(b, "|", nm)]] <- err(apply_struct(e, specs[[nm]])) }
out <- do.call(rbind, res)
gr <- function(e) 50 * sum(abs(e[, 7, ] - kc[, 7, ])) / sum(kc[, 7, ])
print(round(out, 2))
