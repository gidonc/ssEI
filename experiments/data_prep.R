## Data preparation for the Scottish 2007 and New Zealand examples.
##
## prep_scotland(size, n_areas)  -> list(rm, cm, kc, row_names, col_names, bloc, area_idx, actual, s07actual)
## prep_nz(year, size, rows, n_areas) -> same
##
## kc is the array of true cells (areas x rows x columns); it is used only to score the fit.
## Needs the ei.Datasets package (ei_SCO_2007, ei_NZ_<year>) and dplyr / tidyr / purrr.

options(dplyr.summarise.inform = FALSE)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr) })

## ---------------------------------------------------------------------------
## Scotland 2007: constituency vote (rows) by regional list vote (columns)
## ---------------------------------------------------------------------------
## The data do not label the constituency candidates' parties. For each candidate the party is the list-vote party that most of the
## candidate's voters also gave their list vote to; where two candidates in a constituency get the same party, the larger keeps it.
.sco_tables <- function(ei_SCO_2007) {
my_reshape <- function(x){
  x |>
    rename("party_vote" = 1) |>
    pivot_longer(-party_vote, names_to = "cand_vote") |>
    filter(!cand_vote %in% c("Uncertain or Blank",
                             "Voting for too many candidates",
                             "Writing a mark by which the voter could be identified",
                             "Lack of official mark")) |>
    group_by(cand_vote) |>
    arrange(-value) |>
    slice(1)
}

cand_party_proposal <- ei_SCO_2007 |>
  dplyr::select(Number_of_district, District, District_cross_votes) |>
  dplyr::mutate(District_cross_votes = map(District_cross_votes, my_reshape)) |>
  tidyr::unnest(District_cross_votes)


main_parties <- c("Conservative and Unionist Party [The]",
                  "Labour Party [The]",
                  "Liberal Democrats",
                  "Scottish National Party",
                  "Scottish Green Party")

cand_party_proposal2 <- cand_party_proposal |>
  dplyr::group_by(District, party_vote) |>
  dplyr::arrange(Number_of_district, party_vote, -value) |>
  dplyr::mutate(rank = row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(party_vote = ifelse(rank == 1, party_vote, "Other")) |>
  dplyr::mutate(party_vote = ifelse(party_vote %in% main_parties, party_vote, "Other")) |>
  dplyr::select(District, party = party_vote, candidate = cand_vote)


my_cand_reshape <- function(x){
  x |>
    dplyr::rename("party_vote" = 1) |>
    tidyr::pivot_longer(-party_vote, names_to = "candidate", values_to = "votes") |>
    dplyr::mutate(candidate = ifelse(candidate %in% c("Uncertain or Blank",
                                               "Voting for too many candidates",
                                               "Writing a mark by which the voter could be identified",
                                               "Lack of official mark"),
                              "spoilt",
                              candidate)
    )|>
    dplyr::group_by(candidate) |>
    dplyr::summarise(votes = sum(votes, na.rm = TRUE))
}
s07_row_margins_long <- ei_SCO_2007 |>
  dplyr::select(Number_of_district, District, District_cross_votes) |>
  mutate(District_candidate_votes = map(District_cross_votes, my_cand_reshape)) |>
  unnest(District_candidate_votes) |>
  dplyr::select(-District_cross_votes) |>
  left_join(cand_party_proposal2, by = c("candidate", "District")) |>
  mutate(cand_vote = ifelse(is.na(party), "spoilt", party )) |>
  group_by(Number_of_district, District, cand_vote) |>
  summarise(votes = sum(votes, na.rm=TRUE))

s07_row_margins_wide <- s07_row_margins_long |>
  dplyr::select(Number_of_district, District, cand_vote, votes) |>
  pivot_wider(names_from = "cand_vote", values_from="votes", values_fill = 0) |>
  ungroup()


my_party_reshape <- function(x){
  x |>
    rename("party_vote" = 1) |>
    pivot_longer(-party_vote, names_to = "candidate", values_to = "votes") |>
    mutate(party_vote = ifelse(party_vote %in% c("Uncertain or Blank",
                                                 "Voting for too many candidates",
                                                 "Writing a mark by which the voter could be identified",
                                                 "Lack of official mark"),
                               "spoilt",
                               party_vote)
    )|>
    group_by(party_vote) |>
    summarise(votes = sum(votes, na.rm = TRUE))
}

threshold <- .02
s07_col_margins_long <- ei_SCO_2007 |>
  dplyr::select(Number_of_district, District, District_cross_votes) |>
  mutate(District_candidate_votes = map(District_cross_votes, my_party_reshape)) |>
  unnest(District_candidate_votes) |>
  dplyr::select(-District_cross_votes) |>
  mutate(total_votes = sum(votes))|>
  group_by(party_vote) |>
  mutate(n = n(),
         sum_vote = sum(votes, na.rm=TRUE),
         prop_nat_vote = sum_vote/total_votes) |>
  ungroup()|>
  mutate(party_vote = ifelse(prop_nat_vote > threshold,
                             party_vote,
                             "Other")) |>
  group_by(Number_of_district, District, party_vote) |>
  summarise(votes = sum(votes, na.rm=TRUE)) |>
  ungroup()

s07_col_margins_link <- ei_SCO_2007 |>
  dplyr::select(Number_of_district, District, District_cross_votes) |>
  mutate(District_candidate_votes = map(District_cross_votes, my_party_reshape)) |>
  unnest(District_candidate_votes) |>
  dplyr::select(-District_cross_votes) |>
  mutate(total_votes = sum(votes))|>
  group_by(party_vote) |>
  mutate(n = n(),
         sum_vote = sum(votes, na.rm=TRUE),
         prop_nat_vote = sum_vote/total_votes) |>
  ungroup() |>
  mutate(raw_party_vote = party_vote,
         party_vote = ifelse(prop_nat_vote > threshold,
                             party_vote,
                             "Other"))

s07_col_margins_long  <- s07_col_margins_link |>
  group_by(Number_of_district, District, party_vote) |>
  summarise(votes = sum(votes, na.rm=TRUE)) |>
  ungroup()



s07_col_margins_wide <- s07_col_margins_long|> pivot_wider(names_from = "party_vote", values_from = "votes", values_fill = 0) |> ungroup()


s07rm <- s07_row_margins_wide |> dplyr::select(-Number_of_district, -District)
s07cm <- s07_col_margins_wide |> dplyr::select(-Number_of_district, -District)


my_actual_reshape <- function(x, rows_long, cols_long){
  x |>
    rename("party_vote" = 1) |>
    pivot_longer(-party_vote, names_to = "cand_vote", values_to = "votes") |>
    mutate(cand_vote = ifelse(cand_vote %in% c("Uncertain or Blank",
                                               "Voting for too many candidates",
                                               "Writing a mark by which the voter could be identified",
                                               "Lack of official mark"),
                              "spoilt",
                              cand_vote)
    )


}

c_no <- data.frame(col_name = names(s07cm)) |>
  mutate(col_no = row_number())
r_no <- data.frame(row_name = names(s07rm)) |>
  mutate(row_no = row_number())


s07actual <- ei_SCO_2007 |>
  dplyr::select(Number_of_district, District, District_cross_votes) |>
  dplyr::mutate(District_candidate_votes = map(District_cross_votes, my_actual_reshape)) |>
  unnest(District_candidate_votes) |>
  dplyr::select(-District_cross_votes) |>
  left_join(cand_party_proposal2, by = c("cand_vote"="candidate", "District"))  |>
  dplyr::mutate(cand_party = ifelse(is.na(party), "spoilt", party )) |>
  dplyr::mutate(party_vote = ifelse(party_vote %in% c("Uncertain or Blank",
                                               "Voting for too many candidates",
                                               "Writing a mark by which the voter could be identified",
                                               "Lack of official mark"),
                             "spoilt",
                             party_vote)) |>
  left_join(s07_col_margins_link |>
              dplyr::select(District, party_vote_group=party_vote, party_vote=raw_party_vote), by=c("party_vote", "District")) |>
  dplyr::select(Number_of_district, District, col_name = party_vote_group, row_name=cand_party, votes) |>
  group_by(Number_of_district, District, col_name, row_name) |>
  summarise(actual_cell_value = sum(votes)) |>
  ungroup() |>
  mutate(area_no = Number_of_district) |>
  pivot_wider(names_from="row_name", values_from=actual_cell_value, values_fill = 0) |>
  pivot_longer(cols = -c(Number_of_district, District, col_name, area_no),
               names_to = "row_name", values_to = "actual_cell_value") |>
  left_join(c_no) |>
  left_join(r_no) |>
  group_by(area_no, row_no) |>
  mutate(actual_row_margin = sum(actual_cell_value),
         actual_row_rate = actual_cell_value/sum(actual_cell_value))





  list(s07rm = s07rm, s07cm = s07cm, s07actual = s07actual)
}

## 3x3: Lab, SNP, Other (everything else pooled)
collapse_3 <- c(
  "Conservative and Unionist Party [The]" = "Other",
  "Labour Party [The]"                    = "Labour Party [The]",
  "Liberal Democrats"                     = "Other",
  "Scottish National Party"               = "Scottish National Party",
  "Other"                                  = "Other",
  "Scottish Green Party"                   = "Other",
  "spoilt"                                 = "Other"
)

## 5x5: Con, Lab, LibDem, SNP, Other(=Green+Other+spoilt)
collapse_5 <- c(
  "Conservative and Unionist Party [The]" = "Conservative and Unionist Party [The]",
  "Labour Party [The]"                    = "Labour Party [The]",
  "Liberal Democrats"                     = "Liberal Democrats",
  "Scottish National Party"               = "Scottish National Party",
  "Other"                                  = "Other",
  "Scottish Green Party"                   = "Other",
  "spoilt"                                 = "Other"
)


## bloc membership for the collapsed categories (needed by build_diag_first/row_priority)
bloc_3 <- c("Labour Party [The]" = 1, "Scottish National Party" = 2, "Other" = 3)
bloc_5 <- c("Conservative and Unionist Party [The]" = 1, "Labour Party [The]" = 2,
            "Liberal Democrats" = 2, "Scottish National Party" = 3, "Other" = 4)


build_collapsed_margins <- function(rm_full, cm_full, actual_full, collapse_map,
                                    area_idx = NULL) {      # explicit area ids, no hidden sampling
  ## --- 0. optionally subset areas first, so everything downstream is consistent
  if (!is.null(area_idx)) {
    rm_full     <- rm_full[area_idx, , drop = FALSE]
    cm_full     <- cm_full[area_idx, , drop = FALSE]
    ## remap area numbers to 1:n_areas_sub, and keep only the selected areas'
    ## rows in actual_full -- area_no in s07actual must match rownames(rm_full)
    old_to_new  <- setNames(seq_along(area_idx), area_idx)
    actual_full <- actual_full |>
      dplyr::filter(area_no %in% area_idx) |>
      dplyr::mutate(area_no = old_to_new[as.character(area_no)])
  }

  ## --- 1. collapsed margins, as before ---------------------------------------
  collapse_df <- function(df) {
    df2 <- as.data.frame(df)
    names(df2) <- collapse_map[names(df2)]
    groups <- unique(collapse_map)
    out <- sapply(groups, function(g) rowSums(df2[, names(df2) == g, drop = FALSE]))
    as.data.frame(out)
  }
  rm_c <- collapse_df(rm_full)
  cm_c <- collapse_df(cm_full)

  row_names_c <- names(rm_c)
  col_names_c <- names(cm_c)
  r_no_c <- setNames(seq_along(row_names_c), row_names_c)
  c_no_c <- setNames(seq_along(col_names_c), col_names_c)

  ## --- 2. collapse actual_values_wide, by summing within the new groups -----
  actual_wide_c <- actual_full |>
    dplyr::mutate(
      row_name_c = collapse_map[row_name],
      col_name_c = collapse_map[col_name]
    ) |>
    dplyr::group_by(area_no, row_name_c, col_name_c) |>
    dplyr::summarise(actual_cell_value = sum(actual_cell_value, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(row_name = row_name_c, col_name = col_name_c) |>
    dplyr::mutate(
      row_no = r_no_c[row_name],
      col_no = c_no_c[col_name]
    ) |>
    dplyr::group_by(area_no, row_no) |>
    dplyr::mutate(
      actual_row_margin = sum(actual_cell_value),
      actual_row_rate    = actual_cell_value / actual_row_margin
    ) |>
    dplyr::ungroup()

  ## --- 3. rebuild known_cell_values_sco array, same construction as before --
  n_areas_c <- nrow(rm_c)
  known_cell_values_c <- array(0L, dim = c(n_areas_c, length(row_names_c), length(col_names_c)))
  for (i in seq_len(nrow(actual_wide_c))) {
    known_cell_values_c[actual_wide_c$area_no[i],
                        actual_wide_c$row_no[i],
                        actual_wide_c$col_no[i]] <- actual_wide_c$actual_cell_value[i]
  }

  list(
    rm = rm_c, cm = cm_c,
    actual_wide = actual_wide_c,
    known_cell_values = known_cell_values_c,
    row_names = row_names_c, col_names = col_names_c
  )
}
bloc_7 <- c("Conservative and Unionist Party [The]" = 1, "Labour Party [The]" = 2, "Liberal Democrats" = 2, "Other" = 3,
            "spoilt" = 3, "Scottish Green Party" = 4, "Scottish National Party" = 4)

prep_scotland <- function(size = 5, n_areas = 30, seed = 1234) {
  stopifnot(size %in% c(3, 5, 7))
  ei_SCO_2007 <- NULL
  if (requireNamespace("ei.Datasets", quietly = TRUE)) utils::data("ei_SCO_2007", package = "ei.Datasets", envir = environment())
  if (is.null(ei_SCO_2007)) stop("ei_SCO_2007 not found: install the ei.Datasets package")
  tb <- .sco_tables(ei_SCO_2007)
  map <- switch(as.character(size), "3" = collapse_3, "5" = collapse_5, "7" = setNames(names(tb$s07rm), names(tb$s07rm)))   # no pooling
  bloc <- switch(as.character(size), "3" = bloc_3, "5" = bloc_5, "7" = bloc_7)
  set.seed(seed)
  area_idx <- if (n_areas < nrow(tb$s07rm)) sample(seq_len(nrow(tb$s07rm)), n_areas) else seq_len(nrow(tb$s07rm))
  s <- build_collapsed_margins(tb$s07rm, tb$s07cm, tb$s07actual, map, area_idx = area_idx)
  list(rm = as.matrix(s$rm), cm = as.matrix(s$cm), kc = s$known_cell_values,
       row_names = s$row_names, col_names = s$col_names, bloc = bloc, area_idx = area_idx,
       actual = s$actual_wide,          # the true cells in long form, collapsed to this size, areas renumbered 1..n_areas
       s07actual = tb$s07actual)        # the same for all 73 areas and 7 x 7 categories, with the original district numbers
}

## ---------------------------------------------------------------------------
## New Zealand: candidate vote x party vote, by electorate
## ---------------------------------------------------------------------------
## size 6: Labour, National, Green, NZFirst, ACT, Other
## size 4: Labour, National, Green, Other (NZFirst and ACT pooled into Other); only electorates with a candidate in every category
## rows = "candidate" puts the candidate vote on the rows (the party vote is the column); rows = "party" swaps them.
.nz_cat_of <- function(x, cand) {
  p <- if (cand) ifelse(grepl("\\(", x), sub("^.*\\((.*)\\)\\s*$", "\\1", x), x) else x
  p <- sub("\\)$", "", p)
  dplyr::case_when(grepl("Labour", p) ~ "Labour", grepl("^National|New Zealand National", p) ~ "National", grepl("^Green", p) ~ "Green",
                   grepl("New Zealand First|NZ First", p) ~ "NZFirst", grepl("^ACT", p) ~ "ACT", TRUE ~ "Other")
}
prep_nz <- function(year = 2020, size = 6, rows = c("candidate", "party"), n_areas = 30, seed = 1234) {
  rows <- match.arg(rows); stopifnot(size %in% c(4, 6))
  nm <- sprintf("ei_NZ_%d", year)
  env <- new.env(); utils::data(list = nm, package = "ei.Datasets", envir = env); dat <- get(nm, env)
  cats6 <- c("Labour", "National", "Green", "NZFirst", "ACT", "Other")
  to4 <- c(1, 2, 3, 4, 4, 4); cats4 <- c("Labour", "National", "Green", "Other")
  cats <- if (size == 6) cats6 else cats4; K <- length(cats)
  J <- nrow(dat); kc <- array(0, c(J, K, K))                              # [area, candidate vote, party vote]
  for (j in seq_len(J)) {
    x <- as.data.frame(dat$District_cross_votes[[j]])
    pv <- match(.nz_cat_of(as.character(x[[1]]), FALSE), cats6); cv <- match(.nz_cat_of(names(x)[-1], TRUE), cats6)
    if (size == 4) { pv <- to4[pv]; cv <- to4[cv] }
    m <- as.matrix(x[, -1]); m[is.na(m)] <- 0
    for (a in seq_along(pv)) for (b in seq_along(cv)) kc[j, cv[b], pv[a]] <- kc[j, cv[b], pv[a]] + m[a, b]
  }
  ok <- if (size == 4) which(apply(apply(kc, c(1, 2), sum) > 0, 1, all) & apply(apply(kc, c(1, 3), sum) > 0, 1, all)) else seq_len(J)
  set.seed(seed); idx <- if (n_areas < length(ok)) sort(sample(ok, n_areas)) else ok
  kc <- kc[idx, , , drop = FALSE]
  if (rows == "party") kc <- aperm(kc, c(1, 3, 2))
  bloc <- if (size == 6) c(Labour = 1, Green = 1, National = 2, ACT = 2, NZFirst = 3, Other = 3) else c(Labour = 1, National = 2, Green = 1, Other = 3)
  rmm <- apply(kc, c(1, 2), sum); cmm <- apply(kc, c(1, 3), sum); colnames(rmm) <- colnames(cmm) <- cats
  list(rm = rmm, cm = cmm, kc = kc, row_names = cats, col_names = cats, bloc = bloc, area_idx = idx)
}
