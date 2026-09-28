
## --- defaults recovering EXISTING behaviour --------------------------------

## E_rc: old E_rc_hier FALSE -> one group per dimension, all "fixed" at 0
##       old E_rc_hier TRUE  -> one shared group, partial_cp (today's centred hierarchy)
default_E_rc_groups <- function(n_agg_free, E_rc_hier = FALSE, ncp = FALSE) {
  if (!E_rc_hier) {
    list(group_id = seq_len(n_agg_free),
         group_mode = rep("fixed", n_agg_free),      # each dim its own singleton, prior-only
         group_fixed_value = rep(0, n_agg_free))
  } else {
    list(group_id = rep(1L, n_agg_free),
         group_mode = if (ncp) "partial_ncp" else "partial_cp",
         group_fixed_value = 0)
  }
}

## sigma_jrc: old vary_sd FALSE   -> one shared group ("shared")
##            old vary_sd TRUE    -> one group per dim, all "fixed"...
##            -- NO: "vary_sd TRUE" (free/independent) means each dim has its
##            own estimated value with an ORDINARY prior, not a fixed constant.
##            That's the gap flagged earlier: "fixed" bakes in a KNOWN value,
##            but old vary_sd=TRUE wants independent estimation. Use singleton
##            partial_ncp groups for this (a size-1 partial group is
##            auto-collapsed to "shared" by the resolver below, which is
##            EXACTLY "independent, own ordinary prior" for a group of one).
##            old vary_sd "partial" -> one shared pooled group, partial_cp/ncp
default_sigma_groups <- function(Dm1_model, vary_sd = FALSE, ncp = FALSE) {
  if (isFALSE(vary_sd)) {
    list(group_id = rep(1L, Dm1_model), group_mode = "shared", group_fixed_value = 0)
  } else if (isTRUE(vary_sd)) {
    list(group_id = seq_len(Dm1_model),
         group_mode = rep("shared", Dm1_model),       # singleton "shared" == independent own value
         group_fixed_value = rep(0, Dm1_model))
  } else if (identical(vary_sd, "partial")) {
    list(group_id = rep(1L, Dm1_model),
         group_mode = if (ncp) "partial_ncp" else "partial_cp",
         group_fixed_value = 0)
  } else {
    stop("vary_sd must be FALSE, TRUE, or 'partial'")
  }
}

GROUP_MODES <- c(fixed = 0L, shared = 1L, partial_cp = 2L, partial_ncp = 3L, partial_fixed = 4L)

resolve_groups <- function(group_id, group_mode, fixed_value,
                           prior_a, prior_b, tau_a, tau_b, n_dims, label) {
  group_id <- as.integer(group_id)
  stopifnot(length(group_id) == n_dims)
  G <- max(group_id)
  stopifnot(setequal(unique(group_id), seq_len(G)))            # ids must be 1..G, no gaps
  rec <- function(x, nm) {
    if (length(x) == 1) x <- rep(x, G)
    if (length(x) != G) stop(label, ": ", nm, " must have length 1 or ", G)
    x
  }
  group_mode  <- rec(group_mode, "group_mode")
  fixed_value <- rec(fixed_value, "fixed_value")
  prior_a <- rec(prior_a, "prior_a"); prior_b <- rec(prior_b, "prior_b")
  tau_a   <- rec(tau_a, "tau_a");     tau_b   <- rec(tau_b, "tau_b")

  sizes <- tabulate(group_id, G)
  single_partial <- sizes == 1 & group_mode %in% c("partial_cp", "partial_ncp", "partial_fixed")
  if (any(single_partial)) {
    warning(label, ": singleton partial group(s) ", paste(which(single_partial), collapse = ", "),
            " set to SHARED (a singleton cannot be pooled).")
    group_mode[single_partial] <- "shared"
  }
  bad <- group_mode == "partial_fixed" & fixed_value <= 0
  if (any(bad)) stop(label, ": partial_fixed groups need a positive fixed_value (the pooling scale)")

  list(group_id = array(group_id, dim = n_dims), n_groups = G,
       group_mode = array(unname(GROUP_MODES[group_mode]), dim = G),
       group_fixed_value = array(as.numeric(fixed_value), dim = G),
       group_prior_a = array(as.numeric(prior_a), dim = G),
       group_prior_b = array(as.numeric(prior_b), dim = G),
       group_tau_a = array(as.numeric(tau_a), dim = G),
       group_tau_b = array(as.numeric(tau_b), dim = G))
}

## inside build_row_sign_matrix / the diagonal-first builder: record tier per split
## (requires tagging add_split calls with a tier number as they're made -- e.g.
## a global `split_tier` vector paralleling `split_col`)

E_rc_group_id_from_tree <- function(free_dims, split_tier) {
  ## free_dims: which original columns survived reduction (e.g. mr$agg_free_dim)
  ## split_tier: tier label for EVERY column of the original V_ilr, indexed
  ##             the same way; subset to free_dims and renumber consecutively
  raw_tier <- split_tier[free_dims]
  match(raw_tier, sort(unique(raw_tier)))   # consecutive group ids, 1..n_groups
}
