## Symmetric (no margin given priority) version of the sequential margin model: rows are parameters, the hierarchy is on the full table's
## coordinates (row-margin balances, column-margin balances, interaction), b_j ~ N(E_rc, diag(sigma^2)). Rows and columns both enter the
## likelihood (Poisson on the row margins, soft multinomial on the column margins). Uses rw_* helpers for the allocation and starting values.
source("experiments/rowwise_helpers.R")

## margin_reduction list for a full basis (nothing removed): the same fields as ssEI::build_margin_reduction
sym_reduction <- function(V, rm, R, C, eps = 0.5) {
  J <- nrow(rm); D <- R * C; Dm1 <- ncol(V); stopifnot(Dm1 == D - 1)
  sm <- as.matrix(rm) + eps
  out <- list(n_agg_derived = 0L, agg_derived_dim = integer(0), agg_free_dim = as.array(seq_len(Dm1)), Dm1_model = Dm1, n_agg_free = Dm1,
              agg_row_prop = sm / rowSums(sm), V_ilr_full = V, V_ilr_model = V,
              ROT = matrix(0, 0, 0), lflag_rot_agg = 0L, ROT_agg = diag(Dm1), lflag_noncentred_mat = matrix(1L, J, Dm1))
  modifyList(out, ssEI:::margin_param_defaults(Dm1, R, C))
}

## margin parameterisation with rows as parameters: mp_B = d(full coordinates)/d(interior logits) at the Goodman pattern raked to each
## area's margins; whitening of both margins' parameters by the covariance of their log-ratios under the counts
sym_setup <- function(mr, rm, cm, row_order = order(colSums(rm)), rem_col = rep(1L, ncol(rm)), q_bar = rw_goodman(rm, cm)) {
  rm <- as.matrix(rm); cm <- as.matrix(cm); J <- nrow(rm); R <- ncol(rm); C <- ncol(cm); Dm1 <- mr$Dm1_model; K <- (R - 1) * (C - 1)
  V <- mr$V_ilr_model; W <- (rm + 0.01) / (rowSums(rm) + R * 0.01); M <- (cm + 0.01) / (rowSums(cm) + C * 0.01)
  b_of <- function(X) as.vector(crossprod(V, log(as.vector(t(X)))))      # joint table: no division by the row shares
  mp_B <- array(0, c(J, Dm1, K)); h <- 1e-4
  for (j in seq_len(J)) {
    w <- W[j, ]; m <- M[j, ]; X <- w * q_bar
    for (i in 1:1000) { X <- sweep(X, 2, m / colSums(X), "*"); X <- X * (w / rowSums(X)) }
    l0 <- rw_inv(X, row_order, rem_col)
    for (k in seq_len(K)) { e <- numeric(K); e[k] <- h
      mp_B[j, , k] <- (b_of(rw_alloc(w, m, l0 + e, row_order, rem_col)) - b_of(rw_alloc(w, m, l0 - e, row_order, rem_col))) / (2 * h) }
  }
  wh <- function(n) t(chol(diag(1 / n[-length(n)], length(n) - 1) + 1 / n[length(n)]))
  mp_T <- array(0, c(J, C - 1, C - 1)); mp_Tr <- array(0, c(J, R - 1, R - 1))
  for (j in seq_len(J)) { mp_T[j, , ] <- wh(cm[j, ] + 0.01); mp_Tr[j, , ] <- wh(rm[j, ] + 0.01) }
  modifyList(mr, list(
    lflag_margin_param = 1L, lflag_mp_seq = 1L, lflag_mp_seq_anchor = 1L, n_pin = R + C - 2L,
    mp_b_ref = matrix(0, J, Dm1), mp_G = array(0, c(J, R + C - 2, Dm1)), mp_P = array(0, c(J, Dm1, R + C - 2)), mp_N = array(0, c(J, Dm1, K)),
    lflag_seq_or_expected = 3L, mp_row_order = as.array(as.integer(row_order)), mp_rem_col = as.array(as.integer(rem_col)), mp_newton_iters = 10L,
    lflag_mp_beta_whiten = 1L, mp_T = mp_T, lflag_mp_row_whiten = 1L, mp_Tr = mp_Tr,
    lflag_mp_vol_scale = 1L, lflag_mp_seq_scale = 2L, mp_B = mp_B, mp_nc_w = rep(1, Dm1), mp_sigma0 = rep(0.5, Dm1), lflag_mp_scale_fast = 0L))
}

## Fit. V: full orthonormal basis, e.g. rw_symmetric_basis(Vr, Vc) (columns: R - 1 row-margin, C - 1 column-margin, then interaction).
##   sigma_groups: one integer per coordinate (default: one sigma per margin coordinate, one shared by the interaction)
##   E_rc prior: margin coordinates centred on their cross-area mean in the data (sd margin_sd); interaction centred on 0 (independence; sd int_sd)
##   sigma priors (lognormal): margins sigma_prior_margin, interaction sigma_prior_int
sym_fit <- function(kc, V, sigma_groups = NULL, margin_sd = 2, int_sd = 1, sigma_prior_margin = c(log(1), 1), sigma_prior_int = c(log(0.3), 0.5),
                    chains = 2, iter = 1000, warmup = 500, seed = 1234, cores = chains, refresh = max(iter %/% 10, 1), ...) {
  J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]; Dm1 <- R * C - 1; nm <- R + C - 2
  rm <- apply(kc, c(1, 2), sum); cm <- apply(kc, c(1, 3), sum)
  mr <- sym_setup(sym_reduction(V, rm, R, C), rm, cm)
  if (is.null(sigma_groups)) sigma_groups <- c(seq_len(nm), rep(nm + 1L, Dm1 - nm))
  G <- max(sigma_groups); is_m <- seq_len(Dm1) <= nm
  gm <- tapply(is_m, sigma_groups, all)                                  # TRUE = a margin group
  ## margin coordinates from the margins only (independence table of each area): never from the cells
  bt <- t(sapply(seq_len(J), function(j) crossprod(V, log(as.vector(t(outer(rm[j, ] + 0.5, cm[j, ] + 0.5) / sum(rm[j, ])))))))
  E_mean <- ifelse(is_m, colMeans(bt), 0); E_sd <- ifelse(is_m, margin_sd, int_sd)
  pa <- ifelse(gm, sigma_prior_margin[1], sigma_prior_int[1]); pb <- ifelse(gm, sigma_prior_margin[2], sigma_prior_int[2])
  init <- rw_init(J, nm, (R - 1) * (C - 1), E_mean, G, NULL, sigma_prior_int[1])
  init0 <- function() { l <- init(); l$sigma_group_mu <- as.array(ifelse(gm, log(pmax(tapply(apply(bt, 2, sd), sigma_groups, mean), 0.05)), sigma_prior_int[1])); l }
  fit <- ssEI::ei_estimate(
    rm, cm, E_rc_fixed = rep(0, Dm1), sigma_jrc_fixed = rep(.5, Dm1), fix_E_rc = 0, fix_sigma_jrc = 0, E_rc_prior = rep(0, Dm1),
    known_cell_values = kc, use_known_cells = 0, V_ilr = V, n_ilr_rows = R, fit_type = "soft multinom", row_decompose = FALSE,
    rotate_llrep = TRUE, margin_reduction = mr, ROT_E_rc = diag(Dm1), rotate_E_rc = FALSE, link_E_rc = NULL, rotate_lambda = "none",
    sigma_group_id = sigma_groups, sigma_group_mode = rep("shared", G), sigma_group_prior_a = unname(pa), sigma_group_prior_b = unname(pb), sigma_ncp = TRUE,
    E_rc_group_id = 1:Dm1, E_rc_group_mode = rep("shared", Dm1), E_rc_group_prior_a = unname(E_mean), E_rc_group_prior_b = unname(E_sd),
    E_rc_node_logit = FALSE, neutral_logit = "llrep", lambda_centred = FALSE, sigma_floor = 0, prior_sigma_c_scale = 2, prior_lambda_raw_scale = 12,
    prior_gamma_shape = 2, prior_gamma_rate = .5, raw_seq_cell_weights = TRUE,
    chains = chains, cores = cores, refresh = refresh, iter = iter, warmup = warmup, init = init0, seed = seed, ...)
  if (isTRUE(list(...)$return_data)) return(list(data = fit, init = init0))
  invisible(fit)
}
