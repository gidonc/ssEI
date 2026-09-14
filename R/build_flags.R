#' @export
#'
build_flags <- function(Dm1, n_areas,
                        default_agg = 0, default_dev = 0,
                        agg_overrides = integer(0),      # s indices set to 1
                        dev_s_overrides = integer(0),    # s indices: ALL k set to 1
                        dev_ks_overrides = NULL) {        # data.frame(k=, s=) for individual pairs

  lflag_agg <- rep(default_agg, Dm1)
  lflag_agg[agg_overrides] <- 1

  lflag_dev <- matrix(default_dev, nrow = n_areas - 1, ncol = Dm1)
  lflag_noncent_mat <- matrix(default_dev, nrow = n_areas, ncol = Dm1)
  if (length(dev_s_overrides) > 0) lflag_dev[, dev_s_overrides] <- 1
  if (!is.null(dev_ks_overrides) && nrow(dev_ks_overrides) > 0) {
    lflag_dev[cbind(dev_ks_overrides$k, dev_ks_overrides$s)] <- 1
  }
  lflag_noncent_mat[1, 1:Dm1] <- lflag_agg
  lflag_noncent_mat[2:n_areas, 1:Dm1] <- lflag_dev
  list(lflag_agg = lflag_agg, lflag_dev = lflag_dev, lflag_noncent_mat = lflag_noncent_mat)
}
