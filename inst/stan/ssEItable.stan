//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// with various models of row to column rate

functions{
  #include include/allocationfuns.stan
  #include include/realpdf.stan


  vector simplex_constrain_softmax_lp(vector v) {
     int K = size(v) + 1;
     vector[K] v0 = append_row(0, v);
     return softmax(v0);
  }

  // margin param (seq): sequential allocation of an R x C proportion table with row margins w,
  // column margins m and interior logits lam ((R-1)*(C-1), row-major). Same scheme as ss_assign_cvals_lp.
  // Returns (R+1) x C: rows 1:R = table, element [R+1, 1] = log|d interior cells / d lam|,
  // [R+1, 2] = largest share of a cell's true range lost to kink smoothing, [R+1, 3] = number of cells losing more than 1%.
  matrix mp_seq_alloc(vector w, vector m, vector lam, real delta_floor, real delta_min, real delta_rel, int use_or) {
    int R = rows(w);
    int C = rows(m);
    vector[R] sr = w;
    vector[C] sc = m;
    matrix[R + 1, C] out = rep_matrix(0, R + 1, C);
    real lj = 0;
    for (r in 1:(R - 1)) {
      for (c in 1:(C - 1)) {
        real x = lam[(r - 1) * (C - 1) + c];
        if (use_or == 1) {
          real rest = sum(sc[(c + 1):C]);
          out[r, c] = seq_or_cell(sr[r], sc[c], rest, x);
          lj += seq_or_logjac(out[r, c], sr[r], sc[c], rest);
        } else {
          vector[3] bw = seq_cell_bounds(sr[r], sc[c], sum(sc[(c + 1):C]), delta_rel, delta_floor, delta_min);
          real lo = bw[1];
          real width = bw[2];
          real loss = seq_cell_range_loss(sr[r], sc[c], sum(sc[(c + 1):C]), width);
          if (C >= 3) { out[R + 1, 2] = fmax(out[R + 1, 2], loss); out[R + 1, 3] += (loss > 0.01); }
          out[r, c] = lo + inv_logit(x) * width;
          lj += log(width) + log_inv_logit(x) + log1m_inv_logit(x);
        }
        sc[c] -= out[r, c];
        sr[r] -= out[r, c];
      }
      out[r, C] = sr[r];
      sc[C] -= sr[r];
    }
    for (c in 1:C) out[R, c] = sc[c];
    out[R + 1, 1] = lj;
    return out;
  }

  // margin param: implied column log-ratios s (C-1) of a row-conditioned table with row weights w,
  // and their gradient wrt b (V_ilr_model coords). Returns (C-1) x (Dm+1): col 1 = s, cols 2:(Dm+1) = ds/db
  matrix mp_colshare(vector b, vector w, matrix V, int R, int C) {
    int Dm = cols(V);
    vector[C] S = rep_vector(0, C);
    matrix[C, Dm] dS = rep_matrix(0, C, Dm);
    matrix[C - 1, Dm + 1] out;
    for (r in 1:R) {
      matrix[C, Dm] Vr = V[((r - 1) * C + 1):(r * C), ];
      vector[C] q = softmax(Vr * b);
      row_vector[Dm] qV = q' * Vr;
      S += w[r] * q;
      for (c in 1:C) dS[c] += w[r] * q[c] * (Vr[c] - qV);
    }
    for (c in 1:(C - 1)) {
      out[c, 1] = log(S[c]) - log(S[C]);
      out[c, 2:(Dm + 1)] = dS[c] / S[c] - dS[C] / S[C];
    }
    return out;
  }

  vector sum_to_zero(vector y) {
    int N = num_elements(y);
    vector[N + 1] x = zeros_vector(N + 1);
    real sum_w = 0;
    for (n in 1:N) {
        int i = N - n + 1;
        real w = y[i] * inv_sqrt(i * (i + 1));
        sum_w += w;
        x[i] += sum_w;
        x[i + 1] -= w * i;
    }
    return x;
}

}
data{
 int<lower=0> n_areas;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
 matrix<lower=0>[n_areas, C] col_margins; // the column margins in each area
  int<lower = 0> n_agg_derived;
  array[n_agg_derived] int<lower = 1, upper=R*C - 1> agg_derived_dim;
  int<lower = 0> n_agg_free;
  array[n_agg_free] int<lower = 1, upper=R * C - 1> agg_free_dim;
  int<lower=1> Dm1_model;
  matrix[(R * C), (R * C) - 1] V_ilr_full; // basis matrix for ILR transformation
  matrix[(R * C), Dm1_model] V_ilr_model;
  matrix[n_areas*Dm1_model + n_areas - 1, n_areas*Dm1_model + n_areas - 1] ROT_red; //rebuilt for reduced dimensions
  matrix[n_agg_free, n_agg_free] ROT_E_rc; //ilr rotation from raw to E_rc in V_ilr basis

 // vector<lower=0>[(R * C) - 1] sigma_llrep;
 // vector[(R * C) - 1] mu_llrep;
 int<lower=0, upper=1> structural_zeros[n_areas, R, C];  // an array indicating any structural zeros in the data (may include whole rows, whole columns and/or individual cells)
 int<lower=0, upper=2> lflag_dist; // flag indicating whether to use poisson (0), multinomial (1) or negative binomial (2)paramertization
 int<lower=0, upper=2> lflag_family; // flag indicating whether scale of LLrep distribution is based on log-normal (0), cauchy(1) or Gamma (2) family
 int<lower=0, upper=3> lflag_area_re; // flag indicating whether the area mean simplex is uniform (0) or varies with area random effects which are normally distributed (1) or varies with area random effects which are multinormally distributed (non centred paramaterisation) (2) or varies with area random effects which are multinormally distributed (non centred LKJ Onion paramaterisation)
 int<lower  =0, upper=2> lflag_vary_sd; // flag indicating whether variance of area_cell parameters is: (0) shared across cells,  (1) varies by cell,  or (2) has a hierarchical model structure
 int<lower = 0, upper=3> lflag_neutral_logit; // flag indicating the neutral_logit for the allocation process when sequential weights are all zero. (0) gives equality across the rows [-log(cols_remaining)] (1) gives the independent table logit(col_margin/sum_remaining_col_margins) (2) lambda implied by LLrep structure (3) lambda implied by E_rc structure.
 int<lower = 0, upper = 1> lflag_llmod_omit_jr; // flag indicating whether log-linear model should omit area * row interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jc; // flag indicating whether log-linear model should omit area * col interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jrc; // flag indicating whether log-linear model should omit area * row * column interaction
 int<lower = 0, upper = 1> lflag_E_rc_node_logit; // (0) balance paramertisation (1) logit parameterisation
  int<lower =0, upper = 1> lflag_pin_row_effect; // flag indicating the model on row effects
  int<lower = 0, upper =1> lflag_predictors_cm; // flag indicating whether to model columns as well as rows
  int<lower=0, upper=1> lflag_row_decompose;
  int<lower =0, upper = 2> lflag_noncentred; //flag indicating whether to use a centred (0) or non-centred parameterization;
  int<lower=0, upper = 1> lflag_noncentred_mat[n_areas, Dm1_model]; //flag indicated whether cell is centred (0) or non-centred (1)
  int<lower = 0, upper = 1> lflag_rawscw; //flag indicating whether to use raw, or decentred version of sequential cell weights (lamdba coeffients)
  int<lower = 0, upper = 4> lflag_ll_rep; // flag indicating which log linear representation of the final tables should be used.
  int<lower = 0, upper = 1> lflag_rot_llrep;
  int<lower=0, upper=1> lflag_rot_agg;
  matrix[Dm1_model, Dm1_model] ROT_agg;          // orthogonal; identity = initial behaviour

  // ---- margin parameterisation (per-area pinned coords beta + NCP interaction z) ----
  // b_j = mp_b_ref[j] + mp_P[j] * beta_j + mp_N[j] * gamma_j      (b_j = llrep_model_j, V_ilr_model coords)
  // mp_G[j] : n_pin x Dm1 gradient of pinned margin log-ratios at b_ref ; G*P = I, G*N = 0, N'N = I
  int<lower=0, upper=1> lflag_margin_param;
  int<lower=0, upper=1> lflag_mp_seq;            // margin param: 1 = expected table built by sequential allocation: mp_beta = column log-ratio deviation from observed, mp_z = interior logits (centred). Works with rows conditioned or free. Overrides lflag_mp_exact / lflag_mp_gamma_centred; mp_G, mp_P, mp_N, mp_b_ref are unused (pass zeros).
  int<lower=0, upper=1> lflag_mp_seq_anchor;     // margin param (seq): 1 = interior logits = (logits that realise E_rc for this area's rows) + mp_z, i.e. mp_z is the area's deviation from E_rc; 0 = mp_z are the logits themselves
  int<lower=0, upper=1> lflag_mp_exact;          // margin param: 0 = mp_beta is the linearised pinned coord; 1 = mp_beta is the EXACT column log-ratio deviation (Newton solve + Jacobian)
  int<lower=1> mp_newton_iters;                  // margin param: Newton iterations when lflag_mp_exact == 1 (5 is plenty)
  int<lower=0, upper=1> lflag_mp_gamma_centred;  // margin param: 0 = interior non-centred (mp_z = z), 1 = centred (mp_z = gamma)
  int<lower=0> n_pin;
  array[R] int<lower=1, upper=R> mp_row_order;  // sequential expected table, modes 2 and 3: allocation order of the rows (last entry = remainder row)
  array[R] int<lower=1, upper=C> mp_rem_col;    // sequential expected table, modes 2 and 3: for each row, the column found by subtraction / used as reference
  int<lower=0, upper=1> lflag_mp_vol_scale;     // margin param: 1 = log_volume[j] = log(N_j) + (volume parameter) / sqrt(N_j), so the volume parameters are unit scale (N_j = area total)
  int<lower=0, upper=2> lflag_mp_scale_fast;   // seq scale 2 with shared/fixed sigma groups and mp_nc_w = 1. 1 = scaling matrix as a weighted sum of one fixed matrix per sigma group (same map as 0). 2 = low-rank form: the largest sigma group gives a fixed factor, the other groups a small update, so only a small Cholesky is needed per area (same density; mp_z is a rotation of the mode 0/1 mp_z)
  int<lower=0, upper=1> lflag_real_rowwise;    // HYBRID: 1 = the REALISED table is allocated row by row on the observed margins (same scheme and order as the expected table; needs lflag_mp_seq = 1, lflag_seq_or_expected = 3, lflag_neutral_logit = 2). Its logits = the expected table's logits + rl_S * lambda_raw. Empty rows are dropped; empty columns are not supported.
  int<lower=0, upper=1> lflag_real_scale;      // HYBRID row by row: 1 = lambda_raw is also scaled by the Poisson sd of each logit at the CURRENT expected table, sqrt(1/x + 1/rest + 1/x_ref + 1/rest_ref) (x = the cell, rest = the later rows' cells in the same column, ref = the row's reference column); rl_S must then be built on logits divided by the same quantity at the reference table. Removes the dependence of lambda_raw's scale on the expected cell sizes.
  int<lower=0, upper=1> lflag_cpois_norm;      // HYBRID row by row: 1 = divide each realised cell's continuous Poisson density by its integral over x >= 0, so it is a proper density whatever the mean (cells in non-empty rows only; cells of empty rows are fixed at 0 and keep the Poisson probability of 0)
  int<lower=0, upper=1> lflag_soft_smallcell;  // SOFT (fit type 2): 1 = add the log integral of the continuous Poisson density for every expected cell in a non-empty row, i.e. the same implicit prior against very small expected cells that HYBRID has when lflag_cpois_norm = 0 (uses the cz_ grid)
  real<lower=0> smallcell_scale;               // multiplier on the SOFT small-cell term (1 = the HYBRID-equivalent strength)
  int<lower=0, upper=2> lflag_mp_interior;     // sequential margin param: what areas may do beyond moving whole rows and columns. 0 = free interior, the hierarchy b_j ~ N(E_rc, diag(sigma^2)) as before. 1 = interior fixed at E_rc's: each area's table is E_rc's table raked to the area's margins (one move per column, and per row when rows are parameters); mp_z has size 0; the same normal on b_j, restricted to those moves; needs lflag_seq_or_expected = 4 and lflag_mp_seq_anchor = 1. 2 = as 0 but the part of an area's deviation that is NOT such a move has sd kappa * sigma (kappa = 1 is mode 0 exactly, kappa -> 0 is mode 1); any allocation mode
  real<lower=0> mp_kappa_fixed;                // lflag_mp_interior = 2: 0 = kappa is a parameter with a lognormal(prior_mp_kappa_a, prior_mp_kappa_b) prior; > 0 = kappa fixed at this value
  real prior_mp_kappa_a;
  real<lower=0> prior_mp_kappa_b;
  int<lower=0, upper=1> lflag_E_cov;           // sequential margin param: 1 = the hierarchy mean varies by area, b_j ~ N(E_rc + E_cov_gamma .* E_cov_x[j], sigma). One covariate value per area and coordinate (e.g. the centred log share of the coordinate's own row), one coefficient per coordinate.
  int<lower=0, upper=4> lflag_col_eff;         // sequential margin param: area-level column effects shared by every row: in area j, col_eff[j, c] ~ normal(0, col_tau[c]) is added to the log share of column c in every row. 1 = explicit, non-centred; 2 = explicit, centred; 3 = integrated out: b_j ~ multi_normal(E_rc, diag(sigma^2) + mp_Wc * diag(col_tau^2) * mp_Wc'), no per-area effect parameters (needs lflag_mp_seq_scale = 2, lflag_mp_scale_fast = 0, mp_nc_w = 1)
  int<lower=0, upper=1> lflag_col_tau_shared;  // 1 = one col_tau for every column, 0 = one per column
  int<lower=0> n_col_groups;                   // explicit column effects (lflag_col_eff 1 or 2): 0 = one set of effects per area; > 0 = one set per group of areas (e.g. polling places within an electorate share the electorate's candidate effects)
  array[n_col_groups > 0 ? n_areas : 0] int<lower=1, upper=n_col_groups> col_group;
  int<lower=0, upper=1> lflag_col_cov;         // sequential margin param: 1 = a column shift proportional to an observed area-by-column covariate, shared by every row: hierarchy mean = E_rc + mp_Wc * (col_kappa .* col_x[j]). Can be combined with lflag_col_eff.
  int<lower=0, upper=1> lflag_col_kappa_shared; // 1 = one coefficient for every column, 0 = one per column
  matrix[lflag_col_cov ? n_areas : 0, C] col_x; // e.g. log(column share / share of the matching row), centred
  real<lower=0> prior_col_kappa_scale;         // prior sd of col_kappa (normal, mean 0)
  matrix[lflag_col_eff == 4 ? n_areas : 0, C] col_eff_known;   // lflag_col_eff = 4: the column effects are supplied as data (for diagnosis with known cells)
  real<lower=0> prior_col_tau_scale;           // prior sd of col_tau (half-normal)
  matrix[lflag_E_cov ? n_areas : 0, Dm1_model] E_cov_x;
  real<lower=0> prior_E_cov_scale;             // prior sd of the coefficients (normal, mean 0)
  int<lower=0> cz_n;                           // grid for the log integral: cz_n points, t = log(mean) from cz_t0 in steps of cz_h
  real cz_t0;
  real<lower=0> cz_h;
  vector[cz_n] cz_v;                           // log integral at the grid points
  vector[cz_n] cz_d;                           // its derivative with respect to t
  array[lflag_real_rowwise ? n_areas : 0] matrix[(R - 1) * (C - 1), (R - 1) * (C - 1)] rl_S;   // fixed scaling of lambda_raw per area (top-left block of size (non-empty rows - 1) * (C - 1) is used): S S' = inverse Poisson precision of the logits at a reference table, so lambda_raw is roughly unit scale. Any invertible matrix gives the same density.
  int<lower=0, upper=1> lflag_mp_beta_whiten;   // margin param (seq): 1 = the column log-ratio deviations are mp_T[j] * mp_beta[j], so mp_beta is roughly uncorrelated
  array[lflag_mp_beta_whiten ? n_areas : 0] matrix[C - 1, C - 1] mp_T;   // lower-triangular factor of the approximate covariance of the column log-ratios (from the counts)
  int<lower=0, upper=1> lflag_mp_row_whiten;    // margin param (seq), rows as parameters: 1 = the row log-ratio deviations are mp_Tr[j] * (row part of mp_beta[j])
  array[lflag_mp_row_whiten ? n_areas : 0] matrix[R - 1, R - 1] mp_Tr;  // lower-triangular factor of the approximate covariance of the row log-ratios (from the counts)
  int<lower=0, upper=4> lflag_seq_or_expected;  // sequential EXPECTED table: 1 = place each cell by its 2x2 log odds ratio (kink free); 0 = by its position between the bounds; 2 / 3 = ordered cell-wise / row-wise; 4 = adjusted table (raking style: the interior is the log odds ratios against the last row and column)
  real<lower=0> kink_delta_expected;            // kink smoothing in the sequential EXPECTED table (lflag_mp_seq): each bound kink is rounded over delta x (the cell's range); at most 1.4 x delta of any cell's range is lost. 0 = sharp (previous behaviour)
  real<lower=0> kink_delta_realised;            // kink smoothing in the REALISED table (hybrid) and its anchors, same meaning. 0 = sharp (previous behaviour: hinge_delta_* and slack_tol)
  int<lower=0, upper=2> lflag_mp_seq_scale;      // margin param (seq): 1 = interior logits = anchor + s(sigma) .* mp_z, with s^2 = square(mp_A) * sigma^2 (mp_z roughly unit scale); density unchanged (exact), so mp_A is only a preconditioner
  array[lflag_mp_seq_scale == 1 ? n_areas : 0] matrix[Dm1_model - n_pin, Dm1_model] mp_A;
  array[lflag_mp_seq_scale == 2 ? n_areas : 0] matrix[Dm1_model, Dm1_model - n_pin] mp_B;   // d b / d(interior logits) at a reference table (seq scale 2: generalised-least-squares anchor and scaling)   // d(interior logits)/d b at a reference table for each area (pseudo-inverse of d b / d logits)
  vector<lower=0, upper=1>[lflag_mp_seq_scale == 2 ? Dm1_model : 0] mp_nc_w;   // seq scale 2: per-coordinate non-centring weight. 1 = scale mp_z fully by sigma, 0 = scale by the fixed mp_sigma0 only
  vector<lower=0>[lflag_mp_seq_scale == 2 ? Dm1_model : 0] mp_sigma0;         // seq scale 2: fixed reference sigma used where the weight is below 1
  array[lflag_margin_param ? n_areas : 0] vector[Dm1_model] mp_b_ref;
  array[lflag_margin_param ? n_areas : 0] matrix[n_pin, Dm1_model] mp_G;
  array[lflag_margin_param ? n_areas : 0] matrix[Dm1_model, n_pin] mp_P;
  array[lflag_margin_param ? n_areas : 0] matrix[Dm1_model, Dm1_model - n_pin] mp_N;
  int<lower = 0, upper = 2> lflag_rot_lambda; // 0 - none, 1 - per area, 2 - per area plus cross area
  int<lower = 0> n_rot_lam;
  matrix[n_rot_lam, n_rot_lam] ROT_lambda;
  int<lower=0, upper = 3> lflag_fit_type; // flag indicating whether model is (0) HYBRID of sampling and ILR with cell resolution (Poisson) (1) ILR margin resolution (Poisson) (2) ILR margin resolution (multinomial) (3) ILR cell resolution (Poisson)
  int<lower=0, upper = 1> lflag_rot_E_rc;
  // 0 = Additive Log-Ratio 1 (C - 1) log-ratios representing the composition of the (R - 1) the free rows of the matrix
  // 3 = Log Odds Ratios of the
  int<lower=0, upper=1> lflag_fix_E_rc;        // 1 = use fixed values, 0 = estimate
  vector[(R * C) - 1] E_rc_fixed;         // fixed values, ignored if fix_E_rc=0
  int<lower = 0, upper = 1> lflag_lambda_centred; // lambda = lambda_vec (1) or lambda = lambda_vec + neutral_logit (with different prior structures for each case)
  real<lower=0> sigma_floor; //minimum value for all sigma
 real<lower=0> prior_mu_re_scale; // prior for scale of mu_re (mean row effect)
 real<lower=0> prior_mu_ce_scale; // prior for scale of col_effect (mean column effect)
 real<lower=0> prior_sigma_c_scale; //prior for scale of sigma_c (or sigma_c_sigma if lflag_vary_sd == 2)
 real<lower=0> prior_sigma_c_mu_scale; //prior for scale of sigma_c_mu (only if lflag_vary_sd == 2)
 real prior_sigma_mu; //prior for mean of sigma_jrc if lflag_vary_sd == 1 or 0
 real<lower=0> prior_sigma_ce_scale; //prior of scale for sigma_ce
 real<lower=0> prior_sigma_re_scale; //prior of scale for sigma_re
 real<lower=0> prior_cell_effect_scale; //prior of scale for average cell effects
 real<lower=0> prior_lambda_raw_scale; // prior of scale of the raw lambda_raw (greed) parameters (e.g. logit of the consumption of available mass in each cell)
 real<lower=0> prior_gamma_shape;
 real<lower=0> prior_gamma_rate;
 vector[(R * C) - 1] E_rc_prior; // empirically informed prior centres for E_rc
 int<lower=0> known_cell_values[n_areas, R, C]; // for testing purposes
 int<lower=0, upper=1> use_known_cells; // for testing purposes
 real<lower = 0> hinge_delta_floor;
 real<lower = 0> hinge_delta_min;
 real<lower=0.0> slack_tol;
  // E_rc grouping
  array[n_agg_free] int<lower=1> E_rc_group_id;  // which group each dim belongs to


  int<lower=1> E_rc_n_groups;
  vector[E_rc_n_groups] E_rc_group_fixed_value;
  array[E_rc_n_groups] int<lower=0,upper=4> E_rc_group_mode;  // 0 FIXED,1 SHARED,2 CP,3 NCP,4 PARTIAL_FIXED_SCALE
  vector[E_rc_n_groups] E_rc_group_prior_a;
  vector<lower=0>[E_rc_n_groups] E_rc_group_prior_b;
  vector<lower=0>[E_rc_n_groups] E_rc_group_tau_a;
  vector<lower=0>[E_rc_n_groups] E_rc_group_tau_b;

  int<lower=1> sigma_n_groups;
  array[Dm1_model] int<lower=1> sigma_group_id;
  array[sigma_n_groups] int<lower=0,upper=4> sigma_group_mode;
  vector[sigma_n_groups] sigma_group_prior_a;
  vector<lower=0>[sigma_n_groups] sigma_group_prior_b;
  vector<lower=0>[sigma_n_groups] sigma_group_tau_a;
  vector<lower=0>[sigma_n_groups] sigma_group_tau_b;
  vector[sigma_n_groups] sigma_group_fixed_value;
}
transformed data{
  int K;
  int K_j;
  int K_t;
  int K_no_rm;
  int K_c;
  int K_ame;
  int K_jr_start;
  int K_jc_start;
  int K_jrc_start;
  int K_jrc_rstart[R - 1];
  // int R_ll; // rows in the log-linear representation
  // int C_ll; // cols in the log-linear representation
  int has_area_cell_effects;
  int free_R[n_areas];
  int free_C[n_areas];
  int non0_rm;
  int non0_cm;
  int n_poss_cells=0;
  int n_structural_zeros=0;
  int structural_zero_rows[n_areas, R];
  int n_poss_rows=0;
  int structural_zero_cols[n_areas, C];
  int n_poss_cols=0;
  int n_zero_cells;
  int zero_cell_map[n_areas, R, C];
  int zero_cell_map_global[1, R, C];
  int structural_zeros_global[1, R, C];  // an array indicating any structural zeros in the global sums (shouldn't be any)
  int n_free_alr_cells = 0; // free cells in (R-1)*(C-1) submatrix
  int n_ilr_params= 0; // number of semi-free ilr components
  int has_theta;
  int has_area_re;
  int has_area_col_effects;
  int has_area_row_effects;
  int has_L;
  int has_onion;
  int has_L_ame;
  // int mu_re_ce_in_cell_effects  = 1;
  int param_map[n_areas, R - 1, C - 1];
  matrix[n_areas, R - 1] row_margins_lr;
  matrix[n_areas, C - 1] col_margins_lr;
  int row_margins_flat[n_areas * R];
  int col_margins_flat[n_areas * C];
  matrix[n_areas, R] rm_prop;
  matrix[n_areas, R] rm_log;
  matrix[n_areas, C] cm_prop;
  matrix[n_areas, C] cm_log;
  matrix[1, R] global_rm;
  matrix[1, C] global_cm;
  vector[R] mean_rm_prop;
  vector[R] global_rm_prop;
  vector[n_areas] tot_log;
  real<lower=0> prior_phi_scale;
  int n_margin_sigmas;
  int n_jrc_sigmas;
  int n_table_sigmas;
  int is_ilr = 1;
  real sigma_constrain = .001;
  int n_free_areas_rc[R - 1, C - 1]; // count of free areas for each (r, c) for decentred lambdas
  int dev_start_rc[R-1, C-1]; // start for deviation parameters
  int dev_idx = (R-1)*(C-1); // index for deviation parameters in (r, c) order
  // real hinge_delta_floor = 1e-10;
  // real hinge_delta_min = 1e-10;
  prior_phi_scale = 100;
  array[Dm1_model] int<lower=0, upper =1> lflag_agg_noncent = lflag_noncentred_mat[1, 1:(Dm1_model)];
  array[n_areas - 1, Dm1_model] int<lower=0, upper =1> lflag_dev_noncent = lflag_noncentred_mat[2:n_areas, 1:(Dm1_model)];

  if (lflag_rot_agg == 1) {
    if (lflag_rot_llrep != 1) reject("lflag_rot_agg requires lflag_rot_llrep == 1");
    for (s in 1:Dm1_model)
      if (lflag_agg_noncent[s] == 1) reject("lflag_rot_agg requires all agg dims centred");
    {
      matrix[Dm1_model, Dm1_model] chk = ROT_agg' * ROT_agg;
      for (a in 1:Dm1_model) for (b in 1:Dm1_model)
        if (abs(chk[a, b] - (a == b)) > 1e-8) reject("ROT_agg must be orthogonal");
     }
  }

  // margin parameterisation: size of the old stacked vector (0 when margin mode is on)
  int N_llpv = (lflag_margin_param == 1) ? 0 : n_areas * Dm1_model + n_areas - 1;
  if (lflag_margin_param == 1) {
    if (n_pin < 1 || n_pin >= Dm1_model) reject("margin_param needs 1 <= n_pin < Dm1_model; n_pin = ", n_pin);
    if (lflag_rot_agg == 1) reject("margin_param is incompatible with lflag_rot_agg == 1");
    if (lflag_rot_llrep == 1) print("note: lflag_margin_param = 1 overrides lflag_rot_llrep / ROT_red / lflag_noncentred_mat");
    if (lflag_mp_seq == 0) for (j in 1:n_areas) {
      matrix[n_pin, n_pin] gp = mp_G[j] * mp_P[j];
      matrix[n_pin, Dm1_model - n_pin] gn = mp_G[j] * mp_N[j];
      matrix[Dm1_model - n_pin, Dm1_model - n_pin] nn = mp_N[j]' * mp_N[j];
      for (a in 1:n_pin) for (b in 1:n_pin)
        if (abs(gp[a, b] - (a == b)) > 1e-6) reject("mp: G*P != I in area ", j);
      for (a in 1:n_pin) for (b in 1:(Dm1_model - n_pin))
        if (abs(gn[a, b]) > 1e-6) reject("mp: G*N != 0 in area ", j);
      for (a in 1:(Dm1_model - n_pin)) for (b in 1:(Dm1_model - n_pin))
        if (abs(nn[a, b] - (a == b)) > 1e-6) reject("mp: N'N != I in area ", j);
    }
  }

  int n_active_cells = 0;
  for (j in 1:n_areas)
    for (r in 1:R)
      for (c in 1:C)
        if (structural_zeros[j, r, c] == 0)
          n_active_cells += 1;   // only structural zeros excluded now — every other cell active

  array[n_active_cells] int active_j;
  array[n_active_cells] int active_r;
  array[n_active_cells] int active_c;
  {
    int idx = 0;
    for (j in 1:n_areas)
      for (r in 1:R)
        for (c in 1:C)
          if (structural_zeros[j, r, c] == 0) {
            idx += 1;
            active_j[idx] = j; active_r[idx] = r; active_c[idx] = c;
          }
  }

  if(lflag_ll_rep == 0 || lflag_ll_rep == 3){
    // ALR case (0)
    // LOR case (3)
    // R_ll = R - 1;
    // C_ll = C - 1;
    is_ilr = 0;
  } else if(lflag_ll_rep == 1||lflag_ll_rep == 2){
    // Preprogrammed ILR cases
    // R_ll = R;
    // C_ll = C - 1;
    is_ilr = 1;
  } else if(lflag_ll_rep == 4){
    // User defined ILR cases
    // R_ll = n_ilr_rows;
    // C_ll = C - 1;
    is_ilr = 1;
  }

  if(lflag_dist==2){
    has_theta = 1;
  } else{
    has_theta = 0;
  }

  if(lflag_area_re == 0){
    has_area_re = 0;
  } else {
    has_area_re = 1;
  }

  if(lflag_area_re == 2){
    has_L = 1;
  } else {
    has_L = 0;
  }

  if(lflag_area_re == 3){
    has_onion = 1;
  } else {
    has_onion = 0;
  }

  if(lflag_llmod_omit_jrc ==1){
    has_area_cell_effects = 0;
  } else {
    has_area_cell_effects = 1;
  }

  if(lflag_llmod_omit_jc == 1){
    has_area_col_effects = 0;
  } else {
    has_area_col_effects = 1;
  }
  if(lflag_llmod_omit_jr == 1){
    has_area_row_effects = 0;
  } else {
    has_area_row_effects = 1;
  }

  K_j = lflag_predictors_cm ==1 ? R - 1 + C - 1 + (R - 1)*(C - 1) : C - 1 + (R - 1)*(C - 1) ;
  K_jr_start = 0;
  K_jc_start = lflag_predictors_cm == 1? R - 1: 0;
  K_jrc_start = K_jc_start + C - 1;

  for(r in 1:R - 1){
    K_jrc_rstart[r] = K_jrc_start + (r - 1)*(C - 1);
  }



  K_t = (R - 1)* (C - 1);
  K_ame = has_area_row_effects * (R - 1) + has_area_col_effects * (C - 1); // number of area margin effects to estimate (area row effects + area col effects)

  K_c = C * (R - 1);
  K = R * (C - 1);
  K_no_rm = R * (C - 1);

  // if(lflag_inc_rm == 1){
  //   K = (R * C) - 1;
  // } else {
  //   K = R * (C - 1);
  // }
  // K_no_rm = R * (C - 1);


// following code to distinguish structural zeros from sampling zeros in the log-linear model - identifying rows and columns, and number of structural zeros in the data

  for(j in 1:n_areas){
    for(r in 1:R){
      structural_zero_rows[j, r] = sum(structural_zeros[j, r, 1:C])==C ? 1 :0;
      n_poss_rows += sum(structural_zeros[j, r, 1:C])==C ? 0 :1;
      n_structural_zeros += sum(structural_zeros[j, r, 1:C]);
    }
    for(c in 1:C){
      structural_zero_cols[j, c] = sum(structural_zeros[j, 1:R, c])==R ? 1 : 0;
      n_poss_cols += sum(structural_zeros[j, 1:R, c])==R ? 0 :1;
    }
  }
  for(r in 1:R){
    for(c in 1:C){
      structural_zeros_global[1, r, c] = 0;
    }
  }
  n_poss_cells = n_areas*R*C - n_structural_zeros;

  n_zero_cells = 0;
for(j in 1:n_areas){
    for(r in 1:R){
        for(c in 1:C){
            if(structural_zeros[j,r,c] == 0 &&
               (row_margins[j,r] == 0 || col_margins[j,c] == 0)){
                n_zero_cells += 1;
                zero_cell_map[j,r,c] = n_zero_cells;
            } else {
                zero_cell_map[j,r,c] = 0;
            }
        }
    }
}
for(r in 1:R){
  for(c in 1:C){
    zero_cell_map_global[1, r, c] = 0;
  }
}

for(j in 1:n_areas){
    for(r in 1:(R-1)){
        for(c in 1:(C-1)){
            if(structural_zeros[j,r,c] == 0){
                n_free_alr_cells += 1;
            }
        }
        for(c in 1:C){
          if(structural_zeros[j, r, c] == 0){
            n_ilr_params += 1;
          }
        }
    }
}


  // to deal with zeros in the sequential sampling: calculate the number of free parameters (zero row and columns do not need a parameter to allocated cell value of 0)

  int n_param = 0;
  int param_count_from[n_areas];

  for(j in 1:n_areas){
    free_R[j] = 0;
    free_C[j] = 0;
    for(r in 1:R){
      if(row_margins[j, r]>0){
        free_R[j] += 1;
      }
    }
    for(c in 1:C){
      if(col_margins[j, c] > 0){
        free_C[j] += 1;
      }
    }
    param_map[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j] - 1)){
      for (c in 1:(free_C[j] - 1)){
        param_map[j, r, c] = n_param + ((r - 1) * (free_C[j] -1)) + c;
       }
     }
    param_count_from[j] = n_param;
    n_param += max(0, (free_R[j] - 1) * (free_C[j] - 1));
  }
  non0_rm = sum(free_R);
  non0_cm = sum(free_C);




  // Parameter Trackers for cell, row, column and table sequential cell weights
  int n_param_alpha = 0;
  int n_param_beta = 0;
  int n_param_gamma = 0;

  int alpha_start[n_areas];
  int beta_start[n_areas];
  int gamma_start[n_areas];

  for(j in 1:n_areas){
    free_R[j] = 0;
    free_C[j] = 0;
    for(r in 1:R){
      if(row_margins[j, r] > 0){
        free_R[j] += 1;
      }
    }
    for(c in 1:C){
      if(col_margins[j, c] > 0){
        free_C[j] += 1;
      }
    }



    // Store the flat array 1-based start indices for this area
    alpha_start[j] = n_param_alpha + 1;
    beta_start[j] = n_param_beta + 1;
    gamma_start[j] = n_param_gamma + 1;

    int fr = free_R[j] - 1;
    int fc = free_C[j] - 1;

    // Rows scale needs exactly (fr - 1) unconstrained parameters
    if (fr > 1) {
      n_param_alpha += (fr - 1);
    }

    // Columns scale needs exactly (fc - 1) unconstrained parameters
    if (fc > 1) {
      n_param_beta += (fc - 1);
    }
    // 2D Matrix interaction needs exactly (fr - 1) * (fc - 1) parameters
    if (fr > 1 && fc > 1) {
      n_param_gamma += (fr - 1) * (fc - 1);
    }
  }

  array[n_areas, R] int active_row_map = rep_array(0, n_areas, R);
  array[n_areas, C] int active_col_map = rep_array(0, n_areas, C);

  for(j in 1:n_areas){
    int t = 0;
    for (r in 1:R){
      if (row_margins[j,r]>0) {
        t += 1;
         active_row_map[j, t] = r;
        }
    }
    int s = 0;
    for (c in 1:C){
      if (col_margins[j,c]>0){
        s += 1;
        active_col_map[j, s] = c;
        }
      }
    }


  // Param offset tracker
  array[n_areas, R - 1, C - 1] real neutral_logit_array = rep_array(0.0, n_areas, R - 1, C - 1);

  {
  int count = 0;
  for(j in 1:n_areas){
    for(r in 1:(free_R[j] - 1)){
      for(c in 1:(free_C[j] - 1)){
        count += 1;
        if(lflag_neutral_logit == 1){
          real col_sum_active = col_margins[j, active_col_map[j, c]];
          real remaining_sum = 0;
          for (k in c:free_C[j]) remaining_sum += col_margins[j, active_col_map[j, k]];
          real this_share = fmin(fmax(col_sum_active/fmax(remaining_sum, slack_tol), slack_tol), 1 - slack_tol);
          neutral_logit_array[j, r, c] = logit(this_share);
        } else{
          neutral_logit_array[j, r, c] = -log(free_C[j] - c);
        }
      }
    }
  }
  }



  // Parameter Trackers for decentred sequential cell weights with lambda_mu and then sum to zero lambda deviations from lambda_mu in each (r,c)

  // Count of free areas per (r, c)
  for(r in 1:(R-1)){
      for(c in 1:(C-1)){
          n_free_areas_rc[r,c] = 0;
          for(j in 1:n_areas){
              if(param_map[j,r,c] > 0){
                  n_free_areas_rc[r,c] += 1;
              }
          }
      }
  }
  // Assign start indices for lambda deviations
  // First (R-1)*(C-1) slots = lambda_mu parameters
  // Remaining slots = deviation parameters in (r,c) order
    for(r in 1:(R-1)){
        for(c in 1:(C-1)){
            dev_start_rc[r,c] = dev_idx + 1;
            dev_idx += max(0, n_free_areas_rc[r,c] - 1);
        }
    }
  for(j in 1:n_areas){
    for(c in 1:C){
      col_margins_flat[(j - 1)*C + c] = to_int(col_margins[j, c]);
      }
  for(r in 1:R){
    row_margins_flat[(j - 1)*R + r] = to_int(row_margins[j, r]);
    }
  }


  for (j in 1:n_areas){
    tot_log[j] = log(sum(row_margins[j, 1:R]));
    for(r in 1:(R - 1)){
      row_margins_lr[j, r] = log((row_margins[j, r]+ .1)/(row_margins[j, R] + .1));
    }
    for(r in 1:R){
      rm_prop[j, r] = (row_margins[j, r] + .01)/(sum(row_margins[j, 1:R]) + R*.01);
      if(row_margins[j, r] == 0){
        rm_log[j, r] = log(.001);
      } else{
        rm_log[j, r] = log(row_margins[j, r]);
      }
    }
    for(c in 1:(C - 1)){
      col_margins_lr[j, c] = log((col_margins[j, c]+ .1)/(col_margins[j, C] + .1));
    }
    for(c in 1:C){
      cm_prop[j, c] = (col_margins[j, c] + .01)/(sum(col_margins[j, 1:C]) + C*.01);
      if(col_margins[j, c] == 0){
        cm_log[j, c] = log(.001);
      } else{
        cm_log[j, c] = log(col_margins[j, c]);
      }

    }

  }
  for (r in 1:R) mean_rm_prop[r] = mean(rm_prop[1:n_areas, r]);

  for(r in 1:R) global_rm[1, r] = sum(row_margins[1:n_areas, r]);
  for(c in 1:C) global_cm[1, c] = sum(col_margins[1:n_areas, c]);

    int n_full_lam = 0;
  for (j in 1:n_areas)
    if (free_R[j] == R && free_C[j] == C) n_full_lam += (R - 1) * (C - 1);
  array[n_full_lam] int full_lam_idx;
  array[n_full_lam] int full_lam_r;
  array[n_full_lam] int full_lam_c;
  array[n_param - n_full_lam] int part_lam_idx;
  {
    int a = 0; int b = 0;
    for (j in 1:n_areas)
      for (r in 1:(free_R[j] - 1))
        for (c in 1:(free_C[j] - 1)) {
          int idx = param_count_from[j] + (r - 1) * (free_C[j] - 1) + c;
          if (free_R[j] == R && free_C[j] == C) {
            a += 1; full_lam_idx[a] = idx; full_lam_r[a] = r; full_lam_c[a] = c;
          } else {
            b += 1; part_lam_idx[b] = idx;
          }
        }
  }



if(n_param != (n_param_gamma + n_param_alpha + n_param_beta + n_areas)){
  print("paramater map mismatch");
  print(n_param);
  print(n_param_gamma + n_param_alpha + n_param_beta + n_areas);
  print(n_param_gamma);
  print(n_param_alpha);
  print(n_param_beta);
  print(n_areas);

}



  has_L_ame = 0;

  n_margin_sigmas = K_ame;
  int n_row_sigmas = has_area_row_effects * (R - 1);
  int n_col_sigmas = has_area_col_effects * (C - 1);
  n_jrc_sigmas = has_area_cell_effects * ((lflag_vary_sd ==0) ? 1 : (R - 1) * (C - 1));
  n_table_sigmas = n_margin_sigmas + n_jrc_sigmas;

  int K_all = n_table_sigmas;


  int D = R * C;
  int Dm1 = D - 1;
  int Dtot = n_areas * D;
  int Dtot_m1 = Dtot - 1;
  int D_reduced = Dm1_model + 1;
  int Dtot_m1_reduced = n_areas*D_reduced - 1;

  matrix[n_areas, n_areas - 1] V_area = make_helmert_basis(n_areas);

  // --- E_rc: compact hyperparameter index (only SHARED/CP/NCP groups need mu; only CP/NCP need sigma) ---
  int E_rc_n_mu = 0;
  int E_rc_n_sigma = 0;
  array[E_rc_n_groups] int E_rc_mu_pos = rep_array(0, E_rc_n_groups);     // 0 = "no mu param"
  array[E_rc_n_groups] int E_rc_sigma_pos = rep_array(0, E_rc_n_groups); // 0 = "no sigma param"

   for (g in 1:E_rc_n_groups) {
    if (E_rc_group_mode[g] >= 1) { E_rc_n_mu += 1; E_rc_mu_pos[g] = E_rc_n_mu; }
    if (E_rc_group_mode[g] == 2 || E_rc_group_mode[g] == 3) { E_rc_n_sigma += 1; E_rc_sigma_pos[g] = E_rc_n_sigma; }
  }
  // --- E_rc: compact leaf index (only CP/NCP dims need their own raw leaf) ---
  int E_rc_n_leaf = 0;
  for (i in 1:n_agg_free) if (E_rc_group_mode[E_rc_group_id[i]] >= 2) E_rc_n_leaf += 1;
  array[E_rc_n_leaf] int E_rc_leaf_idx;
  {
    int k = 0;
    for (i in 1:n_agg_free)
      if (E_rc_group_mode[E_rc_group_id[i]] >= 2) { k += 1; E_rc_leaf_idx[k] = i; }
  }

  // --- sigma_jrc: identical pattern ---
  int sigma_n_mu = 0;
  int sigma_n_sigma = 0;
  array[sigma_n_groups] int sigma_mu_pos = rep_array(0, sigma_n_groups);
  array[sigma_n_groups] int sigma_sigma_pos = rep_array(0, sigma_n_groups);
  for (g in 1:sigma_n_groups) {
    if (sigma_group_mode[g] >= 1) { sigma_n_mu += 1; sigma_mu_pos[g] = sigma_n_mu; }
    if (sigma_group_mode[g] >= 2) { sigma_n_sigma += 1; sigma_sigma_pos[g] = sigma_n_sigma; }
  }
  int sigma_n_leaf = 0;
  for (i in 1:Dm1_model) if (sigma_group_mode[sigma_group_id[i]] >= 2) sigma_n_leaf += 1;
  array[sigma_n_leaf] int sigma_leaf_idx;
  {
    int k = 0;
    for (i in 1:Dm1_model)
      if (sigma_group_mode[sigma_group_id[i]] >= 2) { k += 1; sigma_leaf_idx[k] = i; }
  }

  // margin param (seq): observed margin log-ratios (rm_prop / cm_prop are the smoothed observed shares)
  //   rows conditioned (lflag_row_decompose == 1): n_pin = C - 1,          Dm1_model = R * (C - 1)
  //   rows free        (lflag_row_decompose == 0): n_pin = R - 1 + C - 1,  Dm1_model = R * C - 1
  //   in both cases mp_z has (R - 1) * (C - 1) interior logits; with log_volume that is R x C per area when rows are free.
  array[(lflag_margin_param == 1 && lflag_mp_seq == 1) ? n_areas : 0] vector[C - 1] mp_s_obs;
  vector[lflag_mp_vol_scale ? n_areas : 0] mp_vol_loc;      // log of the area totals
  vector[lflag_mp_vol_scale ? n_areas : 0] mp_vol_sc;       // 1 / sqrt(area total): posterior sd of a log Poisson volume
  array[lflag_mp_scale_fast == 1 ? n_areas : 0, lflag_mp_scale_fast == 1 ? sigma_n_groups : 0] matrix[Dm1_model - n_pin, Dm1_model - n_pin] mp_C;   // B_j' D_g B_j, D_g = selector of the coordinates in sigma group g
  array[lflag_mp_scale_fast >= 1 ? sigma_n_groups : 0] int mp_sig_first;   // one coordinate from each sigma group
  array[(lflag_margin_param == 1 && lflag_mp_seq == 1 && lflag_row_decompose == 0) ? n_areas : 0] vector[R - 1] mp_r_obs;
  if (lflag_mp_vol_scale == 1) {
    for (j in 1:n_areas) {
      real Nj = sum(row_margins[j, 1:R]);
      mp_vol_loc[j] = log(Nj);
      mp_vol_sc[j] = inv_sqrt(Nj);
    }
  }
  if (lflag_mp_scale_fast >= 1) {
    if (lflag_mp_seq_scale != 2) reject("lflag_mp_scale_fast needs lflag_mp_seq_scale = 2");
    for (g in 1:sigma_n_groups) {
      if (sigma_group_mode[g] > 1) reject("lflag_mp_scale_fast needs shared or fixed sigma groups");
      mp_sig_first[g] = 0;
    }
    for (i in 1:Dm1_model) {
      if (mp_nc_w[i] != 1) reject("lflag_mp_scale_fast needs mp_nc_w = 1 everywhere");
      if (mp_sig_first[sigma_group_id[i]] == 0) mp_sig_first[sigma_group_id[i]] = i;
    }
    if (lflag_mp_scale_fast == 1) for (j in 1:n_areas) {
      for (g in 1:sigma_n_groups) {
        matrix[Dm1_model, Dm1_model - n_pin] Bg = mp_B[j];
        for (i in 1:Dm1_model) if (sigma_group_id[i] != g) Bg[i] = rep_row_vector(0, Dm1_model - n_pin);
        mp_C[j, g] = crossprod(Bg);
      }
    }
  }
  // low-rank form of the scaling (lflag_mp_scale_fast == 2).
  //   precision of the logits: Mi = B' diag(1/sigma^2) B = w_b L0 (I + U W U') L0',  w_b = 1/sigma^2 of the largest sigma group,
  //   L0 = chol(B'B), U = L0^-1 B_u' (B_u = rows of B outside the largest group), W = diag(sigma_b^2 / sigma_u^2 - 1).
  //   With U = Q Rq (thin QR) and S S' = I + Rq W Rq',  F = sqrt(w_b) L0 (I + Q (S - I) Q') satisfies F F' = Mi.
  int mp_Kq = Dm1_model - n_pin;
  int mp_gb = 1;                      // the largest sigma group
  int mp_nu = 0;                      // number of coordinates outside it = size of the Cholesky needed per area
  if (lflag_mp_scale_fast == 2) {
    array[sigma_n_groups] int cnt = rep_array(0, sigma_n_groups);
    for (i in 1:Dm1_model) cnt[sigma_group_id[i]] += 1;
    for (g in 1:sigma_n_groups) if (cnt[g] > cnt[mp_gb]) mp_gb = g;
    mp_nu = Dm1_model - cnt[mp_gb];
  }
  array[lflag_mp_scale_fast == 2 ? n_areas : 0] matrix[mp_Kq, mp_Kq] mp_L0i;            // L0^-1
  array[lflag_mp_scale_fast == 2 ? n_areas : 0] matrix[mp_Kq, Dm1_model] mp_LB;         // L0^-1 B'
  array[lflag_mp_scale_fast == 2 ? n_areas : 0] matrix[mp_Kq, mp_nu] mp_Q;
  array[lflag_mp_scale_fast == 2 ? n_areas : 0] matrix[mp_nu * mp_nu, sigma_n_groups] mp_RG;   // column g = vec(Rq_g Rq_g'), Rq_g = the columns of Rq in group g
  vector[lflag_mp_scale_fast == 2 ? n_areas : 0] mp_ld0;                                // log det L0
  if (lflag_mp_scale_fast == 2) {
    for (j in 1:n_areas) {
      matrix[mp_Kq, mp_Kq] L0 = cholesky_decompose(crossprod(mp_B[j]));
      mp_L0i[j] = mdivide_left_tri_low(L0, diag_matrix(rep_vector(1, mp_Kq)));
      mp_LB[j] = mp_L0i[j] * mp_B[j]';
      mp_ld0[j] = sum(log(diagonal(L0)));
      if (mp_nu > 0) {
        matrix[mp_Kq, mp_nu] U;
        matrix[mp_nu, mp_nu] Rq;
        int k = 0;
        for (i in 1:Dm1_model) if (sigma_group_id[i] != mp_gb) { k += 1; U[, k] = mp_LB[j][, i]; }
        mp_Q[j] = qr_thin_Q(U);
        Rq = qr_thin_R(U);
        for (g in 1:sigma_n_groups) {
          matrix[mp_nu, mp_nu] Rg = rep_matrix(0, mp_nu, mp_nu);
          k = 0;
          for (i in 1:Dm1_model) if (sigma_group_id[i] != mp_gb) { k += 1; if (sigma_group_id[i] == g) Rg[, k] = Rq[, k]; }
          mp_RG[j][, g] = to_vector(tcrossprod(Rg));
        }
      }
    }
  }

  // HYBRID, realised table row by row: the non-empty rows of each area in allocation order
  array[lflag_real_rowwise ? n_areas : 0] int rl_n;          // rows allocated = non-empty rows - 1
  array[lflag_real_rowwise ? n_areas : 0, R] int rl_row;     // table row of the k-th non-empty row in allocation order (entry rl_n + 1 is the remainder row)
  array[lflag_real_rowwise ? n_areas : 0, R] int rl_step;    // its allocation step in the expected table (position in mp_row_order)
  if (lflag_real_rowwise == 1) {
    if (lflag_fit_type != 0 || lflag_margin_param != 1 || lflag_mp_seq != 1 || lflag_seq_or_expected != 3 || lflag_neutral_logit != 2 || lflag_rot_lambda != 0 || lflag_lambda_centred != 0)
      reject("lflag_real_rowwise needs fit type 0 (HYBRID), lflag_margin_param = 1, lflag_mp_seq = 1, lflag_seq_or_expected = 3, lflag_neutral_logit = 2, no lambda rotation, non-centred lambda");
    for (j in 1:n_areas) {
      int t = 0;
      rl_row[j] = rep_array(0, R);
      rl_step[j] = rep_array(0, R);
      if (free_C[j] != C) reject("lflag_real_rowwise: empty columns not supported (area ", j, ")");
      for (k in 1:R) {
        int r = mp_row_order[k];
        if (row_margins[j, r] > 0) { t += 1; rl_row[j, t] = r; rl_step[j, t] = k; }
      }
      rl_n[j] = t - 1;
    }
  }

  // column effects: direction in the model coordinates of adding 1 to the log share of column c in every row
  if (lflag_col_eff == 3) {
    if (lflag_mp_seq != 1 || lflag_mp_seq_scale != 2 || lflag_mp_scale_fast != 0) reject("lflag_col_eff = 3 needs lflag_mp_seq = 1, lflag_mp_seq_scale = 2 and lflag_mp_scale_fast = 0");
    for (i in 1:Dm1_model) if (mp_nc_w[i] != 1) reject("lflag_col_eff = 3 needs mp_nc_w = 1 everywhere");
  }
  matrix[(lflag_col_eff > 0 || lflag_col_cov == 1) ? Dm1_model : 0, (lflag_col_eff > 0 || lflag_col_cov == 1) ? C : 0] mp_Wc;
  if (lflag_col_eff > 0 || lflag_col_cov == 1) {
    for (c in 1:C) {
      vector[R * C] ind = rep_vector(0, R * C);
      for (r in 1:R) ind[(r - 1) * C + c] = 1;
      mp_Wc[, c] = V_ilr_model' * ind;
    }
  }
  // interior fixed or scaled (lflag_mp_interior > 0): the moves an area can make without changing its interaction.
  //   mp_Hc, mp_Hr: orthonormal contrasts (Helmert) for the column and row moves on the log scale
  //   mp_Wm: Dm1_model x n_pin, the direction in the model coordinates of each move (row moves first when rows are parameters)
  int mp_n_mv = lflag_mp_interior > 0 ? n_pin : 0;
  int mp_kappa_free = (lflag_mp_interior == 2 && mp_kappa_fixed == 0) ? 1 : 0;
  matrix[C, C - 1] mp_Hc = rep_matrix(0, C, C - 1);
  matrix[R, R - 1] mp_Hr = rep_matrix(0, R, R - 1);
  matrix[lflag_mp_interior > 0 ? Dm1_model : 0, mp_n_mv] mp_Wm;
  for (k in 1:(C - 1)) { for (i in 1:k) mp_Hc[i, k] = inv_sqrt(k * (k + 1.0)); mp_Hc[k + 1, k] = -k * inv_sqrt(k * (k + 1.0)); }
  for (k in 1:(R - 1)) { for (i in 1:k) mp_Hr[i, k] = inv_sqrt(k * (k + 1.0)); mp_Hr[k + 1, k] = -k * inv_sqrt(k * (k + 1.0)); }
  if (lflag_mp_interior > 0) {
    matrix[R * C, mp_n_mv] Kf = rep_matrix(0, R * C, mp_n_mv);
    int offm = lflag_row_decompose == 1 ? 0 : R - 1;
    if (lflag_margin_param != 1 || lflag_mp_seq != 1) reject("lflag_mp_interior needs lflag_margin_param = 1 and lflag_mp_seq = 1");
    if (lflag_col_eff != 0) reject("lflag_mp_interior cannot be combined with column effects (lflag_col_eff)");
    if (lflag_mp_interior == 1 && (lflag_seq_or_expected != 4 || lflag_mp_seq_anchor != 1)) reject("lflag_mp_interior = 1 needs lflag_seq_or_expected = 4 (adjusted table) and lflag_mp_seq_anchor = 1");
    if (lflag_mp_interior == 2 && lflag_mp_scale_fast != 0) reject("lflag_mp_interior = 2 needs lflag_mp_scale_fast = 0");
    if (lflag_mp_interior == 2 && lflag_mp_seq_scale == 2) for (i in 1:Dm1_model) if (mp_nc_w[i] != 1) reject("lflag_mp_interior = 2 needs mp_nc_w = 1 everywhere");
    if (mp_n_mv != (lflag_row_decompose == 1 ? C - 1 : R - 1 + C - 1)) reject("lflag_mp_interior: n_pin does not match the number of moves");
    for (r in 1:R) for (c in 1:C) {
      if (lflag_row_decompose == 0) for (k in 1:(R - 1)) Kf[(r - 1) * C + c, k] = mp_Hr[r, k];
      for (k in 1:(C - 1)) Kf[(r - 1) * C + c, offm + k] = mp_Hc[c, k];
    }
    mp_Wm = V_ilr_model' * Kf;
  }
  if (lflag_margin_param == 1 && lflag_mp_seq == 1) {
    if (lflag_row_decompose == 1) {
      if (n_pin != C - 1 || Dm1_model != R * (C - 1)) reject("lflag_mp_seq (rows conditioned) needs n_pin == C - 1 and Dm1_model == R * (C - 1)");
    } else {
      if (n_pin != R - 1 + C - 1 || Dm1_model != R * C - 1) reject("lflag_mp_seq (rows free) needs n_pin == R - 1 + C - 1 and Dm1_model == R * C - 1");
    }
    for (j in 1:n_areas) {
      // empty ROWS are fine with the ordered allocations (modes 2, 3): rm_prop is smoothed, so an empty row is a
      // very small row whose cells are tiny fractions of the column capacities and whose log-ratios are prior-led.
      if (free_C[j] != C) reject("lflag_mp_seq: empty columns not supported yet (area ", j, ")");
      if (free_R[j] != R && lflag_seq_or_expected < 2) reject("lflag_mp_seq: empty rows need the ordered allocations or the adjusted table (lflag_seq_or_expected 2, 3 or 4), area ", j);
      for (c in 1:(C - 1)) mp_s_obs[j][c] = log(cm_prop[j, c]) - log(cm_prop[j, C]);
      if (lflag_row_decompose == 0)
        for (r in 1:(R - 1)) mp_r_obs[j][r] = log(rm_prop[j, r]) - log(rm_prop[j, R]);
    }
  }

  // margin param (exact): column log-ratios at each area's reference point, with the model's own row weights
  array[(lflag_margin_param == 1 && lflag_mp_exact == 1 && lflag_mp_seq == 0) ? n_areas : 0] vector[n_pin] mp_s_ref;
  if (lflag_margin_param == 1 && lflag_mp_exact == 1 && lflag_mp_seq == 0) {
    if (lflag_row_decompose != 1) reject("lflag_mp_exact requires lflag_row_decompose == 1");
    if (n_pin != C - 1) reject("lflag_mp_exact requires n_pin == C - 1");
    for (j in 1:n_areas)
      mp_s_ref[j] = mp_colshare(mp_b_ref[j], rm_prop[j]', V_ilr_model, R, C)[, 1];
  }

  matrix[R*C, R*C - 1] S_plus_full;
  matrix[R*C, R*C - 1] S_minus_full;
  for (k in 1:(R*C - 1)) {
    for (i in 1:(R*C)) {
      S_plus_full[i, k]  = V_ilr_full[i, k] >  1e-10 ? 1 : 0;
      S_minus_full[i, k] = V_ilr_full[i, k] < -1e-10 ? 1 : 0;
    }
  }
}
parameters{
  vector[lflag_fit_type==0 ? n_param: 0] lambda_raw;


  // real LLrep_raw[n_areas, (R * C) - 1];
  // vector[n_areas] log_volume;
  vector[N_llpv] LLrep_plus_log_volume_raw;     // size 0 in margin mode
  real log_volume_raw;
  // margin parameterisation
  array[lflag_margin_param ? n_areas : 0] vector[n_pin] mp_beta;              // pinned margin coords (centred)
  array[lflag_margin_param ? n_areas : 0] vector[lflag_mp_interior == 1 ? 0 : Dm1_model - n_pin] mp_z;     // interaction coords: z (NCP | beta) or gamma itself if lflag_mp_gamma_centred; none when the interior is fixed
  vector<lower=0>[mp_kappa_free] mp_kappa;     // lflag_mp_interior = 2: sd of the non-move part of an area's deviation, as a multiple of sigma
  vector[lflag_margin_param ? n_areas - 1 : 0] log_volume_rest;              // log_volume[2:n_areas]
  vector[E_rc_n_leaf] E_rc_raw_leaf;
  vector[E_rc_n_mu] E_rc_group_mu;
  vector<lower=0>[E_rc_n_sigma] E_rc_group_sigma;

  vector[sigma_n_leaf] sigma_raw_leaf;
  vector[sigma_n_mu] sigma_group_mu;
  vector[lflag_E_cov ? Dm1_model : 0] E_cov_gamma;   // coefficients of the area covariate in the hierarchy mean
  vector[lflag_col_cov ? (lflag_col_kappa_shared ? 1 : C) : 0] col_kappa;   // coefficient of the column covariate
  array[(lflag_col_eff == 1 || lflag_col_eff == 2) ? (n_col_groups > 0 ? n_col_groups : n_areas) : 0] vector[C] col_eff_raw;      // area-level column effects when explicit (raw scale if non-centred)
  vector<lower=0>[(lflag_col_eff >= 1 && lflag_col_eff <= 3) ? (lflag_col_tau_shared ? 1 : C) : 0] col_tau;   // their sd across areas (one, or one per column)
  vector<lower=0>[sigma_n_sigma] sigma_group_sigma;


}
transformed parameters{
  vector[lflag_fit_type==0 ? n_param : 0] lambda_vec;
  vector[lflag_fit_type==0 ? n_param : 0] neutral_logit_flat;
  real lambda[lflag_fit_type==0||lflag_fit_type==3 ? n_areas : 0, R - 1, C -1]; // sequential cell weights
  vector[N_llpv] LLrep_plus_log_volume;
  matrix[n_areas, Dm1] LLrep_jrc;
  real log_grand_volume;
  vector[n_areas] log_volume;
  vector[N_llpv] LLrep_plus_log_volume_xform = LLrep_plus_log_volume_raw;
  real mp_kink_loss_max = 0;   // seq expected table: largest share of any cell's true range lost to kink smoothing (this draw)
  real mp_kink_cells = 0;      // seq expected table: number of interior cells losing more than 1% of their range (this draw)
  real rl_lj = 0;        // HYBRID row by row: log Jacobian of the realised tables
  real mp_col_lp = 0;    // column effects integrated out: log density of the tables under the full covariance minus that under diag(sigma^2)
  array[lflag_col_eff == 3 ? n_areas : 0] vector[C] mp_col_eff_mean;   // column effects integrated out: conditional mean of each area's effects given its table
  real rl_norm = 0;      // HYBRID row by row: sum of the log integrals of the continuous Poisson densities (lflag_cpois_norm)
  real mp_beta_lp = 0;   // margin param: log prior of mp_beta (and of gamma if centred) given (E_rc, sigma_jrc)

  real<lower=0> cell_values[n_areas, R, C];
  real log_expected_cell_values[n_areas, R, C];
  // real log_cv[n_areas, R, C] ;
  vector<lower=0>[Dm1_model] sigma_jrc;
  real composition_arr[n_areas, R, C] = rep_array(0.0, n_areas, R, C);
  // vector[(R * C) - 1] E_rc;
  vector[Dm1_model] E_rc;
  matrix[R - 1, C - 1] anchor_lambda_E_rc = rep_matrix(0.0, R - 1, C - 1);
  // real det_J_raw = 0;
  real sigma_scale[lflag_rawscw==0?n_areas:0, R - 1, C - 1];
  real sigma_sq_row_out[lflag_rawscw==0?n_areas:0, (R * C) - 1];
  real sens_sq_out[lflag_rawscw==0?n_areas:0, (R * C) - 1];
  real tmp_J_mu_out[lflag_rawscw==0?n_areas:0, (R * C) - 1];
  array[(lflag_neutral_logit==2||lflag_neutral_logit==3)?n_areas:0, R - 1, C - 1] real anchor_lambda = rep_array(0.0, (lflag_neutral_logit==2||lflag_neutral_logit==3)?n_areas:0, R - 1, C - 1);

vector[n_agg_free] E_rc_raw;
  {
    vector[n_agg_free] u;
    int k = 0;
    for (i in 1:n_agg_free) {
      int g = E_rc_group_id[i];
      int m = E_rc_group_mode[g];
      if (m == 0) {
        u[i] = E_rc_group_fixed_value[g];
      } else if (m == 1) {
        u[i] = E_rc_group_mu[E_rc_mu_pos[g]];
      } else if (m == 2) {
        k += 1;
        u[i] = E_rc_raw_leaf[k];
      } else if (m == 3) {
        k += 1;
        u[i] = E_rc_group_mu[E_rc_mu_pos[g]] + E_rc_group_sigma[E_rc_sigma_pos[g]] * E_rc_raw_leaf[k];
      } else {   // m == 4, partial_fixed
        k += 1;
        u[i] = E_rc_group_mu[E_rc_mu_pos[g]] + E_rc_group_fixed_value[g] * E_rc_raw_leaf[k];
      }
    }

  if (lflag_E_rc_node_logit == 1) {
      vector[R*C - 1] u_full = rep_vector(0, R*C - 1);
      u_full[agg_free_dim] = u;           // derived positions stay at 0 -- annihilated by projection
      vector[R*C] log_p = S_plus_full * log_inv_logit(u_full) + S_minus_full * log1m_inv_logit(u_full);
      E_rc_raw = V_ilr_model' * log_p;
    } else {
      E_rc_raw = u;
    }


  }
  vector[Dm1_model] sigma_jrc_raw;   // log-scale, feeds existing family branch below
  {
    int k = 0;
    for (i in 1:Dm1_model) {
      int g = sigma_group_id[i];
      int m = sigma_group_mode[g];
      if (m == 0) {
        sigma_jrc_raw[i] = sigma_group_fixed_value[g];
      } else if (m == 1) {
        sigma_jrc_raw[i] = sigma_group_mu[sigma_mu_pos[g]];
      } else if (m == 2) {
        k += 1;
        sigma_jrc_raw[i] = sigma_raw_leaf[k];
      } else {
        k += 1;
        sigma_jrc_raw[i] = sigma_group_mu[sigma_mu_pos[g]] + sigma_group_sigma[sigma_sigma_pos[g]] * sigma_raw_leaf[k];
      }
    }
  }

    sigma_jrc = exp(sigma_jrc_raw);

  if (lflag_fix_E_rc == 1) {
    E_rc = E_rc_fixed;
  } else if (lflag_rot_E_rc == 1) {
    E_rc = ROT_E_rc * E_rc_raw;
  } else {
    E_rc = E_rc_raw;
  }

 if(lflag_neutral_logit == 3){
    row_vector[R*C] clr_model_E_rc = to_row_vector(E_rc) * V_ilr_model';

    vector[R*C] log_p_E_rc;
    if (Dm1_model < (R*C - 1)) {
      // reduced dimension: row split isn't in E_rc, use the cross-area average instead
      int idx = 0;
      for (r in 1:R) {
        vector[C] log_q_r;
        for (c in 1:C) { idx += 1; log_q_r[c] = clr_model_E_rc[idx]; }
        real log_norm = log_sum_exp(log_q_r);
        for (c in 1:C) log_p_E_rc[(r-1)*C + c] = log(mean_rm_prop[r]) + log_q_r[c] - log_norm;
      }
    } else {
      // full dimension: E_rc already carries row-split information, plain closure suffices
      log_p_E_rc = log_softmax(to_vector(clr_model_E_rc));
    }

    matrix[R, C] comp_E_rc;
    {
      int idx = 0;
      for (r in 1:R) for (c in 1:C) { idx += 1; comp_E_rc[r, c] = exp(log_p_E_rc[idx]); }
    }
    anchor_lambda_E_rc = inverse_alloc(comp_E_rc, hinge_delta_floor, hinge_delta_min, kink_delta_realised, 0);
}


  {
  vector[n_areas - 1] vol_coeffs;
  vector[n_areas] vol_dev = rep_vector(0, n_areas);
  array[lflag_margin_param ? n_areas : 0] row_vector[Dm1_model] b_mp;
  array[(lflag_margin_param == 1 && lflag_mp_seq == 1) ? n_areas : 0, R - 1, C - 1] real mp_lambda_E;   // expected-table logits, used as the realised table's anchor

  if (lflag_margin_param == 1) {
    // ---- margin parameterisation ----
    // Prior is the existing hierarchy b_j ~ N(E_rc, diag(sigma_jrc^2)), written in the
    // per-area coordinates (beta_j, gamma_j) = (G_j (b_j - b_ref_j), N_j' (b_j - b_ref_j)).
    // beta_j is centred (data-pinned); gamma_j | beta_j is non-centred via z_j. Map is linear
    // with constant Jacobian, so no adjustment is needed.
    if (lflag_mp_vol_scale == 1) log_volume = mp_vol_loc + mp_vol_sc .* append_row(log_volume_raw, log_volume_rest);   // linear: constant Jacobian
    else log_volume = append_row(log_volume_raw, log_volume_rest);
    log_grand_volume = log_sum_exp(log_volume);
    {
      vector[Dm1_model] s2 = square(sigma_jrc);
      matrix[R, C] mp_qE;                                                      // E_rc's table: row-conditional shares (rows conditioned) or joint shares
      vector[lflag_mp_scale_fast == 2 ? sigma_n_groups : 0] mp_wr;             // sigma_b^2 / sigma_g^2 - 1
      if (lflag_mp_seq == 1 && lflag_mp_seq_anchor == 1) {
        vector[R * C] clrE = V_ilr_model * E_rc;
        if (lflag_row_decompose == 1) {
          for (r in 1:R) mp_qE[r] = softmax(clrE[((r - 1) * C + 1):(r * C)])';
        } else {
          vector[R * C] pE = softmax(clrE);
          for (r in 1:R) for (c in 1:C) mp_qE[r, c] = pE[(r - 1) * C + c];
        }
      }
      if (lflag_mp_scale_fast == 2)
        for (g in 1:sigma_n_groups) mp_wr[g] = square(sigma_jrc[mp_sig_first[mp_gb]] / sigma_jrc[mp_sig_first[g]]) - 1;
      // column effects: Sigma = D + Wt Wt', D = diag(sigma^2), Wt = mp_Wc * diag(col_tau).
      //   Sigma^-1 = D^-1 - DiW DiW',  DiW = D^-1 Wt L_M'^-1,  L_M L_M' = I + Wt' D^-1 Wt;  log det Sigma = log det D + 2 * sum(log(diag(L_M)))
      vector[(lflag_col_eff >= 1 && lflag_col_eff <= 3) ? C : 0] mp_tauv;
      matrix[lflag_col_eff == 3 ? Dm1_model : 0, lflag_col_eff == 3 ? C : 0] mp_Wt;
      matrix[lflag_col_eff == 3 ? Dm1_model : 0, lflag_col_eff == 3 ? C : 0] mp_DiW;
      real mp_ldM = 0;
      // moves (lflag_mp_interior > 0): P = Wm' D^-1 Wm is the precision of the move coordinates under the normal on b;
      //   mp_LP = its Cholesky factor, mp_DiW2 = D^-1 Wm LP'^-1, so that D^-1 Wm P^-1 Wm' D^-1 = DiW2 DiW2'.
      //   lflag_mp_interior = 2: Sigma = kappa^2 D + (1 - kappa^2) Wm P^-1 Wm',  Sigma^-1 = kappa^-2 D^-1 + (1 - kappa^-2) DiW2 DiW2',
      //   log det Sigma = log det D + 2 (Dm1_model - n_pin) log kappa.
      matrix[mp_n_mv, mp_n_mv] mp_LP;
      matrix[lflag_mp_interior == 2 ? Dm1_model : 0, lflag_mp_interior == 2 ? mp_n_mv : 0] mp_DiW2;
      real mp_kap = 1;
      real mp_ik2 = 1;                                                         // kappa^-2
      if (lflag_mp_interior > 0) {
        matrix[Dm1_model, mp_n_mv] DW = diag_pre_multiply(inv(s2), mp_Wm);
        matrix[mp_n_mv, mp_n_mv] PP = mp_Wm' * DW;
        mp_LP = cholesky_decompose(0.5 * (PP + PP'));
        if (lflag_mp_interior == 2) {
          mp_DiW2 = mdivide_left_tri_low(mp_LP, DW')';
          mp_kap = mp_kappa_free == 1 ? mp_kappa[1] : mp_kappa_fixed;
          mp_ik2 = inv_square(mp_kap);
        }
      }
      if (lflag_col_eff >= 1 && lflag_col_eff <= 3) mp_tauv = lflag_col_tau_shared == 1 ? rep_vector(col_tau[1], C) : col_tau;
      if (lflag_col_eff == 3) {
        matrix[C, C] LM;
        mp_Wt = diag_post_multiply(mp_Wc, mp_tauv);
        LM = cholesky_decompose(add_diag(crossprod(diag_pre_multiply(inv(sigma_jrc), mp_Wt)), 1));
        mp_DiW = mdivide_left_tri_low(LM, diag_pre_multiply(inv(s2), mp_Wt)')';
        mp_ldM = sum(log(diagonal(LM)));
      }
      for (j in 1:n_areas) {
       if (lflag_mp_seq == 1) {
        // ---- sequential expected table: margins and interior logits are the coordinates ----
        vector[Dm1_model] Ej = E_rc;                                           // hierarchy mean for this area
        matrix[R, C] qEj = mp_qE;                                              // its table
        if (lflag_E_cov == 1 || lflag_col_eff == 1 || lflag_col_eff == 2 || lflag_col_eff == 4 || lflag_col_cov == 1) {
          if (lflag_col_cov == 1) Ej += mp_Wc * ((lflag_col_kappa_shared == 1 ? rep_vector(col_kappa[1], C) : col_kappa) .* col_x[j]');
          if (lflag_E_cov == 1) Ej += E_cov_gamma .* E_cov_x[j]';
          if (lflag_col_eff == 1) Ej += mp_Wc * (mp_tauv .* col_eff_raw[n_col_groups > 0 ? col_group[j] : j]);
          if (lflag_col_eff == 2) Ej += mp_Wc * col_eff_raw[n_col_groups > 0 ? col_group[j] : j];
          if (lflag_col_eff == 4) Ej += mp_Wc * col_eff_known[j]';
          if (lflag_mp_seq_anchor == 1) {
            vector[R * C] clrEj = V_ilr_model * Ej;
            if (lflag_row_decompose == 1) {
              for (r in 1:R) qEj[r] = softmax(clrEj[((r - 1) * C + 1):(r * C)])';
            } else {
              vector[R * C] pEj = softmax(clrEj);
              for (r in 1:R) for (c in 1:C) qEj[r, c] = pEj[(r - 1) * C + c];
            }
          }
        }
        vector[R] wv;                                                          // expected row shares
        vector[C] m;                                                           // expected column shares
        matrix[R + 1, C] al;
        vector[R * C] logq;
        vector[Dm1_model] bb;
        real ljm = 0;
        int off = 0;
        if (lflag_row_decompose == 1) {
          wv = rm_prop[j]';                                                    // rows conditioned on the data
        } else {
          if (lflag_mp_row_whiten == 1) wv = softmax(append_row(mp_r_obs[j] + mp_Tr[j] * mp_beta[j][1:(R - 1)], 0));    // rows are parameters; linear change of variable, constant Jacobian
          else wv = softmax(append_row(mp_r_obs[j] + mp_beta[j][1:(R - 1)], 0));
          ljm += sum(log(wv));
          off = R - 1;
        }
        if (lflag_mp_beta_whiten == 1) m = softmax(append_row(mp_s_obs[j] + mp_T[j] * mp_beta[j][(off + 1):(off + C - 1)], 0));   // linear change of variable: constant Jacobian
        else m = softmax(append_row(mp_s_obs[j] + mp_beta[j][(off + 1):(off + C - 1)], 0));
        ljm += sum(log(m));
        if (lflag_mp_interior == 1) {
          // ---- interior fixed at E_rc's: the table is E_rc's table raked to this area's margins ----
          // b_j - E_rc lies in the span of mp_Wm; the move coordinates have precision P under the normal on b_j.
          // Density of the margin parameters = that normal of the moves x |d moves / d margin log-ratios|
          // (the margin log-ratios are linear in mp_beta, so nothing else is needed).
          matrix[R, C] TE;
          vector[(R - 1) * (C - 1)] lamE;
          vector[Dm1_model] dd;
          if (lflag_row_decompose == 1) TE = diag_pre_multiply(wv, qEj);
          else TE = qEj;
          lamE = alloc_table_inv(TE);
          al = alloc_table(wv, m, lamE, 0);
          for (r in 1:R) for (c in 1:C)
            logq[(r - 1) * C + c] = log(al[r, c]) - (lflag_row_decompose == 1 ? log(wv[r]) : 0);
          bb = V_ilr_model' * logq;
          b_mp[j] = bb';
          dd = bb - Ej;
          mp_beta_lp += -0.5 * dot_product(dd, dd ./ s2) + sum(log(diagonal(mp_LP)))
                        + mp_move_logjac(al[1:R, 1:C], mp_Hr, mp_Hc, lflag_row_decompose);
          for (r in 1:(R - 1)) for (c in 1:(C - 1)) mp_lambda_E[j, r, c] = lamE[(r - 1) * (C - 1) + c];
        } else {
        vector[(R - 1) * (C - 1)] lam = mp_z[j];
        if (lflag_mp_seq_scale == 2) lam = rep_vector(0, (R - 1) * (C - 1));   // built below, after the anchor
        if (lflag_mp_seq_scale == 1) {
          vector[(R - 1) * (C - 1)] sc_z = sqrt(square(mp_A[j]) * s2);     // prior scale of each logit implied by sigma at the reference table
          lam = sc_z .* mp_z[j];
          ljm += sum(log(sc_z));                                            // Jacobian of the rescaling
        }
        if (lflag_mp_seq_anchor == 1) {
          // table implied by E_rc alone (with this area's rows if rows are conditioned), and its logits.
          // A pure location shift of mp_z, so the Jacobian is unchanged.
          matrix[R, C] TE;
          matrix[R - 1, C - 1] aE;
          if (lflag_row_decompose == 1) TE = diag_pre_multiply(wv, qEj);
          else TE = qEj;
          if (lflag_seq_or_expected >= 2) {
            lam += seq_inv_ordered(TE, lflag_seq_or_expected, mp_row_order, mp_rem_col);
          } else {
            aE = inverse_alloc(TE, hinge_delta_floor, hinge_delta_min, kink_delta_expected, lflag_seq_or_expected);
            for (r in 1:(R - 1)) for (c in 1:(C - 1)) lam[(r - 1) * (C - 1) + c] += aE[r, c];
          }
        }
        if (lflag_mp_seq_scale == 2) {
          // lam currently = logits that realise E_rc's table (or 0 if no anchor). Table they give on THIS area's margins:
          int Kq = (R - 1) * (C - 1);
          matrix[R + 1, C] al0;
          if (lflag_seq_or_expected >= 2) al0 = seq_alloc_ordered(wv, m, lam, lflag_seq_or_expected, mp_row_order, mp_rem_col);
          else al0 = mp_seq_alloc(wv, m, lam, hinge_delta_floor, hinge_delta_min, kink_delta_expected, lflag_seq_or_expected);
          vector[R * C] logq0;
          vector[Dm1_model] b0;
          for (r in 1:R) for (c in 1:C)
            logq0[(r - 1) * C + c] = log(al0[r, c]) - (lflag_row_decompose == 1 ? log(wv[r]) : 0);
          b0 = V_ilr_model' * logq0;
          if (lflag_mp_scale_fast == 2) {
            // low-rank form: F F' = Mi with F = L0 (I + Q (S - I) Q') / sigma_b
            real sb = sigma_jrc[mp_sig_first[mp_gb]];
            vector[Kq] t;
            vector[Kq] y;
            t = mp_LB[j] * ((Ej - b0) ./ s2);
            if (mp_nu > 0) {
              matrix[mp_nu, mp_nu] S;
              vector[mp_nu] q;
              vector[mp_nu] u;
              S = cholesky_decompose(add_diag(to_matrix(mp_RG[j] * mp_wr, mp_nu, mp_nu), 1));
              q = mp_Q[j]' * t;
              y = sb * (t + mp_Q[j] * (mdivide_left_tri_low(S, q) - q)) + mp_z[j];
              u = mp_Q[j]' * y;
              lam += sb * (mp_L0i[j]' * (y + mp_Q[j] * (mdivide_right_tri_low(u', S)' - u)));
              ljm += -sum(log(diagonal(S)));
            } else {
              lam += sb * (mp_L0i[j]' * (sb * t + mp_z[j]));
            }
            ljm += Kq * log(sb) - mp_ld0[j];
          } else {
            // scale used for mp_z: sigma^w * sigma0^(1-w) per coordinate (w = 1 everywhere gives full non-centring)
            vector[Dm1_model] s2e = exp(2 * (mp_nc_w .* log(sigma_jrc) + (1 - mp_nc_w) .* log(mp_sigma0)));
            matrix[Kq, Kq] Mi;                                                       // precision of the logits used for the scaling
            matrix[Kq, Kq] Lp;
            if (lflag_mp_scale_fast == 1) {
              Mi = inv(s2[mp_sig_first[1]]) * mp_C[j, 1];
              for (g in 2:sigma_n_groups) Mi += inv(s2[mp_sig_first[g]]) * mp_C[j, g];
            } else {
              Mi = mp_B[j]' * diag_pre_multiply(inv(s2e), mp_B[j]);
            }
            if (lflag_col_eff == 3) Mi -= tcrossprod(mp_B[j]' * mp_DiW);          // precision of the logits under the full covariance
            if (lflag_mp_interior == 2) Mi = mp_ik2 * Mi + (1 - mp_ik2) * tcrossprod(mp_B[j]' * mp_DiW2);
            Lp = cholesky_decompose(0.5 * (Mi + Mi'));
            // conditional mean of the logits given the margins (one GLS step from the anchor) + non-centred deviation
            if (lflag_col_eff == 3) {
              vector[Dm1_model] vv = Ej - b0;
              lam += Lp' \ (mdivide_left_tri_low(Lp, mp_B[j]' * (vv ./ s2 - mp_DiW * (mp_DiW' * vv))) + mp_z[j]);
            } else if (lflag_mp_interior == 2) {
              vector[Dm1_model] vv = Ej - b0;
              lam += Lp' \ (mdivide_left_tri_low(Lp, mp_B[j]' * (mp_ik2 * (vv ./ s2) + (1 - mp_ik2) * (mp_DiW2 * (mp_DiW2' * vv)))) + mp_z[j]);
            } else {
              lam += Lp' \ (mdivide_left_tri_low(Lp, mp_B[j]' * ((Ej - b0) ./ s2e)) + mp_z[j]);
            }
            ljm += -sum(log(diagonal(Lp)));
          }
        }
        if (lflag_seq_or_expected >= 2) al = seq_alloc_ordered(wv, m, lam, lflag_seq_or_expected, mp_row_order, mp_rem_col);
        else al = mp_seq_alloc(wv, m, lam, hinge_delta_floor, hinge_delta_min, kink_delta_expected, lflag_seq_or_expected);
        for (r in 1:R) for (c in 1:C)
          logq[(r - 1) * C + c] = log(al[r, c]) - (lflag_row_decompose == 1 ? log(wv[r]) : 0);
        bb = V_ilr_model' * logq;
        b_mp[j] = bb';
        // hierarchy prior on b + log|d b / d (mp_beta, mp_z)|  (closed form, up to a data-only constant)
        if (lflag_mp_interior == 2) {
          // normal on b with the non-move part of the deviation scaled by kappa (the same constants as normal_lpdf, so kappa = 1 is identical)
          vector[Dm1_model] dd = bb - Ej;
          mp_beta_lp += -0.5 * (mp_ik2 * dot_product(dd, dd ./ s2) + (1 - mp_ik2) * dot_self(mp_DiW2' * dd))
                        - sum(log(sigma_jrc)) - (Dm1_model - mp_n_mv) * log(mp_kap) - 0.5 * Dm1_model * log(2 * pi())
                        + ljm + al[R + 1, 1] - sum(logq);
        } else {
          mp_beta_lp += normal_lpdf(bb | Ej, sigma_jrc) + ljm + al[R + 1, 1] - sum(logq);
        }
        if (lflag_col_eff == 3) {
          vector[Dm1_model] dd = bb - Ej;
          vector[C] uu = mp_DiW' * dd;
          real inc = 0.5 * dot_self(uu) - mp_ldM;
          mp_beta_lp += inc;
          mp_col_lp += inc;
          mp_col_eff_mean[j] = mp_tauv .* (mp_Wt' * (dd ./ s2 - mp_DiW * uu));
        }
        if (C >= 3) { mp_kink_loss_max = fmax(mp_kink_loss_max, al[R + 1, 2]); mp_kink_cells += al[R + 1, 3]; }
        for (r in 1:(R - 1)) for (c in 1:(C - 1)) mp_lambda_E[j, r, c] = lam[(r - 1) * (C - 1) + c];
        }
       } else {
        vector[Dm1_model] d0 = E_rc - mp_b_ref[j];
        matrix[n_pin, Dm1_model] GS = diag_post_multiply(mp_G[j], s2);                 // G Sigma
        matrix[Dm1_model - n_pin, Dm1_model] NS = diag_post_multiply(mp_N[j]', s2);    // N' Sigma
        matrix[n_pin, n_pin] Sbb = GS * mp_G[j]';
        matrix[Dm1_model - n_pin, n_pin] Sgb = NS * mp_G[j]';
        matrix[Dm1_model - n_pin, Dm1_model - n_pin] Sgg = NS * mp_N[j];
        matrix[Dm1_model - n_pin, n_pin] Kg;
        matrix[Dm1_model - n_pin, Dm1_model - n_pin] Cc;
        vector[n_pin] mu_b = mp_G[j] * d0;
        vector[Dm1_model - n_pin] gam;
        Sbb = 0.5 * (Sbb + Sbb');
        Kg = mdivide_right_spd(Sgb, Sbb);           // Sigma_gb Sigma_bb^-1
        Cc = Sgg - Kg * Sgb';                       // Sigma_gg|b
        Cc = 0.5 * (Cc + Cc');
        {
          vector[Dm1_model - n_pin] gmean = mp_N[j]' * d0 + Kg * (mp_beta[j] - mu_b);
          matrix[Dm1_model - n_pin, Dm1_model - n_pin] Lc = cholesky_decompose(Cc);
          if (lflag_mp_gamma_centred == 1) gam = mp_z[j];                    // centred interior
          else                             gam = gmean + Lc * mp_z[j];       // non-centred interior

          if (lflag_mp_exact == 0) {
            // linearised: mp_beta = G (b - b_ref)
            b_mp[j] = (mp_b_ref[j] + mp_P[j] * mp_beta[j] + mp_N[j] * gam)';
            mp_beta_lp += multi_normal_cholesky_lpdf(mp_beta[j] | mu_b, cholesky_decompose(Sbb));
            if (lflag_mp_gamma_centred == 1)
              mp_beta_lp += multi_normal_cholesky_lpdf(mp_z[j] | gmean, Lc);
          } else {
            // exact: mp_beta = s(b) - s(b_ref), the actual column log-ratio deviation.
            // Solve for the P-coordinate bet such that s(b_ref + P bet + N gam) = s_ref + mp_beta.
            vector[n_pin] bet = mp_beta[j];                                  // linear solution as start
            vector[Dm1_model] bb = mp_b_ref[j] + mp_P[j] * bet + mp_N[j] * gam;
            matrix[n_pin, Dm1_model + 1] sg = mp_colshare(bb, rm_prop[j]', V_ilr_model, R, C);
            matrix[n_pin, n_pin] A = sg[, 2:(Dm1_model + 1)] * mp_P[j];
            for (it in 1:mp_newton_iters) {
              bet -= A \ (sg[, 1] - mp_s_ref[j] - mp_beta[j]);
              bb = mp_b_ref[j] + mp_P[j] * bet + mp_N[j] * gam;
              sg = mp_colshare(bb, rm_prop[j]', V_ilr_model, R, C);
              A  = sg[, 2:(Dm1_model + 1)] * mp_P[j];
            }
            if (max(fabs(sg[, 1] - mp_s_ref[j] - mp_beta[j])) > 1e-6)
              reject("margin param: Newton solve did not converge in area ", j);
            b_mp[j] = bb';
            // prior on b itself + Jacobian of (mp_beta, mp_z) -> b :  |det| = det(Lc) / |det(A)|  (x const)
            mp_beta_lp += normal_lpdf(bb | E_rc, sigma_jrc) - log_determinant(A);
            if (lflag_mp_gamma_centred == 0) mp_beta_lp += sum(log(diagonal(Lc)));
          }
        }
       }
      }
    }

  } else if (lflag_rot_llrep == 1) {

    if (lflag_rot_agg == 1){
        LLrep_plus_log_volume_xform[1:Dm1_model] = ROT_agg * LLrep_plus_log_volume_raw[1:Dm1_model];
    }

    // --- centring/noncentring transform on raw vector ---
    // agg block: indices 1:Dm1_model
    for (s in 1:Dm1_model) {
      if (lflag_agg_noncent[s] == 1) {
        LLrep_plus_log_volume_xform[s] =
          sqrt(n_areas) * E_rc[s] + sigma_jrc[s] * LLrep_plus_log_volume_raw[s];
      }
      // centred: xform already initialised to raw (identity), nothing to do
    }

    // dev block: indices after vol block
    {
      int dev_start = Dm1_model + (n_areas - 1);
      for (k in 1:(n_areas - 1)) {
        int base = dev_start + (k - 1) * Dm1_model;
        for (s in 1:Dm1_model) {
          if (lflag_dev_noncent[k, s] == 1) {
            int idx = base + s;
            LLrep_plus_log_volume_xform[idx] =
              sigma_jrc[s] * LLrep_plus_log_volume_raw[idx];
          }
        }
      }
    }

    // single rotation
    LLrep_plus_log_volume = ROT_red * LLrep_plus_log_volume_xform;

    log_grand_volume = log_volume_raw;
    vol_coeffs = LLrep_plus_log_volume[(n_areas * Dm1_model + 1):Dtot_m1_reduced];
    vol_dev    = V_area * vol_coeffs;

  } else {

    vol_coeffs = LLrep_plus_log_volume_raw[(n_areas * Dm1_model + 1):Dtot_m1_reduced];
    log_volume = append_row(log_volume_raw, vol_coeffs);
    log_grand_volume = log_sum_exp(log_volume);

  }

  // --- per-area reconstruction ---
  for (j in 1:n_areas) {

    // 1. get free (within-row) ILR coordinates for this area
    row_vector[Dm1_model] llrep_model_j;
    if (lflag_margin_param == 1) {
      llrep_model_j = b_mp[j];                 // margin param: log_volume already set above
    } else if (lflag_rot_llrep == 1) {
      llrep_model_j = to_row_vector(
        LLrep_plus_log_volume[((j-1)*Dm1_model + 1):(j*Dm1_model)]);
      log_volume[j] = log_grand_volume + vol_dev[j];
    } else {
      llrep_model_j = to_row_vector(
        LLrep_plus_log_volume_raw[((j-1)*Dm1_model + 1):(j*Dm1_model)]);
    }

    // 2. reconstruct full log-composition via p[r,c] = agg_row_prop[j,r] * q[r,c]
    row_vector[R*C] clr_model = llrep_model_j * V_ilr_model';
    vector[R*C] log_p_vec;

    if (lflag_row_decompose == 1) {
      // existing senc-style row-block loop, using rm_prop, unchanged
      int idx = 0;
      for (r in 1:R) {
        vector[C] log_q_r;
        for (c in 1:C) { idx += 1; log_q_r[c] = clr_model[idx]; }
        real log_norm = log_sum_exp(log_q_r);
        for (c in 1:C) log_p_vec[(r-1)*C + c] = log(rm_prop[j, r]) + log_q_r[c] - log_norm;
        }
      } else {
        // Scotland-style: whole-table softmax, no row split, rm_prop unused
        vector[R*C] lp = log_softmax(to_vector(clr_model));
        for (k in 1:(R*C)) log_p_vec[k] = lp[k];
      }

    // 3. project to full Dm1 ILR coordinates (derived dims populated automatically)
    LLrep_jrc[j] = to_row_vector(V_ilr_full' * log_p_vec);

    // 4. downstream — unchanged
    row_vector[R*C] clr_table    = LLrep_jrc[j] * V_ilr_full';
    vector[R*C]     log_composition = log_softmax(to_vector(clr_table));
    {
      int idx = 0;
      for (r in 1:R) {
        for (c in 1:C) {
          idx += 1;
          composition_arr[j, r, c]         = exp(log_composition[idx]);
          log_expected_cell_values[j, r, c] =
            log_volume[j] + log_composition[idx];
        }
      }
    }
    if (lflag_neutral_logit == 2) {

        matrix[R, C] comp_j;
        for (r in 1:R) for (c in 1:C) comp_j[r, c] = composition_arr[j, r, c];
        matrix[R - 1, C - 1] full_anchor;
        if (lflag_margin_param == 1 && lflag_mp_seq == 1) {
          for (r in 1:(R - 1)) for (c in 1:(C - 1)) full_anchor[r, c] = mp_lambda_E[j, r, c];   // margin param (seq): expected logits directly
        } else {
          full_anchor = inverse_alloc(comp_j, hinge_delta_floor, hinge_delta_min, kink_delta_realised, 0);
        }
        for (r in 1:R - 1) {
          for (c in 1:C - 1) {
            anchor_lambda[j, r, c] = full_anchor[r, c];
          }
        }
      } else if(lflag_neutral_logit ==3){
        for(r in 1:R - 1){
          for(c in 1:C - 1){
            anchor_lambda[j, r, c] = anchor_lambda_E_rc[r, c];
          }
        }
      }

if (lflag_fit_type == 0 && lflag_real_rowwise == 1) {
    // realised table row by row on the observed margins, non-empty rows only.
    // logits = the expected table's logits for the same rows + fixed scaling * lambda_raw (Poisson noise, roughly unit scale)
    int nk = rl_n[j];
    int Kr = nk * (C - 1);
    lambda[j] = rep_array(0.0, R - 1, C - 1);
    for (r in 1:R) for (c in 1:C) cell_values[j, r, c] = 0;
    if (nk >= 1) {
      vector[Kr] anc;
      vector[Kr] dev = rl_S[j][1:Kr, 1:Kr] * lambda_raw[(param_count_from[j] + 1):(param_count_from[j] + Kr)];
      vector[Kr] lr;
      if (lflag_real_scale == 1) {
        vector[Kr] sdv;
        vector[C] rest;                                        // expected cells of the later rows, by column
        for (c in 1:C) rest[c] = exp(log_expected_cell_values[j, rl_row[j, nk + 1], c]);
        for (kk in 1:nk) {
          int k = nk + 1 - kk;
          int last = mp_rem_col[rl_row[j, k]];
          int i = 0;
          vector[C] x;
          real vref;
          for (c in 1:C) x[c] = exp(log_expected_cell_values[j, rl_row[j, k], c]);
          vref = inv(x[last]) + inv(rest[last]);
          for (c in 1:C) if (c != last) { i += 1; sdv[(k - 1) * (C - 1) + i] = sqrt(inv(x[c]) + inv(rest[c]) + vref); }
          rest += x;
        }
        dev = sdv .* dev;
        rl_lj += sum(log(sdv));                                // Jacobian of the state-dependent scaling
      }
      vector[nk + 1] wr;
      array[nk + 1] int ord;
      array[nk + 1] int rem;
      matrix[nk + 2, C] alr;
      for (k in 1:nk) for (i in 1:(C - 1)) anc[(k - 1) * (C - 1) + i] = mp_lambda_E[j, rl_step[j, k], i];
      lr = anc + dev;
      for (k in 1:(nk + 1)) { wr[k] = row_margins[j, rl_row[j, k]]; ord[k] = k; rem[k] = mp_rem_col[rl_row[j, k]]; }
      alr = seq_alloc_ordered(wr, col_margins[j]', lr, 3, ord, rem);
      for (k in 1:(nk + 1)) for (c in 1:C) cell_values[j, rl_row[j, k], c] = alr[k, c];
      rl_lj += alr[nk + 2, 1];
      if (lflag_cpois_norm == 1)
        for (k in 1:(nk + 1)) for (c in 1:C)
          rl_norm += cpois_lognorm(log_expected_cell_values[j, rl_row[j, k], c], cz_t0, cz_h, cz_v, cz_d);
      for (i in 1:Kr) {
        lambda_vec[param_count_from[j] + i] = dev[i];
        neutral_logit_flat[param_count_from[j] + i] = anc[i];
      }
      for (k in 1:nk) for (i in 1:(C - 1)) lambda[j, k, i] = lr[(k - 1) * (C - 1) + i];
    } else {
      for (c in 1:C) cell_values[j, rl_row[j, 1], c] = col_margins[j, c];      // one non-empty row: the table is its column margins
    }
} else if(lflag_fit_type == 0) {
    lambda_vec = (lflag_rot_lambda > 0) ? ROT_lambda * lambda_raw : lambda_raw;
    lambda[j] = rep_array(0.0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        int idx =param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c;
        real this_nl;
        if(lflag_neutral_logit==0||lflag_neutral_logit==1){
          this_nl = neutral_logit_array[j, r, c];
        } else if(lflag_neutral_logit == 2||lflag_neutral_logit == 3){
          this_nl = anchor_lambda[j, active_row_map[j, r], active_col_map[j, c]];
        }
        neutral_logit_flat[idx] = this_nl;
        if(lflag_lambda_centred==1){
          lambda[j, r, c] = lambda_vec[idx];
        } else {
          lambda[j, r, c] = this_nl + lambda_vec[idx];
        }
        }
      }
    }


      if(lflag_fit_type==3) {
          lambda[j] = rep_array(0.0, R - 1, C - 1);
        }

  } // end per-area loop
}
if((lflag_fit_type == 0 && lflag_real_rowwise == 0)||lflag_fit_type==3){
   cell_values = ss_assign_cvals_lp(n_areas, R, C, row_margins, col_margins, lambda, hinge_delta_floor, hinge_delta_min, slack_tol, kink_delta_realised);

}



  for(j in 1:n_areas){
    for (r in 1:R){
      for(c in 1:C){
        if(lflag_fit_type == 1||lflag_fit_type==2){
          cell_values[j, r, c] = exp(log_expected_cell_values[j, r, c]);
        }

        // log_cv[j, r, c] = log(robust_hinge_floor_zero(cell_values[j, r, c], hinge_delta_min));
      }
    }
  }



}
model{
// matrix[R, C] overall_values;
// matrix[R, C] global_prop;


if(lflag_fit_type == 0||lflag_fit_type == 3) {
    row_vector[n_active_cells] obs_flat;
    vector[n_active_cells] log_rate_flat;
    for (idx in 1:n_active_cells) {
      obs_flat[idx]  = cell_values[active_j[idx], active_r[idx], active_c[idx]];
      log_rate_flat[idx] = log_expected_cell_values[active_j[idx], active_r[idx], active_c[idx]];
    }
    target += realpoisson_lograte_lpdf(obs_flat | log_rate_flat);
} else {
  vector[n_areas*C] expected_col_margins;
  vector[n_areas*R] expected_row_margins;

for(j in 1:n_areas){
  for(c in 1:C){
      expected_col_margins[(j - 1)*C + c] = sum(cell_values[j, 1:R, c]);
  }
  for(r in 1:R){
      expected_row_margins[(j - 1)*R + r] = sum(cell_values[j, r, 1:C]);
  }
}

if(lflag_fit_type == 1){
  col_margins_flat ~ poisson(expected_col_margins);
  row_margins_flat ~ poisson(expected_row_margins);
} else if(lflag_fit_type == 2){
  row_margins_flat ~ poisson(expected_row_margins);   // scale anchor only
  for(j in 1:n_areas){
  vector[C] mu_j = rep_vector(0, C);
  matrix[C, C] Sigma_j = rep_matrix(0, C, C);
  for(r in 1:R){
    real R_jr = row_margins[j, r];
    if (R_jr > 0){
      vector[C] q_jr;
      for (c in 1:C){
        q_jr[c] = composition_arr[j, r, c]/rm_prop[j, r];
      }
      if (lflag_row_decompose == 0) q_jr /= sum(q_jr);     // rows are parameters: use the model's own row-conditional shares
      mu_j += R_jr *q_jr;
      Sigma_j += R_jr * (diag_matrix(q_jr) - q_jr * q_jr');

    }
  }

  // drop the last column: totals sum to the known area total, so the
  // full C-dim covariance is singular by exactly one dimension
  col_margins[j, 1:(C-1)]' ~ multi_normal(mu_j[1:(C-1)], Sigma_j[1:(C-1), 1:(C-1)]);


  }
}

}

// for(r in 1:R){
//   for(c in 1:C){
//     overall_values[r, c] = sum(cell_values[1:n_areas, r, c]);
//     global_prop[r,c] = overall_values[r,c] / fmax(global_rm[1, r], 1e-10);
//   }
// }

  log_volume ~ normal(0, 10);

if (lflag_margin_param == 1) {
  // hierarchy b_j ~ N(E_rc, diag sigma^2) in (beta, z) coordinates
  target += mp_beta_lp;
  if (lflag_soft_smallcell == 1 && lflag_fit_type == 2)
    for (j in 1:n_areas) for (r in 1:R) if (row_margins[j, r] > 0) for (c in 1:C)
      target += smallcell_scale * cpois_lognorm(log_expected_cell_values[j, r, c], cz_t0, cz_h, cz_v, cz_d);
  if (lflag_real_rowwise == 1) target += rl_lj;      // Jacobian of the realised tables
  if (lflag_real_rowwise == 1 && lflag_cpois_norm == 1) target += -rl_norm;
  if (lflag_mp_gamma_centred == 0 && lflag_mp_exact == 0 && lflag_mp_seq == 0) for (j in 1:n_areas) mp_z[j] ~ std_normal();
} else if(lflag_rot_llrep == 0){
  for(j in 1:n_areas){
    LLrep_jrc[j, agg_free_dim] ~ normal(E_rc, sigma_jrc);
  }
} else if(lflag_rot_llrep == 1){
  for (s in 1:Dm1_model){
    if(lflag_agg_noncent[s]==1){
      LLrep_plus_log_volume_raw[s] ~ std_normal();
    } else {
      LLrep_plus_log_volume_xform[s] ~ normal(sqrt(n_areas)*E_rc[s], sigma_jrc[s]);
    }
  }
  {
    int dev_start = Dm1_model + (n_areas - 1);
    for (k in 1:(n_areas - 1)){
      int base = dev_start + (k - 1) * Dm1_model;
      for (s in  1:Dm1_model){
        int idx = base + s;
        if (lflag_dev_noncent[k, s] == 1){
          LLrep_plus_log_volume_raw[idx] ~ std_normal();
        } else {
          LLrep_plus_log_volume_raw[idx] ~ normal(0, sigma_jrc[s]);
        }
      }
    }
  }

}

  if (mp_kappa_free == 1) mp_kappa ~ lognormal(prior_mp_kappa_a, prior_mp_kappa_b);
  if (lflag_E_cov == 1) E_cov_gamma ~ normal(0, prior_E_cov_scale);
  if (lflag_col_cov == 1) col_kappa ~ normal(0, prior_col_kappa_scale);
  if (lflag_col_eff >= 1 && lflag_col_eff <= 3) {
    col_tau ~ normal(0, prior_col_tau_scale);
    if (lflag_col_eff == 1 || lflag_col_eff == 2)
    for (j in 1:(n_col_groups > 0 ? n_col_groups : n_areas)) {
      if (lflag_col_eff == 1) col_eff_raw[j] ~ std_normal();
      else col_eff_raw[j] ~ normal(0, lflag_col_tau_shared == 1 ? rep_vector(col_tau[1], C) : col_tau);
    }
  }
  // E_rc
  for (g in 1:E_rc_n_groups) {
    int m = E_rc_group_mode[g];
    if (m >= 1) E_rc_group_mu[E_rc_mu_pos[g]] ~ normal(E_rc_group_prior_a[g], E_rc_group_prior_b[g]);
    if (m == 2 || m == 3)
      E_rc_group_sigma[E_rc_sigma_pos[g]] ~ gamma(E_rc_group_tau_a[g], E_rc_group_tau_b[g]);
  }
  for (k in 1:E_rc_n_leaf) {
    int g = E_rc_group_id[E_rc_leaf_idx[k]];
    if (E_rc_group_mode[g] == 2)
      E_rc_raw_leaf[k] ~ normal(E_rc_group_mu[E_rc_mu_pos[g]], E_rc_group_sigma[E_rc_sigma_pos[g]]);
    else
      E_rc_raw_leaf[k] ~ std_normal();                       // modes 3 and 4
  }

  // sigma_jrc: family decides how the centre prior is expressed
  for (g in 1:sigma_n_groups) {
    int m = sigma_group_mode[g];
    if (m >= 1) {
      real lsig = sigma_group_mu[sigma_mu_pos[g]];
      if (lflag_family == 2)
        target += gamma_lpdf(exp(lsig) | sigma_group_prior_a[g], sigma_group_prior_b[g]) + lsig;
      else
        lsig ~ normal(sigma_group_prior_a[g], sigma_group_prior_b[g]);
    }
    if (m == 2 || m == 3)
      sigma_group_sigma[sigma_sigma_pos[g]] ~ gamma(sigma_group_tau_a[g], sigma_group_tau_b[g]);
  }
  for (k in 1:sigma_n_leaf) {
    int g = sigma_group_id[sigma_leaf_idx[k]];
    if (sigma_group_mode[g] == 2)
      sigma_raw_leaf[k] ~ normal(sigma_group_mu[sigma_mu_pos[g]], sigma_group_sigma[sigma_sigma_pos[g]]);
    else
      sigma_raw_leaf[k] ~ std_normal();
  }

if(lflag_fit_type == 0){

  if (lflag_lambda_centred == 1) {
    lambda_vec ~ normal(neutral_logit_flat, prior_lambda_raw_scale);
    } else {
    lambda_vec ~ normal(0, prior_lambda_raw_scale);
  }
}



    if(use_known_cells == 1){
      array[n_poss_cells] int known_cell_values_row_vector;
      vector[n_poss_cells] log_cv_row_vector;
      int counter_cell = 0;
      for(j in 1:n_areas){
        for(r in 1:R){
          for(c in 1:C){
            if(structural_zeros[j,r,c]==0){
              counter_cell += 1;
              log_cv_row_vector[counter_cell] = log_expected_cell_values[j, r, c];
              known_cell_values_row_vector[counter_cell] = known_cell_values[j, r, c];
                }
              }
            }
          }
      target += poisson_lpmf(known_cell_values_row_vector | exp(log_cv_row_vector));
    }
}
generated quantities{
  vector[(R * C) - 1] ilr_mean;
  vector[(R * C) - 1] ilr_var;
  vector[(R * C) - 1] ilr_n;

if(lflag_rawscw == 1||lflag_rawscw==0){
  for(k in 1:((R * C) - 1)){
        real s = 0;
        real s2 = 0;
        real n = 0;
        for(j in 1:n_areas){
            real ilr_val = LLrep_jrc[j, k];
            s += ilr_val;
            s2 += square(ilr_val);
            n += 1;
          }
          ilr_n[k] = n;
          ilr_mean[k] = (n > 0) ? s/n : 0;
          ilr_var[k] = (n > 1) ? s2/n - square(ilr_mean[k]) : 0;
        }
      }
      #include include/generateratesandsummaries.stan

}
