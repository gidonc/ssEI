
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
  vector[n_areas*Dm1_model + n_areas - 1] LLrep_plus_log_volume_raw;
  real log_volume_raw;
  vector[E_rc_n_leaf] E_rc_raw_leaf;
  vector[E_rc_n_mu] E_rc_group_mu;
  vector<lower=0>[E_rc_n_sigma] E_rc_group_sigma;

  vector[sigma_n_leaf] sigma_raw_leaf;
  vector[sigma_n_mu] sigma_group_mu;
  vector<lower=0>[sigma_n_sigma] sigma_group_sigma;


}
transformed parameters{
  vector[lflag_fit_type==0 ? n_param : 0] lambda_vec;
  vector[lflag_fit_type==0 ? n_param : 0] neutral_logit_flat;
  real lambda[lflag_fit_type==0||lflag_fit_type==3 ? n_areas : 0, R - 1, C -1]; // sequential cell weights
  vector[n_areas*Dm1_model + n_areas - 1] LLrep_plus_log_volume;
  matrix[n_areas, Dm1] LLrep_jrc;
  real log_grand_volume;
  vector[n_areas] log_volume;
  vector[n_areas*Dm1_model + n_areas - 1] LLrep_plus_log_volume_xform = LLrep_plus_log_volume_raw;

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
    anchor_lambda_E_rc = inverse_alloc(comp_E_rc, hinge_delta_floor, hinge_delta_min);
}


  {
  vector[n_areas - 1] vol_coeffs;
  vector[n_areas] vol_dev = rep_vector(0, n_areas);

  if (lflag_rot_llrep == 1) {

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
    if (lflag_rot_llrep == 1) {
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
        matrix[R - 1, C - 1] full_anchor = inverse_alloc(comp_j, hinge_delta_floor, hinge_delta_min);
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

if(lflag_fit_type == 0) {
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
if(lflag_fit_type == 0||lflag_fit_type==3){
   cell_values = ss_assign_cvals_lp(n_areas, R, C, row_margins, col_margins, lambda, hinge_delta_floor, hinge_delta_min, slack_tol);

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

if(lflag_rot_llrep == 0){
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
