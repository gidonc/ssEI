
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
  matrix[(R * C), (R * C) - 1] V_ilr; // basis matrix for ILR transformation
  matrix[(R * C) - 1, (R * C) - 1] ROT; //ilr rotation from raw to model basis
 // vector<lower=0>[(R * C) - 1] sigma_llrep;
 // vector[(R * C) - 1] mu_llrep;
 int<lower=0, upper=1> structural_zeros[n_areas, R, C];  // an array indicating any structural zeros in the data (may include whole rows, whole columns and/or individual cells)
 int<lower=0, upper=2> lflag_dist; // flag indicating whether to use poisson (0), multinomial (1) or negative binomial (2)paramertization
 int<lower=0, upper=2> lflag_family; // flag indicating whether scale of LLrep distribution is based on log-normal (0), cauchy(1) or Gamma (2) family
 int<lower=0, upper=1> lflag_E_rc_hier; // does E_rc have hyperparameters (1) or is it based on E_rc_prior (0)
 int<lower=0, upper=3> lflag_area_re; // flag indicating whether the area mean simplex is uniform (0) or varies with area random effects which are normally distributed (1) or varies with area random effects which are multinormally distributed (non centred paramaterisation) (2) or varies with area random effects which are multinormally distributed (non centred LKJ Onion paramaterisation)
 int<lower  =0, upper=2> lflag_vary_sd; // flag indicating whether variance of area_cell parameters is: (0) shared across cells,  (1) varies by cell,  or (2) has a hierarchical model structure
 int<lower = 0, upper=2> lflag_neutral_logit; // flag indicating the neutral_logit for the allocation process when sequential weights are all zero. (0) gives equality across the rows [-log(cols_remaining)] (1) gives the independent table logit(col_margin/sum_remaining_col_margins) (2) gives the current composition implied by LLrep structure.
 int<lower = 0, upper = 1> lflag_llmod_omit_jr; // flag indicating whether log-linear model should omit area * row interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jc; // flag indicating whether log-linear model should omit area * col interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jrc; // flag indicating whether log-linear model should omit area * row * column interaction
  int<lower =0, upper = 1> lflag_pin_row_effect; // flag indicating the model on row effects
  int<lower = 0, upper =1> lflag_predictors_cm; // flag indicating whether to model columns as well as rows
  int<lower =0, upper = 1> lflag_noncentred; //flag indicating whether to use a centred (0) or non-centred parameterization;
  int<lower=0, upper = 1> lflag_noncentred_mat[n_areas, (R * C) - 1]; //flag indicated whether cell is centred (0) or non-centred (1)
  int<lower = 0, upper = 1> lflag_rawscw; //flag indicating whether to use raw, or decentred version of sequential cell weights (lamdba coeffients)
  int<lower = 0, upper = 4> lflag_ll_rep; // flag indicating which log linear representation of the final tables should be used.
  int<lower=0, upper = 1> lflag_rot_llrep;
  int<lower=0, upper = 1> lflag_rot_E_rc;
  // 0 = Additive Log-Ratio 1 (C - 1) log-ratios representing the composition of the (R - 1) the free rows of the matrix
  // 3 = Log Odds Ratios of the
  int<lower=0, upper=1> lflag_fix_E_rc;        // 1 = use fixed values, 0 = estimate
  int<lower=0, upper=1> lflag_fix_sigma_jrc;   // 1 = use fixed values, 0 = estimate
  vector[(R * C) - 1] E_rc_fixed;         // fixed values, ignored if fix_E_rc=0
  vector<lower=0>[R*C - 1] sigma_jrc_fixed;  // fixed values, ignored if fix_sigma_jrc=0
  int<lower = 0, upper = 1> lflag_lambda_raw_offset; // prior on lambda_raw (0) or lambda + neutral_logit (1)
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
  int K_sigmas;
  int K_sigma_c_sigma;
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
  matrix[n_areas, R] rm_prop;
  matrix[n_areas, R] rm_log;
  matrix[n_areas, C] cm_prop;
  matrix[n_areas, C] cm_log;
  matrix[1, R] global_rm;
  matrix[1, C] global_cm;
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
  real neutral_logit_flat[n_param];
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
          neutral_logit_flat[count] = logit(this_share);
          neutral_logit_array[j, r, c] = neutral_logit_flat[count];
        } else{
          neutral_logit_flat[count] = -log(free_C[j] - c);
          neutral_logit_array[j, r, c] = neutral_logit_flat[count];
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
  for(r in 1:R) global_rm[1, r] = sum(row_margins[1:n_areas, r]);
  for(c in 1:C) global_cm[1, c] = sum(col_margins[1:n_areas, c]);


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


  if(lflag_fix_sigma_jrc == 1){
    K_sigmas = 0;
    K_sigma_c_sigma = 0;
  } else if(lflag_vary_sd == 0){
    K_sigmas = 1;
    K_sigma_c_sigma = 0;
  } else if(lflag_vary_sd == 1){
    if(lflag_rawscw==1){
      K_sigmas = (R * C) - 1;
    } else{
      K_sigmas = (R * C) - 1;
    }
    K_sigma_c_sigma = 0;
  } else if(lflag_vary_sd == 2){
      if(lflag_rawscw==1){
      K_sigmas = (R * C) - 1;
    } else{
      K_sigmas = (R * C) - 1;
    }
    K_sigma_c_sigma = 1;
  }

  // int D = R * C;
  //
  // // 1. Generate Helmert bases inside Stan
  // matrix[R, R - 1] Vll_helmert_R = make_helmert_basis(R);
  // matrix[C, C - 1] Vll_helmert_C = make_helmert_basis(C);
  //
  // // 2. Column vectors for grand total scaling
  // matrix[R, 1] j_R = rep_matrix(1.0 / sqrt(R), R, 1);
  // matrix[C, 1] j_C = rep_matrix(1.0 / sqrt(C), C, 1);
  //
  // // 3. Kronecker blocks
  // matrix[D, R - 1] V_row_block = kronecker_prod(Vll_helmert_R, j_C);
  // matrix[D, C - 1] V_col_block = kronecker_prod(j_R, Vll_helmert_C);
  // matrix[D, (R - 1) * (C - 1)] V_resid_block = kronecker_prod(Vll_helmert_R, Vll_helmert_C);
  //
  // // 4. Combine into full space basis and build rotation matrix
  // matrix[D, D - 1] V_row_col = append_col(
  //   append_col(V_row_block, V_col_block),
  //   V_resid_block
  // );
  //
  // matrix[D - 1, D - 1] ROT = V_row_col' * V_ilr;

}
parameters{
  vector[n_param] lambda_raw;
  matrix[lflag_rawscw == 0 ? (R - 1):0, (C - 1)] lambda_mu_rc;
  vector[(lflag_family != 2) ? K_sigmas : 0] sigma_jrc_raw;

// gamma (2): direct positive parameter, no log-transform
vector<lower=0>[(lflag_family == 2) ? K_sigmas : 0] sigma_jrc_direct;

// partial-pooling hyperparameters, additive families
real<lower=0> sigma_c_sigma[(lflag_vary_sd == 2 && lflag_family != 2) ? 1 : 0];
vector[(lflag_vary_sd == 2 && lflag_family != 2) ? 1 : 0] sigma_c_mu;

  real E_rc_mu[(lflag_E_rc_hier == 1) ? 1 : 0];       // shared centre across all cells
  real<lower=0> E_rc_sigma[(lflag_E_rc_hier == 1) ? 1 : 0];

// partial-pooling hyperparameters, gamma (mean + shape/concentration parameterization)
real<lower=0> sigma_c_mu_gamma[(lflag_vary_sd == 2 && lflag_family == 2) ? 1 : 0];
real<lower=0> sigma_c_shape[(lflag_vary_sd == 2 && lflag_family == 2) ? 1 : 0];

  real LLrep_raw[n_areas, (R * C) - 1];
  vector[n_areas] log_volume;

  vector[(lflag_fix_E_rc||lflag_rawscw==0) ? 0 : (R * C) - 1] E_rc_raw;
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  matrix[n_areas, (R * C) - 1] LLrep_jrc;
  real<lower=0> cell_values[n_areas, R, C];
  real log_expected_cell_values[n_areas, R, C];
  real log_cv[n_areas, R, C] ;
  vector<lower=0>[(R * C) - 1] sigma_jrc;
  real composition_arr[n_areas, R, C] = rep_array(0.0, n_areas, R, C);
  vector[(R * C) - 1] ilr_mean;
  vector[(R * C) - 1] ilr_var;
  vector[(R * C) - 1] ilr_n;
  vector[(R * C) - 1] E_rc;
  real det_J_raw = 0;
  real sigma_scale[lflag_rawscw==0?n_areas:0, R - 1, C - 1];
  real sigma_sq_row_out[lflag_rawscw==0?n_areas:0, (R * C) - 1];
  real sens_sq_out[lflag_rawscw==0?n_areas:0, (R * C) - 1];
  real tmp_J_mu_out[lflag_rawscw==0?n_areas:0, (R * C) - 1];

if(lflag_fix_sigma_jrc==1){
  sigma_jrc = sigma_jrc_fixed;
} else if(lflag_rawscw==1||lflag_rawscw==0){
  if(lflag_vary_sd == 0){
    // single shared value, broadcast to every cell
    real shared_val = (lflag_family == 2) ? sigma_jrc_direct[1] : exp(sigma_jrc_raw[1]);
    sigma_jrc = rep_vector(shared_val, (R * C) - 1);
  } else if(lflag_family == 2){
    sigma_jrc = sigma_jrc_direct;
  } else {
    sigma_jrc = exp(sigma_jrc_raw);
  }
  }

if(lflag_rawscw == 1||lflag_rawscw==0){
    for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        // lambda[j, r, c] = lambda_raw[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c]*sigma_scale_r[r];
        lambda[j, r, c] = lambda_raw[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
        }
      }
    }
   }

  if(lflag_fix_E_rc==1){
      E_rc = E_rc_fixed;
  } else if (lflag_rot_E_rc == 1){
    E_rc = ROT' * E_rc_raw;
  } else{
    E_rc = E_rc_raw;
  }

  for (j in 1:n_areas) {
    vector[(R * C) - 1] rowcol_val;
    for(s in 1:((R * C) - 1)){
      if(lflag_noncentred_mat[j, s] == 1){
        rowcol_val[s] = E_rc[s] + sigma_jrc[s] * LLrep_raw[j, s];
      } else {
        rowcol_val[s] = LLrep_raw[j, s];
      }
    }
    if(lflag_rot_llrep == 1){
      LLrep_jrc[j] = to_row_vector(ROT' * rowcol_val);
    } else {
      LLrep_jrc[j] = to_row_vector(rowcol_val);
    }

      row_vector[R * C] clr_table = to_row_vector(LLrep_jrc[j]) * V_ilr';
      vector[R * C] log_composition = log_softmax(to_vector(clr_table));
      int idx = 0;
      for(r in 1:R){
        for(c in 1:C){
          idx += 1;
          composition_arr[j, r, c] = exp(log_composition[idx]);
          log_expected_cell_values[j, r, c] = log_volume[j] + log_composition[idx];
        }
      }

  }
    cell_values = ss_assign_cvals_cpanchor_lp(n_areas, R, C, row_margins, col_margins, lambda, composition_arr, neutral_logit_array, hinge_delta_floor, hinge_delta_min, slack_tol, lflag_neutral_logit);

  for(j in 1:n_areas){
    for (r in 1:R){
      for(c in 1:C){
        log_cv[j, r, c] = log(robust_hinge_floor_zero(cell_values[j, r, c], hinge_delta_min));
      }
    }
  }


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

  // for(r in 1:R_ll){
  //   for(c in 1:C_ll){
  //       real s = 0;
  //       real s2 = 0;
  //       real n = 0;
  //       for(j in 1:n_areas){
  //           if(row_margins[j,r] > 0 && structural_zeros[j,r,c] == 0){
  //               real ilr_val = LLrep_jrc[j,r,c];
  //               s  += ilr_val;
  //               s2 += square(ilr_val);
  //               n  += 1;
  //           }
  //       }
  //       ilr_n[r,c]    = n;
  //       ilr_mean[r,c] = (n > 0) ? s / n : 0;
  //       ilr_var[r,c]  = (n > 1) ? s2/n - square(ilr_mean[r,c]) : 0;
  //   }
  // }
// }

}
model{
matrix[R, C] overall_values;
matrix[R, C] global_prop;
matrix[n_areas, C] implied_col_var;


for(r in 1:R){
  for(c in 1:C){
    overall_values[r, c] = sum(cell_values[1:n_areas, r, c]);
    global_prop[r,c] = overall_values[r,c] / fmax(global_rm[1, r], 1e-10);
  }
}
// for(j in 1:n_areas){
//     for(c in 1:C){
//         implied_col_var[j,c] = 0;
//         for(r in 1:R_ll){
//             for(c_ll in 1:C_ll){
//                 real J = row_margins[j,r] * global_prop[r,c] *
//                          (V_ilr[c,c_ll] -
//                           dot_product(global_prop[r,:], col(V_ilr,c_ll)));
//                 implied_col_var[j,c] += square(J) * square(sigma_jrc[r,c_ll]);
//             }
//         }
//     }
// }
//
// if(lflag_predictors_cm){
//   for(j in 1:n_areas){
//     for(c in 1:C){
//         real implied = dot_product(row_margins[j,:], col(global_prop,c));
//         col_margins[j,c] ~ normal(implied, sqrt(fmax(implied_col_var[j,c], 1e-10)));
//     }
//   }
// }
//

{
    row_vector[n_active_cells] obs_flat;
    vector[n_active_cells] log_rate_flat;
    for (idx in 1:n_active_cells) {
      obs_flat[idx]  = cell_values[active_j[idx], active_r[idx], active_c[idx]];
      log_rate_flat[idx] = log_expected_cell_values[active_j[idx], active_r[idx], active_c[idx]];
    }
    target += realpoisson_lograte_lpdf(obs_flat | log_rate_flat);
  }
  log_volume ~ normal(0, 10);
//
// for (r in 1:R_ll){
//   for(k in 1:C_ll){
//     real ss = 0;
//     for (j in 1:n_areas){
//       if(row_margins[j, r]>0 && structural_zeros[j, r, k] == 0){
//         real extra_var = 0;
//         for (c in 1:C){
//           if(col_margins[j, c]>0){
//             extra_var += square(V_ilr[c, k])/fmax(cell_values[j, r, c], 1e-6);
//           }
//         }
//         real total_var = square(sigma_jrc[r, k]) + extra_var;
//         ss += normal_lpdf(LLrep_jrc[j, r, k] | E_rc[r, k], sqrt(total_var));
//       }
//     }
//     target += ss;
//   }
// }
// if(lflag_rawscw == 1||lflag_rawscw == 0){
//   for(r in 1:R_ll){
//       for(c in 1:C_ll){
//         if(ilr_n[r,c] > 1){
//           real n  = ilr_n[r,c];
//           real mu = ilr_mean[r,c];
//           real v  = ilr_var[r,c];
//           // Sufficient statistic normal log likelihood
//           target += -n * log(sigma_jrc[r,c])
//                     - n * v / (2 * square(sigma_jrc[r,c]))
//                     - n * square(mu - E_rc[r,c]) / (2 * square(sigma_jrc[r,c]));
//           }
//         }
//     }
// }
// if(lflag_rawscw == 0){
//   for(r in 1:R_ll){
//     for(c in 1:C_ll){
//       real n = dev_n[r,c];
//       real v = dev_var[r, c];
//       if(dev_n[r, c] > 1){
//         target += -(n - 1) * log(sigma_jrc[r, c])
//                   -(n - 1) * dev_ss[r, c]/(2*square(sigma_jrc[r, c]));
//       }
//     }
//   }
// }
for(j in 1:n_areas){
  for(k in 1:((R * C) - 1)){
      if(lflag_noncentred_mat[j, k]==1){
        LLrep_raw[j, k] ~ std_normal();
      } else {
        LLrep_raw[j,k] ~ normal(E_rc[k], sigma_jrc[k]);
      }
  }
}
if(lflag_rawscw == 1){
  if(lflag_E_rc_hier == 1){
    E_rc_mu ~ normal(0, prior_mu_re_scale);
    E_rc_sigma ~ gamma(prior_gamma_shape, prior_gamma_rate);
    E_rc ~ normal(E_rc_mu[1], E_rc_sigma[1]);
  } else{
    E_rc ~ normal(E_rc_prior, prior_mu_re_scale);


  }
} else if(lflag_rawscw == 0){
  for(r in 1:(R - 1)){
    lambda_mu_rc[r, 1:(C - 1)] ~ normal(0, prior_mu_re_scale);
  }
}
if(lflag_fix_sigma_jrc == 1){
  // nothing
} else if(lflag_vary_sd == 0){
  // single shared parameter, straightforward prior regardless of family
  if(lflag_family == 2){
    sigma_jrc_direct[1] ~ gamma(prior_gamma_shape, prior_gamma_rate);
  } else if(lflag_family == 0){
    sigma_jrc_raw[1] ~ normal(prior_sigma_mu, prior_sigma_c_scale);
  } else {
    sigma_jrc_raw[1] ~ cauchy(prior_sigma_mu, prior_sigma_c_scale);
  }
} else if(lflag_family == 2){
  // --- Gamma family, vary (1) or partial (2) ---
  if(lflag_vary_sd == 2){
    sigma_c_mu_gamma ~ normal(0, prior_sigma_c_mu_scale) T[0,];
    sigma_c_shape ~ normal(0, prior_sigma_c_scale) T[0,];
    for(s in 1:K_sigmas){
      sigma_jrc_direct[s] ~ gamma(sigma_c_shape[1], sigma_c_shape[1] / sigma_c_mu_gamma[1]);
    }
  } else {  // vary_sd == 1
    to_vector(sigma_jrc_direct) ~ gamma(prior_gamma_shape, prior_gamma_rate);
  }
} else {
  // --- lognormal (0) or cauchy (1), vary (1) or partial (2) ---
  if(lflag_vary_sd == 2){
    sigma_c_mu ~ normal(0, prior_sigma_c_mu_scale);
    sigma_c_sigma ~ normal(0, prior_sigma_c_scale);
    for(s in 1:K_sigmas){
      if(lflag_family == 0){
        sigma_jrc_raw[s] ~ normal(sigma_c_mu[1], sigma_c_sigma[1]);
      } else {
        sigma_jrc_raw[s] ~ cauchy(sigma_c_mu[1], sigma_c_sigma[1]);
      }
    }
  } else {  // vary_sd == 1
    if(lflag_family == 0){
      to_vector(sigma_jrc_raw) ~ normal(prior_sigma_mu, prior_sigma_c_scale);
    } else {
      to_vector(sigma_jrc_raw) ~ cauchy(prior_sigma_mu, prior_sigma_c_scale);
    }
  }
}
// if(lflag_fix_sigma_jrc == 1){
//
// }  else if(lflag_vary_sd == 2){
//     sigma_c_mu ~ normal(0, prior_sigma_c_mu_scale);
//     sigma_c_sigma ~ normal(0, prior_sigma_c_scale);
//     // to_vector(sigma_jrc_raw) ~ normal(0, prior_sigma_c_scale);
//     for(s in 1:R_ll*C_ll){
//         sigma_jrc_raw[s] ~ normal(sigma_c_mu[1], sigma_c_sigma[1]);
//     }
// } else {
//     to_vector(sigma_jrc_raw) ~ normal(prior_sigma_mu, prior_sigma_c_scale);
//     // sigma_jrc_raw ~ normal(0, prior_sigma_c_scale);
// }

    // hinge_delta_floor ~ normal(0, .1);
    // hinge_delta_min ~ normal(0, .1);
    if(lflag_lambda_raw_offset==1){
          if(lflag_neutral_logit == 0||lflag_neutral_logit == 1){
      lambda_raw ~ normal(-neutral_logit_flat, prior_lambda_raw_scale);
    } else if(lflag_neutral_logit==2){
      real comp_logit_offset[n_param];
      int counter = 0;
      for(j in 1:n_areas){
        for(r in 1:(free_R[j] - 1)){
          vector[free_C[j]] comp_row;
          for(k in 1:free_C[j]) {
           comp_row[k] = composition_arr[j, active_row_map[j, r], active_col_map[j, k]];
         }

          for(c in 1:(free_C[j] - 1)){
            counter += 1;
           real remaining_comp_sum = sum(comp_row[c:free_C[j]]);
           real this_share = comp_row[c]/remaining_comp_sum;
           // real neutral_logit = -log(cols_remaining);
           comp_logit_offset[counter] = -1 * logit(this_share);
          }
        }
      }
      lambda_raw ~ normal(comp_logit_offset, prior_lambda_raw_scale);
    }

    } else if(lflag_lambda_raw_offset == 0){
      lambda_raw ~ normal(0, prior_lambda_raw_scale);
    }
    // lambda_raw ~ normal(lambda_raw_mu, 1);
    // lambda_raw_sigma ~ normal(0, prior_sigma_c_scale);
    // lambda_raw_mu ~ normal(0, prior_sigma_c_mu_scale);



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
    #include include/generateratesandsummaries.stan

}
