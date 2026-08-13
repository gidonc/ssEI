
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
 int<lower=0> R_ll; // number of rows in log-linear representation
 int<lower=0> C_ll; // number of cols in log-linear representation
 matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
 matrix<lower=0>[n_areas, C] col_margins; // the column margins in each area
 matrix[C, C - 1] V_ilr_data; // basis matrix for ILR transformation (for ILR 3 case)
 int n_ilr_rows; // number of rows in the basic construction passed to Stan (for ILR 3 case)
 int<lower=0, upper=1> structural_zeros[n_areas, R, C];  // an array indicating any structural zeros in the data (may include whole rows, whole columns and/or individual cells)
 int<lower=0, upper=2> lflag_dist; // flag indicating whether to use poisson (0), multinomial (1) or negative binomial (2)paramertization
 int<lower=0, upper=3> lflag_area_re; // flag indicating whether the area mean simplex is uniform (0) or varies with area random effects which are normally distributed (1) or varies with area random effects which are multinormally distributed (non centred paramaterisation) (2) or varies with area random effects which are multinormally distributed (non centred LKJ Onion paramaterisation)
 int<lower  =0, upper=2> lflag_vary_sd; // flag indicating whether variance of area_cell parameters is: (0) shared across cells,  (1) varies by cell,  or (2) has a hierarchical model structure
 int<lower = 0, upper = 1> lflag_llmod_omit_jr; // flag indicating whether log-linear model should omit area * row interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jc; // flag indicating whether log-linear model should omit area * col interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jrc; // flag indicating whether log-linear model should omit area * row * column interaction
  int<lower = 0, upper =1> lflag_predictors_cm; // flag indicating whether to model columns as well as rows
  int<lower =0, upper = 1> lflag_noncentred; //flag indicating whether to use a centred (0) or non-centred parameterization;
  int<lower=0, upper = 1> lflag_noncentred_mat[n_areas, R_ll, C_ll]; //flag indicated whether cell is centred (0) or non-centred (1)
  int<lower = 0, upper = 1> lflag_rawscw; //flag indicating whether to use raw, or decentred version of sequential cell weights (lamdba coeffients)
  int<lower = 0, upper = 4> lflag_ll_rep; // flag indicating which log linear representation of the final tables should be used.
  // 0 = Additive Log-Ratio 1 (C - 1) log-ratios representing the composition of the (R - 1) the free rows of the matrix
  // 3 = Log Odds Ratios of the
  int<lower=0, upper=1> lflag_fix_E_rc;        // 1 = use fixed values, 0 = estimate
  int<lower=0, upper=1> lflag_fix_sigma_jrc;   // 1 = use fixed values, 0 = estimate
  matrix[R_ll, C_ll] E_rc_fixed;         // fixed values, ignored if fix_E_rc=0
  matrix<lower=0>[R_ll, C_ll] sigma_jrc_fixed;  // fixed values, ignored if fix_sigma_jrc=0
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
 matrix[R_ll, C_ll] E_rc_prior; // empirically informed prior centres for E_rc
 int<lower=0> known_cell_values[n_areas, R, C]; // for testing purposes
 int<lower=0, upper=1> use_known_cells; // for testing purposes
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
  int is_ilr = 0;
  real sigma_constrain = .001;
  int n_free_areas_rc[R - 1, C - 1]; // count of free areas for each (r, c) for decentred lambdas
  int dev_start_rc[R-1, C-1]; // start for deviation parameters
  int dev_idx = (R-1)*(C-1); // index for deviation parameters in (r, c) order
  real hinge_delta_floor = 1e-10;
  real hinge_delta_min = 1e-10;
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

  matrix[C, C - 1] V_ilr = rep_matrix(0.0, C, C - 1);

  if(lflag_ll_rep == 1){
    // Build the Orthonormal Basis matrix via Sequential Binary Partitioning
    for (i in 1:(C - 1)) {
      real r_i = i;
      real weight_left = 1.0 / r_i * sqrt(r_i / (r_i + 1.0));
      real weight_right = -sqrt(r_i / (r_i + 1.0));
      for (j in 1:i) {
        V_ilr[j, i] = weight_left;
      }
      V_ilr[i + 1, i] = weight_right;
    }
  } else if(lflag_ll_rep == 2){

  // Build the Orthonormal Basis matrix via CLR contrasts (geometric mean reference)
  // Each coordinate contrasts category c against the geometric mean of all C categories,
  // giving a symmetric basis with no privileged reference category.
    {
    // CLR contrast matrix
    // Column c has (1 - 1/C) in row c and (-1/C) everywhere else,
    // so each column sums to zero (centred log-ratio structure)
      matrix[C, C - 1] A = rep_matrix(-1.0 / C, C, C - 1);
      for (c in 1:(C - 1)) {
        A[c, c] = 1.0 - 1.0 / C;
      }
      // Gram-Schmidt orthonormalisation column by column
      for (c in 1:(C - 1)) {
        vector[C] v = A[, c];
        for (k in 1:(c - 1)) {
          v = v - dot_product(V_ilr[, k], v) * V_ilr[, k];
        }
        V_ilr[, c] = v / sqrt(dot_self(v));
      }
    }
  } else if (lflag_ll_rep == 4) {
    V_ilr = V_ilr_data;
  }

  if(lflag_fix_sigma_jrc == 1){
    K_sigmas = 0;
    K_sigma_c_sigma = 0;
  } else if(lflag_vary_sd == 0){
    K_sigmas = 1;
    K_sigma_c_sigma = 0;
  } else if(lflag_vary_sd == 1){
    if(lflag_rawscw==1){
      K_sigmas = R_ll*C_ll;
    } else{
      K_sigmas = R_ll*C_ll;
      // K_sigmas = (R - 1)*(C - 1);
    }
    K_sigma_c_sigma = 0;
  } else if(lflag_vary_sd == 2){
      if(lflag_rawscw==1){
      K_sigmas = R_ll*C_ll;
    } else{
      K_sigmas = R_ll*C_ll;
      // K_sigmas = (R - 1)*(C - 1);
    }
    K_sigma_c_sigma = 1;
  }



}
parameters{
  // real gamma_raw[n_param_gamma]; // the raw atomic cell deviations (leading to sequential cell weights)
  // vector[n_areas] mu_scale; //Global table volatility on the raw scale
  // real alpha_raw[n_param_alpha]; // Total sum of (free_R[j] - 1) // Row-specific volatility
  // real log_beta_raw[n_param_beta];  // Total sum of (free_C[j] - 2) // Column-specific volatility
  vector[n_param] lambda_raw;
  matrix[lflag_rawscw == 0 ? (R - 1):0, (C - 1)] lambda_mu_rc;
  // array[n_areas] vector[R*C - 1] E_jrc_raw;
  // array[n_areas] matrix[R, C - 1] E_jrc_raw;
  // array[n_areas] vector[R - 1] E_jr_raw;
  // array[n_areas] vector[C - 1] E_jc_raw;
  // vector[n_areas - 1] E_j_raw;
  // vector[C - 1] E_c_raw;
  // vector[R - 1] E_r_raw;
  // matrix[R, C - 1] E_rc_raw;
  // real E_mu;
  // real E_j_mu;
  // matrix[n_areas, C-1] ALR_jrc_R_raw; // probabilistic ALR for last row

  // vector<lower=0, upper= 1> [n_areas*has_theta] theta;
  // vector<lower=0> [has_theta] phi;
  // real<lower=0> sigma_j;
  // matrix<lower=0>[R - 1, C - 1] sigma_jrc;

  real<lower=0> sigma_jrc_raw[K_sigmas];
  real<lower=0> sigma_c_sigma[(lflag_vary_sd == 2) ? 1 : 0];
  vector[(lflag_vary_sd == 2) ? 1 : 0] sigma_c_mu;
  real LLrep_raw[n_areas, R_ll, C_ll];
  vector[n_areas * R_ll] row_effect_raw;
  real mu_re;
  real<lower=0> sigma_re;

  matrix[(lflag_fix_E_rc||lflag_rawscw==0) ? 0 : R_ll, lflag_fix_E_rc? 0 : C_ll] E_rc_raw;
  // vector[R - 1] E_r_raw;

  // real lambda_zero[n_zero_cells];
  // vector<lower=0>[C] sigma_c;

  // vector<lower=0>[K_j] sigma_j_all;
  // real<lower=0, upper = .1> hinge_delta_floor;
  // real<lower=0, upper = .1> hinge_delta_min;
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real LLrep_jrc[n_areas, R_ll, C_ll];
  real lambda_mu_LLrep[lflag_rawscw==0 ? R_ll : 0, C_ll];
  real<lower=0> cell_values[n_areas, R, C];
  real<lower=0> expected_cell_values[n_areas, R, C];
  real log_cv[n_areas, R, C] ;
  matrix<lower=0>[R_ll, C_ll] sigma_jrc;
  real composition_arr[n_areas, R, C] = rep_array(0.0, n_areas, R, C);
  matrix[R_ll, C_ll] ilr_mean;
  matrix[R_ll, C_ll] ilr_var;
  matrix[R_ll, C_ll] ilr_n;
  matrix[R_ll, C_ll] E_rc;
  real det_J_raw = 0;
  real sigma_scale[lflag_rawscw==0?n_areas:0, R - 1, C - 1];
  real sigma_sq_row_out[lflag_rawscw==0?n_areas:0, R_ll, C_ll];
  real sens_sq_out[lflag_rawscw==0?n_areas:0, R_ll, C_ll];
  real tmp_J_mu_out[lflag_rawscw==0?n_areas:0, R_ll, C_ll];
  matrix[n_areas, R_ll] row_effect;
  // vector[C_ll] LLrep_star[n_areas, R_ll];

  for (j in 1:n_areas) {
    for (r in 1:R_ll) {
      row_effect[j, r] = row_effect_raw[(j - 1) * R_ll + r];
    }
  }


  // matrix[(R - 1)*(C - 1), (R - 1) * (C - 1)] Chol_Sigma_lambda = diag_pre_multiply(sigma_lambda_raw, Lcorr_raw);



  // array[lflag_rawscw == 0 ? n_areas: 0, R - 1, C - 1] real lambda_dev = rep_array(0.0, n_areas, R - 1, C - 1);
  // matrix[lflag_rawscw == 0 ? R_ll: 0, C_ll] E_rc_baseline;

if(lflag_fix_sigma_jrc==1){
  sigma_jrc = sigma_jrc_fixed;
    } else if(lflag_rawscw==1||lflag_rawscw==0){
      int s = 0;
      for(r in 1:R_ll){
          for(c in 1:C_ll){
              s += 1;
              if(lflag_vary_sd == 0){
                sigma_jrc = rep_matrix(sigma_jrc_raw[1], R_ll, C_ll);
              } else if(lflag_vary_sd==1){
                sigma_jrc[r,c] = sigma_jrc_raw[s];
              } else if(lflag_vary_sd == 2){
                sigma_jrc[r,c] = sigma_jrc_raw[s];
              }
          }
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
  } else if (lflag_rawscw==1){
      E_rc = E_rc_raw;
  }

  for (j in 1:n_areas) {
    for (r in 1:R_ll) {
      row_effect[j, r] = row_effect_raw[(j - 1) * R_ll + r];
      for (c in 1:C_ll){
        if(lflag_noncentred_mat[j, r, c]==1){
          LLrep_jrc[j, r, c] = E_rc[r, c] + sigma_jrc[r, c] * LLrep_raw[j, r, c];
        } else {
          LLrep_jrc[j, r, c] = LLrep_raw[j, r, c];
        }

      }

      row_vector[C] clr_row = to_row_vector(LLrep_jrc[j, r]) * V_ilr';
      vector[C] composition = softmax(to_vector(clr_row));
      for (c in 1:C) composition_arr[j, r, c] = composition[c];

      for (c in 1:C)
        expected_cell_values[j, r, c] = exp(row_effect[j, r]) * composition[c];
    }
  }
if(lflag_noncentred == 1){
  cell_values = ss_assign_cvals_cpanchor_lp(n_areas, R, C, row_margins, col_margins, lambda, composition_arr, hinge_delta_floor, hinge_delta_min);
} else if(lflag_noncentred == 0){
  cell_values = ss_assign_cvals_wzeros_hinge_lp(n_areas, R, C, row_margins, col_margins, lambda, hinge_delta_floor, hinge_delta_min);
}
  for(j in 1:n_areas){
    for (r in 1:R){
      for(c in 1:C){
        log_cv[j, r, c] = fmax(log(cell_values[j, r, c]), 1e-10);
      }
    }
  }


if(lflag_rawscw == 1||lflag_rawscw==0){
  for(r in 1:R_ll){
    for(c in 1:C_ll){
        real s = 0;
        real s2 = 0;
        real n = 0;
        for(j in 1:n_areas){
            if(row_margins[j,r] > 0 && structural_zeros[j,r,c] == 0){
                real ilr_val = LLrep_jrc[j,r,c];
                s  += ilr_val;
                s2 += square(ilr_val);
                n  += 1;
            }
        }
        ilr_n[r,c]    = n;
        ilr_mean[r,c] = (n > 0) ? s / n : 0;
        ilr_var[r,c]  = (n > 1) ? s2/n - square(ilr_mean[r,c]) : 0;
    }
  }
}

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
for(j in 1:n_areas){
    for(c in 1:C){
        implied_col_var[j,c] = 0;
        for(r in 1:R_ll){
            for(c_ll in 1:C_ll){
                real J = row_margins[j,r] * global_prop[r,c] *
                         (V_ilr[c,c_ll] -
                          dot_product(global_prop[r,:], col(V_ilr,c_ll)));
                implied_col_var[j,c] += square(J) * square(sigma_jrc[r,c_ll]);
            }
        }
    }
}

if(lflag_predictors_cm){
  for(j in 1:n_areas){
    for(c in 1:C){
        real implied = dot_product(row_margins[j,:], col(global_prop,c));
        col_margins[j,c] ~ normal(implied, sqrt(fmax(implied_col_var[j,c], 1e-10)));
    }
  }
}

// margin-generating layer — every row
  for (j in 1:n_areas)
      target += realpoisson_lpdf(row_margins[j] | exp(to_vector(row_effect[j])));
{
    row_vector[n_active_cells] obs_flat;
    vector[n_active_cells] rate_flat;
    for (idx in 1:n_active_cells) {
      obs_flat[idx]  = cell_values[active_j[idx], active_r[idx], active_c[idx]];
      rate_flat[idx] = fmax(expected_cell_values[active_j[idx], active_r[idx], active_c[idx]], 1e-8);
    }
    target += realpoisson_lpdf(obs_flat | rate_flat);
  }
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
row_effect_raw ~ normal(mu_re, sigma_re);
mu_re ~ normal(0, 10);
sigma_re ~ normal(0, 10);
for(j in 1:n_areas){
  for(r in 1:R_ll){
    for(c in 1:C_ll){
      if(lflag_noncentred_mat[j, r, c]==1){
        LLrep_raw[j, r, c] ~ std_normal();
      } else {
        LLrep_raw[j, r, c] ~ normal(E_rc[r, c], sigma_jrc[r, c]);
      }
    }
  }
}
if(lflag_rawscw == 1){
  for(r in 1:R_ll){
      E_rc[r, 1:C_ll] ~ normal(0, prior_mu_re_scale);
    }
} else if(lflag_rawscw == 0){
  for(r in 1:(R - 1)){
    lambda_mu_rc[r, 1:C_ll] ~ normal(0, prior_mu_re_scale);
  }
}

if(lflag_fix_sigma_jrc == 1){

}  else if(lflag_vary_sd == 2){
    sigma_c_mu ~ normal(0, prior_sigma_c_mu_scale);
    sigma_c_sigma ~ normal(0, prior_sigma_c_scale);
    // to_vector(sigma_jrc_raw) ~ normal(0, prior_sigma_c_scale);
    for(s in 1:R_ll*C_ll){
        sigma_jrc_raw[s] ~ normal(sigma_c_mu[1], sigma_c_sigma[1]);
    }
} else {
    to_vector(sigma_jrc_raw) ~ normal(prior_sigma_mu, prior_sigma_c_scale);
    // sigma_jrc_raw ~ normal(0, prior_sigma_c_scale);
}

    // hinge_delta_floor ~ normal(0, .1);
    // hinge_delta_min ~ normal(0, .1);
    if(lflag_rawscw == 1){
      lambda_raw ~ normal(0, prior_lambda_raw_scale);
    } else {
      // lambda_raw ~ std_normal();
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
              log_cv_row_vector[counter_cell] = log_cv[j, r, c];
              if(use_known_cells == 1){
                known_cell_values_row_vector[counter_cell] = known_cell_values[j, r, c];
                }
              }
            }
          }
        }
      target += poisson_lpmf(known_cell_values_row_vector | exp(log_cv_row_vector) + 1e-10);
    }
}
generated quantities{
    #include include/generateratesandsummaries.stan

}
