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

  // vector sum_to_zero(vector y){
  //   int N = size(y);
  //   int n;
  //   vector[N] omega;
  //   vector[N] S;
  //   vector[N + 1] x;
  //
  //   for(n_rev in 1:(N - 1)){
  //     n = N - n_rev;
  //     omega[n + 1] = y[n + 1]/sqrt((n + 1)*(n + 2));
  //   }
  // S[N] = 0;
  // for(n_rev in 1:(N - 1)){
  //   n = N - n_rev;
  //   S[n] = S[n + 1] + omega[n + 1];
  // }
  // for(n_rev in 1:N){
  //   n = N + 1 - n_rev;
  //   x[n + 1] = S[n] - (n * y[n]/sqrt(n * (n + 1)));
  // }
  // x[1] = S[1] + y[1]/sqrt(2);
  // return(x);
  // }

}
data{
 int<lower=0> n_areas;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
 matrix<lower=0>[n_areas, C] col_margins; // the column margins in each area
 int<lower=0, upper=1> structural_zeros[n_areas, R, C];  // an array indicating any structural zeros in the data (may include whole rows, whole columns and/or individual cells)
 int<lower=0, upper=2> lflag_dist; // flag indicating whether to use poisson (0), multinomial (1) or negative binomial (2) paramertization
 int<lower=0, upper=3> lflag_area_re; // flag indicating whether the area mean simplex is uniform (0) or varies with area random effects which are normally distributed (1) or varies with area random effects which are multinormally distributed (non centred paramaterisation) (2) or varies with area random effects which are multinormally distributed (non centred LKJ Onion paramaterisation)
 int<lower  =0, upper=2> lflag_vary_sd; // flag indicating whether variance of area_cell parameters is: (0) shared across cells,  (1) varies by cell,  or (2) has a hierarchical model structure
 int<lower = 0, upper = 1> lflag_llmod_omit_jr; // flag indicating whether log-linear model should omit area * row interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jc; // flag indicating whether log-linear model should omit area * col interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jrc; // flag indicating whether log-linear model should omit area * row * column interaction
  int<lower = 0, upper =1> lflag_predictors_cm; // flag indicating whether to model columns as well as rows
  int<lower =0, upper = 1> lflag_noncentred; //flag indicating whether to use a centred (0) or non-centred parameterization;
  int<lower = 0, upper = 1> lflag_rawscw; //flag indicating whether to use raw, or row, column and table version of sequential cell weights (lamdba coeffients)
 real<lower=0> prior_mu_re_scale; // prior for scale of mu_re (mean row effect)
 real<lower=0> prior_mu_ce_scale; // prior for scale of col_effect (mean column effect)
 real<lower=0> prior_sigma_c_scale; //prior for scale of sigma_c (or sigma_c_sigma if lflag_vary_sd == 2)
 real<lower=0> prior_sigma_c_mu_scale; //prior for scale of sigma_c_mu (only if lflag_vary_sd == 2)
 real<lower=0> prior_sigma_ce_scale; //prior of scale for sigma_ce
 real<lower=0> prior_sigma_re_scale; //prior of scale for sigma_re
 real<lower=0> prior_cell_effect_scale; //prior of scale for average cell effects
 matrix[R - 1, C - 1] E_rc_prior; // empirically informed prior centres for E_rc
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
  vector[n_areas] tot_log;
  real<lower=0> prior_phi_scale;
  int n_margin_sigmas;
  int n_jrc_sigmas;
  int n_table_sigmas;
  real sigma_constrain = .001;

  prior_phi_scale = 100;

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

    // matrix[C, C - 1] V_ilr;
  // Build the Orthonormal Basis matrix via CLR contrasts (geometric mean reference)
  // Each coordinate contrasts category c against the geometric mean of all C categories,
  // giving a symmetric basis with no privileged reference category.
  // {
    // Step 1: CLR contrast matrix
    // Column c has (1 - 1/C) in row c and (-1/C) everywhere else,
    // so each column sums to zero (centred log-ratio structure)
  //   matrix[C, C - 1] A = rep_matrix(-1.0 / C, C, C - 1);
  //   for (c in 1:(C - 1)) {
  //     A[c, c] = 1.0 - 1.0 / C;
  //   }
  //   // Step 2: Gram-Schmidt orthonormalisation column by column
  //   for (c in 1:(C - 1)) {
  //     vector[C] v = A[, c];
  //     for (k in 1:(c - 1)) {
  //       v = v - dot_product(V_ilr[, k], v) * V_ilr[, k];
  //     }
  //     V_ilr[, c] = v / sqrt(dot_self(v));
  //   }
  // }

}
parameters{
  // real gamma_raw[n_param_gamma]; // the raw atomic cell deviations (leading to sequential cell weights)
  // vector[n_areas] mu_scale; //Global table volatility on the raw scale
  // real alpha_raw[n_param_alpha]; // Total sum of (free_R[j] - 1) // Row-specific volatility
  // real log_beta_raw[n_param_beta];  // Total sum of (free_C[j] - 2) // Column-specific volatility
  vector[n_param] lambda_raw;
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
  real<lower=0> sigma_jrc_raw[(lflag_vary_sd == 0) ? 1 : (R - 1)*(C-1)];
  real<lower=0> sigma_c_sigma[(lflag_vary_sd == 2) ? 1 : 0];
  vector[(lflag_vary_sd == 2) ? 1 : 0] sigma_c_mu;
  // real<lower=0> sigma_jr;
  //real<lower=0> sigma_j;
  matrix[R - 1, C - 1] E_rc;
  // vector[R - 1] E_r_raw;

  real lambda_zero[n_zero_cells];
  // vector<lower=0>[C] sigma_c;

  // vector<lower=0>[K_j] sigma_j_all;
  real<lower=0, upper = .1> hinge_delta_floor;
  real<lower=0, upper = .1> hinge_delta_min;
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real LOR_jrc[n_areas, R - 1, C - 1];

  real<lower=0> cell_values[n_areas, R, C];
  real log_cv[n_areas, R, C] ;
  matrix<lower=0>[R - 1, C-1] sigma_jrc;

  vector[n_areas] E_j = rep_vector(0, n_areas);
  // array[n_areas] vector[R] E_jr = rep_array(rep_vector(0, R), n_areas);
  // array[n_areas] vector[C] E_jc = rep_array(rep_vector(0, C), n_areas);
  // array[n_areas] matrix[R, C] log_e_cell_values =  rep_array(rep_matrix(0, R, C), n_areas);
  // array[n_areas] matrix[R, C] log_r_cell_values =  rep_array(rep_matrix(0, R, C), n_areas);

  vector[R] E_r = rep_vector(0, R);
  // vector[C] E_c;


if(lflag_rawscw == 1){
    for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_raw[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
        }
      }
    }
   } else {
   vector[n_areas] mu_scale = lambda_raw[1:n_areas];
   vector[n_param_gamma] gamma_raw = lambda_raw[(n_areas + 1):(n_areas + n_param_gamma)];
   vector[n_param_alpha] log_alpha_raw = lambda_raw[(n_areas + n_param_gamma + 1):(n_areas + n_param_gamma + n_param_alpha)];
   vector[n_param_beta] log_beta_raw = lambda_raw[(n_areas + n_param_gamma + n_param_alpha + 1):(n_areas + n_param_gamma + n_param_alpha + n_param_beta)];


   for(j in 1:n_areas){
     lambda[j] = rep_array(0, R - 1, C - 1);

     int fr = free_R[j] - 1;
     int fc = free_C[j] - 1;

     if (fr > 0 && fc > 0) {
       vector[fr] a_scale = rep_vector(1.0, fr);
       vector[fc] b_scale = rep_vector(1.0, fc);
       matrix[fr, fc] gamma_mat = rep_matrix(0.0, fr, fc);

       if(fr > 1){
         vector[fr - 1] raw_a;
         for (i in 1:(fr - 1)) {
           raw_a[i] = log_alpha_raw[alpha_start[j] + i - 1];
        }
        a_scale = sum_to_zero(raw_a);

       }

      if (fc > 1) {
        vector[fc - 1] raw_b;
        for (i in 1:(fc - 1)) {
          raw_b[i] = log_beta_raw[beta_start[j] + i - 1];
        }
        b_scale = sum_to_zero(raw_b);
      }

      if (fr > 1 && fc > 1) {
        matrix[fr - 1, fc - 1] raw_gamma;
        for (r in 1:(fr - 1)) {
          for (c in 1:(fc - 1)) {
            raw_gamma[r, c] = gamma_raw[gamma_start[j] + (r - 1) * (fc - 1) + c - 1];
          }
        }

        matrix[fr - 1, fc] pass1_gamma;
        for (r in 1:(fr - 1)) {
          pass1_gamma[r] = to_row_vector(sum_to_zero(to_vector(raw_gamma[r])));
        }

        for (c in 1:fc) {
          gamma_mat[1:fr, c] = sum_to_zero(pass1_gamma[1:(fr - 1), c]);
        }
      }

      for (r in 1:fr) {
        for (c in 1:fc) {
          lambda[j, r, c] = a_scale[r] + b_scale[c] + gamma_mat[r, c];
        }
      }
    }
   }


 }

  // cell_values = ss_assign_cvals_wzeros_hinge_lp(n_areas, R, C, row_margins, col_margins, lambda, hinge_delta_floor, hinge_delta_min);
  // ALR_jrc = ss_assign_alr_wzeros_hinge_newparam_lp(n_areas, R, C, row_margins, col_margins, lambda, structural_zeros, hinge_delta_floor, hinge_delta_min);


// 1. Fetch the Log Cell Values from the engine
  log_cv = ss_assign_lor_wzeros_hinge_newparam_lp(
      n_areas, R, C,
      row_margins, col_margins,
      lambda, lambda_zero, zero_cell_map, structural_zeros,
      hinge_delta_floor, hinge_delta_min
  );

  // 2. Extract cell_values and LOR_jrc
  for (j in 1:n_areas) {
    for (r in 1:R) {
      for (c in 1:C) {
        // Exponentiate to get absolute cell counts
        if (log_cv[j, r, c] > -100.0) {
          cell_values[j, r, c] = exp(log_cv[j, r, c]);
        } else {
          cell_values[j, r, c] = 0.0; // Handle structural zeros gracefully
        }
      }
    }

    // Calculate the Difference-in-Differences interaction directly
    for (r in 1:(R - 1)) {
      for (c in 1:(C - 1)) {
         LOR_jrc[j, r, c] = log_cv[j, r, c]
                          - log_cv[j, r, C]
                          - log_cv[j, R, c]
                          + log_cv[j, R, C];
      }
    }
  }


if(lflag_vary_sd == 0){
    sigma_jrc = rep_matrix(sigma_jrc_raw[1], R - 1, C-1);
} else {
    int s = 0;
    for(r in 1:(R - 1)){
        for(c in 1:(C-1)){
            s += 1;
            sigma_jrc[r,c] = sigma_jrc_raw[s];
        }
    }
}
    // E_r[1:(R - 1)] = E_r_raw[1:(R - 1)];
    // if(lflag_noncentred == 1){
    //   E_j[1:(n_areas - 1)] = E_j_mu + sigma_j * E_j_raw;
    // } else {
    //   E_j[1:(n_areas - 1)] = E_j_raw[1:(n_areas - 1)];
    // }

        // E_mu_j - area mean of log cell values
//     int n_free = 0;
//     for(r in 1:R){
//         for(c in 1:C){
//             if(row_margins[j, r] > 0 && structural_zeros[j, r, c] == 0){
//                 E_mu_j += log_cv[j, r, c];
//                 n_free += 1;
//             }
//         }
//     }
//     E_mu_j = E_mu_j / n_free;
//     // E_j - deviation from null model expectation
//     E_j[j] = E_mu_j - (tot_log[j] - log(R*C));
//
//        for(r in 1:R){
//         if(row_margins[j, r] > 0){
//             real row_mean = mean(to_vector(log_cv[j, r, 1:C]));
//             E_jr[j, r] = row_mean - E_mu_j;
//         } else {
//             E_jr[j, r] = -200;
//         }
//     }
// }

// E_rc - global mean of E_jrc across areas
// for(r in 1:R){
//     for(c in 1:C){
//         real sum_E_jrc = 0;
//         int n_areas_rc = 0;
//         for(j in 1:n_areas){
//             if(row_margins[j, r] > 0 && structural_zeros[j, r, c] == 0){
//                 sum_E_jrc += E_jrc[j, r, c];
//                 n_areas_rc += 1;
//             }
//         }
//         E_rc[r, c] = n_areas_rc > 0 ? sum_E_jrc / n_areas_rc : -200;
//     }
// }

// E_r - global mean of E_jr across areas
// for(r in 1:R){
//     real sum_E_jr = 0;
//     int n_areas_r = 0;
//     for(j in 1:n_areas){
//         if(row_margins[j, r] > 0){
//             sum_E_jr += E_jr[j, r];
//             n_areas_r += 1;
//         }
//     }
//     E_r[r] = n_areas_r > 0 ? sum_E_jr / n_areas_r : -200;


    // for(j in 1:n_areas){
    //   for(r in 1:(R - 1)){
    //     for(c in 1:C){
    //       log_e_cell_values[j, r, c] = E_mu + E_j[j] + E_jr[j, r] + ALR_jrc[j, r, c];
    //     }
    //   }
    //   for(c in 1:C){
    //       log_e_cell_values[j, R, c] = E_mu + E_j[j] + E_jr[j, R] + ALR_jrc_R[j, c];
    //
    //   }
    // }


}
model{
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


// for(j in 1:n_areas){
//     for(r in 1:(R-1)){
//         ALR_jrc[j, r, 1:(C-1)] ~ normal(E_rc[r, 1:(C-1)], sigma_jrc[r, 1:(C-1)]);
//     }
// }
    int counter = 0;
    int n_all_cells = n_free_alr_cells + n_areas * (C - 1);
    vector[n_free_alr_cells] lor_vec;
    vector[n_free_alr_cells] e_rc_vec;
    vector[n_free_alr_cells] sigma_vec;
    row_vector[n_areas * C] cm_vec;
    vector[n_areas * C] e_cm_vec;
    for(j in 1:n_areas){
        for(r in 1:(R - 1)){
            for(c in 1:(C-1)){
                if(structural_zeros[j,r,c] == 0){
                    counter += 1;
                    lor_vec[counter] = LOR_jrc[j,r,c];
                    e_rc_vec[counter] = E_rc[r,c];
                    sigma_vec[counter] = sigma_jrc[r,c];
                }
            }
        }
    }


    lor_vec ~ normal(e_rc_vec, sigma_vec);

  for(r in 1:(R - 1)){
      E_rc[r, 1:(C - 1)] ~ normal(0, prior_mu_re_scale);
      // E_rc[r, 1:(C - 1)] ~ normal(to_row_vector(E_rc_prior[r, 1:(C - 1)]), 3);
      // sigma_jrc[r] ~ normal(0, prior_cell_effect_scale);
      // sigma_jrc[r] ~ normal(0, prior_cell_effect_scale);
  }

if(lflag_vary_sd == 2){
    sigma_c_mu ~ normal(0, prior_sigma_c_mu_scale);
    sigma_c_sigma ~ normal(0, prior_sigma_c_scale);
    for(s in 1:(R - 1)*(C-1)){
        sigma_jrc_raw[s] ~ lognormal(sigma_c_mu[1], sigma_c_sigma[1]);
    }
} else {
    sigma_jrc_raw ~ normal(0, prior_sigma_c_scale);
}


  // sigma_r ~ normal(1, .2);
  // sigma_j ~ normal(0.73, 0.1);


    hinge_delta_floor ~ normal(0, .1);
    hinge_delta_min ~ normal(0, .1);

    // counter = 0;
    // for(j in 1:n_areas){
    //   for(c in 1:C){
    //     counter += 1;
    //     cm_vec[counter] = col_margins[j, c];
    //     e_cm_vec[counter] = exp(ALR_jrc_R[j, c] + E_mu + E_j[j]) + sum(cell_values[j, 1:(R - 1), c]);
    //     // col_margins[j,c] ~ poisson(exp(log_sum_exp(to_vector(log_e_cell_values[j,1:R,c]))));
    //     }
    //   }
    //   target += realpoisson_lpdf(cm_vec | e_cm_vec);

    // Prevent scale collapse (Neal's Funnel) and penalize absurd volatility
  // mu_scale ~ std_normal();

  // Standard normal priors on the raw underlying parameters.
  // The sum_to_zero function relies on these being IID N(0,1) to work perfectly.
  // alpha_raw ~ std_normal();
  // log_beta_raw ~ std_normal();
  // gamma_raw ~ std_normal();
  lambda_zero ~ normal(-5.0, 2.0);
  lambda_raw ~ normal(0, 3);



    if(use_known_cells == 1){
      target += poisson_lpmf(known_cell_values_row_vector |
      exp(log_cv_row_vector) + 1e-10);
      }
}
generated quantities{
    #include include/generateratesandsummaries.stan

}
