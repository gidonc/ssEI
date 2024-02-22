//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// with various models of row to column rate

functions{
  #include include/allocationfuns.stan
  #include include/realpdf.stan
  #include include/lkjonionfun.stan

  vector simplex_constrain_softmax_lp(vector v) {
     int K = size(v) + 1;
     vector[K] v0 = append_row(0, v);
     return softmax(v0);
  }

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
 int<lower = 0, upper = 1> lflag_llmod_const; // flag indicating whether log-linear model should be (0) unconstrained (1) constrained to generate correct area totals
 int<lower = 0, upper = 1> lflag_llmod_omit_jr; // flag indicating whether log-linear model should omit area * row interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jc; // flag indicating whether log-linear model should omit area * col interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jrc; // flag indicating whether log-linear model should omit area * row * column interaction
  int<lower = 0, upper =1> lflag_centred_j; // flag indicating whether to use (1) centred or (0) decentred parameterization
  int<lower = 0, upper =1> lflag_centred_r; // flag indicating whether to use (1) centred or (0) decentred parameterization
  int<lower = 0, upper =1> lflag_centred_c; // flag indicating whether to use (1) centred or (0) decentred parameterization
  int<lower = 0, upper =1> lflag_centred_rc; // flag indicating whether to use (1) centred or (0) decentred parameterization
  int<lower = 0, upper =1> lflag_centred_jr; // flag indicating whether to use (1) centred or (0) decentred parameterization
  int<lower = 0, upper =1> lflag_centred_jc; // flag indicating whether to use (1) centred or (0) decentred parameterization
  int<lower = 0, upper =1> lflag_centred_jrc; // flag indicating whether to use (1) centred or (0) decentred parameterization
  int<lower = 0, upper =1> lflag_centred_m; // flag indicating whether to use (1) centred or (0) decentred parameterization

 int<lower = 0, upper =2> lflag_llmod_structure; // flag indicating whether log-linear model should be (0) saturated (1) omit area cell effect  (areaxrowxcolumn interaction) (2) omit area*column interaction (but include areaxrowXcolumn effect) [replaced by other flags]
 real<lower=0> prior_lkj; // lkj param
 real<lower=0> prior_mu_re_scale; // prior for scale of mu_re (mean row effect)
 real<lower=0> prior_mu_ce_scale; // prior for scale of col_effect (mean column effect)
 real<lower=0> prior_sigma_c_scale; //prior for scale of sigma_c (or sigma_c_sigma if lflag_vary_sd == 2)
 real<lower=0> prior_sigma_c_mu_scale; //prior for scale of sigma_c_mu (only if lflag_vary_sd == 2)
 real<lower=0> prior_sigma_ce_scale; //prior of scale for sigma_ce
 real<lower=0> prior_sigma_re_scale; //prior of scale for sigma_re
 real<lower=0> prior_cell_effect_scale; //prior of scale for average cell effects
}
transformed data{
  int K;
  int K_t;
  int K_no_rm;
  int K_c;
  int K_ame;
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
  int has_phi;
  int has_area_re;
  int has_free_area;
  int has_area_col_effects;
  int has_area_row_effects;
  int has_L;
  int has_onion;
  int has_L_ame;
  // int mu_re_ce_in_cell_effects  = 1;
  real param_map[n_areas, R - 1, C - 1];
  matrix[n_areas, R - 1] row_margins_lr;
  real<lower=0> prior_phi_scale;
  int n_margin_sigmas;
  int n_cell_sigmas;
  int n_table_sigmas;
  real sigma_constrain = .001;
  vector[C - 1] mu_col_effects;
  vector[R - 1] mu_row_effects;

  prior_phi_scale =3;

  if(lflag_dist==2){
    has_phi = 1;
  } else{
    has_phi = 0;
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

  if(lflag_llmod_const == 0){
    has_free_area = 1;
  } else {
    has_free_area = 0;
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

  for (j in 1:n_areas){
    for(r in 1:(R - 1)){
      row_margins_lr[j, r] = log((row_margins[j, r]+ .5)/(sum(row_margins[r, (r + 1):R]) + .5));
    }
  }

  has_L_ame = 0;

  n_margin_sigmas = K_ame;
  int n_row_sigmas = has_area_row_effects * (R - 1);
  int n_col_sigmas = has_area_col_effects * (C - 1);
  n_cell_sigmas = has_area_cell_effects * ((lflag_vary_sd ==0) ? 1 : (R - 1) * (C - 1));
  n_table_sigmas = n_margin_sigmas + n_cell_sigmas;
  int K_all = n_table_sigmas;
  array[2] vector[K_all - 1] shapes = create_shapes(K_all, prior_lkj);

  for(c in 1:C - 1){
    mu_row_effects[c] = log(sum(col(col_margins, c))) - log(sum(col(col_margins, C)));
  }

  for(r in 1:R - 1){
    mu_col_effects[r] = log(sum(col(row_margins, r))) - log(sum(col(row_margins, R)));
  }


}
parameters{
 real lambda_unpadded[n_param]; // sequential cell weights
  matrix[R - 1,  C - 1] L_rc_raw;
  array[n_areas]  matrix[R - 1,  C - 1] L_jrc_raw;
  array[n_areas] vector[C - 1] L_jc_raw;
  array[n_areas] vector[R - 1] L_jr_raw;
  vector[n_areas] L_j;
  vector[R * C - 1]log_mu_cv_raw;
  real<lower=0> sigma_lcv_r[(R - 1) * C];
  real<lower=0> sigma_lcv_c[R * (C - 1)];
  real<lower=0> sigma_jrc;
  real<lower=0> sigma_j;
  real mu_j;
 real<lower=0, upper = .001> hinge_delta_floor;
 real<lower=0, upper = .001> hinge_delta_min;
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];
  array[n_areas] matrix[R, C] log_e_cell_value;
  matrix[R, C] log_mu_cv;
  matrix[R, C] log_mu_cv_r;
  matrix[R, C] log_mu_cv_c;
  matrix[R,  C] L_rc = rep_matrix(0, R, C);
  array[n_areas] matrix[R,  C] L_jrc = rep_array(rep_matrix(0, R, C), n_areas);
  array[n_areas] vector[C] L_jc = rep_array(rep_vector(0, C), n_areas);
  array[n_areas] vector[R] L_jr = rep_array(rep_vector(0, R), n_areas);




  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
     }
   }
 }


  cell_values = ss_assign_cvals_wzeros_hinge_lp(n_areas, R, C, row_margins, col_margins, lambda, hinge_delta_floor, hinge_delta_min);

  for(r in 1:R){
    for(c in 1:C){
      if(c == C && r == R){
        log_mu_cv[r, c] = 0;
      } else{
        log_mu_cv[r, c] = log_mu_cv_raw[c +(r - 1) * C];
      }
    }
  }
  for(r in 1:R){
    for(c in 1:C){
      log_mu_cv_r[r, c] = log_mu_cv[r, c] - log_mu_cv[R, c];
      log_mu_cv_c[r, c] = log_mu_cv[r, c] - log_mu_cv[r, C];
    }
  }


  L_rc[1:R - 1, 1: C - 1] = L_rc_raw;
  for(j in 1:n_areas){
    if(lflag_centred_jrc == 1){
      L_jrc[j, 1:R - 1, 1:C - 1] = L_jrc_raw[j];
    } else if(lflag_centred_jrc ==0){
      for(r in 1:R - 1){
        for(c in 1:C - 1){
          L_jrc[j, r, c] = L_rc[r, c] + sigma_jrc * L_jrc_raw[j, r, c];
        }
      }
    }
    L_jc[j, 1: C - 1] = L_jc_raw[j];
    if(lflag_centred_jr == 1){
      L_jr[j, 1: R - 1] = L_jr_raw[j];
    } else if(lflag_centred_jr ==0) {
      for (r in 1:R - 1){
        for(c in 1:C){
          L_jr[j, r] += exp(-log(C)) * (log_mu_cv_r[r, c] - L_jrc[j, r, c] + L_jr_raw[j, r] * sigma_lcv_r[(r - 1)*C + c]);
        }
      }
    }
    if(lflag_centred_jc == 1){
      L_jc[j, 1: C - 1] = L_jc_raw[j];
    } else if(lflag_centred_jc ==0) {
      for (c in 1:C - 1){
        for(r in 1:R){
          L_jc[j, c] += exp(-log(R)) * (log_mu_cv_c[r, c] - L_jrc[j, r, c] + L_jc_raw[j, c] * sigma_lcv_c[(r - 1)*(C - 1) + c]);
        }
      }
    }
  }


  for(j in 1:n_areas){
    for(r in 1:R){
      for(c in 1:C){
        log_e_cell_value[j, r, c] = L_j[j] + L_jr[j, r] + L_jc[j, c] + L_jrc[j, r, c];
      }
    }
  }

}
model{
  vector[n_poss_cells] e_cell;
  row_vector[n_poss_cells] cell_values_row_vector;

  int counter_cell=0;

  for (j in 1:n_areas){
     for (c in 1:C){
       for(r in 1:R){
         if(structural_zeros[j,r,c]==0){
           counter_cell += 1;
           cell_values_row_vector[counter_cell] = cell_values[j,r,c];
           e_cell[counter_cell] = exp(log_e_cell_value[j, r, c]);
         }
       }
     }
  }

    target +=realpoisson_lpdf(cell_values_row_vector| e_cell);

  for(j in 1:n_areas){

    if(lflag_centred_jrc==1){
      to_vector(L_jrc_raw[j])~normal(to_vector(L_rc_raw), sigma_jrc);
    } else if(lflag_centred_jrc == 0){
      to_vector(L_jrc_raw[j]) ~ std_normal();
    }
    if(lflag_centred_jr==1){
      for(r in 1:R - 1){
        vector[C] ljr_comp = rep_vector(-log(C), C);
        for(c in 1:C){
          ljr_comp[c] += normal_lpdf(L_jr_raw[j, r]| log_mu_cv_r[r, c] - L_jrc[j, r, c], sigma_lcv_r[(r - 1)*C + c]);
        }
        target += log_sum_exp(ljr_comp);
        }
    } else if(lflag_centred_jr==0){
        L_jr_raw[j] ~ std_normal();
    }

    if(lflag_centred_jc==1){
      for(c in 1:C - 1){
        vector[R] ljc_comp = rep_vector(-log(R), R);
        for(r in 1:R){
          ljc_comp[r] = normal_lpdf(L_jc_raw[j, c]| log_mu_cv_c[r, c] - L_jrc[j, r, c], sigma_lcv_c[(r - 1)*(C - 1) + c]);
          }
          target += log_sum_exp(ljc_comp);
        }
     } else if(lflag_centred_jc==0){
       L_jc_raw[j] ~ std_normal();
     }


  }

  L_j ~ normal(mu_j, sigma_j);

  mu_j ~ normal(0, 3);
  sigma_j ~ cauchy(0, 3);
  to_vector(L_rc_raw) ~ normal(0, prior_mu_re_scale);
  to_vector(log_mu_cv) ~ normal(0, prior_mu_re_scale);
  sigma_jrc ~ cauchy(0, prior_cell_effect_scale);
  sigma_lcv_r ~ cauchy(0, prior_mu_re_scale);
  sigma_lcv_c ~ cauchy(0, prior_mu_re_scale);


    hinge_delta_floor ~ normal(0, .001);
    hinge_delta_min ~ normal(0, .001);

}
generated quantities{
  #include include/generateratesandsummaries.stan

}
