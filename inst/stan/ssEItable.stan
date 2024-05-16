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
 matrix<lower=0>[R, C] overall_tots; // the overall total target distribution (if known - used only if lflag_has_target==1)
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
  int<lower = 0, upper =1> lflag_predictors_cm; // flag indicating whether to model columns as well as rows

 int<lower = 0, upper =3> lflag_llmod_structure; // flag indicating whether log-linear model should be (0) saturated (1) omit area cell effect  (areaxrowxcolumn interaction) (2) omit area*column interaction (but include areaxrowXcolumn effect) [replaced by other flags]
 int<lower = 0, upper =1> lflag_has_overall_tots; // flag indicating whether the overall target distribution is known;
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
  int has_theta;
  int has_area_re;
  int has_free_E_j;
  int has_area_col_effects;
  int has_area_row_effects;
  int has_L;
  int has_onion;
  int has_L_ame;
  // int mu_re_ce_in_cell_effects  = 1;
  real param_map[n_areas, R - 1, C - 1];
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

  prior_phi_scale =3;

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

  if(lflag_llmod_const == 0){
    has_free_E_j = 1;
  } else {
    has_free_E_j = 0;
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

  K_j = lflag_predictors_cm ==1 ? R*C : C - 1 + (R - 1)*(C - 1) ;
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

  has_L_ame = 0;

  n_margin_sigmas = K_ame;
  int n_row_sigmas = has_area_row_effects * (R - 1);
  int n_col_sigmas = has_area_col_effects * (C - 1);
  n_jrc_sigmas = has_area_cell_effects * ((lflag_vary_sd ==0) ? 1 : (R - 1) * (C - 1));
  n_table_sigmas = n_margin_sigmas + n_jrc_sigmas;
  int K_all = n_table_sigmas;
  array[2] vector[K_j - 1] shapes = create_shapes(K_j, prior_lkj);

}
parameters{
  real lambda_unpadded[n_param]; // sequential cell weights
  array[n_areas] vector[K_j] E_j_all_raw;
  vector<lower=0, upper= 1> [n_areas*has_theta] theta;
  // vector[has_free_E_j*n_areas] E_j_raw;
  array[(lflag_predictors_cm == 0 ? 1: 0) * n_areas] vector[R - 1] E_jr_raw;
  vector[(lflag_predictors_cm == 0 ? 1: 0) * n_areas] E_j_raw;
  vector[K_j] E_mu_all_raw;
  vector<lower=0>[K_j] sigma_j_all;
  real<lower=0> scale_sigma_j_all;
  real sigma_j_all_mu;
  cholesky_factor_corr[has_L * K_j] L_Omega_raw;
  row_vector[has_onion * (choose(K_j, 2) - 1)] l; // do NOT init with 0 for all elements
  vector<lower = 0, upper = 1>[has_onion * (K_j - 1)] R2; // first element is not really a R^2 but is on (0,1)
  real<lower=0, upper = .001> hinge_delta_floor;
  real<lower=0, upper = .001> hinge_delta_min;
  real E_mu_mu;
  real<lower=0> E_mu_sigma;
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];
  matrix[K_j, K_j] L_Omega;



  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
     }
   }
 }


  cell_values = ss_assign_cvals_wzeros_hinge_lp(n_areas, R, C, row_margins, col_margins, lambda, hinge_delta_floor, hinge_delta_min);


  if(has_onion == 1){
    L_Omega = lkj_onion(K_j, l, R2, shapes); // cholesky_factor corr matrix
  } else if(has_L == 1){
    L_Omega = L_Omega_raw; // cholesky_factor corr matrix
  } else {
    L_Omega = identity_matrix(K_j);
  }

  #include include/generateratesandsummaries.stan

}
model{
  array[n_areas] vector[K_j] E_j_all;
  array[n_areas] vector[K_j] E_mu_all_j;
  array[n_areas] vector[R] E_jr = rep_array(rep_vector(0, R), n_areas);
  array[n_areas] vector[C] E_jc = rep_array(rep_vector(0, C), n_areas);
  vector[n_areas] E_j;
  array[n_areas] matrix[R, C] E_jrc = rep_array(rep_matrix(0, R, C), n_areas);
  array[n_areas] matrix[R, C] log_e_cell_values =  rep_array(rep_matrix(0, R, C), n_areas);
  array[n_areas] matrix[R, C] log_r_cell_values =  rep_array(rep_matrix(0, R, C), n_areas);
  real add_mu;
  vector[K_j] E_mu_all;

  vector[n_poss_cells] e_cell;
  vector[R*C] e_overall_tots = rep_vector(0, R*C);
  row_vector[R*C] overall_tots_row_vector = rep_row_vector(0, R*C);

  row_vector[n_poss_cells] cell_values_row_vector;
  vector[n_poss_cells*has_theta] e_theta;
  int counter_cell = 0;

  if(lflag_centred_m==1){
    E_mu_all = E_mu_all_raw;
  } else if(lflag_centred_m==0) {
    E_mu_all = E_mu_mu + E_mu_sigma*E_mu_all_raw;
  }


  for(j in 1:n_areas){

      if(lflag_centred_jrc == 1){
        E_j_all[j] = E_j_all_raw[j];
      } else if(lflag_centred_jrc == 0) {
        E_j_all[j] = E_mu_all + diag_pre_multiply(sigma_j_all, L_Omega) * E_j_all_raw[j];
      }


      for(r in 1:R-1){
        for(c in 1:C){
          E_jrc[j, r, c] = E_j_all[j, R * (r - 1) + c];
        }
      }
      for(c in 1:C - 1){
        E_jrc[j, R, c] = E_j_all[j, R*(R - 1) + c];
      }
      E_jrc[j, R, C] = -sum(E_jrc[j, 1:(R - 1), 1:C]) - sum(E_jrc[j, R, 1:C - 1]);

      E_j[j] = E_j_all[j, R*C];

      E_mu_all_j[j] = E_mu_all;




     for(r in 1:R){
       for(c in 1:C){
         if(structural_zeros[j,r, c] == 0){
           log_r_cell_values[j, r, c] = E_jrc[j, r, c] + E_jc[j, c] + E_jr[j, r];
         } else{
           log_r_cell_values[j, r, c] = -200;
         }
       }
     }
     // if(has_free_E_j == 1){
     //   E_j[j] = E_j_raw[j];
     // } else {
     //   E_j[j] = tot_log[j] -log_sum_exp(to_vector(log_r_cell_values[j]));
     // }
     for(r in 1:R){
       for (c in 1:C){
         log_e_cell_values[j, r, c] = log_r_cell_values[j, r, c] + E_j[j];
       }
     }


  }


  for (j in 1:n_areas){
    for(r in 1:R){
      for(c in 1:C){
        if(structural_zeros[j,r,c]==0){
          counter_cell += 1;
          if(r<R && c<C){
            cell_values_row_vector[counter_cell] = cell_values[j, r, c];
            if(lflag_dist ==0){
                e_cell[counter_cell]  = exp(log_e_cell_values[j, r, c]);
            } else if(lflag_dist==2){
              e_cell[counter_cell]  = exp(log_e_cell_values[j, r, c]) *(1 - theta[j])/theta[j];
              e_theta[counter_cell] = theta[j];
            }
          } else if (r==R && c ==C){
            cell_values_row_vector[counter_cell] = cell_values[j, r, c];
            if(lflag_dist==0){
              e_cell[counter_cell] = exp(log_e_cell_values[j, r, c]);
            } else if(lflag_dist==2){
              e_cell[counter_cell] = exp(log_e_cell_values[j, r, c])*(1 - theta[j])/theta[j];
              e_theta[counter_cell] = theta[j];
            }

          } else if(r == R && c < C){
            cell_values_row_vector[counter_cell] = cell_values[j, r, c];
            if(lflag_dist==0){
              e_cell[counter_cell] = exp(log_e_cell_values[j, r, c]);
            } else if(lflag_dist==2){
              e_cell[counter_cell] = exp(log_e_cell_values[j, r, c])*(1 - theta[j])/theta[j];
              e_theta[counter_cell] = theta[j];
            }
          } else if(c == C && r < R){
            cell_values_row_vector[counter_cell] = cell_values[j, r, c];
            if(lflag_dist==0){
              e_cell[counter_cell] = exp(log_e_cell_values[j, r, c]);
            } else if(lflag_dist==2){
              e_cell[counter_cell] = exp(log_e_cell_values[j, r, c])*(1 - theta[j])/theta[j];
              e_theta[counter_cell] = theta[j];
            }
          }
        }
      }
    }
   }

   if(lflag_dist==0){
     target +=realpoisson_lpdf(cell_values_row_vector| e_cell);
   } else if (lflag_dist ==2){
     target +=realnegbinom3_lpdf(cell_values_row_vector| e_cell, e_theta);
   }

   if(lflag_has_overall_tots==1){
     for(r in 1:R){
       for(c in 1:C){
         e_overall_tots[c + (r -1)* C] = sum(exp(log_e_cell_values[1:n_areas, r, c]));
         overall_tots_row_vector[c + (r - 1)*C] = overall_tots[r, c];
       }
     }
     target +=realpoisson_lpdf(overall_tots_row_vector| e_overall_tots);
   }


  // cell_values_row_vector ~ normal(e_cell, e_sigma);

  for(j in 1:n_areas){


    if(lflag_centred_jrc==1){
      E_j_all_raw[j] ~ multi_normal_cholesky(E_mu_all_j[j], diag_pre_multiply(sigma_j_all, L_Omega));
    } else{
      E_j_all_raw[j] ~ std_normal();
    }
  }


if(lflag_centred_m==1){
    E_mu_all_raw ~ normal(E_mu_mu, E_mu_sigma);
} else{
   E_mu_all_raw ~ std_normal();
}

E_mu_mu ~ cauchy(0, prior_sigma_c_scale);
E_mu_sigma ~ cauchy(0, prior_cell_effect_scale);

  sigma_j_all ~ lognormal(sigma_j_all_mu, scale_sigma_j_all);
  // sigma_j_all ~ cauchy(0, 2);
  // for(k in 1:K_j){
  //
    // if(k < K_jc_start + 1){
    //   sigma_j_all[k] ~ normal(0, prior_cell_effect_scale);
    // } else if(k<K_jrc_start){
    //   sigma_j_all[k] ~ normal(0, prior_cell_effect_scale);
    // } else if(k< R *C){
    //   sigma_j_all[k] ~ lognormal(sigma_j_all_mu, scale_sigma_j_all);
    // } else{
    //   sigma_j_all[k] ~ normal(0, prior_cell_effect_scale);
    // }
  // }
  sigma_j_all_mu ~ normal(0, prior_cell_effect_scale);
  scale_sigma_j_all ~ normal(0, prior_cell_effect_scale);


  if(has_L == 1){
    L_Omega_raw ~ lkj_corr_cholesky(prior_lkj); // implies L*L'~ lkj_corr(prior_lkj);
  }
  if(has_onion == 1){
    l ~ std_normal();
    R2 ~ beta(shapes[1], shapes[2]);
  }
  // if(has_free_E_j == 1){
  //   E_j_raw ~ normal(0, 10);
  // }

  if(lflag_predictors_cm == 0){
    E_j_raw ~ normal(0, 10);
    for(j in 1:n_areas){
      E_jr_raw[j] ~ normal(0, 10);
    }
  }

    hinge_delta_floor ~ normal(0, .001);
    hinge_delta_min ~ normal(0, .001);

}
generated quantities{


}
