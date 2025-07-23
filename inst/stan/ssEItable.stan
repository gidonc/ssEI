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

  vector sum_to_zero(vector y){
    int N = size(y);
    int n;
    vector[N] omega;
    vector[N] S;
    vector[N + 1] x;

    for(n_rev in 1:(N - 1)){
      n = N - n_rev;
      omega[n + 1] = y[n + 1]/sqrt((n + 1)*(n + 2));
    }
  S[N] = 0;
  for(n_rev in 1:(N - 1)){
    n = N - n_rev;
    S[n] = S[n + 1] + omega[n + 1];
  }
  for(n_rev in 1:N){
    n = N + 1 - n_rev;
    x[n + 1] = S[n] - (n * y[n]/sqrt(n * (n + 1)));
  }
  x[1] = S[1] + y[1]/sqrt(2);
  return(x);
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
 int<lower = 0, upper = 1> lflag_llmod_omit_jr; // flag indicating whether log-linear model should omit area * row interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jc; // flag indicating whether log-linear model should omit area * col interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jrc; // flag indicating whether log-linear model should omit area * row * column interaction
  int<lower = 0, upper =1> lflag_predictors_cm; // flag indicating whether to model columns as well as rows
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
}
parameters{
  real lambda_unpadded[n_param]; // sequential cell weights
  array[n_areas] vector[R*C - 1] E_jrc_raw;
  // array[n_areas] vector[R - 1] E_jr_raw;
  // array[n_areas] vector[C - 1] E_jc_raw;
  vector[n_areas - 1] E_j_raw;
  // vector[C - 1] E_c_raw;
  // vector[R - 1] E_r_raw;
  vector[R * C - 1] E_rc_raw;

  vector<lower=0, upper= 1> [n_areas*has_theta] theta;
  real<lower=0> sigma_j;
  vector<lower=0>[R * C] sigma_jrc;
  // vector<lower=0>[R] sigma_r;
  // vector<lower=0>[C] sigma_c;

  // vector<lower=0>[K_j] sigma_j_all;
  real<lower=0, upper = .001> hinge_delta_floor;
  real<lower=0, upper = .001> hinge_delta_min;
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];
  vector[n_areas] E_j;
  // array[n_areas] vector[R] E_jr = rep_array(rep_vector(0, R), n_areas);
  // array[n_areas] vector[C] E_jc = rep_array(rep_vector(0, C), n_areas);
  array[n_areas] matrix[R, C] E_jrc = rep_array(rep_matrix(0, R, C), n_areas);
  array[n_areas] matrix[R, C] log_e_cell_values =  rep_array(rep_matrix(0, R, C), n_areas);
  // array[n_areas] matrix[R, C] log_r_cell_values =  rep_array(rep_matrix(0, R, C), n_areas);
  vector[R * C] E_rc;
  // vector[R] E_r;
  // vector[C] E_c;


  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
     }
   }
 }


  cell_values = ss_assign_cvals_wzeros_hinge_lp(n_areas, R, C, row_margins, col_margins, lambda, hinge_delta_floor, hinge_delta_min);


  // E_j = sum_to_zero(E_j_raw);
  E_j[1:n_areas - 1] = E_j_raw;
  E_j[n_areas] = 0;
  // E_r = sum_to_zero(E_r_raw);
  // E_c = sum_to_zero(E_c_raw);
  E_rc = sum_to_zero(E_rc_raw);

  for(j in 1:n_areas){
    vector[R * C] E_jrc_vector;
    // E_jr[j, 1:R] = sum_to_zero(E_jr_raw[j, 1:R - 1]);
    // E_jc[j, 1:C] = sum_to_zero(E_jc_raw[j, 1:C - 1]);
    E_jrc_vector[1:R * C] = sum_to_zero(E_jrc_raw[j, 1: R * C - 1]);
    for(c in 1:C){
      E_jrc[j, 1:R, c] = E_jrc_vector[C * (c - 1) + 1:(C *c)];
    }

     for(r in 1:R){
       for(c in 1:C){
         if(structural_zeros[j,r, c] == 0){
           log_e_cell_values[j, r, c] = E_jrc[j, r, c] + E_j[j];
         } else{
           log_e_cell_values[j, r, c] = -200;
         }
       }
     }
  }



}
model{
  vector[n_poss_cells] e_cell;
  row_vector[n_poss_cells] cell_values_row_vector;
  vector[n_poss_cells*has_theta] e_theta;
  int counter_cell = 0;

  for (j in 1:n_areas){
    for(r in 1:R){
      for(c in 1:C){
        if(structural_zeros[j,r,c]==0){
          counter_cell += 1;
          cell_values_row_vector[counter_cell] = cell_values[j, r, c];
          if(lflag_dist ==0){
            e_cell[counter_cell]  = exp(log_e_cell_values[j, r, c]);
          } else if(lflag_dist==2){
            e_cell[counter_cell]  = exp(log_e_cell_values[j, r, c]) *(1 - theta[j])/theta[j];
            e_theta[counter_cell] = theta[j];
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


  // cell_values_row_vector ~ normal(e_cell, e_sigma);

  for(j in 1:n_areas){
    to_vector(E_jrc[j]) ~ normal(E_rc, sigma_jrc);
    // E_jr[j] ~ normal(E_r, sigma_r);
    // E_jc[j] ~ normal(E_c, sigma_c);
  }
  E_j_raw ~ normal(0, sigma_j);

  E_rc ~ normal(0, prior_cell_effect_scale);
  // E_r ~ normal(0, prior_cell_effect_scale);
  // E_c ~ normal(0, prior_cell_effect_scale);

  sigma_jrc ~ normal(0, prior_cell_effect_scale);
  // sigma_r ~ normal(0, prior_cell_effect_scale);
  // sigma_c ~ normal(0, prior_cell_effect_scale);


    hinge_delta_floor ~ normal(0, .001);
    hinge_delta_min ~ normal(0, .001);

}
generated quantities{
    #include include/generateratesandsummaries.stan

}
