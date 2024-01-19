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
 real<lower=0> prior_lkj; // lkj param
 int<lower=0, upper=2> lflag_dist; // flag indicating whether to use poisson (0), multinomial (1) or negative binomial (2) paramertization
 int<lower=0, upper=3> lflag_area_re; // flag indicating whether the area mean simplex is uniform (0) or varies with area random effects which are normally distributed (1) or varies with area random effects which are multinormally distributed (non centred paramaterisation) (2) or varies with area random effects which are multinormally distributed (non centred LKJ Onion paramaterisation)
 int<lower  =0, upper=2> lflag_vary_sd; // flag indicating whether variance of area_cell parameters is shared across cells (0) varies by cell (1) or has a hierarchical model structure (2)
}
transformed data{
  int K;
  int K_t;
  int K_no_rm;
  int K_c;

  K_t = (R - 1)* (C - 1);

  K_c = C * (R - 1);
  K = R * (C - 1);
  K_no_rm = R * (C - 1);

  // if(lflag_inc_rm == 1){
  //   K = (R * C) - 1;
  // } else {
  //   K = R * (C - 1);
  // }
  // K_no_rm = R * (C - 1);

  array[2] vector[K -1] shapes = create_shapes(K_t, prior_lkj);
  int free_R[n_areas];
  int free_C[n_areas];
  int non0_rm;
  int non0_cm;
  int phi_length;
  int has_area_re;
  int has_L;
  int has_onion;
  real param_map[n_areas, R - 1, C - 1];
  matrix[n_areas, R - 1] row_margins_lr;

  // calculate the number of free parameters (zero row and columns do not need a parameter to allocated cell value of 0)

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
  if(lflag_dist==2){
    phi_length = 1;
  } else{
    phi_length = 0;
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




}
parameters{
 real lambda_unpadded[n_param]; // sequential cell weights
 // vector[lflag_mod_cols * K_c] mu_c;  // logit probability row mean to construct theta

 // matrix[has_area_re*n_areas, has_area_re * R * (C - 1)] alpha; // area variation from mu (logit probability row mean)
 // array[has_area_re*n_areas] vector[has_area_re*K] alpha;    // logit column probabilities for area j and rows
 array[n_areas] matrix[(R - 1), (C - 1)] area_cell_effect_raw;    // logit column probabilities for area j and cols
 real area_effect[n_areas];
 matrix[n_areas, (R - 1)] area_row_effect_raw;
 matrix[n_areas, (C - 1)] area_col_effect_raw;
 matrix[(R - 1), (C - 1)] cell_effect_raw;
 vector[C - 1] mu_ce;
 vector[R - 1] mu_re;
 vector<lower=0>[C - 1] sigma_ce;
 vector<lower=0>[R - 1] sigma_re;

 // vector<lower=0>[(R - 1) * (C - 1)] sigma_c; //scale of area col variation from mu_c
 real<lower=0> sigma_c_raw[(lflag_vary_sd ==0) ? 1 : (R - 1) * (C - 1)]; //scale of area col variation from mu_c
 real<lower=0> sigma_c_sigma[(lflag_vary_sd ==2) ? 1 : 0]; //hierarchical standard deviatation on standard deviations
 vector[(lflag_vary_sd ==2) ? 1 : 0] sigma_c_mu; //hierarchical mean on standard deviations
 // cholesky_factor_corr[has_L * K] L_a; // for modelling correlation matrix
 // row_vector[has_onion * (choose(K, 2) - 1)] l; // do NOT init with 0 for all elements
 // vector<lower = 0, upper = 1>[has_onion * (K - 1)] R2; // first element is not really a R^2 but is on (0,1)
 // matrix[lflag_predictors_rm * K_no_rm, lflag_predictors_rm * (R - 1)] betas_rm;

}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];
  // array[n_areas, R] simplex[C] theta_area; // prob row vector for each area
  // array[n_areas, C] simplex[R] theta_c_area; // prob col vector for each area
  real log_e_cell_value[n_areas, R, C];
  matrix[n_areas, R] area_row_effect;
  matrix[n_areas, C] area_col_effect;
  matrix[R, C] cell_effect;
  array[n_areas] matrix[R, C] area_cell_effect;    // logit column probabilities for area j and cols
  real<lower=0> sigma_c[(R - 1) * (C - 1)]; //scale of area col variation from mu_c


  // matrix[max(has_onion, has_L) * K, max(has_onion, has_L) * K] L;

  // if(has_onion == 1){
  //   L = lkj_onion(K, l, R2, shapes); // cholesky_factor corr matrix
  // }
  // if(has_L == 1){
  //   L = L_a; // cholesky_factor corr matrix
  // }

  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
     }
   }
 }


  cell_values = ss_assign_cvals_wzeros_lp(n_areas, R, C, row_margins, col_margins, lambda);


  cell_effect = rep_matrix(0.0, R, C);
  cell_effect[1:(R - 1), 1:(C - 1)] = cell_effect_raw;

  for(j in 1:n_areas){
    area_cell_effect[j] = rep_matrix(0.0, R, C);
    area_cell_effect[j, 1:(R - 1), 1:(C - 1)] = area_cell_effect_raw[j];

    area_col_effect[j, 1:(C - 1)] = area_col_effect_raw[j];
    area_col_effect[j, C] = 0;

    area_row_effect[j, 1:(R - 1)] = area_row_effect_raw[j];
    area_row_effect[j, R] = 0;

    for(r in 1:R){
        for(c in 1:C){
          log_e_cell_value[j, r, c] = area_effect[j] + area_row_effect[j, r] + area_col_effect[j, c] + cell_effect[r, c] + area_cell_effect[j, r, c];
        }
      }
    }

    if(lflag_vary_sd == 0){
      sigma_c = rep_array(sigma_c_raw[1], (R - 1)*(C - 1));
    } else {
      sigma_c = sigma_c_raw;
    }



}
model{
  matrix[non0_cm, R] cell_values_matrix_c;
  int counter_c = 1;
  matrix[non0_cm, R] obs_prob_c;

//
//   if(has_L == 1){
//     L ~ lkj_corr_cholesky(prior_lkj); // implies L*L'~ lkj_corr(prior_lkj);
//   }
//   if(has_onion == 1){
//     l ~ std_normal();
//     R2 ~ beta(shapes[1], shapes[2]);
//   }





  for (j in 1:n_areas){
     for (c in 1:C){
       if(col_margins[j, c] > 0){
         for(r in 1:R){
           cell_values_matrix_c[counter_c, r] = cell_values[j, r, c];
           // obs_prob_c[counter_c, r] = theta_c_area[j, c, r] * col_margins[j, c];
           obs_prob_c[counter_c, r] = exp(log_e_cell_value[j, r, c]);
         }
         counter_c += 1;
        }
      }
    }
  // if(lflag_dist == 1){
  //   target += realmultinom_lpdf(cell_values_matrix | obs_prob);
  // } else if(lflag_dist == 0) {
  //   target += realpoisson_lpdf(to_row_vector(cell_values_matrix) | to_vector(obs_prob));
  // } else if(lflag_dist == 2) {
  //   phi ~ cauchy(0, 3);
  //   target += realnegbinom3_lpdf(to_row_vector(cell_values_matrix) | to_vector(obs_prob), to_vector(phi_matrix));
  // }
  // if(lflag_inc_rm){
  //     target += realmultinom_lpdf(row_margins| rm_prob);
  // }
    target +=realpoisson_lpdf(to_row_vector(cell_values_matrix_c)| to_vector(obs_prob_c));
    // sigma_c ~ normal(0, 5);
    if(lflag_vary_sd == 2){
      sigma_c_mu ~ normal(0, 3);
      sigma_c_sigma ~ normal(0, 5);
      for(s in 1:K_t){
        sigma_c_raw[s]~lognormal(sigma_c_mu, sigma_c_sigma);
      }
    } else {
      sigma_c_raw ~ normal(0, 5);
    }
    mu_ce~ normal(0, 5);
    mu_re~ normal(0, 5);
    sigma_ce~normal(0, 5);
    sigma_re~normal(0, 5);
    to_vector(cell_effect_raw) ~ normal(0, 5);
    for (j in 1:n_areas){
      to_vector(area_cell_effect_raw[j]) ~ normal(0, sigma_c);
      area_row_effect_raw[j] ~ normal(mu_re, sigma_re);
      area_col_effect_raw[j] ~ normal(mu_ce, sigma_ce);
    }
}
generated quantities{

  #include include/generateratesandsummaries.stan

}
