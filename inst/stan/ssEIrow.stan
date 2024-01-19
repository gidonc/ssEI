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
 int<lower=0, upper=1> lflag_inc_rm; //flag indicating whether to include row margins log-ratios in correlation matrix
 int<lower=0, upper=1> lflag_predictors_rm; //flag indicating whether to have row margin log-ratios as predictors of row-to-column (unconstrained) simplex
 int<lower=0, upper=2> lflag_vary_sd; // flag indicating whether variance of area_cell parameters is shared across cells (0) varies by cell (1) or has a hierarchical model structure (2)
}
transformed data{
  int K;
  int K_no_rm;

  if(lflag_inc_rm == 1){
    K = (R * C) - 1;
  } else {
    K = R * (C - 1);
  }
  K_no_rm = R * (C - 1);

  array[2] vector[K -1] shapes = create_shapes(K, prior_lkj);
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
 // simplex[C] theta[R]; // average row to column rates
 vector[K] mu;  // logit probability row mean to construct theta
 vector<lower=0, upper = 1>[phi_length* R] phi;

 // matrix[has_area_re*n_areas, has_area_re * R * (C - 1)] alpha; // area variation from mu (logit probability row mean)
 array[has_area_re*n_areas] vector[has_area_re*K] alpha;    // logit column probabilities for area j and rows
 vector<lower=0>[has_area_re * (lflag_vary_sd == 0 ? 1 : K)] sigma_raw; // scale of area row variation from mu
 real<lower=0> sigma_sigma[has_area_re *( (lflag_vary_sd ==2) ? 1 : 0 )]; //hierarchical standard deviatation on standard deviations
 vector[has_area_re * ((lflag_vary_sd ==2) ? 1 : 0)] sigma_mu; //hierarchical mean on standard deviations
 cholesky_factor_corr[has_L * K] L_a; // for modelling correlation matrix
 row_vector[has_onion * (choose(K, 2) - 1)] l; // do NOT init with 0 for all elements
 vector<lower = 0, upper = 1>[has_onion * (K - 1)] R2; // first element is not really a R^2 but is on (0,1)
 matrix[lflag_predictors_rm * K_no_rm, lflag_predictors_rm * (R - 1)] betas_rm;

}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];
  array[n_areas, R] simplex[C] theta_area; // prob row vector for each area
  array[n_areas] vector[K] eta_area;
  matrix[max(has_onion, has_L) * K, max(has_onion, has_L) * K] L;
  array[lflag_inc_rm * n_areas] simplex[lflag_inc_rm * R] theta_rm_area; // prob vector for row margins in each area
  matrix[n_areas, K] mu_area_rm;
  vector<lower=0>[has_area_re * K] sigma;

  if(has_area_re ==1){
    if(lflag_vary_sd == 0) {
      sigma = rep_vector(sigma_raw[1], K);
    } else{
      sigma = sigma_raw;
    }

  }

  if(has_onion == 1){
    L = lkj_onion(K, l, R2, shapes); // cholesky_factor corr matrix
  }
  if(has_L == 1){
    L = L_a; // cholesky_factor corr matrix
  }

  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
     }
   }
 }


  cell_values = ss_assign_cvals_wzeros_lp(n_areas, R, C, row_margins, col_margins, lambda);


  mu_area_rm = rep_matrix(0, n_areas, K);

  if(lflag_predictors_rm == 1){
      for (j in 1:n_areas){
       for (k in 1:K){
         mu_area_rm[j,k] = row_margins_lr[j] * betas_rm[k]';
       }
      }
  }




  for(j in 1:n_areas){
    if(lflag_area_re == 0){
      eta_area[j] = mu + mu_area_rm[j]';
    } else if(lflag_area_re == 1){
        eta_area[j] = mu + mu_area_rm[j]' + sigma .*(alpha[j]);
    } else if(lflag_area_re == 2 || lflag_area_re ==3){
        eta_area[j] = mu + mu_area_rm[j]' + sigma .*(L * alpha[j]);
    }
  }

  for (j in 1:n_areas) {
    for(r in 1:R){
      theta_area[j, r] = simplex_constrain_softmax_lp(eta_area[j,((r-1)*(C -1) + 1):((r-1)*(C -1) + C -1) ]);
    }
    if(lflag_inc_rm == 1){
       theta_rm_area[j] = simplex_constrain_softmax_lp(eta_area[j, (R*(C - 1) + 1):K]);
    }
  }

}
model{
  matrix[non0_rm, C] obs_prob;
  matrix[non0_rm, C] cell_values_matrix;
  matrix[phi_length*non0_rm, phi_length*C] phi_matrix;
  int counter = 1;
  matrix[lflag_inc_rm * n_areas, lflag_inc_rm * R] rm_prob;
  if(has_area_re==1){
    if(lflag_vary_sd == 2){
      sigma_mu ~ normal(0, 3);
      sigma_sigma ~ normal(0, 5);
      for(s in 1:K){
        sigma_raw[s]~lognormal(sigma_mu[1], sigma_sigma);
      }
    } else {
      sigma_raw ~ normal(0, 5);
    }
    sigma ~ normal(0, 3);
  }


  if(has_L == 1){
    L ~ lkj_corr_cholesky(prior_lkj); // implies L*L'~ lkj_corr(prior_lkj);
  }
  if(has_onion == 1){
    l ~ std_normal();
    R2 ~ beta(shapes[1], shapes[2]);
  }



  mu ~ normal(0, 5);      // vectorized, diffuse
  if(has_area_re==1){
      for (j in 1:n_areas){
        alpha[j] ~ std_normal();
      }
  }

  if(lflag_predictors_rm==1){
    for(k in 1:K_no_rm){
      betas_rm[k] ~ normal(0, 3);
    }
  }



  for (j in 1:n_areas){
    if(lflag_inc_rm == 1){
        rm_prob[j] = theta_rm_area[j]';
    }
    for (r in 1:R){
      if(row_margins[j, r] > 0){
        for (c in 1:C){
          cell_values_matrix[counter, c] = cell_values[j, r, c];
          if(lflag_dist ==0) {
            obs_prob[counter, c] = theta_area[j, r, c]*row_margins[j,r];
          }else if(lflag_dist== 1){
            obs_prob[counter, c] = theta_area[j, r, c];
          } else  if (lflag_dist==2){
            phi_matrix[counter, c] = phi[r];
            obs_prob[counter, c] = theta_area[j, r, c]*row_margins[j, r] - (theta_area[j,r,c]*row_margins[j, r] * phi[r]);
          }
      }
      counter += 1;
    }
   }
  }
  if(lflag_dist == 1){
    target += realmultinom_lpdf(cell_values_matrix | obs_prob);
  } else if(lflag_dist == 0) {
    target += realpoisson_lpdf(to_row_vector(cell_values_matrix) | to_vector(obs_prob));
  } else if(lflag_dist == 2) {
    phi ~ cauchy(0, 3);
    target += realnegbinom3_lpdf(to_row_vector(cell_values_matrix) | to_vector(obs_prob), to_vector(phi_matrix));
  }
  if(lflag_inc_rm){
      target += realmultinom_lpdf(row_margins| rm_prob);
  }
}
generated quantities{

  #include include/generateratesandsummaries.stan

}
