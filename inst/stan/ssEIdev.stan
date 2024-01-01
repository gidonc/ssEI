//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// with various models of row to column rate

functions{
  #include include\allocationfuns.stan
  #include include\realpdf.stan
  #include include\lkjonionfun.stan

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
 real<lower=0> lkj_param; // lkj param
 int<lower=0, upper=2> lflag_mn; // flag indicating whether to use poisson (0), multinomial (1) or negative binomial (2) paramertization
 int<lower=0, upper=2> lflag_area_re; // flag indicating whether the area mean simplex is uniform (0) or varies with area random effects which are normally distributed (1) or varies with area random effects which are multinormally distributed (non centred paramaterisation) (2)
}
transformed data{
  int K = R*(C -1);
  array[2] vector[K -1] shapes = create_shapes(K, lkj_param);
  int free_R[n_areas];
  int free_C[n_areas];
  int non0_rm;
  int phi_length;
  int has_area_re;
  int has_L;
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

  for (j in 1:n_areas){
    for(r in 1:(R - 1)){
      row_margins_lr[j, r] = log((row_margins[j, r]+ .5)/(sum(row_margins[r, (r + 1):R]) + .5));
    }
  }
  if(lflag_mn==2){
    phi_length = 1;
  } else{
    phi_length = 0;
  }

  if(lflag_area_re == 0){
    has_area_re = 0;
  } else {
    has_area_re = 1;
  }

  if(lflag_area_re == 3){
    has_L = 0;
  } else {
    has_L = 1;
  }



}
parameters{
 real lambda_unpadded[n_param]; // sequential cell weights
 // simplex[C] theta[R]; // average row to column rates
 vector[R * (C-1)] mu;  // logit probability row mean to construct theta
 vector<lower=0, upper = 1>[phi_length* R] phi;

 // matrix[has_area_re*n_areas, has_area_re * R * (C - 1)] alpha; // area variation from mu (logit probability row mean)
 array[has_area_re*n_areas] vector[has_area_re*K] alpha;    // logit column probabilities for area j and rows
 vector<lower=0>[has_area_re * R * (C - 1)] sigma; // scale of area variation from mu
 cholesky_factor_corr[has_L * R * (C-1)] L; // for modelling correlation matrix
 // matrix[R - 1, K] betas;
 // row_vector[choose(K, 2) - 1] l; // do NOT init with 0 for all elements
 // vector<lower = 0, upper = 1>[K - 1] R2; // first element is not really a R^2 but is on (0,1)
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];
  array[n_areas, R] simplex[C] theta_area; // prob vector for each area
  array[n_areas] vector[K] eta_area;
  // matrix[n_areas, K] mu_area;
  // vector[n_areas, K] mu_area;
  // matrix[K, K] L = lkj_onion(K, l, R2, shapes); // cholesky_factor corr matrix


  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
     }
   }
 }


  cell_values = ss_assign_cvals_wzeros_lp(n_areas, R, C, row_margins, col_margins, lambda);


 // for (j in 1:n_areas){
 //    for (k in 1:K){
 //      mu_area[j,k] = row_margins_lr[j] * col(betas,k);
 //    }
 // }



  // for(j in 1:n_areas){
  //   eta_area[j] = mu + mu_area[j]' + sigma .*(alpha[j]');
  // }
  for(j in 1:n_areas){
    if(lflag_area_re == 0){
      eta_area[j] = mu;
    } else if(lflag_area_re == 1){
        eta_area[j] = mu + sigma .*(alpha[j]);
    } else if(lflag_area_re == 2){
        eta_area[j] = mu + sigma .*(L * alpha[j]);
    }
  }


  for (j in 1:n_areas) {
    for(r in 1:R){
      theta_area[j, r] = simplex_constrain_softmax_lp(eta_area[j,((r-1)*(C -1) + 1):((r-1)*(C -1) + C -1) ]);
    }
  }

}
model{
  matrix[non0_rm, C] obs_prob;
  matrix[non0_rm, C] cell_values_matrix;
  matrix[phi_length*non0_rm, phi_length*C] phi_matrix;
  int counter = 0;
  sigma ~ normal(0, 3);

  if(has_L == 1){
    L ~ lkj_corr_cholesky(lkj_param); // implies L*L'~ lkj_corr(lkj_param);
  }



  mu ~ normal(0, 5);      // vectorized, diffuse
  for (j in 1:n_areas){
    alpha[j] ~ std_normal();
  }
//
//   for(k in 1:K){
//       col(betas, k) ~ normal(0, 3);
//   }
//

  for (j in 1:n_areas){
    for (r in 1:R){
      if(row_margins[j, r] > 0){
        counter += 1;
        for (c in 1:C){
          cell_values_matrix[counter, c] = cell_values[j, r, c];
          if(lflag_mn ==0) {
            obs_prob[counter, c] = theta_area[j, r, c]*row_margins[j,r];
          }else if(lflag_mn== 1){
            obs_prob[counter, c] = theta_area[j, r, c];
          } else  if (lflag_mn==2){
            phi_matrix[counter, c] = phi[r];
            obs_prob[counter, c] = theta_area[j, r, c]*row_margins[j, r] - (theta_area[j,r,c]*row_margins[j, r] * phi[r]);
          }
      }

    }
   }
  }
  if(lflag_mn == 1){
    target += realmultinom_lpdf(cell_values_matrix | obs_prob);
  } else if(lflag_mn == 0) {
    target += realpoisson_lpdf(to_row_vector(cell_values_matrix) | to_vector(obs_prob));
  } else if(lflag_mn == 2) {
    phi ~ cauchy(0, 3);
    target += realnegbinom3_lpdf(to_row_vector(cell_values_matrix) | to_vector(obs_prob), to_vector(phi_matrix));
  }
}
generated quantities{

  #include include\generateratesandsummaries.stan

}
