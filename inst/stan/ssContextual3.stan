//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// and models row to column rates with:

functions{
  #include include\allocationfuns.stan
  #include include\realpdf.stan

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
}
transformed data{
  int K = R*(C -1); // size correlation matrix excluding the row margins
  int K2 = (R * C) - 1; // size correlation matrix including the row margins
  int free_R[n_areas];
  int free_C[n_areas];
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
  for (j in 1:n_areas){
    for(r in 1:(R - 1)){
      row_margins_lr[j, r] = log((row_margins[j, r]+ .5)/(row_margins[r, R] + .5));
    }
  }
}
parameters{
   real lambda_unpadded[n_param]; // sequential cell weights
   array[n_areas] vector[K2] alpha;    // logit column probabilities for area j and rows
   vector[K2] mu;              // logit probability row mean
   // corr_matrix[R * (C-1)] Omega;      // correlation matrix
   cholesky_factor_corr[K2] L; // instead of corr_matrix[R * (C-1)] Omega;
   vector<lower=0>[K2] sigma;  // scales
   matrix[R - 1, K] betas;
   // real<lower=0> tau; //for scales
   // real mu_sigma; // for scales
}

transformed parameters{
 real lambda[n_areas, R - 1, C -1]; // sequential cell weights padded with zeros
 real<lower=0> cell_values[n_areas, R, C];
 array[n_areas, R] simplex[C] theta_area; // prob vector for each area
 array[n_areas] simplex[R] theta_rm_area; // prob vector for row margins in each area
 array[n_areas] vector[K2] eta_area;
 matrix[n_areas, K2] mu_area;
 // corr_matrix[K] Omega;
 // cov_matrix[K] Sigma;       // covariance matrix


  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j] - 1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] -1)) + c];
     }
   }
 }


 cell_values = ss_assign_cvals_wzeros_lp(n_areas, R, C, row_margins, col_margins, lambda);

 mu_area = rep_matrix(0, n_areas, K2);

 for (j in 1:n_areas){
    for (k in 1:K){
      mu_area[j,k] = row_margins_lr[j] * col(betas,k);
    }
 }



  for(j in 1:n_areas){
    eta_area[j] = mu + mu_area[j]' + sigma .*(L * alpha[j]);
  }

  for (j in 1:n_areas) {
    for(r in 1:R){
      theta_area[j, r] = simplex_constrain_softmax_lp(eta_area[j,((r-1)*(C -1) + 1):((r-1)*(C -1) + C -1) ]);
    }
    theta_rm_area[j] = simplex_constrain_softmax_lp(eta_area[j, (R*(C - 1) +1):K2]);
  }
}
model{
  matrix[n_areas*R, C] obs_prob;
  matrix[n_areas * R, C] cell_values_matrix;
  matrix[n_areas, R] rm_prob;

  int counter = 1;

  for (j in 1:n_areas){
    alpha[j] ~ std_normal();
  }

  for(k in 1:K){
      col(betas, k) ~ normal(0, 3);
  }

  mu ~ normal(0, 5);      // vectorized, diffuse
  // mu_sigma ~ normal(0, 3);
  L ~ lkj_corr_cholesky(2.0); // implies L*L'~ lkj_corr(2.0);
  // Omega ~ lkj_corr(2.0);  // regularize to unit correlation
  // sigma ~ cauchy(0, 5);   // half-Cauchy due to constraint
  sigma ~ normal(0, 3);
  // tau ~ normal(0, 3);   // half-normal due to constraint - implemented instead of Cauchy to see if that helps with divergent transitions
  // sigma ~ lognormal(mu_sigma, tau);
    for (j in 1:n_areas){
      rm_prob[j] = theta_rm_area[j]';
      for (r in 1:R){
        for(c in 1:C){
          cell_values_matrix[(j - 1)*R + r, c] = cell_values[j, r, c];
          // obs_prob[counter, c] = theta_area[j, r, c]*row_margins[j, r];
          obs_prob[counter, c] = theta_area[j, r, c];
        }
      counter += 1;
      }
    }

  target += realmultinom_lpdf(cell_values_matrix | obs_prob);
  target += realmultinom_lpdf(row_margins| rm_prob);
}
generated quantities{
  matrix[K2, K2] Omega;  // Correlation matrix
  matrix[K2, K2] Sigma;  // Covariance matrix

  #include include\generateratesandsummaries.stan

  Omega = L * L';
  for (m in 1:K2) {
    Sigma[m, m] = sigma[m] * sigma[m] * Omega[m, m];
  }
  for (m in 1:(K2-1)) {
    for (n in (m+1):K2) {
      Sigma[m, n] = sigma[m] * sigma[n] * Omega[m, n];
      Sigma[n, m] = Sigma[m, n];
    }
  }

}

