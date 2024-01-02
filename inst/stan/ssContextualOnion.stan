//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// and models row to column rates with:
// implementing the onion method to avoid problems with larger correlation matricies
// onion method implementation described here: https://discourse.mc-stan.org/t/lkj-corr-lpdf-correlation-matrix-is-not-positive-definite/31995/5

functions{
  #include include/allocationfuns.stan
  #include include/realpdf.stan

  vector simplex_constrain_softmax_lp(vector v) {
     int K = size(v) + 1;
     vector[K] v0 = append_row(0, v);
     return softmax(v0);
  }

  array[] vector create_shapes(int K, real eta) {
    real alpha = eta + (K - 2) / 2.0;
    array[2] vector[K - 1] shape;

    shape[1, 1] = alpha;
    shape[2, 1] = alpha;
    for (k in 2 : (K - 1)) {
      alpha -= 0.5;
      shape[1, k] = k / 2.0;
      shape[2, k] = alpha;
    }

    return shape;
  }

  matrix lkj_onion(int K, row_vector l, vector R2, data array[] vector shape) {
    matrix[K, K] L = rep_matrix(0, K, K); // cholesky_factor corr matrix
    {
      int start = 1;
      int end = 2;

      L[1, 1] = 1.0;
      L[2, 1] = 2.0 * R2[1] - 1.0;
      L[2, 2] = sqrt(1.0 - square(L[2, 1]));
      for (k in 2 : (K - 1)) {
        int kp1 = k + 1;
        row_vector[k] l_row = segment(l, start, k);
        real scale = sqrt(R2[k] / dot_self(l_row));
        L[kp1, 1 : k] = l_row[1 : k] * scale;
        L[kp1, kp1] = sqrt(1.0 - R2[k]);
        start = end + 1;
        end = start + k - 1;
      }
    }
    return L;
  }

}
data{
 int<lower=0> n_areas;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
 matrix<lower=0>[n_areas, C] col_margins; // the column margins in each area
 real<lower=0> lkj_param; // lkj param
}
transformed data{
  int K = R*(C - 1);
  array[2] vector[K -1] shapes = create_shapes(K, lkj_param);
  int free_R[n_areas];
  int free_C[n_areas];
  real param_map[n_areas, R - 1, C - 1];
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
}
parameters{
   real lambda_unpadded[n_param]; // sequential cell weights
   array[n_areas] vector[K] alpha;    // logit column probabilities for area j and rows
   vector[R * (C-1)] mu;              // logit probability row mean
   // corr_matrix[R * (C-1)] Omega;      // correlation matrix
   // cholesky_factor_corr[K] L; // instead of corr_matrix[R * (C-1)] Omega;
   vector<lower=0>[K] sigma;  // scales
   row_vector[choose(K, 2) - 1] l; // do NOT init with 0 for all elements
   vector<lower = 0, upper = 1>[K - 1] R2; // first element is not really a R^2 but is on (0,1)
}

transformed parameters{
 real lambda[n_areas, R - 1, C -1]; // sequential cell weights padded with zeros
 real<lower=0> cell_values[n_areas, R, C];
 array[n_areas, R] simplex[C] theta_area; // prob vector for each area
 array[n_areas] vector[K] eta_area;
 matrix[K, K] L = lkj_onion(K, l, R2, shapes); // cholesky_factor corr matrix

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

  for(j in 1:n_areas){
    eta_area[j] = mu + sigma .*(L * alpha[j]);
    // eta_area[j] = mu + (L * alpha[j]);
  }

  for (j in 1:n_areas) {
    for(r in 1:R){
      theta_area[j, r] = simplex_constrain_softmax_lp(eta_area[j,((r-1)*(C -1) + 1):((r-1)*(C -1) + C -1) ]);
    }
  }
}
model{
  matrix[n_areas*R, C] obs_prob;
  matrix[n_areas * R, C] cell_values_matrix;

  int counter = 1;

  l ~ std_normal();
  R2 ~ beta(shapes[1], shapes[2]);

  for (j in 1:n_areas){
    alpha[j] ~ std_normal();
  }

  mu ~ normal(0, 5);      // vectorized, diffuse
  // L ~ lkj_corr_cholesky(2.0); // implies L*L'~ lkj_corr(2.0);
  // Omega ~ lkj_corr(2.0);  // regularize to unit correlation
  // sigma ~ cauchy(0, 5);   // half-Cauchy due to constraint
  sigma ~ normal(0, 3);   // half-normal due to constraint - implemented instead of Cauchy to see if that helps with divergent transitions
    for (j in 1:n_areas){
      for (r in 1:R){
        for(c in 1:C){
          cell_values_matrix[(j - 1)*R + r, c] = cell_values[j, r, c];
          obs_prob[counter, c] = theta_area[j, r, c]*row_margins[j, r];
          // obs_prob[counter, c] = theta_area[j, r, c];
        }
      counter += 1;
      }
    }

  target += realdirmultinom_lpdf(cell_values_matrix | obs_prob);
}
generated quantities{
  matrix[K, K] Omega= multiply_lower_tri_self_transpose(L);  // Correlation matrix
  matrix[K, K] Sigma;  // Covariance matrix

  #include include/generateratesandsummaries.stan

  for (m in 1:K) {
    Sigma[m, m] = sigma[m] * sigma[m] * Omega[m, m];
  }
  for (m in 1:(K-1)) {
    for (n in (m+1):K) {
      Sigma[m, n] = sigma[m] * sigma[n] * Omega[m, n];
      Sigma[n, m] = Sigma[m, n];
    }
  }


}

