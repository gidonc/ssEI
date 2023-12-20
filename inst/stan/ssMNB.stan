//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// and models row to column rates with: Multinomial-Dirichlet

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
}
transformed data{
  int K = R*(C -1);
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
 // simplex[C] theta[R]; // average row to column rates
 vector[R * (C-1)] mu;              // logit probability row mean to construct theta
 matrix[n_areas, R * (C - 1)] alpha; // area variation from mu (logit probability row mean)
 vector<lower=0>[R * (C - 1)] sigma; // scale of area variation from mu

 // row_vector[choose(K, 2) - 1] l; // do NOT init with 0 for all elements
 // vector<lower = 0, upper = 1>[K - 1] R2; // first element is not really a R^2 but is on (0,1)
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];
  array[n_areas, R] simplex[C] theta_area; // prob vector for each area
  array[n_areas] vector[K] eta_area;
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

  for(j in 1:n_areas){
    eta_area[j] = mu + sigma .*(alpha[j]');
  }

  for (j in 1:n_areas) {
    for(r in 1:R){
      theta_area[j, r] = simplex_constrain_softmax_lp(eta_area[j,((r-1)*(C -1) + 1):((r-1)*(C -1) + C -1) ]);
    }
  }

}
model{
  matrix[n_areas * R, C] obs_prob;
  matrix[n_areas * R, C] cell_values_matrix;
  sigma ~ normal(0, 3);

  mu ~ normal(0, 5);      // vectorized, diffuse
  to_vector(alpha) ~ std_normal();

  for (j in 1:n_areas){
    for (r in 1:R){
      for (c in 1:C){
       cell_values_matrix[(j - 1)*R + r, c] = cell_values[j, r, c];
       obs_prob[(j - 1)*R + r, c] = theta_area[j, r, c];
     }
   }
  }
  target += realmultinom_lpdf(cell_values_matrix | obs_prob);
}
generated quantities{

  #include include\generateratesandsummaries.stan

}
