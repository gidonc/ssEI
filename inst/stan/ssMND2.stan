//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// and models row to column rates with: Multinomial-Dirichlet

functions{
  #include include/allocationfuns.stan
  #include include/realpdf.stan
  #include include/lkjonionfun.stan

}
data{
 int<lower=0> n_areas;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
 matrix<lower=0>[n_areas, C] col_margins; // the column margins in each area
}
transformed data{
  int K = R*(C -1);
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
 simplex[C] theta[R]; // average row to column rates
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];

  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
     }
   }
 }


  cell_values = ss_assign_cvals_wzeros_lp(n_areas, R, C, row_margins, col_margins, lambda);
}
model{
  matrix[n_areas * R, C] obs_prob;
  matrix[n_areas * R, C] cell_values_matrix;

  for (j in 1:n_areas){
    for (r in 1:R){
      for (c in 1:C){
       cell_values_matrix[(j - 1)*R + r, c] = cell_values[j, r, c];
       obs_prob[(j - 1)*R + r, c] = theta[r, c]*row_margins[j, r];
     }
   }
  }
  // theta ~ dirichlet(theta_prior);
  target += realmultinom_lpdf(cell_values_matrix | obs_prob);
}
generated quantities{

  #include include/generateratesandsummaries.stan

}
