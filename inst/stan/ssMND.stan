//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// and models row to column rates with: Multinomial-Dirichlet

data{
 int<lower=0> n_areas;
 int<lower=0> R;
 int<lower=0> C;
 matrix<lower=0>[n_areas, R] row_margins;
 matrix<lower=0>[n_areas, C] col_margins;
}
parameters{
 real lambda[n_areas, R - 1, C -1]; // sequential cell weights
 simplex[C] theta[R]; // row to column rates
}
transformed parameters{
 matrix<lower=0>[n_areas * R, C] cell_values;
 real cell_values_flat_plus_detJ[n_areas*R*C + 1];

 real log_det_J;

 cell_values_flat_plus_detJ = assign_cvals(n_areas, R, C, row_margins, col_margins, lambda);

 log_det_J = cell_values_flat_plus_detJ[n_areas*R*C + 1];

  for (j in 1:n_areas){
    for (r in 1:R){
      for (c in 1:C){
       cell_values[(j - 1)*R + r, c] = cell_values_flat_plus_detJ[(j - 1)*R*C + (r-1)*C+c];
     }
   }
  }
}
model{
  matrix[n_areas*R, C] obs_prob;
  int counter = 1;

    for (j in 1:n_areas){
      for (r in 1:R){
        for(c in 1:C){
          obs_prob[counter, c] = theta[r, c]*row_margins[j, r];
        }
      counter += 1;
      }
    }

  target += log_det_J;
  target += extdirmultinom_lpdf(cell_values | obs_prob);
}
