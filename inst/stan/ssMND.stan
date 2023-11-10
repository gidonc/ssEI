//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// and models row to column rates with: Multinomial-Dirichlet

functions{
  #include constraintfunctions.stanfunctions
  #include realpdf.stanfunctions
}
data{
 int<lower=0> n_areas;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
 matrix<lower=0>[n_areas, C] col_margins; // the column margins in each area
}
parameters{
 real lambda[n_areas, R - 1, C -1]; // sequential cell weights
 simplex[C] theta[R]; // average row to column rates
}
transformed parameters{
  matrix<lower=0>[n_areas, R, C] cell_values;

  cell_values = ss_assign_cvals_lp(n_areas, R, C, row_margins, col_margins, lambda);
}
model{
  matrix[n_areas*R, C] obs_prob;
  matrix<lower=0>[n_areas * R, C] cell_values_matrix;

  for (j in 1:n_areas){
    for (r in 1:R){
      for (c in 1:C){
       cell_values_matrix[(j - 1)*R + r, c] = cell_values[j, r, c];
       obs_prob[(j - 1)*R + r, c] = theta[r, c]*row_margins[j, r];
     }
   }
  }
  target += realdirmultinom_lpdf(cell_values_matrix | obs_prob);
}
