//
// This Stan Ecological Inference Program
// Similar to Goodman regression it estimates overall relationship between rows and columns
// The regression parameters provide estimates of average row to column rates

functions{

}
data{
 int<lower=0> n_areas;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
 int<lower=0> col_margins[n_areas, C]; // the column margins in each area
}
transformed data{

}
parameters{
  matrix<lower=0, upper=1>[R,C] beta;
}
transformed parameters{
}
model{
  real predC[n_areas*C];
  int obsC[n_areas*C];
  for(j in 1:n_areas){
    for(c in 1:C){
      predC[(j-1)*C + c] =  dot_product(col(beta, c), row_margins[j]);
      obsC[(j-1)*C + c] = col_margins[j, c];
    }
  }
  to_vector(beta) ~ normal(.5, 1);
  obsC ~ poisson(predC);
}
generated quantities{


}
