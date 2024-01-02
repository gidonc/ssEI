// //
// // This Stan Ecological Inference Program
// // Similar to Goodman regression it estimates overall relationship between rows and columns
// // The regression parameters provide estimates of average row to column rates
//
// functions{
//
// }
// data{
//  int<lower=0> n_areas;
//  int<lower=0> R;  // number of rows
//  int<lower=0> C;  // number of columns
//  matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
//  int<lower=0> col_margins[n_areas, C]; // the column margins in each area
// }
// transformed data{
//   real<lower=0> T[n_areas];
//   matrix<lower=0>[n_areas, C] col_prop;
//   matrix<lower=0>[n_areas, R] row_prop;
//
//   for (j in 1:n_areas){
//     T[j] = sum(col_margins[j]);
//     for(c in 1:C){
//       col_prop[j,c] = col_margins[j,c]/T[j];
//     }
//     for(r in 1:R){
//       row_prop[j,r] = row_margins[j,r]/T[j];
//     }
//   }
//
// }
// parameters{
//   matrix<lower=0, upper=1>[R,C] beta;
//   real<lower=0> sigma;
// }
// transformed parameters{
// }
// model{
//   real predC[n_areas*C];
//   real obsC[n_areas*C];
//   for(j in 1:n_areas){
//     for(c in 1:C){
//       predC[(j-1)*C + c] =  dot_product(col(beta, c), row_prop[j]);
//       obsC[(j-1)*C + c] = col_prop[j, c];
//     }
//   }
//   to_vector(beta) ~ normal(.5, 1);
//   sigma ~ normal(0, 3);
//   obsC ~ normal(predC, sigma);
// }
// generated quantities{
//
//
// }
